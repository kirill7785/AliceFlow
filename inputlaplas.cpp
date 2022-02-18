// файл inputlaplas.cpp
// содержит ввод структур данных  
// построенных графическим редактором.

#pragma once
#ifndef _INPUTLAPLAS_CPP_
#define _INPUTLAPLAS_CPP_ 1

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>  // математические функции
#include <fstream> // Для работы с файлом
#include <iostream>
#include <string>

// недостающий функционал 1 октября 2016.
// эта процедура определена выше по коду в mysolverv0_03.c, поэтому здесь её дефениция излишна.
integer my_imin(integer ia, integer ib)
{
	return (ia < ib ? ia : ib);
} // imin

// Возвращает максимум из двух целых чисел.
integer my_imax(integer ia, integer ib) 
{
	return (ib > ia ? ib : ia);
} // my_imax

// недостающий функционал 1 октября 2016.
// эта процедура определена выше по коду в mysolverv0_03.c, поэтому здесь её дефениция излишна.
doublereal my_imin3(doublereal ra, doublereal rb, doublereal rc)
{
	return std::min(ra, std::min(rb, rc));
} // my_imin3

// Возвращает максимум из двух целых чисел.
doublereal my_imax3(doublereal ra, doublereal rb, doublereal rc)
{
	return std::max(ra, std::max(rb, rc));
} // my_imax3

//#include <string.h> // функции обработки строк
#include "power_temperature_depend.cpp" // сплайновая интерполяция таблично заданной функции (зависимость мощности от температуры)

// передаваемый из файла параметр отвечающий за подробность расчётной сетки.
doublereal etalon_max_size_ratio = 2.0;// Более реалистичный подход 10.0
// параметр отвечающий за качество расчётной сетки.
// 30.0 оптимум по документации FlowVision.
doublereal etalon_max_size_ratio2 = 30.0; // Более реалистичный подход в том чтобы его не ограничивать.
// Принебрежимо малая длина для упрощения 13.08.2019
// Для оптимального быстродействия надо очень аккуратно настраивать.
// 13.08.2019 настраивается автоматом. Здесь начальное значение которое нельзя использовать
// и которое настраивает автомат.
// Количество клеток разбивающих интервал.
// Если увеличить более четырёх то быстродействие теряется.
const doublereal rdivision_interval = 4.0; // 4.0 default 10.0
const doublereal rdivision_intervalCAD = 7.0; // 7.0 default 10.0
// Если bFULL_AUTOMATIC то допуски определяются локально на местах,
// что даёт значительно более лучшую настройку и гибкость управления в автоматическом режиме.
const bool bFULL_AUTOMATIC = true;
doublereal shorter_length_for_simplificationX_BASIC = 1.0e30; 
doublereal shorter_length_for_simplificationY_BASIC = 1.0e30;
doublereal shorter_length_for_simplificationZ_BASIC = 1.0e30;
// 0.	none
// 1.	Snap to grid
// 2.	Snap to grid ALICE
// 3.	Snap to grid ++
integer bsnap_TO_global = 1; // snap to grid
// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
// при повторном обращении идет обращение к ячейки хеш-таблицы.
// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
// Хеш-таблицы для automatic
const integer isize_shorter_hash = 10000000;// Сделан значительный запас точности хеширования.
doublereal* shorter_hash_X = nullptr;
doublereal* shorter_hash_Y = nullptr;
doublereal* shorter_hash_Z = nullptr;
bool* bshorter_hash_X = nullptr;
bool* bshorter_hash_Y = nullptr;
bool* bshorter_hash_Z = nullptr;
// Параметры настройки 
// Блазиус 2мкм.
const doublereal dcabinet_delta_multiplyer = 0.005; // 1.0e-6; // 1.0e-12;// 0.005; // 1/200
const doublereal dbody_dist_multiplyer = 12.5;
const doublereal dbody_R_out_multiplyer = 2.6;
const doublereal dbody_Hcyl_multiplyer = 10.0;
const doublereal dbody_LGATE_situation_multiplyer = 5.2;


// плоскости
const unsigned char XY_PLANE = 1;
const unsigned char XZ_PLANE = 2;
const unsigned char YZ_PLANE = 3;

// тип блока
//const unsigned char SOLID = 1;
//const unsigned char HOLLOW = 2;
//const unsigned char FLUID = 3;
enum class PHYSICS_TYPE_IN_BODY {
	SOLID = 1,
	HOLLOW = 2,
	FLUID = 3
};

// Геометрическая форма блока
const unsigned char PRISM = 0;
const unsigned char CYLINDER = 1;
const unsigned char POLYGON = 2;
const unsigned char CAD_STL = 3; // 13.11.2020




// точка в трёхмерном пространстве
typedef struct TPOinteger {
	doublereal x, y, z;
	TPOinteger() {
		x=0.0, y=0.0, z=0.0;
	}
} TOCHKA;


// точка в трёхмерном пространстве
typedef struct TPOinteger_FLOAT {
	float x, y, z;
	TPOinteger_FLOAT() {
		x = 0.0, y = 0.0, z = 0.0;
	}
} TOCHKA_FLOAT;

// CAD STL кубик состоит из набора треугольников.
// Тип одного треугольника 13.11.2020
typedef struct TgCAD_STL {
	TOCHKA_FLOAT n;
	TOCHKA_FLOAT pa, pb, pc;
	TgCAD_STL* next;

	TgCAD_STL() {
		next = nullptr;
	}

} gCAD_STL;

// Для больших .stl файлов бы разделяем stl модель
// на части для ускорения обработки. 17.11.2020 
typedef struct TbigCAD_STL {

	
	// Окаймляющий сектор прямоугольник данной	части stl файла.
	doublereal min1, min2, max1, max2;

	// Треугольники входящие в окаймляющий сектор прямоугольник.
	gCAD_STL* root_CAD_STL;

	// Односвязный список указатель на следующего.
	TbigCAD_STL* next;

	// индекс для связи с послединим элементом в списке, для того чтобы не проматываеть список каждый раз.
	// или -1 если список пуст.
	int id_endl;

	TbigCAD_STL() {

		root_CAD_STL = nullptr;
		next = nullptr;

		id_endl = -1;
	}

} bigCAD_STL;

// только для enumerate_volume_improved
// и для uniformsimplemeshgen.cpp
// См. также заголовочный файл constr_struct.cpp
typedef struct TBlock_indexes {
	// Именно int64_t так как число контрольных объёмов может быть большим.
	integer iL, iR, jL, jR, kL, kR;

	TBlock_indexes() {
		iL=-1, iR=-2, jL=-1, jR=-2, kL=-1, kR=-2;
	}
} Block_indexes;

// 16.11.2020
double dabs(double d) {
	if (d < 0.0) return -d;
	return d;
}

// Пересечение отрезка и треугольника: 
// 13.11.2020
bool trnLineFacet(TOCHKA p1, TOCHKA p2, TOCHKA pa, TOCHKA pb, TOCHKA pc, TOCHKA& p)
{
	double d;
	double a1, a2, a3;
	double total, denom, mu;
	TOCHKA n, pa1, pa2, pa3;
	const doublereal EPS1 = 0.01; //0.01; // 1.0e-7;
	const doublereal EPS = 5.0e-4; //1.0e-7;//  1.0e-9;

	/* Calculate the parameters for the plane */
	n.x = (pb.y - pa.y) * (pc.z - pa.z) - (pb.z - pa.z) * (pc.y - pa.y);
	n.y = (pb.z - pa.z) * (pc.x - pa.x) - (pb.x - pa.x) * (pc.z - pa.z);
	n.z = (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
	// Normalize vector n.
	doublereal length_n = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
	n.x = n.x / length_n;
	n.y = n.y / length_n;
	n.z = n.z / length_n;
	d = -n.x * pa.x - n.y * pa.y - n.z * pa.z;

	/* Calculate the position on the line that intersects the plane */
	denom = n.x * (p2.x - p1.x) + n.y * (p2.y - p1.y) + n.z * (p2.z - p1.z);
	if (dabs(denom) < EPS1)         /* Line and plane don't intersect */
		return(false);
	mu = -(d + n.x * p1.x + n.y * p1.y + n.z * p1.z) / denom;
	p.x = static_cast <float>(p1.x + mu * (p2.x - p1.x));
	p.y = static_cast <float>(p1.y + mu * (p2.y - p1.y));
	p.z = static_cast <float>(p1.z + mu * (p2.z - p1.z));
	if (mu < 0.0 || mu > 1.0)   /* Intersection not along line segment */
		return(false);

	/* Determine whether or not the intersection point is bounded by pa,pb,pc */
	pa1.x = pa.x - p.x;
	pa1.y = pa.y - p.y;
	pa1.z = pa.z - p.z;
	doublereal length_pa1 = sqrt(pa1.x * pa1.x + pa1.y * pa1.y + pa1.z * pa1.z);
	pa1.x = pa1.x / length_pa1;
	pa1.y = pa1.y / length_pa1;
	pa1.z = pa1.z / length_pa1;

	pa2.x = pb.x - p.x;
	pa2.y = pb.y - p.y;
	pa2.z = pb.z - p.z;
	doublereal length_pa2 = sqrt(pa2.x * pa2.x + pa2.y * pa2.y + pa2.z * pa2.z);
	pa2.x = pa2.x / length_pa2;
	pa2.y = pa2.y / length_pa2;
	pa2.z = pa2.z / length_pa2;

	pa3.x = pc.x - p.x;
	pa3.y = pc.y - p.y;
	pa3.z = pc.z - p.z;
	doublereal length_pa3 = sqrt(pa3.x * pa3.x + pa3.y * pa3.y + pa3.z * pa3.z);
	pa3.x = pa3.x / length_pa3;
	pa3.y = pa3.y / length_pa3;
	pa3.z = pa3.z / length_pa3;

	a1 = pa1.x * pa2.x + pa1.y * pa2.y + pa1.z * pa2.z;
	a2 = pa2.x * pa3.x + pa2.y * pa3.y + pa2.z * pa3.z;
	a3 = pa3.x * pa1.x + pa3.y * pa1.y + pa3.z * pa1.z;
	doublereal RTOD = 180.0 / 3.14159265358979323846;
	total = (acos(a1) + acos(a2) + acos(a3)) * RTOD;
	if (dabs(total - 360.0) > EPS)
		return(false);

	return(true);
}



// Пересечение отрезка и треугольника: 
// 13.11.2020
bool trnLineFacet(TOCHKA_FLOAT p1, TOCHKA_FLOAT p2, TOCHKA_FLOAT pa, TOCHKA_FLOAT pb, TOCHKA_FLOAT pc, TOCHKA_FLOAT& p)
{
	double d;
	double a1, a2, a3;
	double total, denom, mu;
	TOCHKA n, pa1, pa2, pa3;
	// 1.0e-7 Проверено на BSK_Dmitrii
	const doublereal EPS1 =  0.01; //0.01; // 1.0e-7;
	const doublereal EPS = 1.0e-3;// 1.0e-3; 5.0e-4; //1.0e-7;//  1.0e-9;

	/* Calculate the parameters for the plane */
	n.x = (pb.y - pa.y) * (pc.z - pa.z) - (pb.z - pa.z) * (pc.y - pa.y);
	n.y = (pb.z - pa.z) * (pc.x - pa.x) - (pb.x - pa.x) * (pc.z - pa.z);
	n.z = (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
	// Normalize vector n.
	doublereal length_n = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
	n.x = n.x / length_n;
	n.y = n.y / length_n;
	n.z = n.z / length_n;
	d = -n.x * pa.x - n.y * pa.y - n.z * pa.z;

	/* Calculate the position on the line that intersects the plane */
	denom = n.x * (p2.x - p1.x) + n.y * (p2.y - p1.y) + n.z * (p2.z - p1.z);
	if (dabs(denom / sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z))) < EPS1) {         /* Line and plane don't intersect */
		//printf("1 %e\n", dabs(denom / sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z))));
		//system("pause");
		return(false);
	}
	//std::wcout << "|denom|<EPS |denom|=" << fabs(denom) << std::endl;


	mu = -(d + n.x * p1.x + n.y * p1.y + n.z * p1.z) / denom;
	p.x = static_cast <float>(p1.x + mu * (p2.x - p1.x));
	p.y = static_cast <float>(p1.y + mu * (p2.y - p1.y));
	p.z = static_cast <float>(p1.z + mu * (p2.z - p1.z));
	if (mu < 0.0 || mu > 1.0)   /* Intersection not along line segment */
		return(false);

	//std::wcout << "mu (0..1)=" << mu << std::endl;

	/* Determine whether or not the intersection point is bounded by pa,pb,pc */
	pa1.x = pa.x - p.x;
	pa1.y = pa.y - p.y;
	pa1.z = pa.z - p.z;
	doublereal length_pa1 = sqrt(pa1.x * pa1.x + pa1.y * pa1.y + pa1.z * pa1.z);
	pa1.x = pa1.x / length_pa1;
	pa1.y = pa1.y / length_pa1;
	pa1.z = pa1.z / length_pa1;

	pa2.x = pb.x - p.x;
	pa2.y = pb.y - p.y;
	pa2.z = pb.z - p.z;
	doublereal length_pa2 = sqrt(pa2.x * pa2.x + pa2.y * pa2.y + pa2.z * pa2.z);
	pa2.x = pa2.x / length_pa2;
	pa2.y = pa2.y / length_pa2;
	pa2.z = pa2.z / length_pa2;

	pa3.x = pc.x - p.x;
	pa3.y = pc.y - p.y;
	pa3.z = pc.z - p.z;
	doublereal length_pa3 = sqrt(pa3.x * pa3.x + pa3.y * pa3.y + pa3.z * pa3.z);
	pa3.x = pa3.x / length_pa3;
	pa3.y = pa3.y / length_pa3;
	pa3.z = pa3.z / length_pa3;

	a1 = pa1.x * pa2.x + pa1.y * pa2.y + pa1.z * pa2.z;
	a2 = pa2.x * pa3.x + pa2.y * pa3.y + pa2.z * pa3.z;
	a3 = pa3.x * pa1.x + pa3.y * pa1.y + pa3.z * pa1.z;
	doublereal RTOD = 180.0 / 3.14159265358979323846;
	total = (acos(a1) + acos(a2) + acos(a3)) * RTOD;

	//std::wcout << "360 dopusk=" << dabs(total - 360.0) << std::endl;

	if (dabs(total - 360.0) > EPS) {
		//printf("2 %e\n",dabs(total-360));
		//system("pause");
		return(false);
	}
	

	//system("PAUSE");

	return(true);
}

// геометрическое описание
typedef struct TGEOM {
	

	integer itypegeom; // 0 - Prism, 1 - Cylinder, 2 - Polygon, 3 - CAD_STL
	// Prism
	doublereal xS, yS, zS; // координаты начала объекта
	doublereal xE, yE, zE; // координаты конца объекта
	// Cylinder
	integer iPlane; // плоскость в которой лежит нижнее основание цилиндра.
	doublereal xC, yC, zC, Hcyl, R_out_cyl, R_in_cyl;
	// Polygon
	integer iPlane_obj2; // плоскость в которой лежит нижнее основание полигона.
	integer nsizei; // Количество опорных точек образующих выпуклый полигон.
	doublereal *hi, *xi, *yi, *zi;
	// CAD_STL
	bool bminmax, bgetminedge, bvolcad;
	doublereal mem_get_min_edge, volcadstl;
	TOCHKA_FLOAT pmin_mem, pmax_mem;
	gCAD_STL* root_CAD_STL;

	// 17.11.2020
	// В режиме большой модели большая stl модель 
	// делится на небольшие части которые обрабатываются независимо. 
	bool bbigCADmodel; // Режим большой модели.
	const int size_pattern = 64;//16
	bigCAD_STL*** root_big_CAD_STL_model_x_hash_table;
	bigCAD_STL*** root_big_CAD_STL_model_y_hash_table;
	bigCAD_STL*** root_big_CAD_STL_model_z_hash_table;
	bigCAD_STL* root_big_CAD_STL_model_x;
	bigCAD_STL* root_big_CAD_STL_model_y;
	bigCAD_STL* root_big_CAD_STL_model_z;
	int inumber_triangles_for_CAD_STL_model; 

	char name[80]; // имя объёмного геометрического объекта body.

	TGEOM() {
		name[0] = '\0';

		
		itypegeom = PRISM; // 0 - Prism, 1 - Cylinder, 2 - Polygon, 3 - CAD_STL
		// Prism
		xS=0.0, yS=0.0, zS=0.0; // координаты начала объекта
		xE=0.0, yE=0.0, zE=0.0; // координаты конца объекта
		// Cylinder
		iPlane=-2; // плоскость в которой лежит нижнее основание цилиндра.
		xC=0.0, yC=0.0, zC=0.0, Hcyl=0.0, R_out_cyl=0.0, R_in_cyl=0.0;
		// Polygon
		iPlane_obj2=-2; // плоскость в которой лежит нижнее основание полигона.
		nsizei=-2; // Количество опорных точек образующих выпуклый полигон.
		hi=nullptr; xi=nullptr; yi=nullptr; zi=nullptr;
		// CAD_STL
		root_CAD_STL = nullptr;
		bminmax = false;
		bgetminedge=false;
		bvolcad = false;
		volcadstl = 0.0;

		bbigCADmodel = false;
		root_big_CAD_STL_model_x_hash_table = nullptr;
		root_big_CAD_STL_model_y_hash_table = nullptr;
		root_big_CAD_STL_model_z_hash_table = nullptr;

		root_big_CAD_STL_model_x = nullptr;
		root_big_CAD_STL_model_y = nullptr;
		root_big_CAD_STL_model_z = nullptr;
		inumber_triangles_for_CAD_STL_model = 12;// Кубик.

		mem_get_min_edge = 0.0;// длина минимального ребра в cad модели из stl файла.

	}

	void clear_CAD_STL() {
		gCAD_STL* tmp = root_CAD_STL;
		while (root_CAD_STL != nullptr) {
			root_CAD_STL = root_CAD_STL->next;
			tmp->next = nullptr;
			delete tmp;
			tmp = root_CAD_STL;
		}

		inumber_triangles_for_CAD_STL_model = 12;

		if (root_big_CAD_STL_model_x_hash_table != nullptr) {
			for (int i_1 = 0; i_1 < size_pattern; ++i_1) {
				for (int i_2 = 0; i_2 < size_pattern; ++i_2) {
					root_big_CAD_STL_model_x_hash_table[i_1][i_2] = nullptr;
				}
			}
		}

		if (root_big_CAD_STL_model_y_hash_table != nullptr) {
			for (int i_1 = 0; i_1 < size_pattern; ++i_1) {
				for (int i_2 = 0; i_2 < size_pattern; ++i_2) {
					root_big_CAD_STL_model_y_hash_table[i_1][i_2] = nullptr;
				}
			}
		}

		if (root_big_CAD_STL_model_z_hash_table != nullptr) {
			for (int i_1 = 0; i_1 < size_pattern; ++i_1) {
				for (int i_2 = 0; i_2 < size_pattern; ++i_2) {
					root_big_CAD_STL_model_z_hash_table[i_1][i_2] = nullptr;
				}
			}
		}

		bigCAD_STL* tmp1 = root_big_CAD_STL_model_x;
		while (root_big_CAD_STL_model_x != nullptr) {

			tmp = root_big_CAD_STL_model_x->root_CAD_STL;
			while (root_big_CAD_STL_model_x->root_CAD_STL != nullptr) {
				root_big_CAD_STL_model_x->root_CAD_STL = root_big_CAD_STL_model_x->root_CAD_STL->next;
				tmp->next = nullptr;
				delete tmp;
				tmp = root_big_CAD_STL_model_x->root_CAD_STL;
			}


			root_big_CAD_STL_model_x = root_big_CAD_STL_model_x->next;
			tmp1->next = nullptr;
			delete tmp1;
			tmp1 = root_big_CAD_STL_model_x;
		}

		tmp1 = root_big_CAD_STL_model_y;
		while (root_big_CAD_STL_model_y != nullptr) {

			tmp = root_big_CAD_STL_model_y->root_CAD_STL;
			while (root_big_CAD_STL_model_y->root_CAD_STL != nullptr) {
				root_big_CAD_STL_model_y->root_CAD_STL = root_big_CAD_STL_model_y->root_CAD_STL->next;
				tmp->next = nullptr;
				delete tmp;
				tmp = root_big_CAD_STL_model_y->root_CAD_STL;
			}


			root_big_CAD_STL_model_y = root_big_CAD_STL_model_y->next;
			tmp1->next = nullptr;
			delete tmp1;
			tmp1 = root_big_CAD_STL_model_y;
		}


		tmp1 = root_big_CAD_STL_model_z;
		while (root_big_CAD_STL_model_z != nullptr) {

			tmp = root_big_CAD_STL_model_z->root_CAD_STL;
			while (root_big_CAD_STL_model_z->root_CAD_STL != nullptr) {
				root_big_CAD_STL_model_z->root_CAD_STL = root_big_CAD_STL_model_z->root_CAD_STL->next;
				tmp->next = nullptr;
				delete tmp;
				tmp = root_big_CAD_STL_model_z->root_CAD_STL;
			}


			root_big_CAD_STL_model_z = root_big_CAD_STL_model_z->next;
			tmp1->next = nullptr;
			delete tmp1;
			tmp1 = root_big_CAD_STL_model_z;
		}

	}

	// Вычисляет нормаль к треугольнику.
	void calculateNormal_for_triangle(TOCHKA &n, TOCHKA pa, TOCHKA pb, TOCHKA pc) {

		/* Calculate the parameters for the plane */
		n.x = (pb.y - pa.y) * (pc.z - pa.z) - (pb.z - pa.z) * (pc.y - pa.y);
		n.y = (pb.z - pa.z) * (pc.x - pa.x) - (pb.x - pa.x) * (pc.z - pa.z);
		n.z = (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
		// Normalize vector n.
		doublereal length_n = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
		n.x = n.x / length_n;
		n.y = n.y / length_n;
		n.z = n.z / length_n;
	}

	// Вычисляет нормаль к треугольнику.
	void calculateNormal_for_triangle(TOCHKA_FLOAT& n, TOCHKA_FLOAT pa, TOCHKA_FLOAT pb, TOCHKA_FLOAT pc) {

		/* Calculate the parameters for the plane */
		n.x = (pb.y - pa.y) * (pc.z - pa.z) - (pb.z - pa.z) * (pc.y - pa.y);
		n.y = (pb.z - pa.z) * (pc.x - pa.x) - (pb.x - pa.x) * (pc.z - pa.z);
		n.z = (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
		// Normalize vector n.
		float length_n = sqrtf(n.x * n.x + n.y * n.y + n.z * n.z);
		n.x = n.x / length_n;
		n.y = n.y / length_n;
		n.z = n.z / length_n;
	}


	void CYLINDER2CAD_STL() {

		if ((itypegeom == CYLINDER)&&(fabs(R_in_cyl)<1.0e-30)&&(Hcyl>0.0)) {

			// освобождение оперативной памяти.
			clear_CAD_STL();

			inumber_triangles_for_CAD_STL_model = 200;

			// Модель обычная, небольшая.
			bbigCADmodel = false;
			root_big_CAD_STL_model_x = nullptr;
			root_big_CAD_STL_model_y = nullptr;
			root_big_CAD_STL_model_z = nullptr;


			// Внешняя поверхность цилиндра в плоскости XZ,
			// Состоящая из 200 треугольников. Каждый треугольник 
			// нормаль и три вершины, 12 ещественных чисел.
			// Цилидр с центром в начале координат, единичного радиуса и единичной высоты.
			float cylinder_face[2400] = {
				// 1
				static_cast <float>(8.370229e-002), static_cast <float>(0.000000e+000), -9.964908e-001,
				static_cast <float>(1.000000e+000), static_cast <float>(0.000000e+000), 0.000000e+000,
				static_cast <float>(1.125333e+000), static_cast <float>(1.000000e+000), 7.885392e-003,
				static_cast <float>(1.125333e+000), static_cast <float>(0.000000e+000), 7.885392e-003,
				// 2
				static_cast <float>(1.667446e-001), static_cast <float>(0.000000e+000), -9.860001e-001,
				static_cast <float>(1.125333e+000), static_cast <float>(0.000000e+000), 7.885392e-003,
				static_cast <float>(1.125333e+000), static_cast <float>(1.000000e+000), 7.885392e-003,
				static_cast <float>(1.248690e+000), static_cast <float>(1.000000e+000), 3.141689e-002,
				// 3
				static_cast <float>(2.079357e-001), static_cast <float>(0.000000e+000), -9.781425e-001,
				static_cast <float>(1.125333e+000), static_cast <float>(0.000000e+000), 7.885392e-003,
				static_cast <float>(1.248690e+000), static_cast <float>(1.000000e+000), 3.141689e-002,
				static_cast <float>(1.248690e+000), static_cast <float>(0.000000e+000), 3.141689e-002,
				// 4
				 static_cast <float>(2.890083e-001), static_cast <float>(0.000000e+000), -9.573266e-001,
				 static_cast <float>(1.248690e+000), static_cast <float>(0.000000e+000), 3.141689e-002,
				 static_cast <float>(1.248690e+000), static_cast <float>(1.000000e+000), 3.141689e-002,
				 static_cast <float>(1.368125e+000), static_cast <float>(1.000000e+000), 7.022358e-002,
				//5
				static_cast <float>(3.288898e-001), static_cast <float>(0.000000e+000), -9.443683e-001,
				static_cast <float>(1.248690e+000), static_cast <float>(0.000000e+000), 3.141689e-002,
				static_cast <float>(1.368125e+000), static_cast <float>(1.000000e+000), 7.022358e-002,
				static_cast <float>(1.368125e+000), static_cast <float>(0.000000e+000), 7.022358e-002,
				//6
				static_cast <float>(4.067142e-001), 0.000000e+000, static_cast <float>(-9.135554e-001),
				static_cast <float>(1.368125e+000), 0.000000e+000, static_cast <float>(7.022358e-002),
				static_cast <float>(1.368125e+000), 1.000000e+000, static_cast <float>(7.022358e-002),
				static_cast <float>(1.481754e+000), 1.000000e+000, static_cast <float>(1.236934e-001),
				//7
				static_cast <float>(4.446571e-001), 0.000000e+000, -8.957008e-001,
				static_cast <float>(1.368125e+000), 0.000000e+000, 7.022358e-002,
				static_cast <float>(1.481754e+000), 1.000000e+000, 1.236934e-001,
				static_cast <float>(1.481754e+000), 0.000000e+000, 1.236934e-001,
				//8
				static_cast <float>(5.180060e-001), 0.000000e+000, -8.553770e-001,
				static_cast <float>(1.481754e+000), 0.000000e+000, 1.236934e-001,
				static_cast <float>(1.481754e+000), 1.000000e+000, 1.236934e-001,
				static_cast <float>(1.587785e+000), 1.000000e+000, 1.909831e-001,
				//9
				 static_cast <float>(5.534120e-001), 0.000000e+000, -8.329076e-001,
				 static_cast <float>(1.481754e+000), 0.000000e+000, 1.236934e-001,
				 static_cast <float>(1.587785e+000), 1.000000e+000, 1.909831e-001,
				 static_cast <float>(1.587785e+000), 0.000000e+000, 1.909831e-001,
				//10
				6.211286e-001, 0.000000e+000, -7.837087e-001,
				1.587785e+000, 0.000000e+000, 1.909831e-001,
		        1.587785e+000, 1.000000e+000, 1.909831e-001,
		        1.684547e+000, 1.000000e+000, 2.710314e-001,
				//11
				 6.534392e-001, 0.000000e+000, -7.569790e-001,
				1.587785e+000, 0.000000e+000, 1.909831e-001,
        		1.684547e+000, 1.000000e+000, 2.710314e-001,
		        1.684547e+000, 0.000000e+000, 2.710314e-001,
				//12
				7.144555e-001, 0.000000e+000, -6.996809e-001,
				1.684547e+000, 0.000000e+000, 2.710314e-001,
		        1.684547e+000, 1.000000e+000, 2.710314e-001,
		        1.770513e+000, 1.000000e+000, 3.625760e-001,
				//13
				7.431613e-001, 0.000000e+000, -6.691124e-001,
				1.684547e+000, 0.000000e+000, 2.710314e-001,
		        1.770513e+000, 1.000000e+000, 3.625760e-001,
		        1.770513e+000, 0.000000e+000, 3.625760e-001,
				//14
				7.965151e-001, 0.000000e+000, -6.046186e-001,
				1.770513e+000, 0.000000e+000, 3.625760e-001,
		        1.770513e+000, 1.000000e+000, 3.625760e-001,
		        1.844328e+000, 1.000000e+000, 4.641733e-001,
				//15
				8.211632e-001, 0.000000e+000, -5.706934e-001,
				1.770513e+000, 0.000000e+000, 3.625760e-001,
		        1.844328e+000, 1.000000e+000, 4.641733e-001,
		        1.844328e+000, 0.000000e+000, 4.641733e-001,
				//16
				8.660132e-001, 0.000000e+000, -5.000213e-001,
				1.844328e+000, 0.000000e+000, 4.641733e-001,
		        1.844328e+000, 1.000000e+000, 4.641733e-001,
		        1.904827e+000, 1.000000e+000, 5.742208e-001,
				//17
				8.862150e-001, 0.000000e+000, -4.632743e-001,
				1.844328e+000, 0.000000e+000, 4.641733e-001,
		        1.904827e+000, 1.000000e+000, 5.742208e-001,
		        1.904827e+000, 0.000000e+000, 5.742208e-001,
				//18
				9.218537e-001, 0.000000e+000, -3.875382e-001,
				1.904827e+000, 0.000000e+000, 5.742208e-001,
		        1.904827e+000, 1.000000e+000, 5.742208e-001,
		        1.951056e+000, 1.000000e+000, 6.909830e-001,
				//19
				9.372905e-001, 0.000000e+000, -3.485490e-001,
				1.904827e+000, 0.000000e+000, 5.742208e-001,
		        1.951056e+000, 1.000000e+000, 6.909830e-001,
		        1.951056e+000, 0.000000e+000, 6.909830e-001,
				//20
				9.631560e-001, 0.000000e+000, -2.689435e-001,
				1.951056e+000, 0.000000e+000, 6.909830e-001,
		        1.951056e+000, 1.000000e+000, 6.909830e-001,
		        1.982287e+000, 1.000000e+000, 8.126187e-001,
				//21
					9.735845e-001, 0.000000e+000, -2.283270e-001,
					1.951056e+000, 0.000000e+000, 6.909830e-001,
					1.982287e+000, 1.000000e+000, 8.126187e-001,
					1.982287e+000, 0.000000e+000, 8.126187e-001,
				//22
					9.892687e-001, 0.000000e+000, -1.461073e-001,
					1.982287e+000, 0.000000e+000, 8.126187e-001,
					1.982287e+000, 1.000000e+000, 8.126187e-001,
					1.998027e+000, 1.000000e+000, 9.372095e-001,
				//23
					9.945244e-001, 0.000000e+000, -1.045041e-001,
					1.982287e+000, 0.000000e+000, 8.126187e-001,
					1.998027e+000, 1.000000e+000, 9.372095e-001,
					1.998027e+000, 0.000000e+000, 9.372095e-001,
				//24
					 9.997802e-001, 0.000000e+000, -2.096695e-002,
					 1.998027e+000, 0.000000e+000, 9.372095e-001,
					 1.998027e+000, 1.000000e+000, 9.372095e-001,
					 1.998027e+000, 1.000000e+000, 1.062791e+000,
				//25
					9.997802e-001, 0.000000e+000, 2.096695e-002,
					1.998027e+000, 0.000000e+000, 9.372095e-001,
					1.998027e+000, 1.000000e+000, 1.062791e+000,
					1.998027e+000, 0.000000e+000, 1.062791e+000,
				//26
					9.945244e-001, 0.000000e+000, 1.045041e-001,
					1.998027e+000, 0.000000e+000, 1.062791e+000,
					1.998027e+000, 1.000000e+000, 1.062791e+000,
					1.982287e+000, 1.000000e+000, 1.187381e+000,
				//27
					9.892687e-001, 0.000000e+000, 1.461073e-001,
					1.998027e+000, 0.000000e+000, 1.062791e+000,
					1.982287e+000, 1.000000e+000, 1.187381e+000,
					1.982287e+000, 0.000000e+000, 1.187381e+000,
				//28
					9.735845e-001, 0.000000e+000, 2.283270e-001,
					1.982287e+000, 0.000000e+000, 1.187381e+000,
					1.982287e+000, 1.000000e+000, 1.187381e+000,
					1.951056e+000, 1.000000e+000, 1.309017e+000,
				//29
					 9.631560e-001, 0.000000e+000, 2.689435e-001,
					 1.982287e+000, 0.000000e+000, 1.187381e+000,
					 1.951056e+000, 1.000000e+000, 1.309017e+000,
					 1.951056e+000, 0.000000e+000, 1.309017e+000,
			    //30
					9.372905e-001, 0.000000e+000, 3.485490e-001,
					1.951056e+000, 0.000000e+000, 1.309017e+000,
					1.951056e+000, 1.000000e+000, 1.309017e+000,
					1.904827e+000, 1.000000e+000, 1.425779e+000,
				//31
					 9.218537e-001, 0.000000e+000, 3.875382e-001,
					1.951056e+000, 0.000000e+000, 1.309017e+000,
					1.904827e+000, 1.000000e+000, 1.425779e+000,
					1.904827e+000, 0.000000e+000, 1.425779e+000,
				//32
					8.862150e-001, 0.000000e+000, 4.632743e-001,
					1.904827e+000, 0.000000e+000, 1.425779e+000,
					1.904827e+000, 1.000000e+000, 1.425779e+000,
					1.844328e+000, 1.000000e+000, 1.535827e+000,
				//33
					8.660132e-001, 0.000000e+000, 5.000213e-001,
					1.904827e+000, 0.000000e+000, 1.425779e+000,
					1.844328e+000, 1.000000e+000, 1.535827e+000,
					1.844328e+000, 0.000000e+000, 1.535827e+000,
				//34
					8.211632e-001, 0.000000e+000, 5.706934e-001,
					1.844328e+000, 0.000000e+000, 1.535827e+000,
					1.844328e+000, 1.000000e+000, 1.535827e+000,
					1.770513e+000, 1.000000e+000, 1.637424e+000,
				//35
					7.965151e-001, 0.000000e+000, 6.046186e-001,
					1.844328e+000, 0.000000e+000, 1.535827e+000,
					1.770513e+000, 1.000000e+000, 1.637424e+000,
					1.770513e+000, 0.000000e+000, 1.637424e+000,
				//36
					7.431613e-001, 0.000000e+000, 6.691124e-001,
					1.770513e+000, 0.000000e+000, 1.637424e+000,
					1.770513e+000, 1.000000e+000, 1.637424e+000,
					1.684547e+000, 1.000000e+000, 1.728969e+000,
				//37
					7.144555e-001, 0.000000e+000, 6.996809e-001,
					1.770513e+000, 0.000000e+000, 1.637424e+000,
					1.684547e+000, 1.000000e+000, 1.728969e+000,
					1.684547e+000, 0.000000e+000, 1.728969e+000,
				//38
					6.534392e-001, 0.000000e+000, 7.569790e-001,
					1.684547e+000, 0.000000e+000, 1.728969e+000,
					1.684547e+000, 1.000000e+000, 1.728969e+000,
					1.587785e+000, 1.000000e+000, 1.809017e+000,
				//39
					6.211286e-001, 0.000000e+000, 7.837087e-001,
					1.684547e+000, 0.000000e+000, 1.728969e+000,
					1.587785e+000, 1.000000e+000, 1.809017e+000,
					1.587785e+000, 0.000000e+000, 1.809017e+000,
				//40
					5.534120e-001, 0.000000e+000, 8.329076e-001,
					1.587785e+000, 0.000000e+000, 1.809017e+000,
					1.587785e+000, 1.000000e+000, 1.809017e+000,
					1.481754e+000, 1.000000e+000, 1.876307e+000,
				//41
						5.180060e-001, 0.000000e+000, 8.553770e-001,
						1.587785e+000, 0.000000e+000, 1.809017e+000,
						1.481754e+000, 1.000000e+000, 1.876307e+000,
						1.481754e+000, 0.000000e+000, 1.876307e+000,
				//42
						4.446571e-001, 0.000000e+000, 8.957008e-001,
						1.481754e+000, 0.000000e+000, 1.876307e+000,
						1.481754e+000, 1.000000e+000, 1.876307e+000,
						1.368125e+000, 1.000000e+000, 1.929777e+000,
				//43
						4.067142e-001, 0.000000e+000, 9.135554e-001,
						1.481754e+000, 0.000000e+000, 1.876307e+000,
						1.368125e+000, 1.000000e+000, 1.929777e+000,
						1.368125e+000, 0.000000e+000, 1.929777e+000,
				//44
						 3.288898e-001, 0.000000e+000, 9.443683e-001,
						 1.368125e+000, 0.000000e+000, 1.929777e+000,
						 1.368125e+000, 1.000000e+000, 1.929777e+000,
						 1.248690e+000, 1.000000e+000, 1.968583e+000,
				//45
						2.890083e-001, 0.000000e+000, 9.573266e-001,
						1.368125e+000, 0.000000e+000, 1.929777e+000,
						1.248690e+000, 1.000000e+000, 1.968583e+000,
						1.248690e+000, 0.000000e+000, 1.968583e+000,
				//46
						2.079357e-001, 0.000000e+000, 9.781425e-001,
						1.248690e+000, 0.000000e+000, 1.968583e+000,
						1.248690e+000, 1.000000e+000, 1.968583e+000,
						1.125333e+000, 1.000000e+000, 1.992115e+000,
				//47
						1.667446e-001, 0.000000e+000, 9.860001e-001,
						1.248690e+000, 0.000000e+000, 1.968583e+000,
						1.125333e+000, 1.000000e+000, 1.992115e+000,
						1.125333e+000, 0.000000e+000, 1.992115e+000,
				//48
						8.370229e-002, 0.000000e+000, 9.964908e-001,
						1.125333e+000, 0.000000e+000, 1.992115e+000,
						1.125333e+000, 1.000000e+000, 1.992115e+000,
						1.000000e+000, 1.000000e+000, 2.000000e+000,
				//49
						 4.185114e-002, 0.000000e+000, 9.991238e-001,
						 1.125333e+000, 0.000000e+000, 1.992115e+000,
						 1.000000e+000, 1.000000e+000, 2.000000e+000,
						 1.000000e+000, 0.000000e+000, 2.000000e+000,
				//50
						-4.185114e-002, 0.000000e+000, 9.991238e-001,
						1.000000e+000, 0.000000e+000, 2.000000e+000,
						1.000000e+000, 1.000000e+000, 2.000000e+000,
						8.746669e-001, 1.000000e+000, 1.992115e+000,
				//51
						 -8.370229e-002, 0.000000e+000, 9.964908e-001,
						1.000000e+000, 0.000000e+000, 2.000000e+000,
						8.746669e-001, 1.000000e+000, 1.992115e+000,
						8.746669e-001, 0.000000e+000, 1.992115e+000,
				//52
						-1.667446e-001, 0.000000e+000, 9.860001e-001,
						8.746669e-001, 0.000000e+000, 1.992115e+000,
						8.746669e-001, 1.000000e+000, 1.992115e+000,
						7.513102e-001, 1.000000e+000, 1.968583e+000,
				//53
						-2.079357e-001, 0.000000e+000, 9.781425e-001,
						8.746669e-001, 0.000000e+000, 1.992115e+000,
						7.513102e-001, 1.000000e+000, 1.968583e+000,
						7.513102e-001, 0.000000e+000, 1.968583e+000,
				//54
						-2.890083e-001, 0.000000e+000, 9.573266e-001,
						7.513102e-001, 0.000000e+000, 1.968583e+000,
						7.513102e-001, 1.000000e+000, 1.968583e+000,
						6.318755e-001, 1.000000e+000, 1.929777e+000,
				//55
						-3.288898e-001, 0.000000e+000, 9.443683e-001,
						7.513102e-001, 0.000000e+000, 1.968583e+000,
						6.318755e-001, 1.000000e+000, 1.929777e+000,
						6.318755e-001, 0.000000e+000, 1.929777e+000,
				//56
						-4.067142e-001, 0.000000e+000, 9.135554e-001,
						6.318755e-001, 0.000000e+000, 1.929777e+000,
						6.318755e-001, 1.000000e+000, 1.929777e+000,
						5.182464e-001, 1.000000e+000, 1.876307e+000,
				//57
						-4.446571e-001, 0.000000e+000, 8.957008e-001,
						6.318755e-001, 0.000000e+000, 1.929777e+000,
						5.182464e-001, 1.000000e+000, 1.876307e+000,
						5.182464e-001, 0.000000e+000, 1.876307e+000,
				//58
						-5.180060e-001, 0.000000e+000, 8.553770e-001,
						5.182464e-001, 0.000000e+000, 1.876307e+000,
						5.182464e-001, 1.000000e+000, 1.876307e+000,
						4.122148e-001, 1.000000e+000, 1.809017e+000,
				//59
						-5.534120e-001, 0.000000e+000, 8.329076e-001,
						5.182464e-001, 0.000000e+000, 1.876307e+000,
						4.122148e-001, 1.000000e+000, 1.809017e+000,
						4.122148e-001, 0.000000e+000, 1.809017e+000,
				//60
						-6.211286e-001, 0.000000e+000, 7.837087e-001,
						4.122148e-001, 0.000000e+000, 1.809017e+000,
						4.122148e-001, 1.000000e+000, 1.809017e+000,
						3.154529e-001, 1.000000e+000, 1.728969e+000,
				//61
							-6.534392e-001, 0.000000e+000, 7.569790e-001,
							4.122148e-001, 0.000000e+000, 1.809017e+000,
							3.154529e-001, 1.000000e+000, 1.728969e+000,
							3.154529e-001, 0.000000e+000, 1.728969e+000,
				//62
							-7.144555e-001, 0.000000e+000, 6.996809e-001,
							3.154529e-001, 0.000000e+000, 1.728969e+000,
							3.154529e-001, 1.000000e+000, 1.728969e+000,
							2.294868e-001, 1.000000e+000, 1.637424e+000,
				//63
							-7.431613e-001, 0.000000e+000, 6.691124e-001,
							3.154529e-001, 0.000000e+000, 1.728969e+000,
							2.294868e-001, 1.000000e+000, 1.637424e+000,
							2.294868e-001, 0.000000e+000, 1.637424e+000,
				//64
							 -7.965151e-001, 0.000000e+000, 6.046186e-001,
							 2.294868e-001, 0.000000e+000, 1.637424e+000,
							 2.294868e-001, 1.000000e+000, 1.637424e+000,
							 1.556721e-001, 1.000000e+000, 1.535827e+000,
				//65
							-8.211632e-001, 0.000000e+000, 5.706934e-001,
							2.294868e-001, 0.000000e+000, 1.637424e+000,
							1.556721e-001, 1.000000e+000, 1.535827e+000,
							1.556721e-001, 0.000000e+000, 1.535827e+000,
				//66
							-8.660132e-001, 0.000000e+000, 5.000213e-001,
							1.556721e-001, 0.000000e+000, 1.535827e+000,
							1.556721e-001, 1.000000e+000, 1.535827e+000,
							9.517302e-002, 1.000000e+000, 1.425779e+000,
				//67
							-8.862150e-001, 0.000000e+000, 4.632743e-001,
							1.556721e-001, 0.000000e+000, 1.535827e+000,
							9.517302e-002, 1.000000e+000, 1.425779e+000,
							9.517302e-002, 0.000000e+000, 1.425779e+000,
				//68
							-9.218537e-001, 0.000000e+000, 3.875382e-001,
							9.517302e-002, 0.000000e+000, 1.425779e+000,
							9.517302e-002, 1.000000e+000, 1.425779e+000,
							4.894350e-002, 1.000000e+000, 1.309017e+000,
				//69
							 -9.372905e-001, 0.000000e+000, 3.485490e-001,
							 9.517302e-002, 0.000000e+000, 1.425779e+000,
							 4.894350e-002, 1.000000e+000, 1.309017e+000,
							 4.894350e-002, 0.000000e+000, 1.309017e+000,
				//70
							-9.631560e-001, 0.000000e+000, 2.689435e-001,
							4.894350e-002, 0.000000e+000, 1.309017e+000,
							4.894350e-002, 1.000000e+000, 1.309017e+000,
							1.771282e-002, 1.000000e+000, 1.187381e+000,
				//71
							 -9.735845e-001, 0.000000e+000, 2.283270e-001,
							4.894350e-002, 0.000000e+000, 1.309017e+000,
							1.771282e-002, 1.000000e+000, 1.187381e+000,
							1.771282e-002, 0.000000e+000, 1.187381e+000,
				//72
							-9.892687e-001, 0.000000e+000, 1.461073e-001,
							1.771282e-002, 0.000000e+000, 1.187381e+000,
							1.771282e-002, 1.000000e+000, 1.187381e+000,
							1.973356e-003, 1.000000e+000, 1.062791e+000,
				//73
							-9.945244e-001, 0.000000e+000, 1.045041e-001,
							1.771282e-002, 0.000000e+000, 1.187381e+000,
							1.973356e-003, 1.000000e+000, 1.062791e+000,
							1.973356e-003, 0.000000e+000, 1.062791e+000,
				//74
							-9.997802e-001, 0.000000e+000, 2.096695e-002,
							1.973356e-003, 0.000000e+000, 1.062791e+000,
							1.973356e-003, 1.000000e+000, 1.062791e+000,
							1.973356e-003, 1.000000e+000, 9.372095e-001,
				//75
							-9.997802e-001, 0.000000e+000, -2.096695e-002,
							1.973356e-003, 0.000000e+000, 1.062791e+000,
							1.973356e-003, 1.000000e+000, 9.372095e-001,
							1.973356e-003, 0.000000e+000, 9.372095e-001,
				//76
							-9.945244e-001, 0.000000e+000, -1.045041e-001,
							1.973356e-003, 0.000000e+000, 9.372095e-001,
							1.973356e-003, 1.000000e+000, 9.372095e-001,
							1.771282e-002, 1.000000e+000, 8.126187e-001,
				//77
							-9.892687e-001, 0.000000e+000, -1.461073e-001,
							1.973356e-003, 0.000000e+000, 9.372095e-001,
							1.771282e-002, 1.000000e+000, 8.126187e-001,
							1.771282e-002, 0.000000e+000, 8.126187e-001,
				//78
							-9.735845e-001, 0.000000e+000, -2.283270e-001,
							1.771282e-002, 0.000000e+000, 8.126187e-001,
							1.771282e-002, 1.000000e+000, 8.126187e-001,
							4.894350e-002, 1.000000e+000, 6.909830e-001,
				//79
							-9.631560e-001, 0.000000e+000, -2.689435e-001,
							1.771282e-002, 0.000000e+000, 8.126187e-001,
							4.894350e-002, 1.000000e+000, 6.909830e-001,
							4.894350e-002, 0.000000e+000, 6.909830e-001,
				//80
							-9.372905e-001, 0.000000e+000, -3.485490e-001,
							4.894350e-002, 0.000000e+000, 6.909830e-001,
							4.894350e-002, 1.000000e+000, 6.909830e-001,
							9.517302e-002, 1.000000e+000, 5.742208e-001,
				//81
								-9.218537e-001, 0.000000e+000, -3.875382e-001,
								4.894350e-002, 0.000000e+000, 6.909830e-001,
								9.517302e-002, 1.000000e+000, 5.742208e-001,
								9.517302e-002, 0.000000e+000, 5.742208e-001,
				//82
								-8.862150e-001, 0.000000e+000, -4.632743e-001,
								9.517302e-002, 0.000000e+000, 5.742208e-001,
								9.517302e-002, 1.000000e+000, 5.742208e-001,
								1.556721e-001, 1.000000e+000, 4.641733e-001,
				//83
								-8.660132e-001, 0.000000e+000, -5.000213e-001,
								9.517302e-002, 0.000000e+000, 5.742208e-001,
								1.556721e-001, 1.000000e+000, 4.641733e-001,
								1.556721e-001, 0.000000e+000, 4.641733e-001,
				//84
								 -8.211632e-001, 0.000000e+000, -5.706934e-001,
								 1.556721e-001, 0.000000e+000, 4.641733e-001,
								 1.556721e-001, 1.000000e+000, 4.641733e-001,
								 2.294868e-001, 1.000000e+000, 3.625760e-001,
				//85
								-7.965151e-001, 0.000000e+000, -6.046186e-001,
								1.556721e-001, 0.000000e+000, 4.641733e-001,
								2.294868e-001, 1.000000e+000, 3.625760e-001,
								2.294868e-001, 0.000000e+000, 3.625760e-001,
				//86
								-7.431613e-001, 0.000000e+000, -6.691124e-001,
								2.294868e-001, 0.000000e+000, 3.625760e-001,
								2.294868e-001, 1.000000e+000, 3.625760e-001,
								3.154529e-001, 1.000000e+000, 2.710314e-001,
				//87
								-7.144555e-001, 0.000000e+000, -6.996809e-001,
								2.294868e-001, 0.000000e+000, 3.625760e-001,
								3.154529e-001, 1.000000e+000, 2.710314e-001,
								3.154529e-001, 0.000000e+000, 2.710314e-001,
				//88
								-6.534392e-001, 0.000000e+000, -7.569790e-001,
								3.154529e-001, 0.000000e+000, 2.710314e-001,
								3.154529e-001, 1.000000e+000, 2.710314e-001,
								4.122148e-001, 1.000000e+000, 1.909831e-001,
				//89
								 -6.211286e-001, 0.000000e+000, -7.837087e-001,
								 3.154529e-001, 0.000000e+000, 2.710314e-001,
								 4.122148e-001, 1.000000e+000, 1.909831e-001,
								 4.122148e-001, 0.000000e+000, 1.909831e-001,
				//90
								-5.534120e-001, 0.000000e+000, -8.329076e-001,
								4.122148e-001, 0.000000e+000, 1.909831e-001,
								4.122148e-001, 1.000000e+000, 1.909831e-001,
								5.182464e-001, 1.000000e+000, 1.236934e-001,
				//91
								 -5.180060e-001, 0.000000e+000, -8.553770e-001,
								4.122148e-001, 0.000000e+000, 1.909831e-001,
								5.182464e-001, 1.000000e+000, 1.236934e-001,
								5.182464e-001, 0.000000e+000, 1.236934e-001,
				//92
								-4.446571e-001, 0.000000e+000, -8.957008e-001,
								5.182464e-001, 0.000000e+000, 1.236934e-001,
								5.182464e-001, 1.000000e+000, 1.236934e-001,
								6.318755e-001, 1.000000e+000, 7.022358e-002,
				//93
								-4.067142e-001, 0.000000e+000, -9.135554e-001,
								5.182464e-001, 0.000000e+000, 1.236934e-001,
								6.318755e-001, 1.000000e+000, 7.022358e-002,
								6.318755e-001, 0.000000e+000, 7.022358e-002,
				//94
								-3.288898e-001, 0.000000e+000, -9.443683e-001,
								6.318755e-001, 0.000000e+000, 7.022358e-002,
								6.318755e-001, 1.000000e+000, 7.022358e-002,
								7.513102e-001, 1.000000e+000, 3.141689e-002,
				//95
								-2.890083e-001, 0.000000e+000, -9.573266e-001,
								6.318755e-001, 0.000000e+000, 7.022358e-002,
								7.513102e-001, 1.000000e+000, 3.141689e-002,
								7.513102e-001, 0.000000e+000, 3.141689e-002,
				//96
								-2.079357e-001, 0.000000e+000, -9.781425e-001,
								7.513102e-001, 0.000000e+000, 3.141689e-002,
								7.513102e-001, 1.000000e+000, 3.141689e-002,
								8.746669e-001, 1.000000e+000, 7.885392e-003,
				//97
								-1.667446e-001, 0.000000e+000, -9.860001e-001,
								7.513102e-001, 0.000000e+000, 3.141689e-002,
								8.746669e-001, 1.000000e+000, 7.885392e-003,
								8.746669e-001, 0.000000e+000, 7.885392e-003,
				//98
								-8.370229e-002, 0.000000e+000, -9.964908e-001,
								8.746669e-001, 0.000000e+000, 7.885392e-003,
								8.746669e-001, 1.000000e+000, 7.885392e-003,
								1.000000e+000, 1.000000e+000, 0.000000e+000,
				//99
								-4.185114e-002, 0.000000e+000, -9.991238e-001,
								8.746669e-001, 0.000000e+000, 7.885392e-003,
								1.000000e+000, 1.000000e+000, 0.000000e+000,
								1.000000e+000, 0.000000e+000, 0.000000e+000,
				//100
								4.185114e-002, 0.000000e+000, -9.991238e-001,
								1.000000e+000, 0.000000e+000, 0.000000e+000,
								1.000000e+000, 1.000000e+000, 0.000000e+000,
								1.125333e+000, 1.000000e+000, 7.885392e-003,
                //101
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 2.000000e+000,
									1.125333e+000, 1.000000e+000, 1.992115e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
                //102
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 2.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									8.746669e-001, 1.000000e+000, 1.992115e+000,
				//103
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									5.182464e-001, 1.000000e+000, 1.876307e+000,
									6.318755e-001, 1.000000e+000, 1.929777e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
                //104
									 0.000000e+000, 1.000000e+000, 0.000000e+000,
									 1.000000e+000, 1.000000e+000, 1.000000e+000,
									 6.318755e-001, 1.000000e+000, 1.929777e+000,
									 7.513102e-001, 1.000000e+000, 1.968583e+000,
                //105
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									7.513102e-001, 1.000000e+000, 1.968583e+000,
									8.746669e-001, 1.000000e+000, 1.992115e+000,
                //106
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									2.294868e-001, 1.000000e+000, 1.637424e+000,
									3.154529e-001, 1.000000e+000, 1.728969e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
                //107
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									3.154529e-001, 1.000000e+000, 1.728969e+000,
									4.122148e-001, 1.000000e+000, 1.809017e+000,
                //108
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									4.122148e-001, 1.000000e+000, 1.809017e+000,
									5.182464e-001, 1.000000e+000, 1.876307e+000,
                //109
									 0.000000e+000, 1.000000e+000, 0.000000e+000,
									 4.894350e-002, 1.000000e+000, 1.309017e+000,
									 9.517302e-002, 1.000000e+000, 1.425779e+000,
									 1.000000e+000, 1.000000e+000, 1.000000e+000,
                //110
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									9.517302e-002, 1.000000e+000, 1.425779e+000,
									1.556721e-001, 1.000000e+000, 1.535827e+000,
                //111
									 0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									1.556721e-001, 1.000000e+000, 1.535827e+000,
									2.294868e-001, 1.000000e+000, 1.637424e+000,
                //112
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.973356e-003, 1.000000e+000, 9.372095e-001,
									1.973356e-003, 1.000000e+000, 1.062791e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
                //113
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									1.973356e-003, 1.000000e+000, 1.062791e+000,
									1.771282e-002, 1.000000e+000, 1.187381e+000,
                //114
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									1.771282e-002, 1.000000e+000, 1.187381e+000,
									4.894350e-002, 1.000000e+000, 1.309017e+000,
                //115
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									9.517302e-002, 1.000000e+000, 5.742208e-001,
									4.894350e-002, 1.000000e+000, 6.909830e-001,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
                //116
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									4.894350e-002, 1.000000e+000, 6.909830e-001,
									1.771282e-002, 1.000000e+000, 8.126187e-001,
                //117
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									1.771282e-002, 1.000000e+000, 8.126187e-001,
									1.973356e-003, 1.000000e+000, 9.372095e-001,
                //118
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									3.154529e-001, 1.000000e+000, 2.710314e-001,
									2.294868e-001, 1.000000e+000, 3.625760e-001,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
                //119
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									2.294868e-001, 1.000000e+000, 3.625760e-001,
									1.556721e-001, 1.000000e+000, 4.641733e-001,
                //120
									0.000000e+000, 1.000000e+000, 0.000000e+000,
									1.000000e+000, 1.000000e+000, 1.000000e+000,
									1.556721e-001, 1.000000e+000, 4.641733e-001,
									9.517302e-002, 1.000000e+000, 5.742208e-001,
                //121
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										6.318755e-001, 1.000000e+000, 7.022358e-002,
										5.182464e-001, 1.000000e+000, 1.236934e-001,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
				//122
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										5.182464e-001, 1.000000e+000, 1.236934e-001,
										4.122148e-001, 1.000000e+000, 1.909831e-001,
				//123
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										4.122148e-001, 1.000000e+000, 1.909831e-001,
										3.154529e-001, 1.000000e+000, 2.710314e-001,
				//124
										 0.000000e+000, 1.000000e+000, 0.000000e+000,
										 1.000000e+000, 1.000000e+000, 0.000000e+000,
										 8.746669e-001, 1.000000e+000, 7.885392e-003,
										 1.000000e+000, 1.000000e+000, 1.000000e+000,
				//125
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										8.746669e-001, 1.000000e+000, 7.885392e-003,
										7.513102e-001, 1.000000e+000, 3.141689e-002,
				//126
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										7.513102e-001, 1.000000e+000, 3.141689e-002,
										6.318755e-001, 1.000000e+000, 7.022358e-002,
				//127
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.368125e+000, 1.000000e+000, 7.022358e-002,
										1.248690e+000, 1.000000e+000, 3.141689e-002,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
				//128
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.248690e+000, 1.000000e+000, 3.141689e-002,
										1.125333e+000, 1.000000e+000, 7.885392e-003,
				//129
										 0.000000e+000, 1.000000e+000, 0.000000e+000,
										 1.000000e+000, 1.000000e+000, 1.000000e+000,
										 1.125333e+000, 1.000000e+000, 7.885392e-003,
										 1.000000e+000, 1.000000e+000, 0.000000e+000,
				//130
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.684547e+000, 1.000000e+000, 2.710314e-001,
										1.587785e+000, 1.000000e+000, 1.909831e-001,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
				//131
										 0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.587785e+000, 1.000000e+000, 1.909831e-001,
										1.481754e+000, 1.000000e+000, 1.236934e-001,
				//132
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.481754e+000, 1.000000e+000, 1.236934e-001,
										1.368125e+000, 1.000000e+000, 7.022358e-002,
				//133
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.904827e+000, 1.000000e+000, 5.742208e-001,
										1.844328e+000, 1.000000e+000, 4.641733e-001,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
				//134
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.844328e+000, 1.000000e+000, 4.641733e-001,
										1.770513e+000, 1.000000e+000, 3.625760e-001,
				//135
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.770513e+000, 1.000000e+000, 3.625760e-001,
										1.684547e+000, 1.000000e+000, 2.710314e-001,
				//136
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.998027e+000, 1.000000e+000, 9.372095e-001,
										1.982287e+000, 1.000000e+000, 8.126187e-001,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
				//137
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.982287e+000, 1.000000e+000, 8.126187e-001,
										1.951056e+000, 1.000000e+000, 6.909830e-001,
				//138
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.951056e+000, 1.000000e+000, 6.909830e-001,
										1.904827e+000, 1.000000e+000, 5.742208e-001,
				//139
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.951056e+000, 1.000000e+000, 1.309017e+000,
										1.982287e+000, 1.000000e+000, 1.187381e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
				//140
										0.000000e+000, 1.000000e+000, 0.000000e+000,
										1.000000e+000, 1.000000e+000, 1.000000e+000,
										1.982287e+000, 1.000000e+000, 1.187381e+000,
										1.998027e+000, 1.000000e+000, 1.062791e+000,
				//141
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
											1.998027e+000, 1.000000e+000, 1.062791e+000,
											1.998027e+000, 1.000000e+000, 9.372095e-001,
				//142
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.770513e+000, 1.000000e+000, 1.637424e+000,
											1.844328e+000, 1.000000e+000, 1.535827e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
				//143
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
											1.844328e+000, 1.000000e+000, 1.535827e+000,
											1.904827e+000, 1.000000e+000, 1.425779e+000,
				//144
											 0.000000e+000, 1.000000e+000, 0.000000e+000,
											 1.000000e+000, 1.000000e+000, 1.000000e+000,
											 1.904827e+000, 1.000000e+000, 1.425779e+000,
											 1.951056e+000, 1.000000e+000, 1.309017e+000,
				//145
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.481754e+000, 1.000000e+000, 1.876307e+000,
											1.587785e+000, 1.000000e+000, 1.809017e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
				//146
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
											1.587785e+000, 1.000000e+000, 1.809017e+000,
											1.684547e+000, 1.000000e+000, 1.728969e+000,
				//147
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
											1.684547e+000, 1.000000e+000, 1.728969e+000,
											1.770513e+000, 1.000000e+000, 1.637424e+000,
				//148
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.125333e+000, 1.000000e+000, 1.992115e+000,
											1.248690e+000, 1.000000e+000, 1.968583e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
				//149
											 0.000000e+000, 1.000000e+000, 0.000000e+000,
											 1.000000e+000, 1.000000e+000, 1.000000e+000,
											 1.248690e+000, 1.000000e+000, 1.968583e+000,
											 1.368125e+000, 1.000000e+000, 1.929777e+000,
				//150
											0.000000e+000, 1.000000e+000, 0.000000e+000,
											1.000000e+000, 1.000000e+000, 1.000000e+000,
											1.368125e+000, 1.000000e+000, 1.929777e+000,
											1.481754e+000, 1.000000e+000, 1.876307e+000,
				//151
											 0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 2.000000e+000,
											8.746669e-001, 0.000000e+000, 1.992115e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
				//152
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 2.000000e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
											1.125333e+000, 0.000000e+000, 1.992115e+000,
				//153
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.481754e+000, 0.000000e+000, 1.876307e+000,
											1.368125e+000, 0.000000e+000, 1.929777e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
				//154
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
											1.368125e+000, 0.000000e+000, 1.929777e+000,
											1.248690e+000, 0.000000e+000, 1.968583e+000,
				//155
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
											1.248690e+000, 0.000000e+000, 1.968583e+000,
											1.125333e+000, 0.000000e+000, 1.992115e+000,
				//156
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.770513e+000, 0.000000e+000, 1.637424e+000,
											1.684547e+000, 0.000000e+000, 1.728969e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
				//157
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
											1.684547e+000, 0.000000e+000, 1.728969e+000,
											1.587785e+000, 0.000000e+000, 1.809017e+000,
				//158
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
											1.587785e+000, 0.000000e+000, 1.809017e+000,
											1.481754e+000, 0.000000e+000, 1.876307e+000,
				//159
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.951056e+000, 0.000000e+000, 1.309017e+000,
											1.904827e+000, 0.000000e+000, 1.425779e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
				//160
											0.000000e+000, -1.000000e+000, 0.000000e+000,
											1.000000e+000, 0.000000e+000, 1.000000e+000,
											1.904827e+000, 0.000000e+000, 1.425779e+000,
											1.844328e+000, 0.000000e+000, 1.535827e+000,
				//161
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.844328e+000, 0.000000e+000, 1.535827e+000,
												1.770513e+000, 0.000000e+000, 1.637424e+000,
				//162
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.998027e+000, 0.000000e+000, 9.372095e-001,
												1.998027e+000, 0.000000e+000, 1.062791e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
				//163
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.998027e+000, 0.000000e+000, 1.062791e+000,
												1.982287e+000, 0.000000e+000, 1.187381e+000,
				//164
												 0.000000e+000, -1.000000e+000, 0.000000e+000,
												 1.000000e+000, 0.000000e+000, 1.000000e+000,
												 1.982287e+000, 0.000000e+000, 1.187381e+000,
												 1.951056e+000, 0.000000e+000, 1.309017e+000,
				//165
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.904827e+000, 0.000000e+000, 5.742208e-001,
												1.951056e+000, 0.000000e+000, 6.909830e-001,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
				//166
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.951056e+000, 0.000000e+000, 6.909830e-001,
												1.982287e+000, 0.000000e+000, 8.126187e-001,
				//167
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.982287e+000, 0.000000e+000, 8.126187e-001,
												1.998027e+000, 0.000000e+000, 9.372095e-001,
				//168
												static_cast <float>(0.000000e+000), static_cast <float>(-1.000000e+000), static_cast <float>(0.000000e+000),
												static_cast <float>(1.684547e+000), static_cast <float>(0.000000e+000), static_cast <float>(2.710314e-001),
												static_cast <float>(1.770513e+000), static_cast <float>(0.000000e+000), static_cast <float>(3.625760e-001),
												static_cast <float>(1.000000e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.000000e+000),
				//169
												static_cast <float>(0.000000e+000), static_cast <float>(-1.000000e+000), static_cast <float>(0.000000e+000),
												static_cast <float>(1.000000e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.000000e+000),
												static_cast <float>(1.770513e+000), static_cast <float>(0.000000e+000), static_cast <float>(3.625760e-001),
												static_cast <float>(1.844328e+000), static_cast <float>(0.000000e+000), static_cast <float>(4.641733e-001),
				//170
												static_cast <float>(0.000000e+000), static_cast <float>(-1.000000e+000), static_cast <float>(0.000000e+000),
												static_cast <float>(1.000000e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.000000e+000),
												static_cast <float>(1.844328e+000), static_cast <float>(0.000000e+000), static_cast <float>(4.641733e-001),
												static_cast <float>(1.904827e+000), static_cast <float>(0.000000e+000), static_cast <float>(5.742208e-001),
				//171
												static_cast <float>(0.000000e+000), static_cast <float>(-1.000000e+000), static_cast <float>(0.000000e+000),
												static_cast <float>(1.368125e+000), static_cast <float>(0.000000e+000), static_cast <float>(7.022358e-002),
												static_cast <float>(1.481754e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.236934e-001),
												static_cast <float>(1.000000e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.000000e+000),
				//172
												static_cast <float>(0.000000e+000), static_cast <float>(-1.000000e+000), static_cast <float>(0.000000e+000),
												static_cast <float>(1.000000e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.000000e+000),
												static_cast <float>(1.481754e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.236934e-001),
												static_cast <float>(1.587785e+000), static_cast <float>(0.000000e+000), static_cast <float>(1.909831e-001),
				//173
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.587785e+000, 0.000000e+000, 1.909831e-001,
												1.684547e+000, 0.000000e+000, 2.710314e-001,
				//174
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 0.000000e+000,
												1.125333e+000, 0.000000e+000, 7.885392e-003,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
				//175
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.125333e+000, 0.000000e+000, 7.885392e-003,
												1.248690e+000, 0.000000e+000, 3.141689e-002,
				//176
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												1.248690e+000, 0.000000e+000, 3.141689e-002,
												1.368125e+000, 0.000000e+000, 7.022358e-002,
				//177
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												6.318755e-001, 0.000000e+000, 7.022358e-002,
												7.513102e-001, 0.000000e+000, 3.141689e-002,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
				//178
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												7.513102e-001, 0.000000e+000, 3.141689e-002,
												8.746669e-001, 0.000000e+000, 7.885392e-003,
				//179
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
												8.746669e-001, 0.000000e+000, 7.885392e-003,
												1.000000e+000, 0.000000e+000, 0.000000e+000,
				//180
												0.000000e+000, -1.000000e+000, 0.000000e+000,
												3.154529e-001, 0.000000e+000, 2.710314e-001,
												4.122148e-001, 0.000000e+000, 1.909831e-001,
												1.000000e+000, 0.000000e+000, 1.000000e+000,
				//181
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													4.122148e-001, 0.000000e+000, 1.909831e-001,
													5.182464e-001, 0.000000e+000, 1.236934e-001,
				//182
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													5.182464e-001, 0.000000e+000, 1.236934e-001,
													6.318755e-001, 0.000000e+000, 7.022358e-002,
				//183
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													9.517302e-002, 0.000000e+000, 5.742208e-001,
													1.556721e-001, 0.000000e+000, 4.641733e-001,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
				//184
													 0.000000e+000, -1.000000e+000, 0.000000e+000,
													 1.000000e+000, 0.000000e+000, 1.000000e+000,
													 1.556721e-001, 0.000000e+000, 4.641733e-001,
													 2.294868e-001, 0.000000e+000, 3.625760e-001,
				//185
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													2.294868e-001, 0.000000e+000, 3.625760e-001,
													3.154529e-001, 0.000000e+000, 2.710314e-001,
				//186
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.973356e-003, 0.000000e+000, 9.372095e-001,
													1.771282e-002, 0.000000e+000, 8.126187e-001,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
				//187
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													1.771282e-002, 0.000000e+000, 8.126187e-001,
													4.894350e-002, 0.000000e+000, 6.909830e-001,
				//188
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													4.894350e-002, 0.000000e+000, 6.909830e-001,
													9.517302e-002, 0.000000e+000, 5.742208e-001,
				//189
													 0.000000e+000, -1.000000e+000, 0.000000e+000,
													 4.894350e-002, 0.000000e+000, 1.309017e+000,
													 1.771282e-002, 0.000000e+000, 1.187381e+000,
													 1.000000e+000, 0.000000e+000, 1.000000e+000,
				//190
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													1.771282e-002, 0.000000e+000, 1.187381e+000,
													1.973356e-003, 0.000000e+000, 1.062791e+000,
				//191
													 0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													1.973356e-003, 0.000000e+000, 1.062791e+000,
													1.973356e-003, 0.000000e+000, 9.372095e-001,
				//192
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													2.294868e-001, 0.000000e+000, 1.637424e+000,
													1.556721e-001, 0.000000e+000, 1.535827e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
				//193
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													1.556721e-001, 0.000000e+000, 1.535827e+000,
													9.517302e-002, 0.000000e+000, 1.425779e+000,
				//194
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													9.517302e-002, 0.000000e+000, 1.425779e+000,
													4.894350e-002, 0.000000e+000, 1.309017e+000,
				//195
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													5.182464e-001, 0.000000e+000, 1.876307e+000,
													4.122148e-001, 0.000000e+000, 1.809017e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
				//196
													0.000000e+000, -1.000000e+000, 0.000000e+000,
													1.000000e+000, 0.000000e+000, 1.000000e+000,
													4.122148e-001, 0.000000e+000, 1.809017e+000,
													3.154529e-001, 0.000000e+000, 1.728969e+000,
				//197
													static_cast <float>(0.000000e+000), -1.000000e+000, 0.000000e+000,
													static_cast <float>(1.000000e+000), 0.000000e+000, 1.000000e+000,
													static_cast <float>(3.154529e-001), 0.000000e+000, 1.728969e+000,
													static_cast <float>(2.294868e-001), 0.000000e+000, 1.637424e+000,
				//198
													static_cast <float>(0.000000e+000), -1.000000e+000, 0.000000e+000,
													static_cast <float>(8.746669e-001), 0.000000e+000, 1.992115e+000,
													static_cast <float>(7.513102e-001), 0.000000e+000, 1.968583e+000,
													static_cast <float>(1.000000e+000), 0.000000e+000, 1.000000e+000,
				//199
													static_cast <float>(0.000000e+000), -1.000000e+000, 0.000000e+000,
													static_cast <float>(1.000000e+000), 0.000000e+000, 1.000000e+000,
													static_cast <float>(7.513102e-001), 0.000000e+000, 1.968583e+000,
													static_cast <float>(6.318755e-001), 0.000000e+000, 1.929777e+000,
				//200
													static_cast <float>(0.000000e+000), -1.000000e+000, 0.000000e+000,
													static_cast <float>(1.000000e+000), 0.000000e+000, 1.000000e+000,
													static_cast <float>(6.318755e-001), 0.000000e+000, 1.929777e+000,
													static_cast <float>(5.182464e-001), 0.000000e+000, 1.876307e+000
			};


			// XZ_PLANE
			root_CAD_STL = new gCAD_STL;
			gCAD_STL* tmp = root_CAD_STL;
			int i61 = 0;
			for (int i60 = 0; i60 < 200; ++i60) {

				if (iPlane == XZ_PLANE) {

					tmp->n.x = cylinder_face[i61 + 0];
					tmp->n.y = cylinder_face[i61 + 1];
					tmp->n.z = cylinder_face[i61 + 2];
					tmp->pa.x = static_cast <float>(xC + R_out_cyl * (cylinder_face[i61 + 3] - 1.0));
					tmp->pa.y = static_cast <float>(yC + Hcyl * cylinder_face[i61 + 4]);
					tmp->pa.z = static_cast <float>(zC + R_out_cyl * (cylinder_face[i61 + 5] - 1.0));
					tmp->pb.x = static_cast <float>(xC + R_out_cyl * (cylinder_face[i61 + 6] - 1.0));
					tmp->pb.y = static_cast <float>(yC + Hcyl * cylinder_face[i61 + 7]);
					tmp->pb.z = static_cast <float>(zC + R_out_cyl * (cylinder_face[i61 + 8] - 1.0));
					tmp->pc.x = static_cast <float>(xC + R_out_cyl * (cylinder_face[i61 + 9] - 1.0));
					tmp->pc.y = static_cast <float>(yC + Hcyl * cylinder_face[i61 + 10]);
					tmp->pc.z = static_cast <float>(zC + R_out_cyl * (cylinder_face[i61 + 11] - 1.0));
				}
				if (iPlane == XY_PLANE) {

					tmp->n.x = cylinder_face[i61 + 0];
					tmp->n.y = cylinder_face[i61 + 2];
					tmp->n.z = cylinder_face[i61 + 1];
					tmp->pa.x = static_cast <float>(xC + R_out_cyl * (cylinder_face[i61 + 3] - 1.0));
					tmp->pa.y = static_cast <float>(yC + R_out_cyl * (cylinder_face[i61 + 5] - 1.0));
					tmp->pa.z = static_cast <float>(zC + Hcyl *  cylinder_face[i61 + 4]);
					tmp->pb.x = static_cast <float>(xC + R_out_cyl * (cylinder_face[i61 + 6] - 1.0));
					tmp->pb.y = static_cast <float>(yC + R_out_cyl * (cylinder_face[i61 + 8] - 1.0));
					tmp->pb.z = static_cast <float>(zC + Hcyl *  cylinder_face[i61 + 7]);
					tmp->pc.x = static_cast <float>(xC + R_out_cyl * (cylinder_face[i61 + 9] - 1.0));
					tmp->pc.y = static_cast <float>(yC + R_out_cyl * (cylinder_face[i61 + 11] - 1.0));
					tmp->pc.z = static_cast <float>(zC + Hcyl *  cylinder_face[i61 + 10]);
				}
				if (iPlane == YZ_PLANE) {

					tmp->n.x = cylinder_face[i61 + 1];
					tmp->n.y = cylinder_face[i61 + 0];
					tmp->n.z = cylinder_face[i61 + 2];
					tmp->pa.x = static_cast <float>(xC + Hcyl *  cylinder_face[i61 + 4]);
					tmp->pa.y = static_cast <float>(yC + R_out_cyl * (cylinder_face[i61 + 3] - 1.0));
					tmp->pa.z = static_cast <float>(zC + R_out_cyl * (cylinder_face[i61 + 5] - 1.0));
					tmp->pb.x = static_cast <float>(xC + Hcyl *  cylinder_face[i61 + 7]);
					tmp->pb.y = static_cast <float>(yC + R_out_cyl * (cylinder_face[i61 + 6] - 1.0));
					tmp->pb.z = static_cast <float>(zC + R_out_cyl * (cylinder_face[i61 + 8] - 1.0));
					tmp->pc.x = static_cast <float>(xC + Hcyl *  cylinder_face[i61 + 10]);
					tmp->pc.y = static_cast <float>(yC + R_out_cyl * (cylinder_face[i61 + 9] - 1.0));
					tmp->pc.z = static_cast <float>(zC + R_out_cyl * (cylinder_face[i61 + 11] - 1.0));
				}
				if (i60 < 199) {
					tmp->next = new gCAD_STL;
				}
					tmp = tmp->next;
					i61 += 12;

			}
			// 2



		}
	}

	// Преобразование призмы в CAD STL.
	// Шесть граней и 12 треугольников.
	void PRISM2CAD_STL() {
		// 13.11.2020
		if (itypegeom == PRISM) {

			// освобождение оперативной памяти.
			clear_CAD_STL();

			inumber_triangles_for_CAD_STL_model = 12;

			// Модель обычная, небольшая.
			bbigCADmodel = false;
			root_big_CAD_STL_model_x = nullptr;
			root_big_CAD_STL_model_y = nullptr;
			root_big_CAD_STL_model_z = nullptr;

			// E
			root_CAD_STL = new gCAD_STL;
			gCAD_STL* tmp = root_CAD_STL;
			tmp->n.x = 1.0;
			tmp->n.y = 0.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xE);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yS);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			tmp->n.x = 1.0;
			tmp->n.y = 0.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xE);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zS);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			// W
			tmp->n.x = -1.0;
			tmp->n.y = 0.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xS);
			tmp->pb.y = static_cast <float>(yS);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xS);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			tmp->n.x = -1.0;
			tmp->n.y = 0.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xS);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xS);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zS);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			// TOP
			tmp->n.x = 0.0;
			tmp->n.y = 0.0;
			tmp->n.z = 1.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zE);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yS);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			tmp->n.x = 0.0;
			tmp->n.y = 0.0;
			tmp->n.z = 1.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zE);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xS);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			// Bottom
			tmp->n.x = 0.0;
			tmp->n.y = 0.0;
			tmp->n.z = -1.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zS);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yS);
			tmp->pc.z = static_cast <float>(zS);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			tmp->n.x = 0.0;
			tmp->n.y = 0.0;
			tmp->n.z = -1.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xS);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zS);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zS);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			// N
			tmp->n.x = 0.0;
			tmp->n.y = 1.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yE);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xS);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			tmp->n.x = 0.0;
			tmp->n.y = 1.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yE);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yE);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yE);
			tmp->pc.z = static_cast <float>(zS);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			// S
			tmp->n.x = 0.0;
			tmp->n.y = -1.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yS);
			tmp->pb.z = static_cast <float>(zE);
			tmp->pc.x = static_cast <float>(xS);
			tmp->pc.y = static_cast <float>(yS);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = new gCAD_STL;
			tmp = tmp->next;
			tmp->n.x = 0.0;
			tmp->n.y = -1.0;
			tmp->n.z = 0.0;
			tmp->pa.x = static_cast <float>(xS);
			tmp->pa.y = static_cast <float>(yS);
			tmp->pa.z = static_cast <float>(zS);
			tmp->pb.x = static_cast <float>(xE);
			tmp->pb.y = static_cast <float>(yS);
			tmp->pb.z = static_cast <float>(zS);
			tmp->pc.x = static_cast <float>(xE);
			tmp->pc.y = static_cast <float>(yS);
			tmp->pc.z = static_cast <float>(zE);
			tmp->next = nullptr;
			tmp = nullptr;
		}
		if ((itypegeom == POLYGON)&&(nsizei == 4)) {

			// освобождение оперативной памяти.
			clear_CAD_STL();

			inumber_triangles_for_CAD_STL_model = 12;

			// Модель обычная, небольшая.
			bbigCADmodel = false;
			root_big_CAD_STL_model_x = nullptr;
			root_big_CAD_STL_model_y = nullptr;
			root_big_CAD_STL_model_z = nullptr;


			switch (iPlane_obj2) {
			case XY_PLANE:
				if (hi[0] > 0.0) {
					// E
					root_CAD_STL = new gCAD_STL;
					gCAD_STL* tmp = root_CAD_STL;
					tmp->pa.x = static_cast <float>(xi[1]); //xE;
					tmp->pa.y = static_cast <float>(yi[1]); //yS;
					tmp->pa.z = static_cast <float>(zi[1]); //zS;
					tmp->pb.x = static_cast <float>(xi[2]); //xE;
					tmp->pb.y = static_cast <float>(yi[2]); //yE;
					tmp->pb.z = static_cast <float>(zi[2] + hi[2]); //zE;
					tmp->pc.x = static_cast <float>(xi[1]); //xE;
					tmp->pc.y = static_cast <float>(yi[1]); //yS;
					tmp->pc.z = static_cast <float>(zi[1] + hi[1]); //zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[1]); // xE;
					tmp->pa.y = static_cast <float>(yi[1]); // yS;
					tmp->pa.z = static_cast <float>(zi[1]); // zS;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2]); // zS;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2] + hi[2]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// W
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[0]); // xS;
					tmp->pb.y = static_cast <float>(yi[0]); // yS;
					tmp->pb.z = static_cast <float>(zi[0] + hi[0]); // zE;
					tmp->pc.x = static_cast <float>(xi[3]); // xS;
					tmp->pc.y = static_cast <float>(yi[3]); // yE;
					tmp->pc.z = static_cast <float>(zi[3] + hi[3]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[3]); // xS;
					tmp->pb.y = static_cast <float>(yi[3]); // yE;
					tmp->pb.z = static_cast <float>(zi[3] + hi[3]); // zE;
					tmp->pc.x = static_cast <float>(xi[3]); // xS;
					tmp->pc.y = static_cast <float>(yi[3]); // yE;
					tmp->pc.z = static_cast <float>(zi[3]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// TOP
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zE;
					tmp->pb.x = static_cast <float>(xi[1]); // xE;
					tmp->pb.y = static_cast <float>(yi[1]); // yS;
					tmp->pb.z = static_cast <float>(zi[1] + hi[1]); // zE;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2] + hi[2]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zE;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2] + hi[2]); // zE;
					tmp->pc.x = static_cast <float>(xi[3]); // xS;
					tmp->pc.y = static_cast <float>(yi[3]); // yE;
					tmp->pc.z = static_cast <float>(zi[3] + hi[3]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// Bottom
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2]); // zS;
					tmp->pc.x = static_cast <float>(xi[1]); // xE;
					tmp->pc.y = static_cast <float>(yi[1]); // yS;
					tmp->pc.z = static_cast <float>(zi[1]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[3]); // xS;
					tmp->pb.y = static_cast <float>(yi[3]); // yE;
					tmp->pb.z = static_cast <float>(zi[3]); // zS;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// N
					tmp->pa.x = static_cast <float>(xi[3]); // xS;
					tmp->pa.y = static_cast <float>(yi[3]); // yE;
					tmp->pa.z = static_cast <float>(zi[3]); // zS;
					tmp->pb.x = static_cast <float>(xi[3]); // xS;
					tmp->pb.y = static_cast <float>(yi[3]); // yE;
					tmp->pb.z = static_cast <float>(zi[3] + hi[3]); // zE;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2] + hi[2]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[3]); // xS;
					tmp->pa.y = static_cast <float>(yi[3]); // yE;
					tmp->pa.z = static_cast <float>(zi[3]); // zS;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2] + hi[2]); // zE;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// S
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[1]); // xE;
					tmp->pb.y = static_cast <float>(yi[1]); // yS;
					tmp->pb.z = static_cast <float>(zi[1] + hi[1]); // zE;
					tmp->pc.x = static_cast <float>(xi[0]); // xS;
					tmp->pc.y = static_cast <float>(yi[0]); // yS;
					tmp->pc.z = static_cast <float>(zi[0] + hi[0]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[1]); // xE;
					tmp->pb.y = static_cast <float>(yi[1]); // yS;
					tmp->pb.z = static_cast <float>(zi[1]); // zS;
					tmp->pc.x = static_cast <float>(xi[1]); // xE;
					tmp->pc.y = static_cast <float>(yi[1]); // yS;
					tmp->pc.z = static_cast <float>(zi[1] + hi[1]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = nullptr;
					tmp = nullptr;
				}
				else {
				    // Высота полигона отрицательна. 
					// E
					root_CAD_STL = new gCAD_STL;
					gCAD_STL* tmp = root_CAD_STL;
					tmp->pa.x = static_cast <float>(xi[1]); //xE;
					tmp->pa.y = static_cast <float>(yi[1]); //yS;
					tmp->pa.z = static_cast <float>(zi[1]+hi[1]); //zS;
					tmp->pb.x = static_cast <float>(xi[2]); //xE;
					tmp->pb.y = static_cast <float>(yi[2]); //yE;
					tmp->pb.z = static_cast <float>(zi[2]); //zE;
					tmp->pc.x = static_cast <float>(xi[1]); //xE;
					tmp->pc.y = static_cast <float>(yi[1]); //yS;
					tmp->pc.z = static_cast <float>(zi[1]); //zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[1]); // xE;
					tmp->pa.y = static_cast <float>(yi[1]); // yS;
					tmp->pa.z = static_cast <float>(zi[1]+hi[1]); // zS;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2]+hi[2]); // zS;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// W
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]+hi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[0]); // xS;
					tmp->pb.y = static_cast <float>(yi[0]); // yS;
					tmp->pb.z = static_cast <float>(zi[0]); // zE;
					tmp->pc.x = static_cast <float>(xi[3]); // xS;
					tmp->pc.y = static_cast <float>(yi[3]); // yE;
					tmp->pc.z = static_cast <float>(zi[3]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[3]); // xS;
					tmp->pb.y = static_cast <float>(yi[3]); // yE;
					tmp->pb.z = static_cast <float>(zi[3]); // zE;
					tmp->pc.x = static_cast <float>(xi[3]); // xS;
					tmp->pc.y = static_cast <float>(yi[3]); // yE;
					tmp->pc.z = static_cast <float>(zi[3] + hi[3]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// TOP
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zE;
					tmp->pb.x = static_cast <float>(xi[1]); // xE;
					tmp->pb.y = static_cast <float>(yi[1]); // yS;
					tmp->pb.z = static_cast <float>(zi[1]); // zE;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0]); // zE;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2]); // zE;
					tmp->pc.x = static_cast <float>(xi[3]); // xS;
					tmp->pc.y = static_cast <float>(yi[3]); // yE;
					tmp->pc.z = static_cast <float>(zi[3]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// Bottom
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2] + hi[2]); // zS;
					tmp->pc.x = static_cast <float>(xi[1]); // xE;
					tmp->pc.y = static_cast <float>(yi[1]); // yS;
					tmp->pc.z = static_cast <float>(zi[1] + hi[1]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[3]); // xS;
					tmp->pb.y = static_cast <float>(yi[3]); // yE;
					tmp->pb.z = static_cast <float>(zi[3] + hi[3]); // zS;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2] + hi[2]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// N
					tmp->pa.x = static_cast <float>(xi[3]); // xS;
					tmp->pa.y = static_cast <float>(yi[3]); // yE;
					tmp->pa.z = static_cast <float>(zi[3] + hi[3]); // zS;
					tmp->pb.x = static_cast <float>(xi[3]); // xS;
					tmp->pb.y = static_cast <float>(yi[3]); // yE;
					tmp->pb.z = static_cast <float>(zi[3]); // zE;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[3]); // xS;
					tmp->pa.y = static_cast <float>(yi[3]); // yE;
					tmp->pa.z = static_cast <float>(zi[3] + hi[3]); // zS;
					tmp->pb.x = static_cast <float>(xi[2]); // xE;
					tmp->pb.y = static_cast <float>(yi[2]); // yE;
					tmp->pb.z = static_cast <float>(zi[2]); // zE;
					tmp->pc.x = static_cast <float>(xi[2]); // xE;
					tmp->pc.y = static_cast <float>(yi[2]); // yE;
					tmp->pc.z = static_cast <float>(zi[2] + hi[2]); // zS;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					// S
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[1]); // xE;
					tmp->pb.y = static_cast <float>(yi[1]); // yS;
					tmp->pb.z = static_cast <float>(zi[1]); // zE;
					tmp->pc.x = static_cast <float>(xi[0]); // xS;
					tmp->pc.y = static_cast <float>(yi[0]); // yS;
					tmp->pc.z = static_cast <float>(zi[0]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = new gCAD_STL;
					tmp = tmp->next;
					tmp->pa.x = static_cast <float>(xi[0]); // xS;
					tmp->pa.y = static_cast <float>(yi[0]); // yS;
					tmp->pa.z = static_cast <float>(zi[0] + hi[0]); // zS;
					tmp->pb.x = static_cast <float>(xi[1]); // xE;
					tmp->pb.y = static_cast <float>(yi[1]); // yS;
					tmp->pb.z = static_cast <float>(zi[1] + hi[1]); // zS;
					tmp->pc.x = static_cast <float>(xi[1]); // xE;
					tmp->pc.y = static_cast <float>(yi[1]); // yS;
					tmp->pc.z = static_cast <float>(zi[1]); // zE;
					calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
					tmp->next = nullptr;
					tmp = nullptr;

				}

				break;
				case XZ_PLANE:
					if (hi[0] > 0.0) {
						// E
						root_CAD_STL = new gCAD_STL;
						gCAD_STL* tmp = root_CAD_STL;
						tmp->pa.x = static_cast <float>(xi[1]); //xE;
						tmp->pa.y = static_cast <float>(yi[1]); //yS;
						tmp->pa.z = static_cast <float>(zi[1]); //zS;
						tmp->pb.x = static_cast <float>(xi[2]); //xE;
						tmp->pb.y = static_cast <float>(yi[2] + hi[2]); //yE;
						tmp->pb.z = static_cast <float>(zi[2]); //zS;
						tmp->pc.x = static_cast <float>(xi[1]); //xE;
						tmp->pc.y = static_cast <float>(yi[1] + hi[1]); //yE;
						tmp->pc.z = static_cast <float>(zi[1]); //zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[1]); // xE;
						tmp->pa.y = static_cast <float>(yi[1]); // yS;
						tmp->pa.z = static_cast <float>(zi[1]); // zS;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2]); // yS;
						tmp->pb.z = static_cast <float>(zi[2]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// W
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[0]); // xS;
						tmp->pb.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pb.z = static_cast <float>(zi[0]); // zE;
						tmp->pc.x = static_cast <float>(xi[3]); // xS;
						tmp->pc.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pc.z = static_cast <float>(zi[3]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[3]); // xS;
						tmp->pb.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pb.z = static_cast <float>(zi[3]); // zE;
						tmp->pc.x = static_cast <float>(xi[3]); // xS;
						tmp->pc.y = static_cast <float>(yi[3]); // yE;
						tmp->pc.z = static_cast <float>(zi[3]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// TOP
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zE;
						tmp->pb.x = static_cast <float>(xi[1]); // xE;
						tmp->pb.y = static_cast <float>(yi[1] + hi[1]); // yS;
						tmp->pb.z = static_cast <float>(zi[1]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zE;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zE;
						tmp->pc.x = static_cast <float>(xi[3]); // xS;
						tmp->pc.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pc.z = static_cast <float>(zi[3]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// Bottom
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zS;
						tmp->pc.x = static_cast <float>(xi[1]); // xE;
						tmp->pc.y = static_cast <float>(yi[1]); // yS;
						tmp->pc.z = static_cast <float>(zi[1]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[3]); // xS;
						tmp->pb.y = static_cast <float>(yi[3]); // yE;
						tmp->pb.z = static_cast <float>(zi[3]); // zS;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// N
						tmp->pa.x = static_cast <float>(xi[3]); // xS;
						tmp->pa.y = static_cast <float>(yi[3]); // yE;
						tmp->pa.z = static_cast <float>(zi[3]); // zS;
						tmp->pb.x = static_cast <float>(xi[3]); // xS;
						tmp->pb.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pb.z = static_cast <float>(zi[3]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[3]); // xS;
						tmp->pa.y = static_cast <float>(yi[3]); // yE;
						tmp->pa.z = static_cast <float>(zi[3]); // zS;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// S
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[1]); // xE;
						tmp->pb.y = static_cast <float>(yi[1] + hi[1]); // yS;
						tmp->pb.z = static_cast <float>(zi[1]); // zE;
						tmp->pc.x = static_cast <float>(xi[0]); // xS;
						tmp->pc.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pc.z = static_cast <float>(zi[0]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[1]); // xE;
						tmp->pb.y = static_cast <float>(yi[1]); // yS;
						tmp->pb.z = static_cast <float>(zi[1]); // zS;
						tmp->pc.x = static_cast <float>(xi[1]); // xE;
						tmp->pc.y = static_cast <float>(yi[1] + hi[1]); // yS;
						tmp->pc.z = static_cast <float>(zi[1]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = nullptr;
						tmp = nullptr;
					}
					else {
						// Высота полигона отрицательна. 
						// E
						root_CAD_STL = new gCAD_STL;
						gCAD_STL* tmp = root_CAD_STL;
						tmp->pa.x = static_cast <float>(xi[1]); //xE;
						tmp->pa.y = static_cast <float>(yi[1] + hi[1]); //yS;
						tmp->pa.z = static_cast <float>(zi[1]); //zS;
						tmp->pb.x = static_cast <float>(xi[2]); //xE;
						tmp->pb.y = static_cast <float>(yi[2]); //yE;
						tmp->pb.z = static_cast <float>(zi[2]); //zE;
						tmp->pc.x = static_cast <float>(xi[1]); //xE;
						tmp->pc.y = static_cast <float>(yi[1]); //yS;
						tmp->pc.z = static_cast <float>(zi[1]); //zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[1]); // xE;
						tmp->pa.y = static_cast <float>(yi[1] + hi[1]); // yS;
						tmp->pa.z = static_cast <float>(zi[1]); // zS;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zS;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// W
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[0]); // xS;
						tmp->pb.y = static_cast <float>(yi[0]); // yS;
						tmp->pb.z = static_cast <float>(zi[0]); // zE;
						tmp->pc.x = static_cast <float>(xi[3]); // xS;
						tmp->pc.y = static_cast <float>(yi[3]); // yE;
						tmp->pc.z = static_cast <float>(zi[3]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[3]); // xS;
						tmp->pb.y = static_cast <float>(yi[3]); // yE;
						tmp->pb.z = static_cast <float>(zi[3]); // zE;
						tmp->pc.x = static_cast <float>(xi[3]); // xS;
						tmp->pc.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pc.z = static_cast <float>(zi[3]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// TOP
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zE;
						tmp->pb.x = static_cast <float>(xi[1]); // xE;
						tmp->pb.y = static_cast <float>(yi[1]); // yS;
						tmp->pb.z = static_cast <float>(zi[1]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zE;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zE;
						tmp->pc.x = static_cast <float>(xi[3]); // xS;
						tmp->pc.y = static_cast <float>(yi[3]); // yE;
						tmp->pc.z = static_cast <float>(zi[3]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// Bottom
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zS;
						tmp->pc.x = static_cast <float>(xi[1]); // xE;
						tmp->pc.y = static_cast <float>(yi[1] + hi[1]); // yS;
						tmp->pc.z = static_cast <float>(zi[1]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[3]); // xS;
						tmp->pb.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pb.z = static_cast <float>(zi[3]); // zS;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// N
						tmp->pa.x = static_cast <float>(xi[3]); // xS;
						tmp->pa.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pa.z = static_cast <float>(zi[3]); // zS;
						tmp->pb.x = static_cast <float>(xi[3]); // xS;
						tmp->pb.y = static_cast <float>(yi[3]); // yE;
						tmp->pb.z = static_cast <float>(zi[3]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[3]); // xS;
						tmp->pa.y = static_cast <float>(yi[3] + hi[3]); // yE;
						tmp->pa.z = static_cast <float>(zi[3]); // zS;
						tmp->pb.x = static_cast <float>(xi[2]); // xE;
						tmp->pb.y = static_cast <float>(yi[2]); // yE;
						tmp->pb.z = static_cast <float>(zi[2]); // zE;
						tmp->pc.x = static_cast <float>(xi[2]); // xE;
						tmp->pc.y = static_cast <float>(yi[2] + hi[2]); // yE;
						tmp->pc.z = static_cast <float>(zi[2]); // zS;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						// S
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[1]); // xE;
						tmp->pb.y = static_cast <float>(yi[1]); // yS;
						tmp->pb.z = static_cast <float>(zi[1]); // zE;
						tmp->pc.x = static_cast <float>(xi[0]); // xS;
						tmp->pc.y = static_cast <float>(yi[0]); // yS;
						tmp->pc.z = static_cast <float>(zi[0]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = new gCAD_STL;
						tmp = tmp->next;
						tmp->pa.x = static_cast <float>(xi[0]); // xS;
						tmp->pa.y = static_cast <float>(yi[0] + hi[0]); // yS;
						tmp->pa.z = static_cast <float>(zi[0]); // zS;
						tmp->pb.x = static_cast <float>(xi[1]); // xE;
						tmp->pb.y = static_cast <float>(yi[1] + hi[1]); // yS;
						tmp->pb.z = static_cast <float>(zi[1]); // zS;
						tmp->pc.x = static_cast <float>(xi[1]); // xE;
						tmp->pc.y = static_cast <float>(yi[1]); // yS;
						tmp->pc.z = static_cast <float>(zi[1]); // zE;
						calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
						tmp->next = nullptr;
						tmp = nullptr;
					}
					break;
					case YZ_PLANE:
						if (hi[0] > 0.0) {
							// E
							root_CAD_STL = new gCAD_STL;
							gCAD_STL* tmp = root_CAD_STL;
							tmp->pa.x = static_cast <float>(xi[1]); //xE;
							tmp->pa.y = static_cast <float>(yi[1]); //yS;
							tmp->pa.z = static_cast <float>(zi[1]); //zS;
							tmp->pb.x = static_cast <float>(xi[2] + hi[2]); //xE;
							tmp->pb.y = static_cast <float>(yi[2]); //yE;
							tmp->pb.z = static_cast <float>(zi[2]); //zS;
							tmp->pc.x = static_cast <float>(xi[1] + hi[1]); //xE;
							tmp->pc.y = static_cast <float>(yi[1]); //yE;
							tmp->pc.z = static_cast <float>(zi[1]); //zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[1]); // xE;
							tmp->pa.y = static_cast <float>(yi[1]); // yS;
							tmp->pa.z = static_cast <float>(zi[1]); // zS;
							tmp->pb.x = static_cast <float>(xi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yS;
							tmp->pb.z = static_cast <float>(zi[2]); // zE;
							tmp->pc.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// W
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pb.y = static_cast <float>(yi[0]); // yS;
							tmp->pb.z = static_cast <float>(zi[0]); // zE;
							tmp->pc.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pc.y = static_cast <float>(yi[3]); // yE;
							tmp->pc.z = static_cast <float>(zi[3]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pb.y = static_cast <float>(yi[3]); // yE;
							tmp->pb.z = static_cast <float>(zi[3]); // zE;
							tmp->pc.x = static_cast <float>(xi[3]); // xS;
							tmp->pc.y = static_cast <float>(yi[3]); // yE;
							tmp->pc.z = static_cast <float>(zi[3]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// TOP
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zE;
							tmp->pb.x = static_cast <float>(xi[1] + hi[1]); // xE;
							tmp->pb.y = static_cast <float>(yi[1]); // yS;
							tmp->pb.z = static_cast <float>(zi[1]); // zE;
							tmp->pc.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zE;
							tmp->pb.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zE;
							tmp->pc.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pc.y = static_cast <float>(yi[3]); // yE;
							tmp->pc.z = static_cast <float>(zi[3]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// Bottom
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zS;
							tmp->pc.x = static_cast <float>(xi[1]); // xE;
							tmp->pc.y = static_cast <float>(yi[1]); // yS;
							tmp->pc.z = static_cast <float>(zi[1]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[3]); // xS;
							tmp->pb.y = static_cast <float>(yi[3]); // yE;
							tmp->pb.z = static_cast <float>(zi[3]); // zS;
							tmp->pc.x = static_cast <float>(xi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// N
							tmp->pa.x = static_cast <float>(xi[3]); // xS;
							tmp->pa.y = static_cast <float>(yi[3]); // yE;
							tmp->pa.z = static_cast <float>(zi[3]); // zS;
							tmp->pb.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pb.y = static_cast <float>(yi[3]); // yE;
							tmp->pb.z = static_cast <float>(zi[3]); // zE;
							tmp->pc.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[3]); // xS;
							tmp->pa.y = static_cast <float>(yi[3]); // yE;
							tmp->pa.z = static_cast <float>(zi[3]); // zS;
							tmp->pb.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zE;
							tmp->pc.x = static_cast <float>(xi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// S
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[1] + hi[1]); // xE;
							tmp->pb.y = static_cast <float>(yi[1]); // yS;
							tmp->pb.z = static_cast <float>(zi[1]); // zE;
							tmp->pc.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pc.y = static_cast <float>(yi[0]); // yS;
							tmp->pc.z = static_cast <float>(zi[0]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[1]); // xE;
							tmp->pb.y = static_cast <float>(yi[1]); // yS;
							tmp->pb.z = static_cast <float>(zi[1]); // zS;
							tmp->pc.x = static_cast <float>(xi[1] + hi[1]); // xE;
							tmp->pc.y = static_cast <float>(yi[1]); // yS;
							tmp->pc.z = static_cast <float>(zi[1]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = nullptr;
							tmp = nullptr;
						}
						else {
							// Высота полигона отрицательна. 
							// E
							root_CAD_STL = new gCAD_STL;
							gCAD_STL* tmp = root_CAD_STL;
							tmp->pa.x = static_cast <float>(xi[1] + hi[1]); //xE;
							tmp->pa.y = static_cast <float>(yi[1]); //yS;
							tmp->pa.z = static_cast <float>(zi[1]); //zS;
							tmp->pb.x = static_cast <float>(xi[2]); //xE;
							tmp->pb.y = static_cast <float>(yi[2]); //yE;
							tmp->pb.z = static_cast <float>(zi[2]); //zE;
							tmp->pc.x = static_cast <float>(xi[1]); //xE;
							tmp->pc.y = static_cast <float>(yi[1]); //yS;
							tmp->pc.z = static_cast <float>(zi[1]); //zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[1] + hi[1]); // xE;
							tmp->pa.y = static_cast <float>(yi[1]); // yS;
							tmp->pa.z = static_cast <float>(zi[1]); // zS;
							tmp->pb.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zS;
							tmp->pc.x = static_cast <float>(xi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// W
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[0]); // xS;
							tmp->pb.y = static_cast <float>(yi[0]); // yS;
							tmp->pb.z = static_cast <float>(zi[0]); // zE;
							tmp->pc.x = static_cast <float>(xi[3]); // xS;
							tmp->pc.y = static_cast <float>(yi[3]); // yE;
							tmp->pc.z = static_cast <float>(zi[3]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[3]); // xS;
							tmp->pb.y = static_cast <float>(yi[3]); // yE;
							tmp->pb.z = static_cast <float>(zi[3]); // zE;
							tmp->pc.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pc.y = static_cast <float>(yi[3]); // yE;
							tmp->pc.z = static_cast <float>(zi[3]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// TOP
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zE;
							tmp->pb.x = static_cast <float>(xi[1]); // xE;
							tmp->pb.y = static_cast <float>(yi[1]); // yS;
							tmp->pb.z = static_cast <float>(zi[1]); // zE;
							tmp->pc.x = static_cast <float>(xi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zE;
							tmp->pb.x = static_cast <float>(xi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zE;
							tmp->pc.x = static_cast <float>(xi[3]); // xS;
							tmp->pc.y = static_cast <float>(yi[3]); // yE;
							tmp->pc.z = static_cast <float>(zi[3]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// Bottom
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zS;
							tmp->pc.x = static_cast <float>(xi[1] + hi[1]); // xE;
							tmp->pc.y = static_cast <float>(yi[1]); // yS;
							tmp->pc.z = static_cast <float>(zi[1]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pb.y = static_cast <float>(yi[3]); // yE;
							tmp->pb.z = static_cast <float>(zi[3]); // zS;
							tmp->pc.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// N
							tmp->pa.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pa.y = static_cast <float>(yi[3]); // yE;
							tmp->pa.z = static_cast <float>(zi[3]); // zS;
							tmp->pb.x = static_cast <float>(xi[3]); // xS;
							tmp->pb.y = static_cast <float>(yi[3]); // yE;
							tmp->pb.z = static_cast <float>(zi[3]); // zE;
							tmp->pc.x = static_cast <float>(xi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[3] + hi[3]); // xS;
							tmp->pa.y = static_cast <float>(yi[3]); // yE;
							tmp->pa.z = static_cast <float>(zi[3]); // zS;
							tmp->pb.x = static_cast <float>(xi[2]); // xE;
							tmp->pb.y = static_cast <float>(yi[2]); // yE;
							tmp->pb.z = static_cast <float>(zi[2]); // zE;
							tmp->pc.x = static_cast <float>(xi[2] + hi[2]); // xE;
							tmp->pc.y = static_cast <float>(yi[2]); // yE;
							tmp->pc.z = static_cast <float>(zi[2]); // zS;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							// S
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[1]); // xE;
							tmp->pb.y = static_cast <float>(yi[1]); // yS;
							tmp->pb.z = static_cast <float>(zi[1]); // zE;
							tmp->pc.x = static_cast <float>(xi[0]); // xS;
							tmp->pc.y = static_cast <float>(yi[0]); // yS;
							tmp->pc.z = static_cast <float>(zi[0]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = new gCAD_STL;
							tmp = tmp->next;
							tmp->pa.x = static_cast <float>(xi[0] + hi[0]); // xS;
							tmp->pa.y = static_cast <float>(yi[0]); // yS;
							tmp->pa.z = static_cast <float>(zi[0]); // zS;
							tmp->pb.x = static_cast <float>(xi[1] + hi[1]); // xE;
							tmp->pb.y = static_cast <float>(yi[1]); // yS;
							tmp->pb.z = static_cast <float>(zi[1]); // zS;
							tmp->pc.x = static_cast <float>(xi[1]); // xE;
							tmp->pc.y = static_cast <float>(yi[1]); // yS;
							tmp->pc.z = static_cast <float>(zi[1]); // zE;
							calculateNormal_for_triangle(tmp->n, tmp->pa, tmp->pb, tmp->pc);
							tmp->next = nullptr;
							tmp = nullptr;
						}
						break;
			}
		}
		if (itypegeom == CYLINDER) {
			CYLINDER2CAD_STL();
		}
}


	

	// Центр призмы с точкой p2 образует направляющий отрезок.
	// Данный метод вычисляет расстояние от центра призмы до точки пересечения отрезка с гранью призмы.
	doublereal distance_center_to_gran_PRISM_along_orientation_vector(TOCHKA_FLOAT p2)
	{
		if (itypegeom == PRISM) {

			TOCHKA_FLOAT pc, p;
			pc.x = static_cast <float>((0.5 * (xS + xE)));
			pc.y = static_cast <float>((0.5 * (yS + yE)));
			pc.z = static_cast <float>((0.5 * (zS + zE)));

			PRISM2CAD_STL();

			gCAD_STL* tmp = root_CAD_STL;

			while (tmp != nullptr) {

				if (trnLineFacet(pc, p2,  tmp->pa, tmp->pb, tmp->pc, p)) {
					tmp = nullptr;
					return sqrt((pc.x-p.x)* (pc.x - p.x)+ (pc.y - p.y) * (pc.y - p.y)+ (pc.z - p.z) * (pc.z - p.z));
				}
					

				tmp = tmp->next;
			}
		}
		return 1.0e20; // Нет пересечения, бесконечное расстояние 21,02,2021
	}

	// Вычисляет окаймляющую призму для CAD STL объекта.
	void minimaxCAD_STL(TOCHKA_FLOAT &pmin, TOCHKA_FLOAT &pmax) {

		if (!bminmax) {

			pmax.x = -1.0e30f;
			pmax.y = -1.0e30f;
			pmax.z = -1.0e30f;

			pmin.x = 1.0e30f;
			pmin.y = 1.0e30f;
			pmin.z = 1.0e30f;


			gCAD_STL* tmp = root_CAD_STL;

			while (tmp != nullptr) {

				// maximum

				if (tmp->pa.x > pmax.x) {
					pmax.x = tmp->pa.x;
				}
				if (tmp->pb.x > pmax.x) {
					pmax.x = tmp->pb.x;
				}
				if (tmp->pc.x > pmax.x) {
					pmax.x = tmp->pc.x;
				}

				if (tmp->pa.y > pmax.y) {
					pmax.y = tmp->pa.y;
				}
				if (tmp->pb.y > pmax.y) {
					pmax.y = tmp->pb.y;
				}
				if (tmp->pc.y > pmax.y) {
					pmax.y = tmp->pc.y;
				}

				if (tmp->pa.z > pmax.z) {
					pmax.z = tmp->pa.z;
				}
				if (tmp->pb.z > pmax.z) {
					pmax.z = tmp->pb.z;
				}
				if (tmp->pc.z > pmax.z) {
					pmax.z = tmp->pc.z;
				}

				// minimum

				if (tmp->pa.x < pmin.x) {
					pmin.x = tmp->pa.x;
				}
				if (tmp->pb.x < pmin.x) {
					pmin.x = tmp->pb.x;
				}
				if (tmp->pc.x < pmin.x) {
					pmin.x = tmp->pc.x;
				}

				if (tmp->pa.y < pmin.y) {
					pmin.y = tmp->pa.y;
				}
				if (tmp->pb.y < pmin.y) {
					pmin.y = tmp->pb.y;
				}
				if (tmp->pc.y < pmin.y) {
					pmin.y = tmp->pc.y;
				}

				if (tmp->pa.z < pmin.z) {
					pmin.z = tmp->pa.z;
				}
				if (tmp->pb.z < pmin.z) {
					pmin.z = tmp->pb.z;
				}
				if (tmp->pc.z < pmin.z) {
					pmin.z = tmp->pc.z;
				}


				tmp = tmp->next;
			}

			pmin_mem = pmin;
			pmax_mem = pmax;
			bminmax = true;

		}
		else {
			pmin = pmin_mem;
			pmax = pmax_mem;
		}
	}

	// Принадлежит ли точка внутренности CAD STL объекта.
	bool in_CAD_STL_check(TOCHKA p, int& k, int ib) {

		// 13.11.2020

		TOCHKA_FLOAT pf;
		pf.x = static_cast <float>(p.x);
		pf.y = static_cast <float>(p.y);
		pf.z = static_cast <float>(p.z);

		TOCHKA_FLOAT pmin;
		TOCHKA_FLOAT pmax;

		minimaxCAD_STL(pmin, pmax);



		if ((pf.x > pmin.x) && (pf.x < pmax.x) && (pf.y > pmin.y) && (pf.y < pmax.y) && (pf.z > pmin.z) && (pf.z < pmax.z)) {


			// решение принимается большинством голосов за три теста.
			bool bout1 = false, bout2 = false, bout3 = false;

			TOCHKA_FLOAT pf2, pintersect;
			pf2.x = static_cast <float>(p.x);
			pf2.y = static_cast <float>(p.y);
			pf2.z = static_cast <float>((pmax.z+ 0.02* (pmax.z - pmin.z)));//+0.02

			int inumberintersect = 0;

			if (!bbigCADmodel) {

				gCAD_STL* tmp = root_CAD_STL;

				while (tmp != nullptr) {

					if (trnLineFacet(pf, pf2, tmp->pa, tmp->pb, tmp->pc, pintersect)) {
						inumberintersect++;
					}

					tmp = tmp->next;
				}

			}
			else {
				// Большая .stl модель.
				//bigCAD_STL* tmp1 = root_big_CAD_STL_model_z;

				gCAD_STL* tmp = nullptr;//, * tmp_loc = nullptr;
				
				/*while (tmp1 != nullptr) {

					if ((p.x>=tmp1->min1)&& (p.x <= tmp1->max1)&& (p.y >= tmp1->min2) && (p.y <= tmp1->max2)) {
						tmp_loc= tmp1->root_CAD_STL;
					}

					tmp1 = tmp1->next;
				}*/
				
				TOCHKA_FLOAT pmin;
				TOCHKA_FLOAT pmax;

				minimaxCAD_STL(pmin, pmax);

				//int key_x = my_imax(0,my_imin((int)(p.x- root_big_CAD_STL_model_z_hash_table[0][0]->min1)/(root_big_CAD_STL_model_z_hash_table[size_pattern - 1][0]->max1 -root_big_CAD_STL_model_z_hash_table[0][0]->min1),size_pattern-1));
				//int key_y = my_imax(0, my_imin((int)(p.y - root_big_CAD_STL_model_z_hash_table[0][0]->min2) / (root_big_CAD_STL_model_z_hash_table[0][size_pattern - 1]->max2 - root_big_CAD_STL_model_z_hash_table[0][0]->min2), size_pattern - 1));
				const double dsize_pattern = 1.0 * size_pattern;

				integer key_x = my_imax(0, my_imin(static_cast<integer>((dsize_pattern)*(p.x - pmin.x) / (pmax.x - pmin.x)), size_pattern - 1));
				integer key_y = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (p.y - pmin.y) / (pmax.y - pmin.y)), size_pattern - 1));
				
				/*if ((tmp_loc!=nullptr)/*&&(root_big_CAD_STL_model_z_hash_table[key_x][key_y]->root_CAD_STL != nullptr)*//*) {
					std::wcout << "p.x=" << p.x << " min=" << pmin.x << " max=" << pmax.x <<  std::endl;
					std::wcout << "p.y=" << p.y << " min=" << pmin.y << " max=" << pmax.y << std::endl;
					std::wcout << "x=" << key_x << " y=" << key_y << std::endl;
				}*/
				//else {
					//printf("ok");
				//}
				//tmp = tmp_loc;
				if (root_big_CAD_STL_model_z_hash_table[key_x][key_y] != nullptr) {
					//bigCAD_STL* tmp1t = root_big_CAD_STL_model_z_hash_table[key_x][key_y];
					//if ((p.x >= tmp1t->min1) && (p.x <= tmp1t->max1) && (p.y >= tmp1t->min2) && (p.y <= tmp1t->max2)) {
						tmp = root_big_CAD_STL_model_z_hash_table[key_x][key_y]->root_CAD_STL;
					//}
				}
				/*
				if (tmp_loc != tmp) {
					std::wcout << "p.x=" << p.x << " min=" << pmin.x << " max=" << pmax.x << std::endl;
					std::wcout << "p.y=" << p.y << " min=" << pmin.y << " max=" << pmax.y << std::endl;
					std::wcout << "x=" << key_x << " y=" << key_y << std::endl;

					//std::wcout << "x=" << static_cast<integer>((dsize_pattern) * (p.x - pmin.x) / (pmax.x - pmin.x)) - 32 << " y=" << static_cast<integer>((dsize_pattern) * (p.y - pmin.y) / (pmax.y - pmin.y)) - 32 << std::endl;


					for (int i_3 = 0; i_3 < size_pattern; ++i_3) {
						for (int i_4 = 0; i_4 < size_pattern; i_4++) {
							if (root_big_CAD_STL_model_z_hash_table[i_3][i_4] != nullptr) {
								if (tmp_loc == root_big_CAD_STL_model_z_hash_table[i_3][i_4]->root_CAD_STL) {
									std::wcout << "i_3=" << i_3 << " i_4=" << i_4 << "nakonec found\n";
									system("pause");
								}
							}
						}
					}

					system("pause");
				}
				*/
				/*
				if ((tmp_loc != nullptr) /*&& (root_big_CAD_STL_model_z_hash_table[key_x][key_y]->root_CAD_STL != nullptr)*//*) {
					std::cout << "tmp " << tmp << " tmp_loc " << tmp_loc << std::endl;


					for (int i_3 = 0; i_3 < size_pattern; ++i_3) {
						for (int i_4 = 0; i_4 < size_pattern; i_4++) {
							if (root_big_CAD_STL_model_z_hash_table[i_3][i_4] != nullptr) {
								if (tmp_loc == root_big_CAD_STL_model_z_hash_table[i_3][i_4]->root_CAD_STL) {
									std::wcout << "i_3=" << i_3 << " i_4=" << i_4 << "nakonec found\n";
									system("pause");
								}
							}
						}
					}

					/*if (root_big_CAD_STL_model_z_hash_table[key_x + 1][key_y] != nullptr) {
						std::cout << "x+1 y0" << root_big_CAD_STL_model_z_hash_table[key_x + 1][key_y]->root_CAD_STL << std::endl;
					}
					if (root_big_CAD_STL_model_z_hash_table[key_x - 1][key_y] != nullptr) {
						std::cout << "x-1 y0" << root_big_CAD_STL_model_z_hash_table[key_x - 1][key_y]->root_CAD_STL << std::endl;
					}
					if (root_big_CAD_STL_model_z_hash_table[key_x][key_y + 1] != nullptr) {
						std::cout << "x0 y+1" << root_big_CAD_STL_model_z_hash_table[key_x][key_y + 1]->root_CAD_STL << std::endl;
					}
					if (root_big_CAD_STL_model_z_hash_table[key_x + 1][key_y + 1] != nullptr) {
						std::cout << "x+1 y+1" << root_big_CAD_STL_model_z_hash_table[key_x + 1][key_y + 1]->root_CAD_STL << std::endl;
					}
					if (root_big_CAD_STL_model_z_hash_table[key_x - 1][key_y + 1] != nullptr) {
						std::cout << "x+1 y+1" << root_big_CAD_STL_model_z_hash_table[key_x - 1][key_y + 1]->root_CAD_STL << std::endl;
					}*//*

					system("pause");
				}*/

				while (tmp != nullptr) {

					if (trnLineFacet(pf, pf2, tmp->pa, tmp->pb, tmp->pc, pintersect)) {
						inumberintersect++;
					}

					tmp = tmp->next;
				}

			}

			if (inumberintersect % 2 == 0) {
				// снаружи
				bout1 = true;
			}

			pf2.x = static_cast <float>((pmax.x+ 0.02*(pmax.x - pmin.x)));//+0.02
			pf2.y = static_cast <float>(p.y);
			pf2.z = static_cast <float>(p.z);

			inumberintersect = 0;
			
			if (!bbigCADmodel) {

				gCAD_STL* tmp = root_CAD_STL;

				while (tmp != nullptr) {

					if (trnLineFacet(pf, pf2, tmp->pa, tmp->pb, tmp->pc, pintersect)) {
						inumberintersect++;
					}

					tmp = tmp->next;
				}

			}
			else {
				// Большая .stl модель.
				//bigCAD_STL* tmp1 = root_big_CAD_STL_model_x;

				gCAD_STL* tmp = nullptr;

				/*while (tmp1 != nullptr) {

					if ((p.y >= tmp1->min1) && (p.y <= tmp1->max1) && (p.z >= tmp1->min2) && (p.z <= tmp1->max2)) {
						tmp = tmp1->root_CAD_STL;
					}

					tmp1 = tmp1->next;
				}*/
				
				TOCHKA_FLOAT pmin;
				TOCHKA_FLOAT pmax;

				minimaxCAD_STL(pmin, pmax);

				//int key_y = my_imax(0, my_imin((int)(p.y - root_big_CAD_STL_model_x_hash_table[0][0]->min1) / (root_big_CAD_STL_model_x_hash_table[size_pattern - 1][0]->max1 - root_big_CAD_STL_model_x_hash_table[0][0]->min1), size_pattern - 1));
				//int key_z = my_imax(0, my_imin((int)(p.z - root_big_CAD_STL_model_x_hash_table[0][0]->min2) / (root_big_CAD_STL_model_x_hash_table[0][size_pattern - 1]->max2 - root_big_CAD_STL_model_x_hash_table[0][0]->min2), size_pattern - 1));
				
				integer key_y = my_imax(0, my_imin(static_cast<integer>(((size_pattern)*(p.y - pmin.y)) / (pmax.y - pmin.y)), size_pattern - 1));
				integer key_z = my_imax(0, my_imin(static_cast<integer>(((size_pattern) * (p.z - pmin.z)) / (pmax.z - pmin.z)), size_pattern - 1));

				//if (root_big_CAD_STL_model_x_hash_table[key_y][key_z]->root_CAD_STL == nullptr) {
				//	std::wcout << "y=" << key_y << " z=" << key_z << std::endl;
				//}
				//else {
				//	printf("ok");
				//}

				if (root_big_CAD_STL_model_x_hash_table[key_y][key_z] != nullptr) {
					//bigCAD_STL* tmp1t = root_big_CAD_STL_model_x_hash_table[key_y][key_z];
					//if ((p.y >= tmp1t->min1) && (p.y <= tmp1t->max1) && (p.z >= tmp1t->min2) && (p.z <= tmp1t->max2)) {
						tmp = root_big_CAD_STL_model_x_hash_table[key_y][key_z]->root_CAD_STL;
					//}
				}
				
				while (tmp != nullptr) {

					if (trnLineFacet(pf, pf2, tmp->pa, tmp->pb, tmp->pc, pintersect)) {
						inumberintersect++;
					}

					tmp = tmp->next;
				}

			}

			if (inumberintersect % 2 == 0) {
				// снаружи
				bout2 = true;
			}

			pf2.x = static_cast <float>(p.x);
			pf2.y = static_cast <float>((pmax.y+ 0.02 * (pmax.y - pmin.y)));//+0.02
			pf2.z = static_cast <float>(p.z);

			inumberintersect = 0;

			if (!bbigCADmodel) {

				gCAD_STL*  tmp = root_CAD_STL;

			    while (tmp != nullptr) {

				   if (trnLineFacet(pf, pf2, tmp->pa, tmp->pb, tmp->pc, pintersect)) {
				      inumberintersect++;
				   }

				   tmp = tmp->next;
			    }

			}
			else {
				// Большая .stl модель.
				//bigCAD_STL* tmp1 = root_big_CAD_STL_model_y;

				gCAD_STL* tmp = nullptr;

				/*while (tmp1 != nullptr) {

					if ((p.x >= tmp1->min1) && (p.x <= tmp1->max1) && (p.z >= tmp1->min2) && (p.z <= tmp1->max2)) {
						tmp = tmp1->root_CAD_STL;
					}

				    tmp1 = tmp1->next;
				}*/
				
				TOCHKA_FLOAT pmin;
				TOCHKA_FLOAT pmax;

				minimaxCAD_STL(pmin, pmax);

				// 
				//int key_x = my_imax(0, my_imin((int)(p.x - root_big_CAD_STL_model_y_hash_table[0][0]->min1) / (root_big_CAD_STL_model_y_hash_table[size_pattern - 1][0]->max1 - root_big_CAD_STL_model_y_hash_table[0][0]->min1), size_pattern - 1));
				//int key_z = my_imax(0, my_imin((int)(p.z - root_big_CAD_STL_model_y_hash_table[0][0]->min2) / (root_big_CAD_STL_model_y_hash_table[0][size_pattern - 1]->max2 - root_big_CAD_STL_model_y_hash_table[0][0]->min2), size_pattern - 1));
				
				integer key_x = my_imax(0, my_imin(static_cast<integer>(1.0*(size_pattern)*(p.x - pmin.x) / (pmax.x - pmin.x)), size_pattern - 1));
				integer key_z = my_imax(0, my_imin(static_cast<integer>(1.0*(size_pattern) * (p.z - pmin.z) / (pmax.z - pmin.z)), size_pattern - 1));

				//if (root_big_CAD_STL_model_y_hash_table[key_x][key_z]->root_CAD_STL == nullptr) {
					//std::wcout << "x=" << key_x << " z=" << key_z << std::endl;
				//}
				//else {
				//	printf("ok");
				//}

				if (root_big_CAD_STL_model_y_hash_table[key_x][key_z] != nullptr) {
					//bigCAD_STL* tmp1t = root_big_CAD_STL_model_y_hash_table[key_x][key_z];
					//if ((p.x >= tmp1t->min1) && (p.x <= tmp1t->max1) && (p.z >= tmp1t->min2) && (p.z <= tmp1t->max2)) {
						tmp = root_big_CAD_STL_model_y_hash_table[key_x][key_z]->root_CAD_STL;
					//}
				}
				
				while (tmp != nullptr) {

					if (trnLineFacet(pf, pf2, tmp->pa, tmp->pb, tmp->pc, pintersect)) {
						inumberintersect++;
					}

					tmp = tmp->next;
				}

			}

			if (inumberintersect % 2 == 0) {
				// снаружи
				bout3 = true;
			}

			int iout = 0;
			if (bout1) iout++;
			if (bout2) iout++;
			if (bout3) iout++;

			
			if (iout >= 2) {
				// Голосование за то что точка снаружи.
				return false;
			}
			else {
				// Точка внутри.
				k = ib;
				return true;
			}
			
			/*
			if (iout == 3) {
				// Голосование за то что точка снаружи.
				return false;
			}
			else {
				// Точка внутри.
				k = ib;
				return true;
			}
			*/
		}
		else {
			return false;
		}

	}

	bool CAD_is_PRISM() {

		bool b1 = false, b2 = false, b3 = false;

		//integer iPlane_obj3 = -1;

		TOCHKA_FLOAT pmin, pmax;

		minimaxCAD_STL(pmin, pmax);

		double epsx = 0.001*(pmax.x - pmin.x);
		double epsy = 0.001*(pmax.y - pmin.y);
		double epsz = 0.001*(pmax.z - pmin.z);


		if (!bbigCADmodel) {

			gCAD_STL* tmp = root_CAD_STL;

			double x1 = 1.0e30, x2 = 1.0e30;
			bool bPlane = true;

			while (tmp != nullptr) {

				if (x1 > 1.0e29) {
					x1 = tmp->pa.x;
				}
				else {
					if (fabs(tmp->pa.x - x1) < epsx) {
					}
					else {
						if (x2 > 1.0e29) {
							x2 = tmp->pa.x;
						}
						else {
							if (fabs(tmp->pa.x - x2) < epsx) {
							}
							else {
								bPlane = false; // YZ не является плоскостью полигона.
							}
						}
					}
				}
				if (x1 > 1.0e29) {
					x1 = tmp->pb.x;
				}
				else {
					if (fabs(tmp->pb.x - x1) < epsx) {
					}
					else {
						if (x2 > 1.0e29) {
							x2 = tmp->pb.x;
						}
						else {
							if (fabs(tmp->pb.x - x2) < epsx) {
							}
							else {
								bPlane = false; // YZ не является плоскостью полигона.
							}
						}
					}
				}
				if (x1 > 1.0e29) {
					x1 = tmp->pc.x;
				}
				else {
					if (fabs(tmp->pc.x - x1) < epsx) {
					}
					else {
						if (x2 > 1.0e29) {
							x2 = tmp->pc.x;
						}
						else {
							if (fabs(tmp->pc.x - x2) < epsx) {
							}
							else {
								bPlane = false; // YZ не является плоскостью полигона.
							}
						}
					}
				}


				tmp = tmp->next;
			}

			if (bPlane) {
				//iPlane_obj3 = YZ_PLANE;
				b1 = true;
			}

			//iPlane_obj3 = -1;
			//if (iPlane_obj3 == -1) 
			{
				tmp = root_CAD_STL;

				double y1 = 1.0e30, y2 = 1.0e30;
				bPlane = true;

				while (tmp != nullptr) {

					if (y1 > 1.0e29) {
						y1 = tmp->pa.y;
					}
					else {
						if (fabs(tmp->pa.y - y1) < epsy) {
						}
						else {
							if (y2 > 1.0e29) {
								y2 = tmp->pa.y;
							}
							else {
								if (fabs(tmp->pa.y - y2) < epsy) {
								}
								else {
									bPlane = false; // YZ не является плоскостью полигона.
								}
							}
						}
					}
					if (y1 > 1.0e29) {
						y1 = tmp->pb.y;
					}
					else {
						if (fabs(tmp->pb.y - y1) < epsy) {
						}
						else {
							if (y2 > 1.0e29) {
								y2 = tmp->pb.y;
							}
							else {
								if (fabs(tmp->pb.y - y2) < epsy) {
								}
								else {
									bPlane = false; // YZ не является плоскостью полигона.
								}
							}
						}
					}
					if (y1 > 1.0e29) {
						y1 = tmp->pc.y;
					}
					else {
						if (fabs(tmp->pc.y - y1) < epsy) {
						}
						else {
							if (y2 > 1.0e29) {
								y2 = tmp->pc.y;
							}
							else {
								if (fabs(tmp->pc.y - y2) < epsy) {
								}
								else {
									bPlane = false; // XZ не является плоскостью полигона.
								}
							}
						}
					}


					tmp = tmp->next;
				}

				if (bPlane) {
					//iPlane_obj3 = XZ_PLANE;
					b2 = true;
				}
			}

			//iPlane_obj3 = -1;
			//if (iPlane_obj3 == -1) 
			{
				tmp = root_CAD_STL;

				double z1 = 1.0e30, z2 = 1.0e30;
				bPlane = true;

				while (tmp != nullptr) {

					if (z1 > 1.0e29) {
						z1 = tmp->pa.z;
					}
					else {
						if (fabs(tmp->pa.z - z1) < epsz) {
						}
						else {
							if (z2 > 1.0e29) {
								z2 = tmp->pa.z;
							}
							else {
								if (fabs(tmp->pa.z - z2) < epsz) {
								}
								else {
									bPlane = false; // XY не является плоскостью полигона.
								}
							}
						}
					}
					if (z1 > 1.0e29) {
						z1 = tmp->pb.z;
					}
					else {
						if (fabs(tmp->pb.z - z1) < epsz) {
						}
						else {
							if (z2 > 1.0e29) {
								z2 = tmp->pb.z;
							}
							else {
								if (fabs(tmp->pb.z - z2) < epsz) {
								}
								else {
									bPlane = false; // XY не является плоскостью полигона.
								}
							}
						}
					}
					if (z1 > 1.0e29) {
						z1 = tmp->pc.z;
					}
					else {
						if (fabs(tmp->pc.z - z1) < epsz) {
						}
						else {
							if (z2 > 1.0e29) {
								z2 = tmp->pc.z;
							}
							else {
								if (fabs(tmp->pc.z - z2) < epsz) {
								}
								else {
									bPlane = false; // XY не является плоскостью полигона.
								}
							}
						}
					}


					tmp = tmp->next;
				}

				if (bPlane) {
					//iPlane_obj3 = XY_PLANE;
					b3 = true;
				}
			}

		}

		return (b1&&b2&&b3);
	}


	// Возвращает объём CAD_STL объекта.
	doublereal volume_CAD_STL() {

		doublereal vol = 0.0;

		TOCHKA_FLOAT pmin;
		TOCHKA_FLOAT pmax;

		minimaxCAD_STL(pmin, pmax);

		double epsx = 0.001*(pmax.x - pmin.x);
		double epsy = 0.001*(pmax.y - pmin.y);
		double epsz = 0.001*(pmax.z - pmin.z);

		if (!bvolcad) {

			// 14.11.2020
			// Работает только для CAD_STL объекта.

			


			integer iPlane_obj3=-1; // плоскость в которой лежит нижнее основание полигона.
			
			if (!bbigCADmodel) {

				gCAD_STL* tmp = root_CAD_STL;

				double x1=1.0e30, x2=1.0e30;
				bool bPlane = true;

				while (tmp != nullptr) {

					if (x1 > 1.0e29) {
						x1=tmp->pa.x;
					}
					else {
						if (fabs(tmp->pa.x - x1) < epsx) {
						}
						else {
							if (x2 > 1.0e29) {
								x2 = tmp->pa.x;
							}
							else {
								if (fabs(tmp->pa.x - x2) < epsx) {
								}
								else {
									bPlane = false; // YZ не является плоскостью полигона.
								}
							}
						}
					}
					if (x1 > 1.0e29) {
						x1 = tmp->pb.x;
					}
					else {
						if (fabs(tmp->pb.x - x1) < epsx) {
						}
						else {
							if (x2 > 1.0e29) {
								x2 = tmp->pb.x;
							}
							else {
								if (fabs(tmp->pb.x - x2) < epsx) {
								}
								else {
									bPlane = false; // YZ не является плоскостью полигона.
								}
							}
						}
					}
					if (x1 > 1.0e29) {
						x1 = tmp->pc.x;
					}
					else {
						if (fabs(tmp->pc.x - x1) < epsx) {
						}
						else {
							if (x2 > 1.0e29) {
								x2 = tmp->pc.x;
							}
							else {
								if (fabs(tmp->pc.x - x2) < epsx) {
								}
								else {
									bPlane = false; // YZ не является плоскостью полигона.
								}
							}
						}
					}


					tmp = tmp->next;
				}

				if (bPlane) {
					iPlane_obj3 = YZ_PLANE;
				}

				if (iPlane_obj3 == -1) {
					tmp = root_CAD_STL;

					double y1 = 1.0e30, y2 = 1.0e30;
					bPlane = true;

					while (tmp != nullptr) {

						if (y1 > 1.0e29) {
							y1 = tmp->pa.y;
						}
						else {
							if (fabs(tmp->pa.y - y1) < epsy) {
							}
							else {
								if (y2 > 1.0e29) {
									y2 = tmp->pa.y;
								}
								else {
									if (fabs(tmp->pa.y - y2) < epsy) {
									}
									else {
										bPlane = false; // YZ не является плоскостью полигона.
									}
								}
							}
						}
						if (y1 > 1.0e29) {
							y1 = tmp->pb.y;
						}
						else {
							if (fabs(tmp->pb.y - y1) < epsy) {
							}
							else {
								if (y2 > 1.0e29) {
									y2 = tmp->pb.y;
								}
								else {
									if (fabs(tmp->pb.y - y2) < epsy) {
									}
									else {
										bPlane = false; // YZ не является плоскостью полигона.
									}
								}
							}
						}
						if (y1 > 1.0e29) {
							y1 = tmp->pc.y;
						}
						else {
							if (fabs(tmp->pc.y - y1) < epsy) {
							}
							else {
								if (y2 > 1.0e29) {
									y2 = tmp->pc.y;
								}
								else {
									if (fabs(tmp->pc.y - y2) < epsy) {
									}
									else {
										bPlane = false; // XZ не является плоскостью полигона.
									}
								}
							}
						}


						tmp = tmp->next;
					}

					if (bPlane) {
						iPlane_obj3 = XZ_PLANE;
					}
				}

				if (iPlane_obj3 == -1) {
					tmp = root_CAD_STL;

					double z1 = 1.0e30, z2 = 1.0e30;
					bPlane = true;

					while (tmp != nullptr) {

						if (z1 > 1.0e29) {
							z1 = tmp->pa.z;
						}
						else {
							if (fabs(tmp->pa.z - z1) < epsz) {
							}
							else {
								if (z2 > 1.0e29) {
									z2 = tmp->pa.z;
								}
								else {
									if (fabs(tmp->pa.z - z2) < epsz) {
									}
									else {
										bPlane = false; // XY не является плоскостью полигона.
									}
								}
							}
						}
						if (z1 > 1.0e29) {
							z1 = tmp->pb.z;
						}
						else {
							if (fabs(tmp->pb.z - z1) < epsz) {
							}
							else {
								if (z2 > 1.0e29) {
									z2 = tmp->pb.z;
								}
								else {
									if (fabs(tmp->pb.z - z2) < epsz) {
									}
									else {
										bPlane = false; // XY не является плоскостью полигона.
									}
								}
							}
						}
						if (z1 > 1.0e29) {
							z1 = tmp->pc.z;
						}
						else {
							if (fabs(tmp->pc.z - z1) < epsz) {
							}
							else {
								if (z2 > 1.0e29) {
									z2 = tmp->pc.z;
								}
								else {
									if (fabs(tmp->pc.y - z2) < epsz) {
									}
									else {
										bPlane = false; // XY не является плоскостью полигона.
									}
								}
							}
						}


						tmp = tmp->next;
					}

					if (bPlane) {
						iPlane_obj3 = XY_PLANE;
					}
				}

			}

			

			float resolution_x = 100.0;
			float resolution_y = 100.0;
			float resolution_z = 100.0;

			switch (iPlane_obj3) {
			case XY_PLANE:
				resolution_x = 100.0;
				resolution_y = 100.0;
				resolution_z = 1.0;
				break;
			case XZ_PLANE:
				resolution_x = 100.0;
				resolution_y = 1.0;
				resolution_z = 100.0;
				break;
			case YZ_PLANE:
				resolution_x = 1.0;
				resolution_y = 100.0;
				resolution_z = 100.0;
				break;
			}

			float hx, hy, hz;

			hx = fabs(pmax.x - pmin.x) / resolution_x;
			hy = fabs(pmax.y - pmin.y) / resolution_y;
			hz = fabs(pmax.z - pmin.z) / resolution_z;

			// Наблюдается сильное замедление работы программы.
			//int i_my_num_core_parallelesation = 1;
#ifdef _OPENMP 
			//i_my_num_core_parallelesation = omp_get_max_threads();
			//int inum_cores = number_cores();
			//omp_set_num_threads(inum_cores); // установка числа потоков
#endif

			integer inx7 = static_cast<integer>(resolution_x);
			integer iny7 = static_cast<integer>(resolution_y);
			integer inz7 = static_cast<integer>(resolution_z);

#pragma omp parallel for reduction(+: vol)
			for (integer i_73 = 0; i_73 <inx7; i_73++) {
				for (integer j_73 = 0; j_73 < iny7; j_73++) {
					for (integer k_73 = 0; k_73 < inz7; k_73++) {

						TOCHKA p;

						p.x = pmin.x + i_73 * hx + 0.5 * hx;
						p.y = pmin.y + j_73 * hy + 0.5 * hy;
						p.z = pmin.z + k_73 * hz + 0.5 * hz;

						int ib, k;
						ib = -1;
						if (in_CAD_STL_check(p, k, ib)) {
							vol += hx * hy * hz;
						}
					}
				}
			}

			volcadstl = vol;
			bvolcad = true;

#ifdef _OPENMP 
			//omp_set_num_threads(i_my_num_core_parallelesation);
#endif

		}
		else {
			vol = volcadstl;
		}

		// Возвращает объём CAD_STL объекта.
		return vol;
	}


	// Возвращает минимальное ребро треугольника.
	doublereal min_size_edge() {

		doublereal dist = 1.0e30;

		if (!bgetminedge) {

			

			if (itypegeom == CAD_STL) {

				gCAD_STL* tmp = root_CAD_STL;

				while (tmp != nullptr) {

					doublereal dist_loc = ((tmp->pa.x - tmp->pb.x) * (tmp->pa.x - tmp->pb.x) + (tmp->pa.y - tmp->pb.y) * (tmp->pa.y - tmp->pb.y) + (tmp->pa.z - tmp->pb.z) * (tmp->pa.z - tmp->pb.z));

					if (dist > dist_loc) {
						dist = dist_loc;
					}

					dist_loc = ((tmp->pa.x - tmp->pc.x) * (tmp->pa.x - tmp->pc.x) + (tmp->pa.y - tmp->pc.y) * (tmp->pa.y - tmp->pc.y) + (tmp->pa.z - tmp->pc.z) * (tmp->pa.z - tmp->pc.z));

					if (dist > dist_loc) {
						dist = dist_loc;
					}

					dist_loc = ((tmp->pb.x - tmp->pc.x) * (tmp->pb.x - tmp->pc.x) + (tmp->pb.y - tmp->pc.y) * (tmp->pb.y - tmp->pc.y) + (tmp->pb.z - tmp->pc.z) * (tmp->pb.z - tmp->pc.z));

					if (dist > dist_loc) {
						dist = dist_loc;
					}


					tmp = tmp->next;

				}

			}

			mem_get_min_edge = dist;
			bgetminedge = true;

		}
		else {
			dist = mem_get_min_edge;
		}

		return sqrt(dist);

	}


	// Преобразование призмы в CAD STL.
	// Шесть граней и 12 треугольников.
	// Данное преобразование нужно для 
	// тестирования методов обрабатывающих CAD_STL.
	void PRISM2CAD_STL_forTEST() {

		// 14.11.2020



		if ((itypegeom == PRISM)||
			((itypegeom == POLYGON) && (nsizei == 4))) 
		{
			PRISM2CAD_STL();

			inumber_triangles_for_CAD_STL_model = 12;

			itypegeom = CAD_STL;

			TOCHKA_FLOAT pmin;
			TOCHKA_FLOAT pmax;

			minimaxCAD_STL(pmin, pmax);

			// Окаймляющая прямоугольная призма для
			// CAD STL объекта.

			xS = pmin.x;
			yS = pmin.y;
			zS = pmin.z;

			xE = pmax.x;
			yE = pmax.y;
			zE = pmax.z;

		}
		if (itypegeom == CYLINDER) {
			PRISM2CAD_STL();

			inumber_triangles_for_CAD_STL_model = 200;

			itypegeom = CAD_STL;

			TOCHKA_FLOAT pmin;
			TOCHKA_FLOAT pmax;

			minimaxCAD_STL(pmin, pmax);

			// Окаймляющая прямоугольная призма для
			// CAD STL объекта.

			xS = pmin.x;
			yS = pmin.y;
			zS = pmin.z;

			xE = pmax.x;
			yE = pmax.y;
			zE = pmax.z;

			//std::wcout << name << std::endl;
			//std::cout << "xS=" << xS << " xE=" << xE << " yS=" << yS << " yE=" << yE << " zS=" << zS << " zE=" << zE << std::endl;
			//std::wcout << "center x="<<(0.5*(xS+xE))<< "  y = "<<(0.5*(yS+yE))<< "   z = "<<(0.5*(zS+zE))<<std::endl;
			//system("pause");
		}
	}

	// Считывает бинарный STL файл с CAD геометрией.
	// 15.11.2020; 17.11.2020 для большой .stl модели.
	void ReadSTL_binary(int lb, TGEOM &gCabinet) {

		std::cout << "please, input *.stl file name...\n";

		std::string FName;

		std::cin >> FName;
		std::cout << "\nFile name '" << FName << "' succsefull reading. \n";


		// Считывание.

		int x = 0;//Переменные
		short int z = 0;
		float y = 0;
		char title[80];
		int inum = 0;

		/*НАЧАЛО РАБОТЫ С ФАЙЛОМ*/
		std::ifstream in(FName, std::ios::binary);

		in.read((char*)&title, sizeof(title));
		in.read((char*)&x, sizeof(x)); //перенос байтов из файла в "х"
		inum = x;
		std::cout << "number triangles " << inum << std::endl;
		system("PAUSE");		

		// освобождение оперативной памяти.
		clear_CAD_STL();

		inumber_triangles_for_CAD_STL_model = inum;

		root_CAD_STL = new gCAD_STL;
		gCAD_STL* tmp = root_CAD_STL;

		for (int i = 0; i < inum; ++i) {
			in.read((char*)&y, sizeof(y));  //перенос байтов из файла в "y"
			tmp->n.x = y;
			in.read((char*)&y, sizeof(y));
			tmp->n.y = y;
			in.read((char*)&y, sizeof(y));
			tmp->n.z = y;

			in.read((char*)&y, sizeof(y));
			tmp->pa.x = y;
			in.read((char*)&y, sizeof(y));
			tmp->pa.y = y;
			in.read((char*)&y, sizeof(y));
			tmp->pa.z = y;

			in.read((char*)&y, sizeof(y));
			tmp->pb.x = y;
			in.read((char*)&y, sizeof(y));
			tmp->pb.y = y;
			in.read((char*)&y, sizeof(y));
			tmp->pb.z = y;

			in.read((char*)&y, sizeof(y));
			tmp->pc.x = y;
			in.read((char*)&y, sizeof(y));
			tmp->pc.y = y;
			in.read((char*)&y, sizeof(y));
			tmp->pc.z = y;

			in.read((char*)&z, sizeof(z)); // это может быть информация о цвете.

			if (i < 2) {
				std::cout << "Test print to the console of the first two read triangles:" << std::endl;
				std::cout << "nx=" << tmp->n.x << " ny=" << tmp->n.y << " nz=" << tmp->n.z << std::endl;
				std::cout << "pa.x=" << tmp->pa.x << " pa.y=" << tmp->pa.y << " pa.z=" << tmp->pa.z << std::endl;
				std::cout << "pb.x=" << tmp->pb.x << " pb.y=" << tmp->pb.y << " pb.z=" << tmp->pb.z << std::endl;
				std::cout << "pc.x=" << tmp->pc.x << " pc.y=" << tmp->pc.y << " pc.z=" << tmp->pc.z << std::endl;
				//system("PAUSE");
			}

			if (i < inum - 1) {
				tmp->next = new gCAD_STL;
			}
			tmp = tmp->next;

		}
		        
		in.close();
		/*КОНЕЦ РАБОТЫ С ФАЙЛОМ*/


		// Вывод статистики.

		

		int number_faces=0;

		tmp = root_CAD_STL;

		while (tmp != nullptr) {

			number_faces++;

			tmp = tmp->next;
		}

		TOCHKA_FLOAT pmin;
		TOCHKA_FLOAT pmax;

		std::wcout << "number triangle faces is equal=" << number_faces << std::endl;

		minimaxCAD_STL(pmin, pmax);

		if (number_faces > 1000) {

			bbigCADmodel = true;

			// oz 256=16*16.

			root_big_CAD_STL_model_z = new bigCAD_STL;
			bigCAD_STL * tmp1 = root_big_CAD_STL_model_z;

			
			const double dsize_pattern = 1.0 * size_pattern;

			root_big_CAD_STL_model_z_hash_table = new bigCAD_STL**[size_pattern];
			for (int i = 0; i < size_pattern; ++i) {
				root_big_CAD_STL_model_z_hash_table[i]= new bigCAD_STL*[size_pattern];
				for (int j = 0; j < size_pattern; ++j) {
					root_big_CAD_STL_model_z_hash_table[i][j] = nullptr;
				}
			}

			// Всего 256 сегментов.
			for (int i = 0; i < size_pattern; ++i) {
				for (int j = 0; j < size_pattern; ++j) {
					tmp1->min1 = pmin.x + i * (pmax.x - pmin.x) / dsize_pattern;
					tmp1->max1 = pmin.x + (i + 1) * (pmax.x - pmin.x) / dsize_pattern;
					tmp1->min2 = pmin.y + j * (pmax.y - pmin.y) / dsize_pattern;
					tmp1->max2 = pmin.y + (j + 1) * (pmax.y - pmin.y) / dsize_pattern;

					tmp1->root_CAD_STL = nullptr;
					tmp1->id_endl = i * size_pattern + j;

					if (!((i == size_pattern-1) && (j == size_pattern-1))) {
						tmp1->next = new bigCAD_STL;
					}
					else {
						tmp1->next = nullptr;
					}

					//root_big_CAD_STL_model_z_hash_table[i][j] = tmp1;
					tmp1 = tmp1->next;
				}
			}

			// Oz распределение по ректанглам.

			gCAD_STL** endl_link = new gCAD_STL * [size_pattern* size_pattern];

			//bool badd = false;

			tmp = root_CAD_STL;

			while (tmp != nullptr) {

				//double x1 = tmp->pa.x;
				//double y1 = tmp->pa.y;

				double tr_xmin = my_imin3(tmp->pa.x, tmp->pb.x, tmp->pc.x);
				double tr_ymin = my_imin3(tmp->pa.y, tmp->pb.y, tmp->pc.y);

				double tr_xmax = my_imax3(tmp->pa.x, tmp->pb.x, tmp->pc.x);
				double tr_ymax = my_imax3(tmp->pa.y, tmp->pb.y, tmp->pc.y);
				
				tmp1 = root_big_CAD_STL_model_z;

				while (tmp1 != nullptr) {

					//if ((x1>= tmp1->min1)&&(x1<= tmp1->max1)&& (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) 
					if ((tr_xmin<= tmp1->max1)&&(tr_xmax>=tmp1->min1)&& (tr_ymin <= tmp1->max2) && (tr_ymax >= tmp1->min2))
					{

						//badd = true;

						// Добавление.
						if (tmp1->root_CAD_STL == nullptr) {
							tmp1->root_CAD_STL = new gCAD_STL;
							tmp1->root_CAD_STL->n = tmp->n;
							tmp1->root_CAD_STL->pa = tmp->pa;
							tmp1->root_CAD_STL->pb = tmp->pb;
							tmp1->root_CAD_STL->pc = tmp->pc;
							tmp1->root_CAD_STL->next=nullptr;
							endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
						}
						else {
							endl_link[tmp1->id_endl]->next= new gCAD_STL;
							endl_link[tmp1->id_endl]->next->n = tmp->n;
							endl_link[tmp1->id_endl]->next->pa = tmp->pa;
							endl_link[tmp1->id_endl]->next->pb = tmp->pb;
							endl_link[tmp1->id_endl]->next->pc = tmp->pc;
							endl_link[tmp1->id_endl]->next->next = nullptr;
							endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
						}
					}

					tmp1 = tmp1->next;
				}

				

				tmp = tmp->next;
			}
			/*
			if (!badd) {
				badd = false;

				tmp = root_CAD_STL;

				while (tmp != nullptr) {

					double x1 = tmp->pb.x;
					double y1 = tmp->pb.y;

					tmp1 = root_big_CAD_STL_model_z;

					while (tmp1 != nullptr) {

						if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) {

							badd = true;

							// Добавление.
							if (tmp1->root_CAD_STL == nullptr) {
								tmp1->root_CAD_STL = new gCAD_STL;
								tmp1->root_CAD_STL->n = tmp->n;
								tmp1->root_CAD_STL->pa = tmp->pa;
								tmp1->root_CAD_STL->pb = tmp->pb;
								tmp1->root_CAD_STL->pc = tmp->pc;
								tmp1->root_CAD_STL->next = nullptr;
								endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
							}
							else {
								endl_link[tmp1->id_endl]->next = new gCAD_STL;
								endl_link[tmp1->id_endl]->next->n = tmp->n;
								endl_link[tmp1->id_endl]->next->pa = tmp->pa;
								endl_link[tmp1->id_endl]->next->pb = tmp->pb;
								endl_link[tmp1->id_endl]->next->pc = tmp->pc;
								endl_link[tmp1->id_endl]->next->next = nullptr;
								endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
							}
						}

						tmp1 = tmp1->next;
					}



					tmp = tmp->next;
				}

			}

			if (!badd) {
				badd = false;

				tmp = root_CAD_STL;

				while (tmp != nullptr) {

					double x1 = tmp->pc.x;
					double y1 = tmp->pc.y;

					tmp1 = root_big_CAD_STL_model_z;

					while (tmp1 != nullptr) {

						if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) {

							badd = true;

							// Добавление.
							if (tmp1->root_CAD_STL == nullptr) {
								tmp1->root_CAD_STL = new gCAD_STL;
								tmp1->root_CAD_STL->n = tmp->n;
								tmp1->root_CAD_STL->pa = tmp->pa;
								tmp1->root_CAD_STL->pb = tmp->pb;
								tmp1->root_CAD_STL->pc = tmp->pc;
								tmp1->root_CAD_STL->next = nullptr;
								endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
							}
							else {
								endl_link[tmp1->id_endl]->next = new gCAD_STL;
								endl_link[tmp1->id_endl]->next->n = tmp->n;
								endl_link[tmp1->id_endl]->next->pa = tmp->pa;
								endl_link[tmp1->id_endl]->next->pb = tmp->pb;
								endl_link[tmp1->id_endl]->next->pc = tmp->pc;
								endl_link[tmp1->id_endl]->next->next = nullptr;
								endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
							}
						}

						tmp1 = tmp1->next;
					}



					tmp = tmp->next;
				}

			}
			*/
			for (int i = 0; i < size_pattern * size_pattern; ++i) {
				endl_link[i] = nullptr;
			}

			bigCAD_STL* tmp1t = root_big_CAD_STL_model_z;

			

			while (tmp1t != nullptr) {

				doublereal dkey_x= 0.5*(tmp1t->min1+ tmp1t->max1);
				doublereal dkey_y = 0.5 * (tmp1t->min2 + tmp1t->max2);

				integer key_x = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (dkey_x - pmin.x) / (pmax.x - pmin.x)), size_pattern - 1));
				integer key_y = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (dkey_y - pmin.y) / (pmax.y - pmin.y)), size_pattern - 1));

				root_big_CAD_STL_model_z_hash_table[key_x][key_y] = tmp1t;

				

				tmp1t = tmp1t->next;
			}



			// ox 256=16*16.

			root_big_CAD_STL_model_x = new bigCAD_STL;
			tmp1 = root_big_CAD_STL_model_x;

			root_big_CAD_STL_model_x_hash_table = new bigCAD_STL**[size_pattern];
			for (int i = 0; i < size_pattern; ++i) {
				root_big_CAD_STL_model_x_hash_table[i] = new bigCAD_STL*[size_pattern];
				for (int j = 0; j < size_pattern; ++j) {
					root_big_CAD_STL_model_x_hash_table[i][j] = nullptr;
				}
			}

			// Всего 256 сегментов.
			for (int i = 0; i < size_pattern; ++i) {
				for (int j = 0; j < size_pattern; ++j) {
					tmp1->min1 = pmin.y + i * (pmax.y - pmin.y) / dsize_pattern;
					tmp1->max1 = pmin.y + (i + 1) * (pmax.y - pmin.y) / dsize_pattern;
					tmp1->min2 = pmin.z + j * (pmax.z - pmin.z) / dsize_pattern;
					tmp1->max2 = pmin.z + (j + 1) * (pmax.z - pmin.z) / dsize_pattern;

					tmp1->root_CAD_STL = nullptr;
					tmp1->id_endl = i * size_pattern + j;

					if (!((i == size_pattern - 1) && (j == size_pattern - 1))) {
						tmp1->next = new bigCAD_STL;
					}
					else {
						tmp1->next = nullptr;
					}
					//root_big_CAD_STL_model_x_hash_table[i][j] = tmp1;
					tmp1 = tmp1->next;
				}
			}

			// Ox распределение по ректанглам.			

			//badd = false;

			tmp = root_CAD_STL;

			while (tmp != nullptr) {

				//double x1 = tmp->pa.y;
				//double y1 = tmp->pa.z;

				double tr_xmin = my_imin3(tmp->pa.y, tmp->pb.y, tmp->pc.y);
				double tr_ymin = my_imin3(tmp->pa.z, tmp->pb.z, tmp->pc.z);

				double tr_xmax = my_imax3(tmp->pa.y, tmp->pb.y, tmp->pc.y);
				double tr_ymax = my_imax3(tmp->pa.z, tmp->pb.z, tmp->pc.z);

				
				tmp1 = root_big_CAD_STL_model_x;

				while (tmp1 != nullptr) {

					//if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) 
						if ((tr_xmin <= tmp1->max1) && (tr_xmax >= tmp1->min1) && (tr_ymin <= tmp1->max2) && (tr_ymax >= tmp1->min2))
					{
						// Добавление.

						//badd = true;

						if (tmp1->root_CAD_STL == nullptr) {
							tmp1->root_CAD_STL = new gCAD_STL;
							tmp1->root_CAD_STL->n = tmp->n;
							tmp1->root_CAD_STL->pa = tmp->pa;
							tmp1->root_CAD_STL->pb = tmp->pb;
							tmp1->root_CAD_STL->pc = tmp->pc;
							tmp1->root_CAD_STL->next = nullptr;
							endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
						}
						else {
							endl_link[tmp1->id_endl]->next = new gCAD_STL;
							endl_link[tmp1->id_endl]->next->n = tmp->n;
							endl_link[tmp1->id_endl]->next->pa = tmp->pa;
							endl_link[tmp1->id_endl]->next->pb = tmp->pb;
							endl_link[tmp1->id_endl]->next->pc = tmp->pc;
							endl_link[tmp1->id_endl]->next->next = nullptr;
							endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
						}
					}

					tmp1 = tmp1->next;
				}



				tmp = tmp->next;
			}

			/*
			if (!badd) {
				badd = false;

				tmp = root_CAD_STL;

				while (tmp != nullptr) {

					double x1 = tmp->pb.y;
					double y1 = tmp->pb.z;

					tmp1 = root_big_CAD_STL_model_x;

					while (tmp1 != nullptr) {

						if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) {
							// Добавление.

							badd = true;

							if (tmp1->root_CAD_STL == nullptr) {
								tmp1->root_CAD_STL = new gCAD_STL;
								tmp1->root_CAD_STL->n = tmp->n;
								tmp1->root_CAD_STL->pa = tmp->pa;
								tmp1->root_CAD_STL->pb = tmp->pb;
								tmp1->root_CAD_STL->pc = tmp->pc;
								tmp1->root_CAD_STL->next = nullptr;
								endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
							}
							else {
								endl_link[tmp1->id_endl]->next = new gCAD_STL;
								endl_link[tmp1->id_endl]->next->n = tmp->n;
								endl_link[tmp1->id_endl]->next->pa = tmp->pa;
								endl_link[tmp1->id_endl]->next->pb = tmp->pb;
								endl_link[tmp1->id_endl]->next->pc = tmp->pc;
								endl_link[tmp1->id_endl]->next->next = nullptr;
								endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
							}
						}

						tmp1 = tmp1->next;
					}



					tmp = tmp->next;
				}
			}

			if (!badd) {
				badd = false;

				tmp = root_CAD_STL;

				while (tmp != nullptr) {

					double x1 = tmp->pc.y;
					double y1 = tmp->pc.z;

					tmp1 = root_big_CAD_STL_model_x;

					while (tmp1 != nullptr) {

						if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) {
							// Добавление.

							badd = true;

							if (tmp1->root_CAD_STL == nullptr) {
								tmp1->root_CAD_STL = new gCAD_STL;
								tmp1->root_CAD_STL->n = tmp->n;
								tmp1->root_CAD_STL->pa = tmp->pa;
								tmp1->root_CAD_STL->pb = tmp->pb;
								tmp1->root_CAD_STL->pc = tmp->pc;
								tmp1->root_CAD_STL->next = nullptr;
								endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
							}
							else {
								endl_link[tmp1->id_endl]->next = new gCAD_STL;
								endl_link[tmp1->id_endl]->next->n = tmp->n;
								endl_link[tmp1->id_endl]->next->pa = tmp->pa;
								endl_link[tmp1->id_endl]->next->pb = tmp->pb;
								endl_link[tmp1->id_endl]->next->pc = tmp->pc;
								endl_link[tmp1->id_endl]->next->next = nullptr;
								endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
							}
						}

						tmp1 = tmp1->next;
					}



					tmp = tmp->next;
				}
			}
			*/
			for (int i = 0; i < size_pattern * size_pattern; ++i) {
				endl_link[i] = nullptr;
			}

			 tmp1t = root_big_CAD_STL_model_x;



			while (tmp1t != nullptr) {

				doublereal dkey_y = 0.5 * (tmp1t->min1 + tmp1t->max1);
				doublereal dkey_z = 0.5 * (tmp1t->min2 + tmp1t->max2);

				integer key_y = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (dkey_y - pmin.y) / (pmax.y - pmin.y)), size_pattern - 1));
				integer key_z = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (dkey_z - pmin.z) / (pmax.z - pmin.z)), size_pattern - 1));

				root_big_CAD_STL_model_x_hash_table[key_y][key_z] = tmp1t;

				

				tmp1t = tmp1t->next;
			}


			// oy 256=16*16.

			root_big_CAD_STL_model_y = new bigCAD_STL;
			tmp1 = root_big_CAD_STL_model_y;

			root_big_CAD_STL_model_y_hash_table = new bigCAD_STL**[size_pattern];
			for (int i = 0; i < size_pattern; ++i) {
				root_big_CAD_STL_model_y_hash_table[i] = new bigCAD_STL*[size_pattern];
				for (int j = 0; j < size_pattern; ++j) {
					root_big_CAD_STL_model_y_hash_table[i][j] = nullptr;
				}
			}

			// Всего 256 сегментов.
			for (int i = 0; i < size_pattern; ++i) {
				for (int j = 0; j < size_pattern; ++j) {
					tmp1->min1 = pmin.x + i * (pmax.x - pmin.x) / dsize_pattern;
					tmp1->max1 = pmin.x + (i + 1) * (pmax.x - pmin.x) / dsize_pattern;
					tmp1->min2 = pmin.z + j * (pmax.z - pmin.z) / dsize_pattern;
					tmp1->max2 = pmin.z + (j + 1) * (pmax.z - pmin.z) / dsize_pattern;

					tmp1->root_CAD_STL = nullptr;
					tmp1->id_endl = i * size_pattern + j;

					if (!((i == size_pattern - 1) && (j == size_pattern - 1))) {
						tmp1->next = new bigCAD_STL;
					}
					else {
						tmp1->next = nullptr;
					}
					//root_big_CAD_STL_model_y_hash_table[i][j] = tmp1;
					tmp1 = tmp1->next;
				}
			}

			// Oz распределение по ректанглам.			

			//badd = false;

			tmp = root_CAD_STL;

			while (tmp != nullptr) {

				//double x1 = tmp->pa.x;
				//double y1 = tmp->pa.z;

				double tr_xmin = my_imin3(tmp->pa.x, tmp->pb.x, tmp->pc.x);
				double tr_ymin = my_imin3(tmp->pa.z, tmp->pb.z, tmp->pc.z);

				double tr_xmax = my_imax3(tmp->pa.x, tmp->pb.x, tmp->pc.x);
				double tr_ymax = my_imax3(tmp->pa.z, tmp->pb.z, tmp->pc.z);
				

				tmp1 = root_big_CAD_STL_model_y;

				while (tmp1 != nullptr) {

					//if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2))
						if ((tr_xmin <= tmp1->max1) && (tr_xmax >= tmp1->min1) && (tr_ymin <= tmp1->max2) && (tr_ymax >= tmp1->min2))
					{
						// Добавление.

						//badd = true;

						if (tmp1->root_CAD_STL == nullptr) {
							tmp1->root_CAD_STL = new gCAD_STL;
							tmp1->root_CAD_STL->n = tmp->n;
							tmp1->root_CAD_STL->pa = tmp->pa;
							tmp1->root_CAD_STL->pb = tmp->pb;
							tmp1->root_CAD_STL->pc = tmp->pc;
							tmp1->root_CAD_STL->next = nullptr;
							endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
						}
						else {
							endl_link[tmp1->id_endl]->next = new gCAD_STL;
							endl_link[tmp1->id_endl]->next->n = tmp->n;
							endl_link[tmp1->id_endl]->next->pa = tmp->pa;
							endl_link[tmp1->id_endl]->next->pb = tmp->pb;
							endl_link[tmp1->id_endl]->next->pc = tmp->pc;
							endl_link[tmp1->id_endl]->next->next = nullptr;
							endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
						}
					}

					tmp1 = tmp1->next;
				}



				tmp = tmp->next;
			}
			/*
			if (!badd) {
				badd = false;

				tmp = root_CAD_STL;

				while (tmp != nullptr) {

					double x1 = tmp->pb.x;
					double y1 = tmp->pb.z;

					tmp1 = root_big_CAD_STL_model_y;

					while (tmp1 != nullptr) {

						if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) {
							// Добавление.

							badd = true;

							if (tmp1->root_CAD_STL == nullptr) {
								tmp1->root_CAD_STL = new gCAD_STL;
								tmp1->root_CAD_STL->n = tmp->n;
								tmp1->root_CAD_STL->pa = tmp->pa;
								tmp1->root_CAD_STL->pb = tmp->pb;
								tmp1->root_CAD_STL->pc = tmp->pc;
								tmp1->root_CAD_STL->next = nullptr;
								endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
							}
							else {
								endl_link[tmp1->id_endl]->next = new gCAD_STL;
								endl_link[tmp1->id_endl]->next->n = tmp->n;
								endl_link[tmp1->id_endl]->next->pa = tmp->pa;
								endl_link[tmp1->id_endl]->next->pb = tmp->pb;
								endl_link[tmp1->id_endl]->next->pc = tmp->pc;
								endl_link[tmp1->id_endl]->next->next = nullptr;
								endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
							}
						}

						tmp1 = tmp1->next;
					}



					tmp = tmp->next;
				}
			}


			if (!badd) {
				badd = false;

				tmp = root_CAD_STL;

				while (tmp != nullptr) {

					double x1 = tmp->pc.x;
					double y1 = tmp->pc.z;

					tmp1 = root_big_CAD_STL_model_y;

					while (tmp1 != nullptr) {

						if ((x1 >= tmp1->min1) && (x1 <= tmp1->max1) && (y1 >= tmp1->min2) && (y1 <= tmp1->max2)) {
							// Добавление.

							badd = true;

							if (tmp1->root_CAD_STL == nullptr) {
								tmp1->root_CAD_STL = new gCAD_STL;
								tmp1->root_CAD_STL->n = tmp->n;
								tmp1->root_CAD_STL->pa = tmp->pa;
								tmp1->root_CAD_STL->pb = tmp->pb;
								tmp1->root_CAD_STL->pc = tmp->pc;
								tmp1->root_CAD_STL->next = nullptr;
								endl_link[tmp1->id_endl] = tmp1->root_CAD_STL;
							}
							else {
								endl_link[tmp1->id_endl]->next = new gCAD_STL;
								endl_link[tmp1->id_endl]->next->n = tmp->n;
								endl_link[tmp1->id_endl]->next->pa = tmp->pa;
								endl_link[tmp1->id_endl]->next->pb = tmp->pb;
								endl_link[tmp1->id_endl]->next->pc = tmp->pc;
								endl_link[tmp1->id_endl]->next->next = nullptr;
								endl_link[tmp1->id_endl] = endl_link[tmp1->id_endl]->next;
							}
						}

						tmp1 = tmp1->next;
					}



					tmp = tmp->next;
				}
			}
			*/
			for (int i = 0; i < size_pattern * size_pattern; ++i) {
				endl_link[i] = nullptr;
			}

			delete[] endl_link;


			 tmp1t = root_big_CAD_STL_model_y;



			while (tmp1t != nullptr) {

				doublereal dkey_x = 0.5 * (tmp1t->min1 + tmp1t->max1);
				doublereal dkey_z = 0.5 * (tmp1t->min2 + tmp1t->max2);

				integer key_x = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (dkey_x - pmin.x) / (pmax.x - pmin.x)), size_pattern - 1));
				integer key_z = my_imax(0, my_imin(static_cast<integer>((dsize_pattern) * (dkey_z - pmin.z) / (pmax.z - pmin.z)), size_pattern - 1));

				root_big_CAD_STL_model_y_hash_table[key_x][key_z] = tmp1t;



				tmp1t = tmp1t->next;
			}
		}

		std::wcout << "import .stl CAD obj dimensions:\n";
		std::wcout << "xS=" << pmin.x << "  xE=" << pmax.x << std::endl;
		std::wcout << "yS=" << pmin.y << "  yE=" << pmax.y << std::endl;
		std::wcout << "zS=" << pmin.z << "  zE=" << pmax.z << std::endl;

		std::wcout << "min_edge=" << min_size_edge() << std::endl;

		//system("PAUSE");

		// Окаймляющая призма.

		itypegeom = CAD_STL;
		xS = pmin.x;
		xE = pmax.x;
		yS = pmin.y;
		yE = pmax.y;
		zS = pmin.z;
		zE = pmax.z;

		/*
		if (lb == 1) {
			std::wcout << "please, enter cabinet xS=\n";
			std::cin >> gCabinet.xS;
			std::wcout << "please, enter cabinet xE=\n";
			std::cin >> gCabinet.xE;

			std::wcout << "please, enter cabinet yS=\n";
			std::cin >> gCabinet.yS;
			std::wcout << "please, enter cabinet yE=\n";
			std::cin >> gCabinet.yE;

			std::wcout << "please, enter cabinet zS=\n";
			std::cin >> gCabinet.zS;
			std::wcout << "please, enter cabinet zE=\n";
			std::cin >> gCabinet.zE;
		}
		*/

		std::wcout << "stl binary file succseffully reading...\n";

		system("PAUSE");

	}

	// Автоматический деструктор.
	~TGEOM() {
		// освобождение оперативной памяти.
		clear_CAD_STL();
	}

} GEOM;

// Мощность рассеиваемая в тепло зависит от
// температуры и смещения стока.
// Информация о зависимости хранится в таблично заданной 
// функция. Таблица взята из результатов расчёта в программе
// Г.З. Гарбера.
typedef struct TTEMP_DEP_POWER {

	integer intemp; // количество различных дискретных значений температуры.
	integer inoffset_drain; // количество различных дискретных значений смещения стока.

	// таблица: левый столбец температуры,
	// верхняя строка без первого значения - смещения стока.
	// Остальные значения в таблице соответствуют рассеиваемой мощности.
	doublereal* rtemp; // значения температуры
	doublereal* roffset_drain; // значения смещения стока
	doublereal** rpower_table; // таблица мощностей.

	char* sname; // имя текстового файла содержащего табличную функцию.

	TTEMP_DEP_POWER() {
		 sname = nullptr; // имя текстового файла содержащего табличную функцию.

		intemp = 0; // количество различных дискретных значений температуры.
		inoffset_drain = 0; // количество различных дискретных значений смещения стока.

		// таблица: левый столбец температуры,
		// верхняя строка без первого значения - смещения стока.
		// Остальные значения в таблице соответствуют рассеиваемой мощности.
		rtemp = nullptr; // значения температуры
		roffset_drain = nullptr; // значения смещения стока
		rpower_table = nullptr; // таблица мощностей.
	}
} TEMP_DEP_POWER;

// источник тепла
typedef struct TSOURCE {
	

	// мощность источника тепла
	// в случае если мощность зависит от максимальной температуры и смещения стока то
	// это значение автоматически пересчитывается перед сборкой матрицы, 
	// так что в коде сборки граничного условия ничего менять не нужно.
	doublereal power;
	
	
	doublereal roperation_offset_drain; // рабочее значение смещения стока.
	doublereal power_multiplyer; // на эту константу домножается реальная мощность, это удобно для корректирования
	// табличной мощности например когда нужно считать четверть транзистора.

	doublereal square; // площадь источника тепла
	integer iPlane; // плоскость в которой лежит источник тепла
	// 1 - XY, 2 - XZ, 3 - YZ.

	// В случае bgarber_depend  , значение power является коэффициентом на который домножается 
	// истинное табличное значение мощности.
	integer igarber_depend; // уникальный номер табличной зависимости. (нумерация начинается с нуля).
	
	GEOM g; // границы объекта.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id;

	// Мощность может зависеть от температуры и смещения стока.
	// Зависимость задаётся таблично на основе расчётов в программе Г.З. Гарбера.
	bool bgarber_depend; // задана ли зависимость рассеиваемой мощности от температуры и смещения стока.

	char name[80]; // имя источника тепла.


	TSOURCE() {
		name[0] = '\0';
		// мощность источника тепла
		// в случае если мощность зависит от максимальной температуры и смещения стока то
		// это значение автоматически пересчитывается перед сборкой матрицы, 
		// так что в коде сборки граничного условия ничего менять не нужно.
		power = 0.0;
		// Мощность может зависеть от температуры и смещения стока.
		// Зависимость задаётся таблично на основе расчётов в программе Г.З. Гарбера.
		 bgarber_depend = false; // задана ли зависимость рассеиваемой мощности от температуры и смещения стока.
		// В случае bgarber_depend  , значение power является коэффициентом на который домножается 
		// истинное табличное значение мощности.
		igarber_depend = -1; // уникальный номер табличной зависимости. (нумерация начинается с нуля).
		roperation_offset_drain = 28.0; // рабочее значение смещения стока.
		power_multiplyer = 1.0; // на эту константу домножается реальная мощность, это удобно для корректирования
		// табличной мощности например когда нужно считать четверть транзистора.

		square = 0.0; // площадь источника тепла
		iPlane = XZ_PLANE; // плоскость в которой лежит источник тепла
		// 1 - XY, 2 - XZ, 3 - YZ.
		//GEOM g; // границы объекта.

		// принадлежность объединению 
		// 0 - не принадлежит (принадлежит кабинету),
		// n > 0 принадлежит объединению с номером n.
		iunion_id = 0;
	}
} SOURCE;

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
enum class THERMAL_STRESS_BOUNDARY_CONDITION {
	FREE=0, X_FIXIT=1, Y_FIXIT=2, Z_FIXIT=3,
	XY_FIXIT=4, XZ_FIXIT=5, YZ_FIXIT=6, ALL_FIXIT=7,
	X_FORCE=8, Y_FORCE=9, Z_FORCE=10
};

// стенка (идеальный теплоотвод)
typedef struct TWALL {
	

	// 1 - Дирихле, 2 - однородное условие Неймана,
	// 3 - Ньютон-Рихман, 4 - Стефан-Больцман.
	WALL_BOUNDARY_CONDITION ifamily; // род краевого условия по температуре
	doublereal Tamb, hf; // значение температуры на идеальном теплоотводе и тепловой поток.
	doublereal emissivity, film_coefficient;
	doublereal ViewFactor; // Фактор видимости.

	doublereal Vx, Vy, Vz, P;
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
	THERMAL_STRESS_BOUNDARY_CONDITION ithermal_Stress_boundary_condition;
	doublereal xForce, yForce, zForce; // Три компоненты силы, приложенные к узлу в Ньютонах.
	integer iPlane; // плоскость в которой лежит стенка
	// 1 - XY, 2 - XZ, 3 - YZ.
	GEOM g; // границы объекта.
	// Фиксированная граница с нулевым смещением для
	// расчёта прочности конструкции.
	//bool bfixboundary; // true фиксированная, false свободная.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id;

	bool bsymmetry, bpressure, bopening;

	char name[80]; // имя стенки.

	TWALL() {
		name[0] = '\0';

		// 1 - Дирихле, 2 - однородное условие Неймана(теплоизоляция),
		// 3 - Ньютон-Рихман(конвекция),
		// 4 - Стефан-Больцман(излучение).
		ifamily = WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; // однородное условие Неймана по температуре (теплоизоляция).
		Tamb = 20.0, hf = 0.0; // значение температуры на идеальном теплоотводе и тепловой поток.
		emissivity = 0.8, film_coefficient = 3.0;
		ViewFactor = 1.0;
		bsymmetry = false, bpressure = false, bopening = false;
		Vx = 0.0, Vy = 0.0, Vz = 0.0, P = 0.0;
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
		ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::ALL_FIXIT;
		// Три компоненты силы, приложенные к узлу в Ньютонах.
		xForce = 0.0, yForce = 0.0, zForce = 0.0; 
		iPlane = XZ_PLANE; // плоскость в которой лежит стенка
		// 1 - XY, 2 - XZ, 3 - YZ.
		//GEOM g; // границы объекта.
		// Фиксированная граница с нулевым смещением для
		// расчёта прочности конструкции.
		//bool bfixboundary; // true фиксированная, false свободная.

		// принадлежность объединению 
		// 0 - не принадлежит (принадлежит кабинету),
		// n > 0 принадлежит объединению с номером n.
		iunion_id = 0;
	}
} WALL;





// Проверка типа геометрии при вводе.
// CHECK_TYPE_GEOM(din);
void CHECK_TYPE_GEOM(integer itypegeom) {
	switch (itypegeom) {
	case PRISM: break;
	case CYLINDER: break;
	case POLYGON: break;
	case CAD_STL: break;
	default: 
		printf("ERROR undefined type geom in inputlaplas.\n");
		system("PAUSE");
		exit(1);
		break;
	}
}// CHECK_TYPE_GEOM

// свойства материалов
typedef struct TPROP {
	// Постоянные значения заданные пользователем:
	// rho - плотность,
	// cp - удельная теплоёмкость при постоянном давлении,
	// lam - теплопроводность,
	// mu - динамическая вязкость,
	// beta_t - коэффициент линейного температурного расширения.
	// doublereal rho, cp, lam;
	float rho;
	// 16_11_2016
	// температура передаётся в градусах Цельсия.
	// температурно - зависимые теплопроводность и теплоёмкость.
	integer n_lam, n_cp, n_beta_t_solid, n_YoungModule, n_Poisson_ratio; // число точек для линейной интерполяции.
	// arr_lam=f(temp_lam);
	// arr_cp=f(temp_cp);
	// beta_t_solid для Механики.
	float *temp_lam, *temp_cp, *temp_beta_t_solid, *arr_lam, *arr_cp, *arr_beta_t_solid;
	// Ортотропность теплопроводности:
	// позволяет моделировать тепловые трубы и материалы с ортотропной теплопроводностью:
	// CuMoCu, SiC и текстолитовые платы.
	// теплопроводность в заданном координатном направлении домножается на этот множитель.
	float orthotropy_multiplyer_x, orthotropy_multiplyer_y, orthotropy_multiplyer_z;
	float mu,  beta_t;
	// Следующие два параметра относятся к внутренней библиотеке 
	// материалов:
	integer blibmat; // 0 - материал пользователя и используются 
	// постоянные параметры предложенные выше. 1 - библиотечный материал.
	integer ilibident; // идентификатор библиотечного материала.
	// идентификатор: 
	// 0 - текучая среда которой нет в библиотеке, 
	// 100 - твёрдое тело которого нет в библиотеке.
	// В случае если идентификатор принимает значения 0 и 100 то
	// используются постоянные параметры заданные пользователем, т.к.
	// в этом случае blibmat - обязательно равен 0.
	// 1-99 - библиотечная текучая среда: жидкости, газы и жидкие металлы.
	// 101-<infinity библиотечное твёрдое тело.
	bool bBussineskApproach; // нужно ли использовать приближение Обербека-Буссинеска

	// Также предусмотрено моделирование неньютоновских
	// вязкопластических сред подчиняющихся некоторым простейшим законам:
	/* номер закона название закона
	*  0            постоянное значение (ньютоновская жидкость) 
	*  1            Оствальд-де Вель (power-law fluid)
	*  2            Кессон
	*  3            Прандтль
	*  4            Carreau
	*  5            Пауэлл-Эйринг
	*  6            Уильямсон
	*/
	integer ilawmu; // номер закона

	float mumin, mumax; // ограничители динамической вязкости
	float Amu, Bmu, Cmu, degreennmu; // константы моделей для зависимости вязкости от напряжения сдвига.

	// 8.4.2017 
	// Параметры прочности материала.
	float *temp_Poisson_ratio, *arr_Poisson_ratio;
	float *temp_Young_Module, *arr_Young_Module;
	float orthotropy_multiplyer_Poisson_ratio_xy, orthotropy_multiplyer_Poisson_ratio_yz, orthotropy_multiplyer_Poisson_ratio_xz;
	float orthotropy_multiplyer_Poisson_ratio_yx, orthotropy_multiplyer_Poisson_ratio_zy, orthotropy_multiplyer_Poisson_ratio_zx;
	float orthotropy_multiplyer_x_Young_Module, orthotropy_multiplyer_y_Young_Module, orthotropy_multiplyer_z_Young_Module;
	// Модуль сдвига
	bool bActive_ShearModule;
	float ShearModule_xy, ShearModule_yz, ShearModule_xz;
	float orthotropy_multiplyer_x_beta_t_solid, orthotropy_multiplyer_y_beta_t_solid, orthotropy_multiplyer_z_beta_t_solid;
	// Коэффициенты Ламе не используются.
	TPROP() {
		// Постоянные значения заданные пользователем:
		// rho - плотность,
		// cp - теплоёмкость,
		// lam - теплопроводность,
		// mu - динамическая вязкость,
		// beta_t - коэффициент линейного температурного расширения.
		//doublereal rho, cp, lam;
		rho = -1.0f;
		// 16_11_2016
		// температура передаётся в градусах Цельсия.
		// температурно - зависимые теплопроводность и теплоёмкость.
		n_lam = -1, n_cp = -1, n_beta_t_solid=-1, n_YoungModule=-1, n_Poisson_ratio=-1; // число точек для линейной интерполяции.
		// arr_lam=f(temp_lam);
		// arr_cp=f(temp_cp);
		temp_lam=nullptr; temp_cp=nullptr; arr_lam=nullptr;
		arr_cp = nullptr; 
		temp_beta_t_solid = nullptr; arr_beta_t_solid = nullptr;// beta_t_solid для Механики.
		// Ортотропность теплопроводности:
		// позволяет моделировать тепловые трубы и материалы с ортотропной теплопроводностью:
		// CuMoCu, SiC и текстолитовые платы.
		// теплопроводность в заданном координатном направлении домножается на этот множитель.
		orthotropy_multiplyer_x=1.0, orthotropy_multiplyer_y = 1.0, orthotropy_multiplyer_z = 1.0;
		// Poisson_ratio в заданном координатном направлении домножается на этот множитель.
		orthotropy_multiplyer_Poisson_ratio_xy = 1.0; orthotropy_multiplyer_Poisson_ratio_yz = 1.0; orthotropy_multiplyer_Poisson_ratio_xz=1.0;
		orthotropy_multiplyer_Poisson_ratio_yx = 1.0; orthotropy_multiplyer_Poisson_ratio_zy = 1.0; orthotropy_multiplyer_Poisson_ratio_zx = 1.0;

		// Young_Module в заданном координатном направлении домножается на этот множитель.
		orthotropy_multiplyer_x_Young_Module = 1.0; orthotropy_multiplyer_y_Young_Module = 1.0; orthotropy_multiplyer_z_Young_Module=1.0;
		// Коэффициент линейного теплового расширения в заданном координатном направлении домножается на этот множитель.
		orthotropy_multiplyer_x_beta_t_solid = 1.0; orthotropy_multiplyer_y_beta_t_solid = 1.0; orthotropy_multiplyer_z_beta_t_solid=1.0;

		mu = -1.0f, beta_t = -1.0f;
		// Следующие два параметра относятся к внутренней библиотеке 
		// материалов:
		blibmat=-3; // 0 - материал пользователя и используются 
		// постоянные параметры предложенные выше. 1 - библиотечный материал.
		ilibident=-3; // идентификатор библиотечного материала.
		// идентификатор: 
		// 0 - текучая среда которой нет в библиотеке, 
		// 100 - твёрдое тело которого нет в библиотеке.
		// В случае если идентификатор принимает значения 0 и 100 то
		// используются постоянные параметры заданные пользователем, т.к.
		// в этом случае blibmat - обязательно равен 0.
		// 1-99 - библиотечная текучая среда: жидкости, газы и жидкие металлы.
		// 101-<infinity библиотечное твёрдое тело.
		bBussineskApproach=true; // нужно ли использовать приближение Обербека-Буссинеска

		// Также предусмотрено моделирование неньютоновских
		// вязкопластических сред подчиняющихся некоторым простейшим законам:
		/* номер закона название закона
		*  0            постоянное значение (ньютоновская жидкость) 
		*  1            Оствальд-де Вель (power-law fluid)
		*  2            Кессон
		*  3            Прандтль
		*  4            Carreau
		*  5            Пауэлл-Эйринг
		*  6            Уильямсон
		*/
		ilawmu=-3; // номер закона

		mumin = -1.0f, mumax = -1.0f; // ограничители динамической вязкости
		Amu = -1.0f, Bmu = -1.0f, Cmu = -1.0f, degreennmu = -1.0f; // константы моделей для зависимости вязкости от напряжения сдвига.

		// 8.4.2017 
		// Параметры прочности материала.
		// Коэффициенты Ламе не используются в данном коде 23.08.2020.
		// Коэффициент Пуассона.
		temp_Poisson_ratio = nullptr;
		arr_Poisson_ratio = nullptr;
		
		// Модуль Юнга.
		temp_Young_Module = nullptr;
		arr_Young_Module = nullptr;

		// Модуль сдвига.
		bActive_ShearModule=false; // не активен
		ShearModule_xy = 200.0; ShearModule_yz = 200.0; ShearModule_xz=200.0;
	}

} PROP;

const integer FULL_F_to_F_INTERPOLATION = 0; // Интерполяция по расширенному шаблону.
const integer AMG1R5_IN_HOUSE = 1;// Собственная реализация интерполяции amg1r5.

enum class AMGCL_ITERATOR_ALG {BiCGStab=0, FGMRes=1};
enum class MY_SORT_ALGORITHM {COUNTING_SORT=0, QUICK_SORT=1, HEAP_SORT=2, TIM_SORT=3};

// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
// default - 3.
enum class RS_COARSENING_KERNEL_DATA_STRUCTURE {AVL_TREE=0, SPLAY_TREE=1, BINARY_HEAP=2, TREAP=3, RED_BLACK_TREE=4, FIBONACCI_HEAP=5, VAN_EMDE_BOAS_TREE=6};

// Алгоритм построения грубой сетки.
// AMG Splitting (coarsening)
// Способ построения C-F разбиения: 0 - standart, 1 - RS2, 2 - ST classical standart (по умолчанию), 3 - RS2 ST.	
// RS2 улучшенная версия построения C-F разбиения содержащая второй проход.
// ST - на основе STRONG Transpose (именно применение ST опции рекомендовано в литературе).
// 8 - PMIS, 9 -HMIS.
// 10 - PMIS применённое к квадрату матрицы начиная с 4 уровня вложености.
// 11 - PMIS применённое к квадрату матрицы на всех уровнях (мой аналог aggressive coarsening). Должно использоваься в сочетании с дальнобойной интерполяцией.
enum class MY_AMG_SPLITTING_COARSENING_ALGORITHM {
	CLASSICAL_ALL_CONNECTION=0,
	RS2_ALL_CONNECTION=1,
	CLASSICAL_ST_ALL_CONNECTION=2, // ST - Strong Transpose.
	RS2_ST_ALL_CONNECTION=3,
	CLASSICAL_NEG_CONNECTION = 4,
	RS2_NEG_CONNECTION = 5,
	CLASSICAL_ST_NEG_CONNECTION = 6,
	RS2_ST_NEG_CONNECTION = 7,
	PMIS=8,
	HMIS=9,
	PMIS2=10,
	PMIS2_full=11
};


// Управление алгебраическим многосеточным методом из интерфейса.
typedef struct TMY_AMG_MANAGER {

	// Алгоритм сортировки используемый в многосеточном методе РУМБА.
	// 0 - Counting Sort, 1 - QUICKSORT, 3 - HEAPSORT, 4 - TimSort.
	MY_SORT_ALGORITHM imySortAlgorithm; // 0 - Counting Sort Default.

	// 0 - не печатать портрет матрицы, 
	// 1 - печатать портрет матрицы.
	bool bTemperatureMatrixPortrait;
	bool bSpeedMatrixPortrait;
	bool bPressureMatrixPortrait;
	bool bStressMatrixPortrait;
	bool bMatrixPortrait;
	
	// Temperature
	doublereal theta_Temperature;
	integer maximum_delete_levels_Temperature;
	integer nFinnest_Temperature, nu1_Temperature, nu2_Temperature;
	integer memory_size_Temperature; //13*size(matrix A)
	integer ilu2_smoother_Temperature; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
	// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
	// default - 3.
	RS_COARSENING_KERNEL_DATA_STRUCTURE iCFalgorithm_and_data_structure_Temperature;
	// Speed
	doublereal theta_Speed;
	integer maximum_delete_levels_Speed;
	integer nFinnest_Speed, nu1_Speed, nu2_Speed;
	integer memory_size_Speed;
	integer ilu2_smoother_Speed; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
	// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
	// default - 3.
	RS_COARSENING_KERNEL_DATA_STRUCTURE iCFalgorithm_and_data_structure_Speed;
	// Pressure
	doublereal theta_Pressure;
	integer maximum_delete_levels_Pressure;
	integer nFinnest_Pressure, nu1_Pressure, nu2_Pressure;
	integer memory_size_Pressure;
	integer ilu2_smoother_Pressure; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
	// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
	// default - 3.
	RS_COARSENING_KERNEL_DATA_STRUCTURE iCFalgorithm_and_data_structure_Pressure;
	// Stress
	doublereal theta_Stress;
	integer maximum_delete_levels_Stress;
	integer nFinnest_Stress, nu1_Stress, nu2_Stress;
	integer memory_size_Stress;
	integer ilu2_smoother_Stress; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
	// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
	// default - 3.
	RS_COARSENING_KERNEL_DATA_STRUCTURE iCFalgorithm_and_data_structure_Stress;
	// global
	bool bCFJacoby; // C/F упорядочивание при сглаживании.
	integer iRunge_Kutta_smoother; // 3 - третьего порядка, 5 - пятого порядка, любое другое число не используется. 
	integer iFinnest_ilu; // 0 не используется, 1 - ilu0. Только на самой подробной сетке.
	// Использование iluk разложения на глубоких уровнях вложенности для которых
	// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
	bool b_ilu_smoothers_in_nnz_n_LE_6;
	real_mix_precision theta; // strength threshold
	//integer maximum_levels; // максимальное количество уровней вложенности (уровни выше редуцируются).
	integer maximum_delete_levels; // Количество уровней отсекаемых снизу в области грубой сетки.
	integer nFinnest, nu1, nu2; // Количества сглаживаний.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
	// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
	// default - 2.
	RS_COARSENING_KERNEL_DATA_STRUCTURE iCFalgorithm_and_data_structure;
	integer memory_size; // В размерах матрицы А.
	// Для метода верхней релаксации в сглаживателе.
	
	integer number_interpolation_procedure; // идентификатор процедуры интерполяции.
	integer number_interpolation_procedure_Temperature;
	integer number_interpolation_procedure_Speed;
	integer number_interpolation_procedure_Pressure;
	integer number_interpolation_procedure_Stress;

	// 6 december 2016.
	// Подлежит удалению Refactoring.
	//integer itypemodifyinterpol=0; // номер модификации интерполяции. // Подлежит удалению Refactoring.
	//integer inumberadaptpass=0; // максимальное количество сканов-проходов с модификациями. // Подлежит удалению Refactoring.
	//integer baglomeration_with_consistency_scaling = 0;
	// Принудительное усиление диагонали
	// в случае обнаружения внедиагональных
	// positive connections.
	integer bdiagonal_dominant;


	// Параметры сглаживателя SOR Роуч.
	doublereal gold_const;
	doublereal gold_const_Temperature;
	doublereal gold_const_Speed;
	doublereal gold_const_Pressure;
	doublereal gold_const_Stress;

	// Пороги отсечек
	// Отсечка же называется theta threshold. 
		
	doublereal magic;
	doublereal F_to_F_Temperature, F_to_F_Speed, F_to_F_Pressure, F_to_F_Stress; // magic
	integer ilu2_smoother; // 0 - не использовать, 1 - использовать.
	// AMG Splitting (coarsening)
	// Способ построения C-F разбиения: 0 - standart, 1 - RS2, 2 - ST classical standart (по умолчанию), 3 - RS2 ST.	
	// RS2 улучшенная версия построения C-F разбиения содержащая второй проход.
	// ST - на основе STRONG Transpose.
	// 8 - PMIS, 9 -HMIS.
	MY_AMG_SPLITTING_COARSENING_ALGORITHM icoarseningTemp, icoarseningSpeed, icoarseningPressure, icoarseningStress;
    MY_AMG_SPLITTING_COARSENING_ALGORITHM icoarseningtype;
	// Stabilization BiCGStab.
	// 8.01.2017
	// предобусловленный алгебраическим многосеточным методом.
	// 0 - используется просто алгебраический многосеточный метод без какого-либо привлечения алгоритмов подпространства Крылова,
	// 1 - Используется алгоритм Хенка ван дер Ворста BiCGStab [1992], предобусловленный алгебраическим многосеточным методом.
	// 2 - Используется алгоритм Юсефа Саада и Мартина Г. Шульца FGMRes [1986], предобусловленный алгебраическим многосеточным методом.
	// 3 - Нелинейный многосеточный метод (обновление правой части на каждом V цикле). Для нелинейных граничных условий.
	integer istabilizationTemp, istabilizationSpeed, istabilizationPressure, istabilizationStress; // 0 - none
	integer istabilization; // 0 - none
	// ipatch - номер патча.
	integer ipatch_number;

	integer iprint_log, iprint_log_Temperature, iprint_log_Speed, iprint_log_Pressure, iprint_log_Stress;

	// truncation for interpolation.
	integer itruncation_interpolation;
	integer itruncation_interpolation_Temperature;
	integer itruncation_interpolation_Speed;
	integer itruncation_interpolation_Pressure;
	integer itruncation_interpolation_Stress;
	doublereal truncation_interpolation;
	doublereal truncation_interpolation_Temperature;
	doublereal truncation_interpolation_Speed;
	doublereal truncation_interpolation_Pressure;
	doublereal truncation_interpolation_Stress;

	// 01.07.2021
	// Задавать ли порог threshold автоматическим образом или всёже вручную.
	bool bthreshold_Pressure_auto, bthreshold_Temperature_auto, bthreshold_Speed_auto, bthreshold_Stress_auto;
	bool bthreshold_auto;
	// Применять ли CF переупорядочивание (порядок обхода) при итерировании, например в методе Якоби.
	bool bcf_reorder_Pressure, bcf_reorder_Temperature, bcf_reorder_Speed, bcf_reorder_Stress;
	bool bcf_reorder;

	// gmres smoother
	// Ю.Саад, Мартин Г. Шульц [1986].
	bool b_gmresTemp, b_gmresSpeed, b_gmresPressure, b_gmresStress;
	bool b_gmres;

	// Chebyshev smoother - сглаживатель Пафнутия Львовича Чебышева.
	// Пафнутий Львович Чебышев  4 (16) мая 1821 - 26 ноября (8 декабря) 1894.
	bool b_ChebyshevSmootherTemp, b_ChebyshevSmootherSpeed, b_ChebyshevSmootherPressure, b_ChebyshevSmootherStress;
	bool b_ChebyshevSmoother;

	// Степень полинома П.Л. Чебышева в сглаживателе Чебышева.
	integer Chebyshev_degree;

	// spai-0 smoother
	bool b_spai0Temp, b_spai0Speed, b_spai0Pressure, b_spai0Stress;
	bool b_spai0;

	// fgmres(m_restart)
	int m_restart;

	// lfil for BiCGStab+ILU2 and fgmres.
	integer lfil;

	// AMGCL parameters
	// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped_jacobi;
	// 4 - spai1; 5 - chebyshev; 6 - iluk, k=1; 7 - iluk, k=2;
	// 8 - iluk, k=4; 9 - iluk, k=6; 10 - iluk, k=8; 11 - iluk, k=10.
	integer amgcl_smoother; 
	integer amgcl_selector; // 0 - Ruge and Stuben (amg1r5 analog); 1 - smoother aggregation.
	AMGCL_ITERATOR_ALG amgcl_iterator; // 0 - BiCGStab; 1 - FGMRes;

	
	

	TMY_AMG_MANAGER() {
		// Алгоритм сортировки используемый в многосеточном методе РУМБА.
		// 0 - Counting Sort, 1 - QUICKSORT, 3 - HEAPSORT, 4 - TimSort.
		imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT; // 0 - Counting Sort Default.

		// 0 - не печатать портрет матрицы, 
		// 1 - печатать портрет матрицы.
		bTemperatureMatrixPortrait = false;
		bSpeedMatrixPortrait = false;
		bPressureMatrixPortrait = false;
		bStressMatrixPortrait = false;
		bMatrixPortrait = false;

	

		// Temperature
		theta_Temperature = 0.24;
		maximum_delete_levels_Temperature = 0;
		nFinnest_Temperature = 2; nu1_Temperature = 1; nu2_Temperature = 2;
		memory_size_Temperature = 13; //13*size(matrix A)
		ilu2_smoother_Temperature = 0; // 0 - не использовать, 1 - использовать.
		// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
		// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
		// default - 3.
		iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
		// Speed
		theta_Speed = 0.24;
		maximum_delete_levels_Speed = 0;
		nFinnest_Speed = 2; nu1_Speed = 1; nu2_Speed = 2;
		memory_size_Speed = 13;
		ilu2_smoother_Speed = 0; // 0 - не использовать, 1 - использовать.
		// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
		// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
		// default - 3.
		iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
		// Pressure
		theta_Pressure = 0.24;
		maximum_delete_levels_Pressure = 0;
		nFinnest_Pressure = 2; nu1_Pressure = 1; nu2_Pressure = 2;
		memory_size_Pressure = 15;
		ilu2_smoother_Pressure = 0; // 0 - не использовать, 1 - использовать.
		// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
		// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
		// default - 3.
		iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
		// Stress
		theta_Stress = 0.24;
		maximum_delete_levels_Stress = 0;
		nFinnest_Stress = 2; nu1_Stress = 1; nu2_Stress = 2;
		memory_size_Stress = 22;
		ilu2_smoother_Stress = 0; // 0 - не использовать, 1 - использовать.
		// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
		// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
		// default - 3.
		iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
		// global
		bCFJacoby = true; // C/F упорядочивание при сглаживании.
		iRunge_Kutta_smoother = 0; // 3 - третьего порядка, 5 - пятого порядка, любое другое число не используется. 
		iFinnest_ilu = 0; // 0 не используется, 1 - ilu0. Только на самой подробной сетке.
		// Использование iluk разложения на глубоких уровнях вложенности для которых
		// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
		b_ilu_smoothers_in_nnz_n_LE_6 = false;
		theta = static_cast<real_mix_precision>(0.24); // strength threshold
		//integer maximum_levels; // максимальное количество уровней вложенности (уровни выше редуцируются).
		maximum_delete_levels = 0; // Количество уровней отсекаемых снизу в области грубой сетки.
		nFinnest = 2; nu1 = 1; nu2 = 2; // Количества сглаживаний.
		// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
		// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
		// default - 2.
		iCFalgorithm_and_data_structure = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
		memory_size = 13; // В размерах матрицы А.
		// Для метода верхней релаксации в сглаживателе.
		//const integer FULL_F_to_F_INTERPOLATION = 0; // Интерполяция по расширенному шаблону.
		//const integer AMG1R5_IN_HOUSE = 1;// Собственная реализация интерполяции amg1r5.
		number_interpolation_procedure = AMG1R5_IN_HOUSE; // идентификатор процедуры интерполяции.
		number_interpolation_procedure_Temperature = AMG1R5_IN_HOUSE;
		number_interpolation_procedure_Speed = AMG1R5_IN_HOUSE;
		number_interpolation_procedure_Pressure = AMG1R5_IN_HOUSE;
		number_interpolation_procedure_Stress = AMG1R5_IN_HOUSE;

		// 6 december 2016.
		// Подлежит удалению Refactoring.
		//integer itypemodifyinterpol=0; // номер модификации интерполяции. // Подлежит удалению Refactoring.
		//integer inumberadaptpass=0; // максимальное количество сканов-проходов с модификациями. // Подлежит удалению Refactoring.
		//integer baglomeration_with_consistency_scaling = 0;
		// Принудительное усиление диагонали
		// в случае обнаружения внедиагональных
		// positive connections.
		bdiagonal_dominant = 1;


		// Параметры сглаживателя SOR Роуч.
		gold_const = 0.24;
		gold_const_Temperature = 0.24;
		gold_const_Speed = 0.24;
		gold_const_Pressure = 0.24;
		gold_const_Stress = 0.24;

		// Пороги отсечек
		// Отсечка же называется theta threshold. 
		
		magic = 0.4;
		F_to_F_Temperature = 0.4; F_to_F_Speed = 0.4;
		F_to_F_Pressure = 0.4; F_to_F_Stress = 0.4; // magic
		ilu2_smoother = 0; // 0 - не использовать, 1 - использовать.
		// AMG Splitting (coarsening)
		// Способ построения C-F разбиения: 0 - standart, 1 - RS2, 2 - ST classical standart (по умолчанию), 3 - RS2 ST.	
		// RS2 улучшенная версия построения C-F разбиения содержащая второй проход.
		// ST - на основе STRONG Transpose.
		// 8 - PMIS, 9 -HMIS.
		icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
		icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
		icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; 
		icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
		icoarseningtype = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
		// Stabilization BiCGStab.
		// 8.01.2017
		// предобусловленный алгебраическим многосеточным методом.
		// 0 - используется просто алгебраический многосеточный метод без какого-либо привлечения алгоритмов подпространства Крылова,
		// 1 - Используется алгоритм Хенка ван дер Ворста BiCGStab [1992], предобусловленный алгебраическим многосеточным методом.
		// 2 - Используется алгоритм Юсефа Саада и Мартина Г. Шульца FGMRes [1986], предобусловленный алгебраическим многосеточным методом.
		// 3 - Нелинейный многосеточный метод (обновление правой части на каждом V цикле). Для нелинейных граничных условий.
		istabilizationTemp = 0; istabilizationSpeed = 0; istabilizationPressure = 0;
		istabilizationStress = 0; // 0 - none
		istabilization = 0; // 0 - none
		// ipatch - номер патча.
		ipatch_number = 7;

		iprint_log = 1, iprint_log_Temperature = 1, iprint_log_Speed = 1, iprint_log_Pressure = 1, iprint_log_Stress = 1;

		// truncation for interpolation.
		itruncation_interpolation = 0;
		itruncation_interpolation_Temperature = 0;
		itruncation_interpolation_Speed = 0;
		itruncation_interpolation_Pressure = 0;
		itruncation_interpolation_Stress = 0;
		truncation_interpolation = 0.2;
		truncation_interpolation_Temperature = 0.2;
		truncation_interpolation_Speed = 0.2;
		truncation_interpolation_Pressure = 0.2;
		truncation_interpolation_Stress = 0.2;

		// gmres smoother
		// Ю.Саад, Мартин Г. Шульц [1986].
		b_gmresTemp = false; b_gmresSpeed = false;
		b_gmresPressure = false; b_gmresStress = false;
		b_gmres = false;

		// bicgstab smoother
		// Хенк ван дер Ворст [1992].
		b_ChebyshevSmootherTemp = false; b_ChebyshevSmootherSpeed = false;
		b_ChebyshevSmootherPressure = false; b_ChebyshevSmootherStress = false;
		b_ChebyshevSmoother = false;

		// Степень полинома П.Л. Чебышева в сглаживателе Чебышева.
		// Могут быть только значения 5; 8; 16; 25; 32; 64; 128.
		Chebyshev_degree = 5;

		// spai-0 smoother
		b_spai0Temp = false; b_spai0Speed = false; b_spai0Pressure = false;
		b_spai0Stress = false;
		b_spai0 = false;

		// 01.07.2021
	    // Задавать ли порог threshold автоматическим образом или всёже вручную.
		bthreshold_Pressure_auto = true; bthreshold_Temperature_auto = true; bthreshold_Speed_auto = true; bthreshold_Stress_auto = true;
		bthreshold_auto = true;
		// Применять ли CF переупорядочивание (порядок обхода) при итерировании, например в методе Якоби.
		bcf_reorder_Pressure = true; bcf_reorder_Temperature = true; bcf_reorder_Speed = true; bcf_reorder_Stress = true;
		bcf_reorder=true;

		// fgmres(m_restart)
		m_restart = 20;

		// lfil for BiCGStab+ILU2 and fgmres.
		lfil = 2;

		// AMGCL parameters
		// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped_jacobi;
		// 4 - spai1; 5 - chebyshev; 6 - iluk, k=1; 7 - iluk, k=2;
		// 8 - iluk, k=4; 9 - iluk, k=6; 10 - iluk, k=8; 11 - iluk, k=10.
		amgcl_smoother = 0; 
		amgcl_selector = 0; // 0 - Ruge and Stuben (amg1r5 analog); 1 - smoother aggregation.
		amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab; // 0 - BiCGStab; 1 - FGMRes;


		
		
	}

} MY_AMG_MANAGER;

// Параметры настройки алгебраического многосеточного метода Румба v.0.14
MY_AMG_MANAGER my_amg_manager;




void my_amg_manager_init() {
	
	// Алгоритм сортировки используемый в 
	// алгебраическом многосеточном методе РУМБА.
	// 0 - COUNTING SORT Х.Г. Сьювард 1954г.
	// 1 - QUICK SORT Чарльз Хоар 1962г.
	// 2 - HEAP SORT Вильямс 1964г.
	// 3 - Tim SORT Тим Петерсом 2002 г.
	my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT; // default value COUNTING SORT

	// Параметры собственного многосеточного метода о умолчанию.
	// Настройки решателя СЛАУ зависят от типа уравнения которое подаётся на вход:
	// симметричность <-> анизотропность, несимметричность, диффузионная задача <-> конвективная задача,
	// степень преобладания конвекции (число Рейнольдса).
	// my_amg_manager.maximum_levels = 20; // максимальное число уровней начиная с которого начинается усечение.
	my_amg_manager.maximum_delete_levels = 0; // Количество уровней отсекаемых в нижней части где грубая сетка.
	my_amg_manager.number_interpolation_procedure = 3; // номер интерполяционной процедуры.
	my_amg_manager.number_interpolation_procedure_Temperature = 3;
	my_amg_manager.number_interpolation_procedure_Speed = 3;
	my_amg_manager.number_interpolation_procedure_Pressure = 3;
	my_amg_manager.number_interpolation_procedure_Stress = 3;

	//my_amg_manager.baglomeration_with_consistency_scaling = 0;
	my_amg_manager.bdiagonal_dominant = 1;

	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap, 4 - Red Black Tree,
	// 5 - Fibonacci Heap, 6 - van Emde Boas Tree.
	my_amg_manager.iCFalgorithm_and_data_structure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP; // 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;// 3-Treap.	

	my_amg_manager.bTemperatureMatrixPortrait = false; // NO_PRINT
	my_amg_manager.bSpeedMatrixPortrait = false; // NO_PRINT
	my_amg_manager.bPressureMatrixPortrait = false; // NO_PRINT
	my_amg_manager.bStressMatrixPortrait = false; // NO_PRINT
	my_amg_manager.bMatrixPortrait = false; // NO_PRINT


	my_amg_manager.nFinnest = 2; // число итераций на подробной сетке.
	my_amg_manager.nu1 = 1; // число предсглаживаний.
	my_amg_manager.nu2 = 2; // число пост сглаживаний.	
	my_amg_manager.memory_size = 9; // количество оперативной памяти в размерностях матрицы А.
	my_amg_manager.gold_const = 0.2; // Параметр верхней релаксации в сглаживателе.
	my_amg_manager.gold_const_Temperature = 0.2;
	my_amg_manager.gold_const_Speed = 0.2;
	my_amg_manager.gold_const_Pressure = 0.2;
	my_amg_manager.gold_const_Stress = 0.2;

	my_amg_manager.bCFJacoby = true; // CF-Jacobi smoothers 12% сокращение числа V циклов. 5.06.2017
	// Runge-Kutt smoother: 3 - третьего порядка, 5 - пятого порядка, любое другое число не используется.
	my_amg_manager.iRunge_Kutta_smoother = 0;

	// ВНИМАНИЕ !!! устаревшие параметры, они больше не используются. 24.09.2019
	// 0 - не используется, 1 используется ILU разложение в качестве предобуславливателя.
	my_amg_manager.iFinnest_ilu = 0; 
	// Использование iluk разложения на глубоких уровнях вложенности для которых
	// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
	my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6 = false;
	// ВНИМАНИЕ !!! устаревшие параметры, они больше не используются. 24.09.2019


	my_amg_manager.theta = static_cast<real_mix_precision>(0.24);
	my_amg_manager.magic = 0.4;
	my_amg_manager.F_to_F_Temperature = 0.4;
	my_amg_manager.F_to_F_Speed = 0.4;
	my_amg_manager.F_to_F_Pressure = 0.4;
	my_amg_manager.F_to_F_Stress = 0.4;
	my_amg_manager.ilu2_smoother = 0; // 0 - не использовать, 1 - использовать.

	// Устаревшие переменные, более не используются и подлежат удалению. 11.05.2019
	//my_amg_manager.itypemodifyinterpol=0; // номер модификации интерполяции.
	//my_amg_manager.inumberadaptpass=0; // максимальное количество сканов-проходов с модификациями.

	my_amg_manager.theta_Temperature = 0.24;
	my_amg_manager.maximum_delete_levels_Temperature = 0;
	my_amg_manager.nFinnest_Temperature = 2;
	my_amg_manager.nu1_Temperature = 1;
	my_amg_manager.nu2_Temperature = 2;
	my_amg_manager.memory_size_Temperature = 9;
	my_amg_manager.ilu2_smoother_Temperature = 0; // 0 - не использовать, 1 - использовать.
	// Speed
	my_amg_manager.theta_Speed = 0.24;
	my_amg_manager.maximum_delete_levels_Speed = 0;
	my_amg_manager.nFinnest_Speed = 2;
	my_amg_manager.nu1_Speed = 1;
	my_amg_manager.nu2_Speed = 2;
	my_amg_manager.memory_size_Speed = 9;
	my_amg_manager.ilu2_smoother_Speed = 0; // 0 - не использовать, 1 - использовать.
	// Pressure
	my_amg_manager.theta_Pressure = 0.24;
	my_amg_manager.maximum_delete_levels_Pressure = 0;
	my_amg_manager.nFinnest_Pressure = 2;
	my_amg_manager.nu1_Pressure = 1;
	my_amg_manager.nu2_Pressure = 2;
	my_amg_manager.memory_size_Pressure = 9;
	my_amg_manager.ilu2_smoother_Pressure = 0; // 0 - не использовать, 1 - использовать.
	// Stress
	my_amg_manager.theta_Stress = 0.24;
	my_amg_manager.maximum_delete_levels_Stress = 0;
	my_amg_manager.nFinnest_Stress = 2;
	my_amg_manager.nu1_Stress = 1;
	my_amg_manager.nu2_Stress = 2;
	my_amg_manager.memory_size_Stress = 9;
	my_amg_manager.ilu2_smoother_Stress = 0; // 0 - не использовать, 1 - использовать.
	// AMG Splitting (coarsening)
	// Способ построения C-F разбиения: 0 - standart, 1 - RS 2, 2 - standart ST, 3 - RS2 ST. ST - strong transpose.
	// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
	my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; // standart strong transpose ST.
	my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; // standart strong transpose ST.
	my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; // standart strong transpose ST.
	my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; // standart strong transpose ST.
	my_amg_manager.icoarseningtype = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; // standart vs RS 2.
	// Stabilization BiCGStab.
	// 8.01.2017 Метод ван дер Ворста BiCGStab 
	// предобусловленный алгебраичесеким многосеточным методом.
	// 0 - используется просто алгебраический многосеточный метод без какого-либо привлечения алгоритмов подпространства Крылова,
	// 1 - Используется алгоритм Х. Ван дер Ворста BiCGStab [1992], предобусловленный алгебраическим многосеточным методом.
	// 2 - Используется алгоритм Саада и Шульца FGMRes [1986], предобусловленный алгебраическим многосеточным методом.
	// 3 - Нелинейный многосеточный метод (обновление правой части на каждом V цикле). Для нелинейных граничных условий.
	my_amg_manager.istabilizationTemp = 0; // none
	my_amg_manager.istabilizationSpeed = 0; // none
	my_amg_manager.istabilizationPressure = 0; // none
	my_amg_manager.istabilizationStress = 0; // none
	my_amg_manager.istabilization = 0; // none

	// номер применяемого патча.
	my_amg_manager.ipatch_number = 0; // 0 - патч не используется.

	// Печать лога на консоль.
	my_amg_manager.iprint_log = 1;
	my_amg_manager.iprint_log_Temperature = 1;
	my_amg_manager.iprint_log_Speed = 1;
	my_amg_manager.iprint_log_Pressure = 1;
	my_amg_manager.iprint_log_Stress = 1;

	// 01.07.2021 begin
	my_amg_manager.bthreshold_Temperature_auto = true;
	my_amg_manager.bthreshold_Speed_auto = true;
	my_amg_manager.bthreshold_Pressure_auto = true;
	my_amg_manager.bthreshold_Stress_auto = true;
	my_amg_manager.bthreshold_auto = true;

	my_amg_manager.bcf_reorder_Temperature = true;
	my_amg_manager.bcf_reorder_Speed = true;
	my_amg_manager.bcf_reorder_Pressure = true;
	my_amg_manager.bcf_reorder_Stress = true;
	my_amg_manager.bcf_reorder = true;
	// 01.07.2021 end

	// truncation for interpolation.
	// По умолчанию усечение интерполяции не используется.
	my_amg_manager.itruncation_interpolation = 0; // 0 - off
	my_amg_manager.itruncation_interpolation_Temperature = 0;
	my_amg_manager.itruncation_interpolation_Speed = 0;
	my_amg_manager.itruncation_interpolation_Pressure = 0;
	my_amg_manager.itruncation_interpolation_Stress = 0;
	// 0.2 recomended Stuben.
	my_amg_manager.truncation_interpolation = 0.2; // 0.2 recomended default value.
	my_amg_manager.truncation_interpolation_Temperature = 0.2;
	my_amg_manager.truncation_interpolation_Speed = 0.2;
	my_amg_manager.truncation_interpolation_Pressure = 0.2;
	my_amg_manager.truncation_interpolation_Stress = 0.2;

	// GMRES smoother.
	my_amg_manager.b_gmresTemp = false;
	my_amg_manager.b_gmresSpeed = false;
	my_amg_manager.b_gmresPressure = false;
	my_amg_manager.b_gmresStress = false;
	my_amg_manager.b_gmres = false;

	// BiCGStab smoother.
	my_amg_manager.b_ChebyshevSmootherTemp = false;
	my_amg_manager.b_ChebyshevSmootherSpeed = false;
	my_amg_manager.b_ChebyshevSmootherPressure = false;
	my_amg_manager.b_ChebyshevSmootherStress = false;
	my_amg_manager.b_ChebyshevSmoother = false;

	my_amg_manager.Chebyshev_degree = 5; // степень полинома Чебышева используемая в сглаживателе Чебышева.


	//spai0 smoother
	my_amg_manager.b_spai0Temp = false;
	my_amg_manager.b_spai0Speed = false;
	my_amg_manager.b_spai0Pressure = false;
	my_amg_manager.b_spai0Stress = false;
	my_amg_manager.b_spai0 = false;

	// 01.07.2021
	// Задавать ли порог threshold автоматическим образом или всёже вручную.
	my_amg_manager.bthreshold_Pressure_auto = true; my_amg_manager.bthreshold_Temperature_auto = true; my_amg_manager.bthreshold_Speed_auto = true; my_amg_manager.bthreshold_Stress_auto = true;
	my_amg_manager.bthreshold_auto = true;
	// Применять ли CF переупорядочивание (порядок обхода) при итерировании, например в методе Якоби.
	my_amg_manager.bcf_reorder_Pressure = true; my_amg_manager.bcf_reorder_Temperature = true; my_amg_manager.bcf_reorder_Speed = true; my_amg_manager.bcf_reorder_Stress = true;
	my_amg_manager.bcf_reorder = true;

	//fgmres(m_restart)
	my_amg_manager.m_restart = 20; // Количество итераций алгоритма fgmres перед перезапуском.

	// amg default settings:
	my_amg_manager.lfil = 2; // default value


	// AMGCL parameters
	// 0 - spai0; 1 - ilu0; 2- gauss-seidel; 3 - damped-jacobi;
	// 4 - spai1; 5 - chebyshev; 6 - iluk, k=1; 7 - iluk, k=2;
	// 8 - iluk, k=4; 9 - iluk, k=6; 10 - iluk, k=8; 11 - iluk, k=10.
	my_amg_manager.amgcl_smoother = 0; // 0 - spai0; default
	my_amg_manager.amgcl_selector = 1; // 0 - Ruge-Stueben (amg1r5 analog); 1 - smoother aggregation.
	my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab; // 0 - BiCGStab; 1 - FGMRes.

}

doublereal return_gold_const(integer irelx) {
	doublereal double_gold_const = -1.0;

	switch (irelx) {
	case 0:// Gauss-Seidel
		double_gold_const = -1.0;
		break;
	case 1:// ILUk, k=lfil; 
		double_gold_const  = -1.0;
		break;
	case 2:// Рунге - Кутта 3 порядка
		double_gold_const = -1.0;
		break;
	case 3:// Рунге - Кутта 5 порядка
		double_gold_const = -1.0;
		break;
	case 4:// damped Jacoby
		double_gold_const = -0.6667;
		break;
	case 5:// П. Роуч SOR
		double_gold_const  = 0.2;
		break;
	default: double_gold_const = -1.0;
		break;
	}
	return double_gold_const;
}

// Печатает логотип.
/*
void printLOGO() {
	printf("     #     #        O     ###   #####		######  #        ###   #               #\n");
	printf("   #  #    #             #   #  #			#       #       #   #   #             # \n");
	printf("   #  #    #        #   #       #====  		#       #      #     #   #     #     #  \n");
	printf("  ######   #        #   #       #====		###     #      #     #    #   # #   #   \n");
	printf(" #     #   #        #    #   #  #			#       #       #   #      # #   # #    \n");
	printf("#       #  #######  #     ###   #####		#       #######  ###        #     #    #\n");
	printf("\n");
	printf("version v0.57 2009 - 2021\n");
	Sleep(3000);
}
*/




// Вычисляет view factor между двумя пластинами составляющими между собой прямой угол
// и контактирующих по общему ребру.
// 2 августа 2016.
doublereal view_factor_perpendicular_rectangles_with_a_common_edge(doublereal x34, doublereal y34, doublereal z34)
{
	// Вычисление View Factors между двумя перпендикулярными пластинами имеющими общее ребро
	// на основе аналитически вычисленного значения интеграла.
	// Heat and Mass Transfer: Fundamentals & Applications. Fourth Edition.
	// Yunus A.Cengel, Afshin J.Ghajar. McGraw-Hill, 2011.
	doublereal h = z34 / x34;
	doublereal w = y34 / x34;
	if ((h > 120.0) || (w > 120.0)) {
		printf("H=%e W=%e. calculate view factor error...\n",h,w);
		//system("PAUSE");
		system("PAUSE");
		exit(1);
	}
	// Isidoro Martinez Radiative View Factors. 1995-2016.
	//doublereal M_PI = 3.1415926;//math.h
	doublereal a = ((1.0 + h*h)*(1.0 + w*w)) / (1.0 + h*h + w*w);
	doublereal b = (w*w*(1.0 + h*h + w*w)) / ((1.0 + w*w)*(h*h + w*w));
	doublereal c = (h*h*(1.0 + h*h + w*w)) / ((1.0 + h*h)*(h*h + w*w));
	doublereal Fi_j = (1.0 / (M_PI*w))*(h*atan(1.0 / h) + w*atan(1.0 / w) -
		sqrt(h*h + w*w)*atan(1.0 / (sqrt(h*h + w*w))) + 
		0.25*log(a*pow(b, w*w)*pow(c, h*h)));
	return Fi_j;
} // view_factor_perpendicular_rectangles_with_a_common_edge

// Вычисляет view factor между двумя параллельными одинаковыми пластинами.
doublereal view_factor_aligned_parallel_rectangles(doublereal x34, doublereal y34, doublereal L34)
{
    // Вычисление View Factors между двумя параллельными пластинами на основе аналитически вычисленного 
	// значения интеграла.
	// Heat and Mass Transfer: Fundamentals & Applications. Fourth Edition.
	// Yunus A. Cengel, Afshin J.Ghajar. McGraw-Hill, 2011.
	doublereal x_tilda = x34 / L34;
	doublereal y_tilda = y34 / L34;
	if ((x_tilda > 120.0) || (y_tilda > 120.0)) {
		printf("X=%e Y=%e. calculate view factor error...\n", x_tilda, y_tilda);
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	// Isidoro Martinez Radiative View Factors. 1995-2016.
	//doublereal M_PI = 3.1415926;// math.h
	doublereal x_tilda1 = sqrt(1.0 + x_tilda*x_tilda);
	doublereal y_tilda1 = sqrt(1.0 + y_tilda*y_tilda);
	doublereal Fi_j = (1.0 / (M_PI*x_tilda*y_tilda))*(log((x_tilda1*x_tilda1*y_tilda1*y_tilda1) / (x_tilda1*x_tilda1 + y_tilda1*y_tilda1 - 1.0)) +
		2.0*x_tilda*(y_tilda1*atan(x_tilda / y_tilda1) - atan(x_tilda)) + 2.0*y_tilda*(x_tilda1*atan(y_tilda / x_tilda1) - atan(y_tilda)));
	return Fi_j;
} //view_factor_aligned_parallel_rectangles

// Для вычисления усреднённой температуры на грани нужен
// список узлов грани (два внутренних соседа для внутренней грани) и 
// (граничный узел и ближайший к нему соседний внутренний для узлов на границе расчётной области).
typedef struct TMY_PAIR {
	integer node1, node21, node22, node23, node24; // если отсутствует то стоит -1.
	doublereal  dS1, dS2, dS3, dS4; // площадь ячейки.
	// 20 сентября 2016.
	// dS должно быть позже полностью удалено.
	//doublereal dS;
	TMY_PAIR() {
		node1=-1; node21=-1; node22=-1; node23=-1; node24=-1; // если отсутствует то стоит -1.
		dS1=0.0; dS2=0.0; dS3=0.0; dS4=0.0; // площадь ячейки.
		// 20 сентября 2016.
		// dS должно быть позже полностью удалено.
		//doublereal dS;
	}
} MY_PAIR;

// Параметр нижней релаксации для вычисления плотностей радиационных потоков.
const doublereal alpha_Radiation_block = 0.2;//0.1 optimum
// Постоянная Стефана-Больцмана для вычисления плотностей радиационных потоков.
// постоянная Стефана-Больцмана для расчёта теплообмена илучением.
const doublereal STEFAN_BOLCMAN_CONST = 5.670367e-8; // W*m!-2*K!-4.
// фиксированное число итераций для нахождения плотностей радиационных потоков.
const integer maxiter_Radiation_block = 1000; // вообще хватает и 200, здесь с запасом.


typedef struct TBLOCKRADIATION {
	// emissivity:
	doublereal emissW, emissE, emissS, emissN, emissB, emissT;

	// View Factors for Prism Object.
	doublereal FWE, FWS, FWN, FWB, FWT;
	doublereal FEW, FES, FEN, FEB, FET;
	doublereal FSW, FSE, FSN, FSB, FST;
	doublereal FNW, FNE, FNS, FNB, FNT;
	doublereal FBW, FBE, FBS, FBN, FBT;
	doublereal FTW, FTE, FTS, FTN, FTB;
	// Температуры в Кельвинах на гранях Prism Object:
	// среднее арифметическое температуры на грани Prism Object.
	doublereal TempW, TempE, TempS, TempN, TempB, TempT;
	// Плотности радиационных тепловых потоков на гранях Prism Object.
	doublereal JW, JE, JS, JN, JB, JT;
	bool binternalRadiation;
	

	// список узлов на каждой из граней и их количество:
	MY_PAIR* nodelistW;
	integer nodelistWsize;
	MY_PAIR* nodelistE;
	integer nodelistEsize ;
	MY_PAIR* nodelistS;
	integer nodelistSsize ;
	MY_PAIR* nodelistN;
	integer nodelistNsize ;
	MY_PAIR* nodelistB;
	integer nodelistBsize ;
	MY_PAIR* nodelistT;
	integer nodelistTsize;

	TBLOCKRADIATION() {
// emissivity:
	 emissW = -1.0e30; emissE = -1.0e30; emissS = -1.0e30;
	 emissN = -1.0e30; emissB = -1.0e30; emissT = -1.0e30;
	 binternalRadiation=false;
	// View Factors for Prism Object.
	 FWE = -1.0e30; FWS = -1.0e30; FWN = -1.0e30; FWB = -1.0e30; FWT = -1.0e30;
	 FEW = -1.0e30; FES = -1.0e30; FEN = -1.0e30; FEB = -1.0e30; FET = -1.0e30;
	 FSW = -1.0e30; FSE = -1.0e30; FSN = -1.0e30; FSB = -1.0e30; FST = -1.0e30;
	 FNW = -1.0e30; FNE = -1.0e30; FNS = -1.0e30; FNB = -1.0e30; FNT = -1.0e30;
	 FBW = -1.0e30; FBE = -1.0e30; FBS = -1.0e30; FBN = -1.0e30; FBT = -1.0e30;
	 FTW = -1.0e30; FTE = -1.0e30; FTS = -1.0e30; FTN = -1.0e30; FTB = -1.0e30;
	// Температуры в Кельвинах на гранях Prism Object:
	// среднее арифметическое температуры на грани Prism Object.
	 TempW = -1.0e30; TempE = -1.0e30; TempS = -1.0e30;
	 TempN = -1.0e30; TempB = -1.0e30; TempT = -1.0e30;
	// Плотности радиационных тепловых потоков на гранях Prism Object.
	 JW=0.0; JE=0.0; JS = 0.0; JN = 0.0; JB = 0.0; JT = 0.0;
	

	// список узлов на каждой из граней и их количество:
	 nodelistW=nullptr;
	 nodelistWsize=-4;
	 nodelistE=nullptr;
	 nodelistEsize = -4;
	 nodelistS=nullptr;
	 nodelistSsize = -4;
	 nodelistN=nullptr;
	 nodelistNsize = -4;
	 nodelistB=nullptr;
	 nodelistBsize = -4;
	 nodelistT=nullptr;
	 nodelistTsize = -4;
	}

} BLOCKRADIATION;

// стиль зависимости мощности тепловыделения от времени.
// 0 - не зависит, 1 - square wave зависимость,
// 2 - square wave 2, 3 - hot cold режим,
// 4 - piecewise const.
enum class POWER_TIME_DEPEND {
	CONST_POWER = 0,
	SQUARE_WAVE=1,
	SQUARE_WAVE2 = 2,
	HOT_COLD=3,
	PIECEWISE_CONST=4
};


// блок
typedef struct TBLOCK {
	char name[80]; // имя объёмного геометрического объекта.

	PHYSICS_TYPE_IN_BODY itype; // тип SOLID, HOLLOW или FLUID.
	GEOM g;
	integer imatid; // идентификатор материала в библиотеке
	//doublereal Sc; // мощность тепловыделения на единицу объёма
	// 19 11 2016 Температурно зависимая мощность тепловыделения.
	integer n_Sc;
	doublereal *arr_Sc, *temp_Sc;
    // стиль зависимости мощности тепловыделения от времени.
	// 0 - не зависит, 1 - square wave зависимость,
	// 2 - square wave 2, 3 - hot cold режим,
	// 4 - piecewise const.
	POWER_TIME_DEPEND ipower_time_depend;
	// Всё относящееся к теплообмену излучением:
	// излучательные способности поверхностей, модель вакуумного промежутка.
	BLOCKRADIATION radiation;
	bool bvisible; // Виден ли блок при экспорте в tecplot.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id;
	// Фиксировать боковую стенку цилиндра ?
	bool CylinderFixed;

	TBLOCK() {
		name[0] = '\0';

		itype = PHYSICS_TYPE_IN_BODY::HOLLOW; // тип SOLID, HOLLOW или FLUID.
		//GEOM g;
		imatid = -4; // идентификатор материала в библиотеке
		//doublereal Sc; // мощность тепловыделения на единицу объёма
		// 19 11 2016 Температурно зависимая мощность тепловыделения.
		n_Sc = -4;
		arr_Sc=nullptr; temp_Sc=nullptr;
		// стиль зависимости мощности тепловыделения от времени.
		// 0 - не зависит, 1 - square wave зависимость,
		// 2 - square wave 2, 3 - hot cold режим,
		// 4 - piecewise const.
		ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
		// Всё относящееся к теплообмену излучением:
		// излучательные способности поверхностей, модель вакуумного промежутка.
		//BLOCKRADIATION radiation;
		bvisible=true; // Виден ли блок при экспорте в tecplot.

		// принадлежность объединению 
		// 0 - не принадлежит (принадлежит кабинету),
		// n > 0 принадлежит объединению с номером n.
		iunion_id = -4;
		// Фиксировать боковую стенку цилиндра ?
		CylinderFixed=false;
	}
} BLOCK;


// 1.09.2017.
doublereal my_sign(doublereal a) {
	if (a >= 0.0) {
		return 1.0;
	}
	else {
		return -1.0;
	}
} // my_sign





  // метод суммирования углов.
  // 1.09.2017.
bool in_polygon_check(TOCHKA p, integer nsizei,
	doublereal* &xi, doublereal* &yi, doublereal* &zi,
	doublereal* &hi, integer iPlane_obj2, int &k,
	int ib)
{

	bool bfound = false;

	doublereal sumphi = 0.0;
	integer i_73 = 0;

	switch (iPlane_obj2) {
	case XY_PLANE:
		for (i_73 = 0; i_73 < nsizei - 1; i_73++) {
			sumphi += acos(((xi[i_73] - p.x)*(xi[i_73 + 1] - p.x) + (yi[i_73] - p.y)*(yi[i_73 + 1] - p.y)) / (sqrt((xi[i_73] - p.x)*(xi[i_73] - p.x) + (yi[i_73] - p.y)*(yi[i_73] - p.y))*sqrt((xi[i_73 + 1] - p.x)*(xi[i_73 + 1] - p.x) + (yi[i_73 + 1] - p.y)*(yi[i_73 + 1] - p.y))))*my_sign((xi[i_73] - p.x)*(yi[i_73 + 1] - p.y) - (xi[i_73 + 1] - p.x)*(yi[i_73] - p.y));
		}
		i_73 = nsizei - 1;
		sumphi += acos(((xi[i_73] - p.x)*(xi[0] - p.x) + (yi[i_73] - p.y)*(yi[0] - p.y)) / (sqrt((xi[i_73] - p.x)*(xi[i_73] - p.x) + (yi[i_73] - p.y)*(yi[i_73] - p.y))*sqrt((xi[0] - p.x)*(xi[0] - p.x) + (yi[0] - p.y)*(yi[0] - p.y))))*my_sign((xi[i_73] - p.x)*(yi[0] - p.y) - (xi[0] - p.x)*(yi[i_73] - p.y));
		break;
	case XZ_PLANE:
		for (i_73 = 0; i_73 < nsizei - 1; i_73++) {
			sumphi += acos(((xi[i_73] - p.x)*(xi[i_73 + 1] - p.x) + (zi[i_73] - p.z)*(zi[i_73 + 1] - p.z)) / (sqrt((xi[i_73] - p.x)*(xi[i_73] - p.x) + (zi[i_73] - p.z)*(zi[i_73] - p.z))*sqrt((xi[i_73 + 1] - p.x)*(xi[i_73 + 1] - p.x) + (zi[i_73 + 1] - p.z)*(zi[i_73 + 1] - p.z))))*my_sign((xi[i_73] - p.x)*(zi[i_73 + 1] - p.z) - (xi[i_73 + 1] - p.x)*(zi[i_73] - p.z));
		}
		i_73 = nsizei - 1;
		sumphi += acos(((xi[i_73] - p.x)*(xi[0] - p.x) + (zi[i_73] - p.z)*(zi[0] - p.z)) / (sqrt((xi[i_73] - p.x)*(xi[i_73] - p.x) + (zi[i_73] - p.z)*(zi[i_73] - p.z))*sqrt((xi[0] - p.x)*(xi[0] - p.x) + (zi[0] - p.z)*(zi[0] - p.z))))*my_sign((xi[i_73] - p.x)*(zi[0] - p.z) - (xi[0] - p.x)*(zi[i_73] - p.z));
		break;
	case YZ_PLANE:
		for (i_73 = 0; i_73 < nsizei - 1; i_73++) {
			sumphi += acos(((yi[i_73] - p.y)*(yi[i_73 + 1] - p.y) + (zi[i_73] - p.z)*(zi[i_73 + 1] - p.z)) / (sqrt((yi[i_73] - p.y)*(yi[i_73] - p.y) + (zi[i_73] - p.z)*(zi[i_73] - p.z))*sqrt((yi[i_73 + 1] - p.y)*(yi[i_73 + 1] - p.y) + (zi[i_73 + 1] - p.z)*(zi[i_73 + 1] - p.z))))*my_sign((yi[i_73] - p.y)*(zi[i_73 + 1] - p.z) - (yi[i_73 + 1] - p.y)*(zi[i_73] - p.z));
		}
		i_73 = nsizei - 1;
		sumphi += acos(((yi[i_73] - p.y)*(yi[0] - p.y) + (zi[i_73] - p.z)*(zi[0] - p.z)) / (sqrt((yi[i_73] - p.y)*(yi[i_73] - p.y) + (zi[i_73] - p.z)*(zi[i_73] - p.z))*sqrt((yi[0] - p.y)*(yi[0] - p.y) + (zi[0] - p.z)*(zi[0] - p.z))))*my_sign((yi[i_73] - p.y)*(zi[0] - p.z) - (yi[0] - p.y)*(zi[i_73] - p.z));
		break;
	}

	if (fabs(sumphi) > 1.0e-7) {
		// точка лежит внутри полигона.
		switch (iPlane_obj2) {
		case XY_PLANE:
			if (hi[0] > 0) {
				if ((p.z >= zi[0]) && (p.z <= zi[0] + hi[0])) {
					k = ib;
					bfound = true;
				}
			}
			else {
				if ((p.z >= zi[0] + hi[0]) && (p.z <= zi[0] )) {
					k = ib;
					bfound = true;
				}
			}
			break;
		case XZ_PLANE:
			if (hi[0] > 0) {
				if ((p.y >= yi[0]) && (p.y <= yi[0] + hi[0])) {
					k = ib;
					bfound = true;
				}
			}
			else {
				if ((p.y >= yi[0] + hi[0]) && (p.y <= yi[0])) {
					k = ib;
					bfound = true;
				}
			}
			break;
		case YZ_PLANE:
			if (hi[0] > 0) {
				if ((p.x >= xi[0]) && (p.x <= xi[0] + hi[0])) {
					k = ib;
					bfound = true;
				}
			}
			else
			{
				if ((p.x >= xi[0] + hi[0]) && (p.x <= xi[0] )) {
					k = ib;
					bfound = true;
				}
			}
			break;
		}
	}

	return bfound;

} // in_polygon_check

bool in_polygon_check_first(TOCHKA p, integer nsizei,
	doublereal* &xi, doublereal* &yi, doublereal* &zi, doublereal* &hi,
	integer iPlane_obj2, integer &k, integer ib) {

	bool bfound = false;

	doublereal dpolygon_tolerance = 1e-3;

	// Polygon
	// точка принадлежит полигону если она принадлежит одному из треугольников его образующих.
	// Задача триангуляции выпуклого полигона состоит из вычисления его центра масс и образования треугольников
	// каждый из которых имеет вершину в этом центре масс и ребро на одной из сторон выпуклого многогранника.
	// Точка принадлежит треугольнику только когда площадь треугольника в точности равна сумме трёх площадей треугольников,
	// образуемых этой точкой и одной из строн первоначального треугольника.

	// Работает только для выпуклого многоугольника и только когда
	// все начальные уровни и все высоты одинаковы.
	doublereal xavg = 0.0, yavg = 0.0, zavg = 0.0;
	for (integer i_73 = 0; i_73 < nsizei; i_73++) {
		xavg += xi[i_73] / nsizei;
		yavg += yi[i_73] / nsizei;
		zavg += zi[i_73] / nsizei;
	}
	doublereal Sgl = 0.0;
	doublereal Sloc = 0.0;
	doublereal minx84, maxx84, miny84, maxy84, minz84, maxz84;
	doublereal epsilon_tri = 1.0e-8;

	switch (iPlane_obj2) {
	case XY_PLANE:
		minx84 = 1.0e36;
		maxx84 = -1.0e36;
		miny84 = 1.0e36;
		maxy84 = -1.0e36;
		minz84 = 1.0e36;
		maxz84 = -1.0e36;
		for (integer i_73 = 0; i_73 < nsizei; i_73++) {
			// minimum
			if (xi[i_73] < minx84) minx84 = xi[i_73];
			if (yi[i_73] < miny84) miny84 = yi[i_73];
			if (zi[i_73] < minz84) minz84 = zi[i_73];
			// maximum
			if (xi[i_73] > maxx84) maxx84 = xi[i_73];
			if (yi[i_73] > maxy84) maxy84 = yi[i_73];
			if (zi[i_73] > maxz84) maxz84 = zi[i_73];
		}
		//printf("epsilon_tri=%e length=%e\n", epsilon_tri, fabs((maxx84 - minx84)*(maxy84 - miny84)));
		//system("pause");
		//epsilon_tri = dpolygon_tolerance*sqrt((maxx84 - minx84)*(maxx84 - minx84) + (maxy84 - miny84)*(maxy84 - miny84));
		epsilon_tri = dpolygon_tolerance*fabs((maxx84 - minx84)*(maxy84 - miny84));

		for (integer i_73 = 0; i_73 < nsizei - 1; i_73++) {
			Sgl = 0.0;
			Sgl = 0.5*((xi[i_73 + 1] - xi[i_73])*(yavg - yi[i_73]) - (yi[i_73 + 1] - yi[i_73])*(xavg - xi[i_73]));
			Sgl = fabs(Sgl);

			Sloc = 0.5*((p.x - xi[i_73])*(yavg - yi[i_73]) - (p.y - yi[i_73])*(xavg - xi[i_73]));
			Sgl -= fabs(Sloc);
			Sloc = 0.5*((xi[i_73 + 1] - p.x)*(yavg - p.y) - (yi[i_73 + 1] - p.y)*(xavg - p.x));
			Sgl -= fabs(Sloc);
			Sloc = 0.5*((xi[i_73 + 1] - xi[i_73])*(p.y - yi[i_73]) - (yi[i_73 + 1] - yi[i_73])*(p.x - xi[i_73]));
			Sgl -= fabs(Sloc);

			if (fabs(Sgl) < epsilon_tri) {
				if ((p.z >= zi[0]) && (p.z <= zi[0] + hi[0])) {
					k = ib;
					bfound = true;
				}
			}
		}
		Sgl = 0.0;
		Sgl = 0.5*((xi[0] - xi[nsizei - 1])*(yavg - yi[nsizei - 1]) - (yi[0] - yi[nsizei - 1])*(xavg - xi[nsizei - 1]));
		Sgl = fabs(Sgl);

		Sloc = 0.5*((p.x - xi[nsizei - 1])*(yavg - yi[nsizei - 1]) - (p.y - yi[nsizei - 1])*(xavg - xi[nsizei - 1]));
		Sgl -= fabs(Sloc);
		Sloc = 0.5*((xi[0] - p.x)*(yavg - p.y) - (yi[0] - p.y)*(xavg - p.x));
		Sgl -= fabs(Sloc);
		Sloc = 0.5*((xi[0] - xi[nsizei - 1])*(p.y - yi[nsizei - 1]) - (yi[0] - yi[nsizei - 1])*(p.x - xi[nsizei - 1]));
		Sgl -= fabs(Sloc);

		if (fabs(Sgl) < epsilon_tri) {
			if ((p.z >= zi[0]) && (p.z <= zi[0] + hi[0])) {
				k = ib;
				bfound = true;
			}
		}
		break;
	case XZ_PLANE:
		minx84 = 1.0e36;
		maxx84 = -1.0e36;
		miny84 = 1.0e36;
		maxy84 = -1.0e36;
		minz84 = 1.0e36;
		maxz84 = -1.0e36;
		for (integer i_73 = 0; i_73 < nsizei; i_73++) {
			// minimum
			if (xi[i_73] < minx84) minx84 = xi[i_73];
			if (yi[i_73] < miny84) miny84 = yi[i_73];
			if (zi[i_73] < minz84) minz84 = zi[i_73];
			// maximum
			if (xi[i_73] > maxx84) maxx84 = xi[i_73];
			if (yi[i_73] > maxy84) maxy84 = yi[i_73];
			if (zi[i_73] > maxz84) maxz84 = zi[i_73];
		}
		//printf("epsilon_tri=%e length=%e\n", epsilon_tri, fabs((maxx84 - minx84)*(maxz84 - minz84)));
		//system("pause");
		//epsilon_tri = dpolygon_tolerance*sqrt((maxx84 - minx84)*(maxx84 - minx84) + (maxz84 - minz84)*(maxz84 - minz84));
		epsilon_tri = dpolygon_tolerance*fabs((maxx84 - minx84)*(maxz84 - minz84));

		for (integer i_73 = 0; i_73 < nsizei - 1; i_73++) {
			Sgl = 0.0;
			Sgl = 0.5*((xi[i_73 + 1] - xi[i_73])*(zavg - zi[i_73]) - (zi[i_73 + 1] - zi[i_73])*(xavg - xi[i_73]));
			Sgl = fabs(Sgl);

			Sloc = 0.5*((p.x - xi[i_73])*(zavg - zi[i_73]) - (p.z - zi[i_73])*(xavg - xi[i_73]));
			Sgl -= fabs(Sloc);
			Sloc = 0.5*((xi[i_73 + 1] - p.x)*(zavg - p.z) - (zi[i_73 + 1] - p.z)*(xavg - p.x));
			Sgl -= fabs(Sloc);
			Sloc = 0.5*((xi[i_73 + 1] - xi[i_73])*(p.z - zi[i_73]) - (zi[i_73 + 1] - zi[i_73])*(p.x - xi[i_73]));
			Sgl -= fabs(Sloc);

			if (fabs(Sgl) < epsilon_tri) {
				if ((p.y >= yi[0]) && (p.y <= yi[0] + hi[0])) {
					k = ib;
					bfound = true;
				}
			}
		}
		Sgl = 0.0;
		Sgl = 0.5*((xi[0] - xi[nsizei - 1])*(zavg - zi[nsizei - 1]) - (zi[0] - zi[nsizei - 1])*(xavg - xi[nsizei - 1]));
		Sgl = fabs(Sgl);

		Sloc = 0.5*((p.x - xi[nsizei - 1])*(zavg - zi[nsizei - 1]) - (p.z - zi[nsizei - 1])*(xavg - xi[nsizei - 1]));
		Sgl -= fabs(Sloc);
		Sloc = 0.5*((xi[0] - p.x)*(zavg - p.z) - (zi[0] - p.z)*(xavg - p.x));
		Sgl -= fabs(Sloc);
		Sloc = 0.5*((xi[0] - xi[nsizei - 1])*(p.z - zi[nsizei - 1]) - (zi[0] - zi[nsizei - 1])*(p.x - xi[nsizei - 1]));
		Sgl -= fabs(Sloc);

		if (fabs(Sgl) < epsilon_tri) {
			if ((p.y >= yi[0]) && (p.y <= yi[0] + hi[0])) {
				k = ib;
				bfound = true;
			}
		}
		break;
	case YZ_PLANE:
		minx84 = 1.0e36;
		maxx84 = -1.0e36;
		miny84 = 1.0e36;
		maxy84 = -1.0e36;
		minz84 = 1.0e36;
		maxz84 = -1.0e36;
		for (integer i_73 = 0; i_73 < nsizei; i_73++) {
			// minimum
			if (xi[i_73] < minx84) minx84 = xi[i_73];
			if (yi[i_73] < miny84) miny84 = yi[i_73];
			if (zi[i_73] < minz84) minz84 = zi[i_73];
			// maximum
			if (xi[i_73] > maxx84) maxx84 = xi[i_73];
			if (yi[i_73] > maxy84) maxy84 = yi[i_73];
			if (zi[i_73] > maxz84) maxz84 = zi[i_73];
		}
		//printf("epsilon_tri=%e length=%e\n", epsilon_tri, fabs((maxz84 - minz84)*(maxy84 - miny84)));
		//system("pause");
		//epsilon_tri = dpolygon_tolerance*sqrt((maxz84 - minz84)*(maxz84 - minz84) + (maxy84 - miny84)*(maxy84 - miny84));
		epsilon_tri = dpolygon_tolerance*fabs((maxz84 - minz84)*(maxy84 - miny84));

		for (integer i_73 = 0; i_73 < nsizei - 1; i_73++) {
			Sgl = 0.0;
			Sgl = 0.5*((zi[i_73 + 1] - zi[i_73])*(yavg - yi[i_73]) - (yi[i_73 + 1] - yi[i_73])*(zavg - zi[i_73]));
			Sgl = fabs(Sgl);

			Sloc = 0.5*((p.z - zi[i_73])*(yavg - yi[i_73]) - (p.y - yi[i_73])*(zavg - zi[i_73]));
			Sgl -= fabs(Sloc);
			Sloc = 0.5*((zi[i_73 + 1] - p.z)*(yavg - p.y) - (yi[i_73 + 1] - p.y)*(zavg - p.z));
			Sgl -= fabs(Sloc);
			Sloc = 0.5*((zi[i_73 + 1] - zi[i_73])*(p.y - yi[i_73]) - (yi[i_73 + 1] - yi[i_73])*(p.z - zi[i_73]));
			Sgl -= fabs(Sloc);

			if (fabs(Sgl) < epsilon_tri) {
				if ((p.x >= xi[0]) && (p.x <= xi[0] + hi[0])) {
					k = ib;
					bfound = true;
				}
			}
		}
		Sgl = 0.0;
		Sgl = 0.5*((zi[0] - zi[nsizei - 1])*(yavg - yi[nsizei - 1]) - (yi[0] - yi[nsizei - 1])*(zavg - zi[nsizei - 1]));
		Sgl = fabs(Sgl);

		Sloc = 0.5*((p.z - zi[nsizei - 1])*(yavg - yi[nsizei - 1]) - (p.y - yi[nsizei - 1])*(zavg - zi[nsizei - 1]));
		Sgl -= fabs(Sloc);
		Sloc = 0.5*((zi[0] - p.z)*(yavg - p.y) - (yi[0] - p.y)*(zavg - p.z));
		Sgl -= fabs(Sloc);
		Sloc = 0.5*((zi[0] - zi[nsizei - 1])*(p.y - yi[nsizei - 1]) - (yi[0] - yi[nsizei - 1])*(p.z - zi[nsizei - 1]));
		Sgl -= fabs(Sloc);

		if (fabs(Sgl) < epsilon_tri) {
			if ((p.x >= xi[0]) && (p.x <= xi[0] + hi[0])) {
				k = ib;
				bfound = true;
			}
		}
		break;
	}

	return bfound;

} // in_polygon_check_first


bool in_polygon(TOCHKA p, integer nsizei, doublereal* &xi, doublereal* &yi, doublereal* &zi, doublereal* &hi, integer iPlane_obj2, int &k, int ib) {

	bool bfound = false;

	// Первоначальная самописная реализация основанная на идее равентва площадей.
	// Сначала мноугольник разбивается на треугольники (триангулируется), а затем
	// сканируются все треугольники по очереди и если точка принадлежит треугольнику то площадь 
	// треугольника равна сумме площадей трёх треугольников образованных исследуемой точкой и вршинами первоначального треугольника. 
	//bfound = in_polygon_check_first(p, nsizei, xi, yi, zi, hi, iPlane_obj2, k, ib);

	// Теоретическиобоснованная версия проверки алгоритм которой найден в википедии.
	bfound = in_polygon_check(p, nsizei, xi, yi, zi, hi, iPlane_obj2, k, ib);

	return bfound;

}

// Площадь треугольника по трем вершинам.
doublereal square_triangle(doublereal x1, doublereal x2, doublereal x3, 
	doublereal y1, doublereal y2, doublereal y3) {
	
	return 0.5 * fabs((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3));
}

// Приближенно вычисляет объем полигона.
// 23.10.2020
doublereal Volume_polygon(integer nsizei,
	doublereal* &xi, doublereal* &yi, doublereal* &zi,
	doublereal* &hi, integer iPlane_obj2)
{

	// точность увеличивается при увеличении величины resolution_x и resolution_y.
	doublereal resolution_x = 100.0;
	doublereal resolution_y = 100.0;
	doublereal resolution_z = 100.0;

	doublereal xS = 1.0e30;
	doublereal xE = -1.0e30;
	doublereal yS = 1.0e30;
	doublereal yE = -1.0e30;
	doublereal zS = 1.0e30;
	doublereal zE = -1.0e30;
	doublereal hx = 0.0, hy = 0.0, hz = 0.0;
	TOCHKA p;
	int k, ib;
	doublereal vol = 0.0;

	switch (iPlane_obj2) {
	case XY_PLANE:
		resolution_x = 100.0;
		resolution_y = 100.0;
		resolution_z = 1.0;

		switch (nsizei) {
		case 3:
			vol = hi[0] * square_triangle(xi[0],xi[1],xi[2], yi[0], yi[1], yi[2]);
			break;
		case 4:
			// Внимание только если четырехугольник выпуклый!!!
			vol = hi[0] * (square_triangle(xi[0], xi[1], xi[2], yi[0], yi[1], yi[2])+ square_triangle(xi[0], xi[3], xi[2], yi[0], yi[3], yi[2]));
			break;
		default:

			for (integer i_73 = 0; i_73 <= nsizei - 1; i_73++) {
				if (xi[i_73] < xS) xS = xi[i_73];
				if (xi[i_73] > xE) xE = xi[i_73];
				if (yi[i_73] < yS) yS = yi[i_73];
				if (yi[i_73] > yE) yE = yi[i_73];
			}
			hx = fabs(xE - xS) / resolution_x;
			hy = fabs(yE - yS) / resolution_y;
#pragma omp parallel for private(k, ib) reduction(+: vol)
			for (integer i_73 = 0; i_73 < static_cast<integer>(resolution_x); i_73++) {
				for (integer j_73 = 0; j_73 < static_cast<integer>(resolution_y); j_73++) {
					p.x = xS + i_73 * hx + 0.5 * hx;
					p.y = yS + j_73 * hy + 0.5 * hy;
					p.z = zi[0] + 0.5 * hi[0];
					ib = -1;
					if (in_polygon(p, nsizei, xi, yi, zi, hi, iPlane_obj2, k, ib)) {
						vol += hi[0] * hx * hy;
					}
				}
			}
			break;
		}
		break;

	case XZ_PLANE:
		resolution_x = 100.0;
		resolution_y = 1.0;
		resolution_z = 100.0;

		switch (nsizei) {
		case 3:
			vol = hi[0] * square_triangle(xi[0], xi[1], xi[2], zi[0], zi[1], zi[2]);
			break;
		case 4:
			// Внимание только если четырехугольник выпуклый!!!
			vol = hi[0] * (square_triangle(xi[0], xi[1], xi[2], zi[0], zi[1], zi[2]) + square_triangle(xi[0], xi[3], xi[2], zi[0], zi[3], zi[2]));
			break;
		default:

			for (integer i_73 = 0; i_73 <= nsizei - 1; i_73++) {
				if (xi[i_73] < xS) xS = xi[i_73];
				if (xi[i_73] > xE) xE = xi[i_73];
				if (zi[i_73] < zS) zS = zi[i_73];
				if (zi[i_73] > zE) zE = zi[i_73];
			}
			hx = fabs(xE - xS) / resolution_x;
			hz = fabs(zE - zS) / resolution_z;
#pragma omp parallel for private(k, ib) reduction(+: vol)
			for (integer i_73 = 0; i_73 < static_cast<integer>(resolution_x); i_73++) {
				for (integer j_73 = 0; j_73 < static_cast<integer>(resolution_z); j_73++) {
					p.x = xS + i_73 * hx + 0.5 * hx;
					p.z = zS + j_73 * hz + 0.5 * hz;
					p.y = yi[0] + 0.5 * hi[0];
					ib = -1;
					if (in_polygon(p, nsizei, xi, yi, zi, hi, iPlane_obj2, k, ib)) {
						vol += hi[0] * hx * hz;
					}
				}
			}
			break;
		}
		break;

	case YZ_PLANE:
		resolution_x = 1.0;
		resolution_y = 100.0;
		resolution_z = 100.0;


		switch (nsizei) {
		case 3:
			vol = hi[0] * square_triangle(yi[0], yi[1], yi[2], zi[0], zi[1], zi[2]);
			break;
		case 4:
			// Внимание только если четырехугольник выпуклый!!!
			vol = hi[0] * (square_triangle(yi[0], yi[1], yi[2], zi[0], zi[1], zi[2]) + square_triangle(yi[0], yi[3], yi[2], zi[0], zi[3], zi[2]));
			break;
		default:

			for (integer i_73 = 0; i_73 <= nsizei - 1; i_73++) {
				if (yi[i_73] < yS) yS = yi[i_73];
				if (yi[i_73] > yE) yE = yi[i_73];
				if (zi[i_73] < zS) zS = zi[i_73];
				if (zi[i_73] > zE) zE = zi[i_73];
			}
			hy = fabs(yE - yS) / resolution_y;
			hz = fabs(zE - zS) / resolution_z;
#pragma omp parallel for private(k, ib) reduction(+: vol)
			for (integer i_73 = 0; i_73 < static_cast<integer>(resolution_y); i_73++) {
				for (integer j_73 = 0; j_73 < static_cast<integer>(resolution_z); j_73++) {
					p.y = yS + i_73 * hy + 0.5 * hy;
					p.z = zS + j_73 * hz + 0.5 * hz;
					p.x = xi[0] + 0.5 * hi[0];
					ib = -1;
					if (in_polygon(p, nsizei, xi, yi, zi, hi, iPlane_obj2, k, ib)) {
						vol += hi[0] * hy * hz;
					}
				}
			}
			break;
		}

		break;
	}

	return vol;
}




integer ihash_key_shorter(doublereal g, doublereal min_g, doublereal max_g )
{
	return (static_cast<integer>(isize_shorter_hash*((g - min_g) / (1.05*(max_g - min_g)))));

} // ihash_key_73


  //здесь будет управление допуском на местах.
doublereal shorter_length_for_simplificationXold(doublereal g, BLOCK* &b, integer lb,
	WALL* &w, integer lw, SOURCE* &s, integer ls) {
	if (!bFULL_AUTOMATIC) {
		return shorter_length_for_simplificationX_BASIC;
	}
	else {

		doublereal db = 1.0e30;


		//integer key = static_cast<integer>(isize_shorter_hash*((g - b[0].g.xS) / (1.05*(b[0].g.xE - b[0].g.xS))));
		integer key = ihash_key_shorter(g, b[0].g.xS, b[0].g.xE);

		if (!bshorter_hash_X[key]) {

			for (integer i = 0; i < lb; ++i) {
				if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == POLYGON)) {
					doublereal dmult = dbody_R_out_multiplyer;
					if (((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.yE - b[i].g.yS)) && ((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.zE - b[i].g.zS))) {
						// Ситуация слой термопасты или клея.
						dmult = dbody_dist_multiplyer;
					}
					if (((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.yE - b[i].g.yS)) || ((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.zE - b[i].g.zS))) {
						// Ситуация палец затвора.
						dmult = dbody_LGATE_situation_multiplyer;
					}
					if ((g >= b[i].g.xS - dmult *(b[i].g.xE - b[i].g.xS)) && (g <= b[i].g.xE + dmult *(b[i].g.xE - b[i].g.xS))) {
						// g в зоне влияния блока с номером i.
						if (((1.0 / rdivision_interval)*(b[i].g.xE - b[i].g.xS)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.xE - b[i].g.xS));
					}
				}
				if (b[i].g.itypegeom == CYLINDER) {
					switch (b[i].g.iPlane) {
					case XY_PLANE: case XZ_PLANE:
						if (b[i].g.R_in_cyl < 1.0e-36) {
							if ((g >= b[i].g.xC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.xC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl));
							}
						}
						else {
							if ((g >= b[i].g.xC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.xC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval)*(b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_in_cyl));
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl));
							}
						}
						break;
					case YZ_PLANE:
						if ((g >= b[i].g.xC - dbody_Hcyl_multiplyer *b[i].g.Hcyl) && (g <= b[i].g.xC + b[i].g.Hcyl + dbody_Hcyl_multiplyer *b[i].g.Hcyl)) {
							// g в зоне влияния блока с номером i.
							if (((1.0 / rdivision_interval)*(b[i].g.Hcyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.Hcyl));
						}
						break;
					}
				}
			}
			if (lb == 1) {
				db = 1.0e-10;
			}
			else {
				// 1/200
				if (db > dcabinet_delta_multiplyer *(fabs(b[0].g.xE - b[0].g.xS))) db = dcabinet_delta_multiplyer *(fabs(b[0].g.xE - b[0].g.xS));
			}
			for (integer i = 0; i < lw; ++i) {
				doublereal dmult = dbody_R_out_multiplyer;
				switch (w[i].g.iPlane) {
				case XY_PLANE: case XZ_PLANE:
					if ((g >= w[i].g.xS - dmult * (w[i].g.xE - w[i].g.xS)) && (g <= w[i].g.xE + dmult * (w[i].g.xE - w[i].g.xS))) {
						// g в зоне влияния стенки с номером i.
						if (((1.0 / rdivision_interval)*(w[i].g.xE - w[i].g.xS)) < db) db = ((1.0 / rdivision_interval)*(w[i].g.xE - w[i].g.xS));
					}
					break;
				}
			}

			for (integer i = 0; i < ls; ++i) {
				doublereal dmult = dbody_R_out_multiplyer;
				switch (s[i].g.iPlane) {
				case XY_PLANE: case XZ_PLANE:
					if ((g >= s[i].g.xS - dmult * (s[i].g.xE - s[i].g.xS)) && (g <= s[i].g.xE + dmult * (s[i].g.xE - s[i].g.xS))) {
						// g в зоне влияния источника с номером i.
						if (((1.0 / rdivision_interval)*(s[i].g.xE - s[i].g.xS)) < db) db = ((1.0 / rdivision_interval)*(s[i].g.xE - s[i].g.xS));
					}
					break;
				}
			}

			bshorter_hash_X[key] = true;
			shorter_hash_X[key] = db;
		}
		else {
			db = shorter_hash_X[key];
		}

		return db;
	}
}

//здесь будет управление допуском на местах.
// 27.10.2020 функция модифицирована (eps уменьшено) для полигонов повернутых относительно 
// координатных линий сетки под углом, полигонов
// в которых выделяется тепловая мощность. Тестирование произведено, точность удовлетворительна.
doublereal shorter_length_for_simplificationX(doublereal g, BLOCK* &b, integer lb, 
	WALL* &w, integer lw, SOURCE* &s, integer ls) {

	if (bBlasiusIM) {
		// Закоментировал для задачи Блазиуса. Потом требуется раскоментировать.

		return 1.0e-30;
	}
	else {


		if (!bFULL_AUTOMATIC) {
			return shorter_length_for_simplificationX_BASIC;
		}
		else {

			doublereal db = 1.0e30;


			//integer key = static_cast<integer>(isize_shorter_hash*((g - b[0].g.xS) / (1.05*(b[0].g.xE - b[0].g.xS))));
			integer key = ihash_key_shorter(g, b[0].g.xS, b[0].g.xE);

			if (key >= isize_shorter_hash) {
				printf("error in X : shorter_length_for_simplificationX CAD obj outside the cabinet...\n");
				std::cout << "g= " << g << " b[0].g.xS=" << b[0].g.xS << "  b[0].g.xE=" << b[0].g.xE << " key= " << key << std::endl;
				system("PAUSE");
				exit(1);
			}

			if (!bshorter_hash_X[key]) {

				for (integer i = 0; i < lb; ++i) {
					if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == POLYGON)) {
						doublereal dmult = dbody_R_out_multiplyer;
						if (((5.0 * (b[i].g.xE - b[i].g.xS)) < (b[i].g.yE - b[i].g.yS)) && ((5.0 * (b[i].g.xE - b[i].g.xS)) < (b[i].g.zE - b[i].g.zS))) {
							// Ситуация слой термопасты или клея.
							dmult = dbody_dist_multiplyer;
						}
						if (((5.0 * (b[i].g.xE - b[i].g.xS)) < (b[i].g.yE - b[i].g.yS)) || ((5.0 * (b[i].g.xE - b[i].g.xS)) < (b[i].g.zE - b[i].g.zS))) {
							// Ситуация палец затвора.
							dmult = dbody_LGATE_situation_multiplyer;
						}
						if ((g >= b[i].g.xS - dmult * (b[i].g.xE - b[i].g.xS)) && (g <= b[i].g.xE + dmult * (b[i].g.xE - b[i].g.xS))) {


							if ((fabs(dmult - dbody_R_out_multiplyer) < 0.1) && (b[i].g.itypegeom == POLYGON) && (b[i].arr_Sc[0] > 1.0e-20)) {

								// Синопсис: Это источник тепла в полигоне и полигон не яляется копией призмы преобразованной в полигон которая паралельна осям координат.
								// Это сложная полигональная геометрическая форма расположенная под углом к сеточным линиям к тому же выделяющая тепловую мощность.

								doublereal minimum_edge_length = fabs((b[i].g.xE - b[i].g.xS));
								if ((b[i].g.iPlane_obj2 == XY_PLANE)) {
									for (integer i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
										doublereal dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
									{
										doublereal dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
								}
								if ((b[i].g.iPlane_obj2 == XZ_PLANE)) {
									for (integer i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
										doublereal dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]) * (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
									{
										doublereal dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
								}
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval) * (minimum_edge_length)) < db) db = ((1.0 / rdivision_interval) * (minimum_edge_length));
							}
							else {
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval) * (b[i].g.xE - b[i].g.xS)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.xE - b[i].g.xS));
							}
						}
					}
					if (b[i].g.itypegeom == CYLINDER) {
						switch (b[i].g.iPlane) {
						case XY_PLANE: case XZ_PLANE:
							if (b[i].g.R_in_cyl < 1.0e-36) {
								if ((g >= b[i].g.xC - dbody_R_out_multiplyer * b[i].g.R_out_cyl) && (g <= b[i].g.xC + dbody_R_out_multiplyer * b[i].g.R_out_cyl)) {
									// g в зоне влияния блока с номером i.
									if (((1.0 / rdivision_interval) * (b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.R_out_cyl));
								}
							}
							else {
								if ((g >= b[i].g.xC - dbody_R_out_multiplyer * b[i].g.R_out_cyl) && (g <= b[i].g.xC + dbody_R_out_multiplyer * b[i].g.R_out_cyl)) {
									// g в зоне влияния блока с номером i.
									if (((1.0 / rdivision_interval) * (b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.R_in_cyl));
									// g в зоне влияния блока с номером i.
									if (((1.0 / rdivision_interval) * (b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.R_out_cyl - b[i].g.R_in_cyl));
								}
							}
							break;
						case YZ_PLANE:
							if ((g >= b[i].g.xC - dbody_Hcyl_multiplyer * b[i].g.Hcyl) && (g <= b[i].g.xC + b[i].g.Hcyl + dbody_Hcyl_multiplyer * b[i].g.Hcyl)) {
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval) * (b[i].g.Hcyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.Hcyl));
							}
							break;
						}
					}
					if (b[i].g.itypegeom == CAD_STL) {

						doublereal dmult = dbody_R_out_multiplyer;

						if (b[i].g.inumber_triangles_for_CAD_STL_model > 140) {
							dmult = 0.0001;
						}

						TOCHKA_FLOAT pmin, pmax;

						b[i].g.minimaxCAD_STL(pmin, pmax);

						if ((g >= pmin.x - dmult * (pmax.x - pmin.x)) && (g <= pmax.x + dmult * (pmax.x - pmin.x))) {

							TOCHKA_FLOAT pmin, pmax;

							b[i].g.minimaxCAD_STL(pmin, pmax);

							doublereal dist_loc = fmax((1.0 / 200.0) * (pmax.x - pmin.x), b[i].g.min_size_edge());

							// Поправка на очень тонкие слои.
							if (((1.0 / 200.0) * (pmax.x - pmin.x)) < 0.2 * b[i].g.min_size_edge()) {
								dist_loc = ((1.0 / rdivision_intervalCAD) * b[i].g.min_size_edge());
							}

							if (dist_loc < db) db = dist_loc;
						}

					}
				}
				if (lb == 1) {
					db = 1.0e-10;
				}
				else {
					// 1/200
					if (db > dcabinet_delta_multiplyer * (fabs(b[0].g.xE - b[0].g.xS))) db = dcabinet_delta_multiplyer * (fabs(b[0].g.xE - b[0].g.xS));
				}
				for (integer i = 0; i < lw; ++i) {
					doublereal dmult = dbody_R_out_multiplyer;
					switch (w[i].g.iPlane) {
					case XY_PLANE: case XZ_PLANE:
						if ((g >= w[i].g.xS - dmult * (w[i].g.xE - w[i].g.xS)) && (g <= w[i].g.xE + dmult * (w[i].g.xE - w[i].g.xS))) {
							// g в зоне влияния стенки с номером i.
							if (((1.0 / rdivision_interval) * (w[i].g.xE - w[i].g.xS)) < db) db = ((1.0 / rdivision_interval) * (w[i].g.xE - w[i].g.xS));
						}
						break;
					}
				}

				for (integer i = 0; i < ls; ++i) {
					doublereal dmult = dbody_R_out_multiplyer;
					switch (s[i].g.iPlane) {
					case XY_PLANE: case XZ_PLANE:
						if ((g >= s[i].g.xS - dmult * (s[i].g.xE - s[i].g.xS)) && (g <= s[i].g.xE + dmult * (s[i].g.xE - s[i].g.xS))) {
							// g в зоне влияния источника с номером i.
							if (((1.0 / rdivision_interval) * (s[i].g.xE - s[i].g.xS)) < db) db = ((1.0 / rdivision_interval) * (s[i].g.xE - s[i].g.xS));
						}
						break;
					}
				}

				bshorter_hash_X[key] = true;
				shorter_hash_X[key] = db;
			}
			else {
				db = shorter_hash_X[key];
			}

			return db;
		}
	}
}

// Здесь будет управление допуском на местах.
// 27.10.2020 функция модифицирована (eps уменьшено) для полигонов повернутых относительно 
// координатных линий сетки под углом, полигонов
// в которых выделяется тепловая мощность. Тестирование произведено, точность удовлетворительна.
doublereal shorter_length_for_simplificationY(doublereal g, BLOCK* &b, integer lb,
	WALL* &w, integer lw, SOURCE* &s, integer ls) {

	if (bBlasiusIM) {
		// Закоментировал для задачи Блазиуса. Потом требуется раскоментировать.

		return 1.0e-30;
	}
	else {

		if (!bFULL_AUTOMATIC) {
			return shorter_length_for_simplificationY_BASIC;
		}
		else {
			doublereal db = 1.0e30;

			//integer key = static_cast<integer>(isize_shorter_hash*((g - b[0].g.yS) / (1.05*(b[0].g.yE - b[0].g.yS))));
			integer key = ihash_key_shorter(g, b[0].g.yS, b[0].g.yE);

			if (key >= isize_shorter_hash) {
				printf("error in Y : shorter_length_for_simplificationY CAD obj outside the cabinet...\n");
				system("PAUSE");
				exit(1);
			}

			if (!bshorter_hash_Y[key]) {



				for (integer i = 0; i < lb; ++i) {
					if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == POLYGON)) {
						doublereal dmult = dbody_R_out_multiplyer;
						if (((5.0 * (b[i].g.yE - b[i].g.yS)) < (b[i].g.xE - b[i].g.xS)) && ((5.0 * (b[i].g.yE - b[i].g.yS)) < (b[i].g.zE - b[i].g.zS))) {
							// Ситуация слой термопасты или клея.
							dmult = dbody_dist_multiplyer;
						}
						if (((5.0 * (b[i].g.yE - b[i].g.yS)) < (b[i].g.xE - b[i].g.xS)) || ((5.0 * (b[i].g.yE - b[i].g.yS)) < (b[i].g.zE - b[i].g.zS))) {
							// Ситуация палец затвора.
							dmult = dbody_LGATE_situation_multiplyer;
						}

						if ((g >= b[i].g.yS - dmult * (b[i].g.yE - b[i].g.yS)) && (g <= b[i].g.yE + dmult * (b[i].g.yE - b[i].g.yS))) {


							if ((fabs(dmult - dbody_R_out_multiplyer) < 0.1) && (b[i].g.itypegeom == POLYGON) && (b[i].arr_Sc[0] > 1.0e-20)) {

								// Синопсис: Это источник тепла в полигоне и полигон не яляется копией призмы преобразованной в полигон которая паралельна осям координат.
								// Это сложная полигональная геометрическая форма расположенная под углом к сеточным линиям к тому же выделяющая тепловую мощность.

								doublereal minimum_edge_length = fabs((b[i].g.yE - b[i].g.yS));
								if ((b[i].g.iPlane_obj2 == XY_PLANE)) {
									for (integer i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
										doublereal dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
									{
										doublereal dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
								}
								if ((b[i].g.iPlane_obj2 == YZ_PLANE)) {
									for (integer i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
										doublereal dist = sqrt((b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]) * (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
									{
										doublereal dist = sqrt((b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
										if (dist < minimum_edge_length) {
											minimum_edge_length = dist;
										}
									}
								}
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval) * (minimum_edge_length)) < db) db = ((1.0 / rdivision_interval) * (minimum_edge_length));
							}
							else {
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval) * (b[i].g.yE - b[i].g.yS)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.yE - b[i].g.yS));
							}
						}
					}
					if (b[i].g.itypegeom == CYLINDER) {
						switch (b[i].g.iPlane) {
						case XY_PLANE: case YZ_PLANE:
							if (b[i].g.R_in_cyl < 1.0e-36) {
								if ((g >= b[i].g.yC - dbody_R_out_multiplyer * b[i].g.R_out_cyl) && (g <= b[i].g.yC + dbody_R_out_multiplyer * b[i].g.R_out_cyl)) {
									// g в зоне влияния блока с номером i.
									if (((1.0 / rdivision_interval) * (b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.R_out_cyl));
								}
							}
							else {
								if ((g >= b[i].g.yC - dbody_R_out_multiplyer * b[i].g.R_out_cyl) && (g <= b[i].g.yC + dbody_R_out_multiplyer * b[i].g.R_out_cyl)) {
									// g в зоне влияния блока с номером i.
									if (((1.0 / rdivision_interval) * (b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.R_in_cyl));
									// g в зоне влияния блока с номером i.
									if (((1.0 / rdivision_interval) * (b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.R_out_cyl - b[i].g.R_in_cyl));
								}
							}
							break;
						case XZ_PLANE:
							if ((g >= b[i].g.yC - dbody_Hcyl_multiplyer * b[i].g.Hcyl) && (g <= b[i].g.yC + b[i].g.Hcyl + dbody_Hcyl_multiplyer * b[i].g.Hcyl)) {
								// g в зоне влияния блока с номером i.
								if (((1.0 / rdivision_interval) * (b[i].g.Hcyl)) < db) db = ((1.0 / rdivision_interval) * (b[i].g.Hcyl));
							}
							break;
						}
					}
					if (b[i].g.itypegeom == CAD_STL) {

						TOCHKA_FLOAT pmin, pmax;
						doublereal dmult = dbody_R_out_multiplyer;
						if (b[i].g.inumber_triangles_for_CAD_STL_model > 140) {
							dmult = 0.0001;
						}

						b[i].g.minimaxCAD_STL(pmin, pmax);

						if ((g >= pmin.y - dmult * (pmax.y - pmin.y)) && (g <= pmax.y + dmult * (pmax.y - pmin.y))) {

							TOCHKA_FLOAT pmin, pmax;

							b[i].g.minimaxCAD_STL(pmin, pmax);

							doublereal dist_loc = fmax((1.0 / 200.0) * (pmax.y - pmin.y), b[i].g.min_size_edge());

							// Поправка на очень тонкие слои.
							if (((1.0 / 200.0) * (pmax.y - pmin.y)) < 0.2 * b[i].g.min_size_edge()) {
								dist_loc = ((1.0 / rdivision_intervalCAD) * b[i].g.min_size_edge());
							}

							if (dist_loc < db) db = dist_loc;
						}

					}
				}
				if (lb == 1) {
					db = 1.0e-10;
				}
				else {
					// 1/200				
					if (db > dcabinet_delta_multiplyer * (fabs(b[0].g.yE - b[0].g.yS))) db = dcabinet_delta_multiplyer * (fabs(b[0].g.yE - b[0].g.yS));
				}

				for (integer i = 0; i < lw; ++i) {
					doublereal dmult = dbody_R_out_multiplyer;
					switch (w[i].g.iPlane) {
					case XY_PLANE: case YZ_PLANE:
						if ((g >= w[i].g.yS - dmult * (w[i].g.yE - w[i].g.yS)) && (g <= w[i].g.yE + dmult * (w[i].g.yE - w[i].g.yS))) {
							// g в зоне влияния стенки с номером i.
							if (((1.0 / rdivision_interval) * (w[i].g.yE - w[i].g.yS)) < db) db = ((1.0 / rdivision_interval) * (w[i].g.yE - w[i].g.yS));
						}
						break;
					}
				}

				for (integer i = 0; i < ls; ++i) {
					doublereal dmult = dbody_R_out_multiplyer;
					switch (s[i].g.iPlane) {
					case XY_PLANE: case YZ_PLANE:
						if ((g >= s[i].g.yS - dmult * (s[i].g.yE - s[i].g.yS)) && (g <= s[i].g.yE + dmult * (s[i].g.yE - s[i].g.yS))) {
							// g в зоне влияния источника с номером i.
							if (((1.0 / rdivision_interval) * (s[i].g.yE - s[i].g.yS)) < db) db = ((1.0 / rdivision_interval) * (s[i].g.yE - s[i].g.yS));
						}
						break;
					}
				}
				bshorter_hash_Y[key] = true;
				shorter_hash_Y[key] = db;
			}
			else {
				db = shorter_hash_Y[key];
			}
			return db;
		}
	}
}

// Здесь будет управление допуском на местах.
// 27.10.2020 функция модифицирована (eps уменьшено) для полигонов повернутых относительно 
// координатных линий сетки под углом, полигонов
// в которых выделяется тепловая мощность. Тестирование произведено, точность удовлетворительна.
doublereal shorter_length_for_simplificationZ(doublereal g, BLOCK* &b, integer lb,
	WALL* &w, integer lw, SOURCE* &s, integer ls) {

	if (!bFULL_AUTOMATIC) {
		return shorter_length_for_simplificationZ_BASIC;
	}
	else {
		doublereal db = 1.0e30;

		//integer key = static_cast<integer>(isize_shorter_hash*((g - b[0].g.zS) / (1.05*(b[0].g.zE - b[0].g.zS))));
		integer key = ihash_key_shorter(g, b[0].g.zS, b[0].g.zE);

		if (key >= isize_shorter_hash) {
			printf("error in Z : shorter_length_for_simplificationZ CAD obj outside the cabinet...\n");
			system("PAUSE");
			exit(1);
		}

		if (!bshorter_hash_Z[key]) {



		for (integer i = 0; i < lb; ++i) {
			if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == POLYGON)) {

				doublereal dmult = dbody_R_out_multiplyer;
				if (((5.0*(b[i].g.zE - b[i].g.zS)) < (b[i].g.yE - b[i].g.yS)) && ((5.0*(b[i].g.zE - b[i].g.zS)) < (b[i].g.xE - b[i].g.xS))) {
					// Ситуация слой термопасты или клея.
					dmult = dbody_dist_multiplyer;
				}
				if (((5.0*(b[i].g.zE - b[i].g.zS)) < (b[i].g.yE - b[i].g.yS)) || ((5.0*(b[i].g.zE - b[i].g.zS)) < (b[i].g.xE - b[i].g.xS))) {
					// Ситуация палец затвора.
					dmult = dbody_LGATE_situation_multiplyer;
				}


				if ((g >= b[i].g.zS - dmult*(b[i].g.zE - b[i].g.zS)) && (g <= b[i].g.zE + dmult*(b[i].g.zE - b[i].g.zS))) {

					if ((fabs(dmult - dbody_R_out_multiplyer)<0.1) && (b[i].g.itypegeom == POLYGON) && (b[i].arr_Sc[0]>1.0e-20)) {

						// Синопсис: Это источник тепла в полигоне и полигон не яляется копией призмы преобразованной в полигон которая паралельна осям координат.
						// Это сложная полигональная геометрическая форма расположенная под углом к сеточным линиям к тому же выделяющая тепловую мощность.

						doublereal minimum_edge_length = fabs((b[i].g.zE - b[i].g.zS));
						if ((b[i].g.iPlane_obj2 == XZ_PLANE)) {
							for (integer i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
								doublereal dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7])*(b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7])*(b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
								if (dist < minimum_edge_length) {
									minimum_edge_length = dist;
								}
							}
							{
								doublereal dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
								if (dist < minimum_edge_length) {
									minimum_edge_length = dist;
								}
							}
						}
						if ((b[i].g.iPlane_obj2 == YZ_PLANE)) {
							for (integer i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
								doublereal dist = sqrt((b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7])*(b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7])*(b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
								if (dist < minimum_edge_length) {
									minimum_edge_length = dist;
								}
							}
							{
								doublereal dist = sqrt((b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
								if (dist < minimum_edge_length) {
									minimum_edge_length = dist;
								}
							}
						}
						// g в зоне влияния блока с номером i.
						if (((1.0 / rdivision_interval)*(minimum_edge_length)) < db) db = ((1.0 / rdivision_interval)*(minimum_edge_length));
					}
					else {
						// g в зоне влияния блока с номером i.
						if (((1.0 / rdivision_interval)*(b[i].g.zE - b[i].g.zS)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.zE - b[i].g.zS));
					}
				}
			}
			if (b[i].g.itypegeom == CYLINDER) {
				switch (b[i].g.iPlane) {
				case XZ_PLANE: case YZ_PLANE:
					if (b[i].g.R_in_cyl < 1.0e-36) {
						if ((g >= b[i].g.zC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.zC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
							// g в зоне влияния блока с номером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl));
						}
					}
					else {
						if ((g >= b[i].g.zC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.zC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
							// g в зоне влияния блока с номером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_in_cyl));
							// g в зоне влияния блока с номером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl));
						}
					}
					break;
				case XY_PLANE:
					if ((g >= b[i].g.zC- dbody_Hcyl_multiplyer *b[i].g.Hcyl) && (g <= b[i].g.zC + b[i].g.Hcyl+ dbody_Hcyl_multiplyer *b[i].g.Hcyl)) {
						// g в зоне влияния блока с номером i.
						if (((1.0 / rdivision_interval)*(b[i].g.Hcyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.Hcyl));
					}
					break;
				}
			}
			if (b[i].g.itypegeom == CAD_STL) {

				doublereal dmult = dbody_R_out_multiplyer;

				if (b[i].g.inumber_triangles_for_CAD_STL_model > 140) {
					dmult = 0.0001;
				}

				TOCHKA_FLOAT pmin, pmax;

				b[i].g.minimaxCAD_STL(pmin, pmax);

				if ((g >= pmin.z - dmult * (pmax.z - pmin.z)) && (g <= pmax.z + dmult * (pmax.z - pmin.z))) {

					TOCHKA_FLOAT pmin, pmax;

					b[i].g.minimaxCAD_STL(pmin, pmax);

					doublereal dist_loc = fmax((1.0 / 200.0) * (pmax.z - pmin.z), b[i].g.min_size_edge());

					// Поправка на очень тонкие слои.
					if (((1.0 / 200.0) * (pmax.z - pmin.z)) < 0.2 * b[i].g.min_size_edge()) {
						dist_loc = ((1.0 / rdivision_intervalCAD) * b[i].g.min_size_edge());
					}

					if (dist_loc < db) db = dist_loc;
				}

			}
		}
		if (lb == 1) {
			db = 1.0e-10;
		}
		else {
			// 1/200
			if (db > dcabinet_delta_multiplyer *(fabs(b[0].g.zE - b[0].g.zS))) db = dcabinet_delta_multiplyer *(fabs(b[0].g.zE - b[0].g.zS));
		}

		for (integer i = 0; i < lw; ++i) {
			doublereal dmult = dbody_R_out_multiplyer;
			switch (w[i].g.iPlane) {				
			case YZ_PLANE: case XZ_PLANE:
				if ((g >= w[i].g.zS - dmult * (w[i].g.zE - w[i].g.zS)) && (g <= w[i].g.zE + dmult * (w[i].g.zE - w[i].g.zS))) {
					// g в зоне влияния стенки с номером i.
					if (((1.0 / rdivision_interval)*(w[i].g.zE - w[i].g.zS)) < db) db = ((1.0 / rdivision_interval)*(w[i].g.zE - w[i].g.zS));
				}
				break;
			}
		}

		for (integer i = 0; i < ls; ++i) {
			doublereal dmult = dbody_R_out_multiplyer;
			switch (s[i].g.iPlane) {				
			case YZ_PLANE: case XZ_PLANE:
				if ((g >= s[i].g.zS - dmult * (s[i].g.zE - s[i].g.zS)) && (g <= s[i].g.zE + dmult * (s[i].g.zE - s[i].g.zS))) {
					// g в зоне влияния источника с номером i.
					if (((1.0 / rdivision_interval)*(s[i].g.zE - s[i].g.zS)) < db) db = ((1.0 / rdivision_interval)*(s[i].g.zE - s[i].g.zS));
				}
				break;
			}
		}

		bshorter_hash_Z[key] = true;
		shorter_hash_Z[key] = db;
		}
		else {
			db = shorter_hash_Z[key];
		}

		return db;
	}
}


// вычисляет 30 view факторов внутри вакуумного промежутка.
// 3 августа 2016.
void calculate_view_factors(BLOCK &b)
{
	doublereal x34, y34, z34, s34_W, s34_E, s34_S, s34_N, s34_B, s34_T;
	// min X (W)
	x34 = fabs(b.g.yE - b.g.yS); // x
	y34 = fabs(b.g.zE - b.g.zS); // y
	z34 = fabs(b.g.xE - b.g.xS); // z
	b.radiation.FWE = view_factor_aligned_parallel_rectangles(x34, y34, z34);
	b.radiation.FWS = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FWN = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FWB = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FWT = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	s34_W = b.radiation.FWE + 2.0*(b.radiation.FWS + b.radiation.FWB); // проверочная сумма должна быть равна 1.0.
	// max X (E)
	x34 = fabs(b.g.yE - b.g.yS);
	y34 = fabs(b.g.zE - b.g.zS);
	z34 = fabs(b.g.xE - b.g.xS);
	b.radiation.FEW = view_factor_aligned_parallel_rectangles(y34,x34,z34);
	b.radiation.FES = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FEN = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FEB = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FET = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	s34_E = b.radiation.FEW + 2.0*(b.radiation.FES + b.radiation.FEB); // проверочная сумма должна быть равна 1.0.
	// min Y (S)
	z34 = fabs(b.g.yE - b.g.yS);
	x34 = fabs(b.g.zE - b.g.zS);
	y34 = fabs(b.g.xE - b.g.xS);
	b.radiation.FSW = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FSE = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FSN = view_factor_aligned_parallel_rectangles(x34,y34,z34);
	b.radiation.FSB = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FST = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	s34_S = b.radiation.FSN + 2.0*(b.radiation.FSW + b.radiation.FSB);
    // max Y (N)
	z34 = fabs(b.g.yE - b.g.yS);
	x34 = fabs(b.g.zE - b.g.zS);
	y34 = fabs(b.g.xE - b.g.xS);
	b.radiation.FNW = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FNE = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FNS = view_factor_aligned_parallel_rectangles(x34,y34,z34);
	b.radiation.FNB = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FNT = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	s34_N = b.radiation.FNS + 2.0*(b.radiation.FNW + b.radiation.FNB);
	// min Z (B)
	x34 = fabs(b.g.yE - b.g.yS);
	z34 = fabs(b.g.zE - b.g.zS);
	y34 = fabs(b.g.xE - b.g.xS);
	b.radiation.FBW = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FBE = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FBS = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FBN = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FBT = view_factor_aligned_parallel_rectangles(x34,y34,z34);
	s34_B = b.radiation.FBT + 2.0*(b.radiation.FBW + b.radiation.FBS);
	// max Z (T)
	x34 = fabs(b.g.yE - b.g.yS);
	z34 = fabs(b.g.zE - b.g.zS);
	y34 = fabs(b.g.xE - b.g.xS);
	b.radiation.FTW = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FTE = view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
	b.radiation.FTS = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FTN = view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
	b.radiation.FTB = view_factor_aligned_parallel_rectangles(x34,y34,z34);
	s34_T = b.radiation.FTB + 2.0*(b.radiation.FTW + b.radiation.FTS);
	// log:
	printf("view factors:\n");
	printf("FWE=%1.3f FWS=%1.3f FWN=%1.3f FWB=%1.3f FWT=%1.3f sum=%1.3f\n", b.radiation.FWE, b.radiation.FWS, b.radiation.FWN, b.radiation.FWB, b.radiation.FWT, s34_W);
	printf("FEW=%1.3f FES=%1.3f FEN=%1.3f FEB=%1.3f FET=%1.3f sum=%1.3f\n", b.radiation.FEW, b.radiation.FES, b.radiation.FEN, b.radiation.FEB, b.radiation.FET, s34_E);
	printf("FSW=%1.3f FSE=%1.3f FSN=%1.3f FSB=%1.3f FST=%1.3f sum=%1.3f\n", b.radiation.FSW, b.radiation.FSE, b.radiation.FSN, b.radiation.FSB, b.radiation.FST, s34_S);
	printf("FNW=%1.3f FNE=%1.3f FNS=%1.3f FNB=%1.3f FNT=%1.3f sum=%1.3f\n", b.radiation.FNW, b.radiation.FNE, b.radiation.FNS, b.radiation.FNB, b.radiation.FNT, s34_N);
	printf("FBW=%1.3f FBE=%1.3f FBS=%1.3f FBN=%1.3f FBT=%1.3f sum=%1.3f\n", b.radiation.FBW, b.radiation.FBE, b.radiation.FBS, b.radiation.FBN, b.radiation.FBT, s34_B);
	printf("FTW=%1.3f FTE=%1.3f FTS=%1.3f FTN=%1.3f FTB=%1.3f sum=%1.3f\n", b.radiation.FTW, b.radiation.FTE, b.radiation.FTS, b.radiation.FTN, b.radiation.FTB, s34_T);
	//system("PAUSE");
	
} // calculate_view_factors


// Возвращает значение теплопроводности при температуре t_C в градусах.
float get_lam(integer nsize_lam, float* &temp_lam, float* &arr_lam, doublereal t_C) {
	if (nsize_lam == 1) {
		return arr_lam[0];
	}
	else if (nsize_lam>1) {
		if (t_C < temp_lam[0]) {
			return arr_lam[0];
		}
		else if (t_C > temp_lam[nsize_lam - 1]) {
			return arr_lam[nsize_lam - 1];
		}
		else {
			integer if_3 = -1;
			for (integer i_sc = 0; i_sc < nsize_lam - 1; i_sc++) {
				if ((t_C >= temp_lam[i_sc]) && (t_C <= temp_lam[i_sc + 1])) {
					if_3 = i_sc;
					break;
				}
			}
			if (if_3 > -1) {
				doublereal af = (arr_lam[if_3+1]-arr_lam[if_3]) / (temp_lam[if_3 + 1] - temp_lam[if_3]);
				doublereal bf = (arr_lam[if_3] - af*temp_lam[if_3]);
				return static_cast <float>((af*t_C+bf));
			}
			else {
				printf("error in get_lam function if_3==-1\n");
				system("pause");
				exit(1);
			}
		}
	}
	else {
		printf("error in get_lam function\n");
		system("pause");
		exit(1);
	}
	return 0.0;
} // get_lam

// 15.08.2020
// Возвращает значение коэффициента линейного теплового расширения при температуре t_C в градусах.
float get_beta_t_solid(integer nsize_beta_t_solid, float* &temp_beta_t_solid, float* &arr_beta_t_solid, doublereal t_C) {
	if (nsize_beta_t_solid == 1) {
		return arr_beta_t_solid[0];
	}
	else if (nsize_beta_t_solid > 1) {
		if (t_C < temp_beta_t_solid[0]) {
			return arr_beta_t_solid[0];
		}
		else if (t_C > temp_beta_t_solid[nsize_beta_t_solid - 1]) {
			return arr_beta_t_solid[nsize_beta_t_solid - 1];
		}
		else {
			integer if_3 = -1;
			for (integer i_sc = 0; i_sc < nsize_beta_t_solid - 1; i_sc++) {
				if ((t_C >= temp_beta_t_solid[i_sc]) && (t_C <= temp_beta_t_solid[i_sc + 1])) {
					if_3 = i_sc;
					break;
				}
			}
			if (if_3 > -1) {
				doublereal af = (arr_beta_t_solid[if_3 + 1] - arr_beta_t_solid[if_3]) / (temp_beta_t_solid[if_3 + 1] - temp_beta_t_solid[if_3]);
				doublereal bf = (arr_beta_t_solid[if_3] - af * temp_beta_t_solid[if_3]);
				return static_cast <float>((af*t_C + bf));
			}
			else {
				printf("error in get_beta_t_solid function if_3==-1\n");
				system("pause");
				exit(1);
			}
		}
	}
	else {
		printf("error in get_beta_t_solid function\n");
		system("pause");
		exit(1);
	}
	return 0.0;
} // get_beta_t_solid


// 23.08.2020
// Возвращает значение модуля Юнга при температуре t_C в градусах.
float get_Young_Module(integer n_YoungModule, float* &temp_Young_Module, float* &arr_Young_Module, doublereal t_C) {
	if (n_YoungModule == 1) {
		return arr_Young_Module[0];
	}
	else if (n_YoungModule > 1) {
		if (t_C < temp_Young_Module[0]) {
			return arr_Young_Module[0];
		}
		else if (t_C > temp_Young_Module[n_YoungModule - 1]) {
			return arr_Young_Module[n_YoungModule - 1];
		}
		else {
			integer if_3 = -1;
			for (integer i_sc = 0; i_sc < n_YoungModule - 1; i_sc++) {
				if ((t_C >= temp_Young_Module[i_sc]) && (t_C <= temp_Young_Module[i_sc + 1])) {
					if_3 = i_sc;
					break;
				}
			}
			if (if_3 > -1) {
				doublereal af = (arr_Young_Module[if_3 + 1] - arr_Young_Module[if_3]) / (temp_Young_Module[if_3 + 1] - temp_Young_Module[if_3]);
				doublereal bf = (arr_Young_Module[if_3] - af * temp_Young_Module[if_3]);
				return static_cast <float>((af*t_C + bf));
			}
			else {
				printf("error in get_Young_Module function if_3==-1\n");
				system("pause");
				exit(1);
			}
		}
	}
	else {
		printf("error in get_Young_Module function\n");
		system("pause");
		exit(1);
	}
	return 0.0;
} // get_Young_Module

// 23.08.2020
// Возвращает значение коэффициента Пуассона при температуре t_C в градусах.
float get_Poisson_ratio(integer n_Poisson_ratio, float* &temp_Poisson_ratio, float* &arr_Poisson_ratio, doublereal t_C) {
	if (n_Poisson_ratio == 1) {
		return arr_Poisson_ratio[0];
	}
	else if (n_Poisson_ratio > 1) {
		if (t_C < temp_Poisson_ratio[0]) {
			return arr_Poisson_ratio[0];
		}
		else if (t_C > temp_Poisson_ratio[n_Poisson_ratio - 1]) {
			return arr_Poisson_ratio[n_Poisson_ratio - 1];
		}
		else {
			integer if_3 = -1;
			for (integer i_sc = 0; i_sc < n_Poisson_ratio - 1; i_sc++) {
				if ((t_C >= temp_Poisson_ratio[i_sc]) && (t_C <= temp_Poisson_ratio[i_sc + 1])) {
					if_3 = i_sc;
					break;
				}
			}
			if (if_3 > -1) {
				doublereal af = (arr_Poisson_ratio[if_3 + 1] - arr_Poisson_ratio[if_3]) / (temp_Poisson_ratio[if_3 + 1] - temp_Poisson_ratio[if_3]);
				doublereal bf = (arr_Poisson_ratio[if_3] - af * temp_Poisson_ratio[if_3]);
				return static_cast <float>((af*t_C + bf));
			}
			else {
				printf("error in get_Poisson_ratio function if_3==-1\n");
				system("pause");
				exit(1);
			}
		}
	}
	else {
		printf("error in get_Poisson_ratio function\n");
		system("pause");
		exit(1);
	}
	return 0.0;
} // get_Poisson_ratio

  // Возвращает значение теплоёмкости при температуре t_C в градусах.
float get_cp(integer nsize_cp, float* &temp_cp, float* &arr_cp, doublereal t_C) {
	if (nsize_cp == 1) {
		return arr_cp[0];
	}
	else if (nsize_cp>1) {
		if (t_C < temp_cp[0]) {
			return arr_cp[0];
		}
		else if (t_C > temp_cp[nsize_cp - 1]) {
			return arr_cp[nsize_cp - 1];
		}
		else {
			integer if_3 = -1;
			for (integer i_sc = 0; i_sc < nsize_cp - 1; i_sc++) {
				if ((t_C >= temp_cp[i_sc]) && (t_C <= temp_cp[i_sc + 1])) {
					if_3 = i_sc;
					break;
				}
			}
			if (if_3 > -1) {
				doublereal af = (arr_cp[if_3 + 1] - arr_cp[if_3]) / (temp_cp[if_3 + 1] - temp_cp[if_3]);
				doublereal bf = (arr_cp[if_3] - af*temp_cp[if_3]);
					return static_cast <float>((af*t_C + bf));
			}
			else {
				printf("error in get_cp function if_3==-1\n");
				system("pause");
				exit(1);
			}
		}
	}
	else {
		printf("error in get_cp function\n");
		system("pause");
		exit(1);
	}
	return 0.0;
} // get_cp

// 19 november 2016.
// Возвращает значение объёмной плотности мощности тепловыделения при температуре t_C в градусах.
doublereal get_power(integer nsize_Sc, doublereal* &temp_Sc, doublereal* &arr_Sc, doublereal t_C) {
	if (nsize_Sc == 1) {
		return arr_Sc[0];
	}
	else if (nsize_Sc>1) {
		if (t_C < temp_Sc[0]) {
			return arr_Sc[0];
		}
		else if (t_C > temp_Sc[nsize_Sc - 1]) {
			return arr_Sc[nsize_Sc - 1];
		}
		else {
			integer if_3 = -1;
			for (integer i_sc = 0; i_sc < nsize_Sc - 1; i_sc++) {
				if ((t_C >= temp_Sc[i_sc]) && (t_C <= temp_Sc[i_sc + 1])) {
					if_3 = i_sc;
					break;
				}
			}
			if (if_3 > -1) {
				doublereal af = (arr_Sc[if_3 + 1] - arr_Sc[if_3]) / (temp_Sc[if_3 + 1] - temp_Sc[if_3]);
				doublereal bf = (arr_Sc[if_3] - af*temp_Sc[if_3]);
				return (af*t_C + bf);
			}
			else {
				printf("error in get_power function if_3==-1\n");
				system("pause");
				exit(1);
			}
		}
	}
	else {
		printf("error in get_power function\n");
		system("pause");
		exit(1);
	}
	return 0.0;
} // get_power

// Вычисляет плотность теплового потока излучения на гранях вакуумного промежутка
// по заданным значениям:
// 1. излучательных способностей на гранях промежутка: emissivity,
// 2. Усреднённым температурам на гранях промежутка.
// 3. массиву View Factors на гранях промежутка.
void calculation_density_radiation_heat_flux(BLOCK &b)
{
	// Только если внутри Prism Object включено излучение.
	if (b.radiation.binternalRadiation)
	{
		
		// инициализация:
		b.radiation.JW = 0.0;
		b.radiation.JE = 0.0;
		b.radiation.JS = 0.0;
		b.radiation.JN = 0.0;
		b.radiation.JB = 0.0;
		b.radiation.JT = 0.0;
		// диагональные коэффициенты.
		doublereal apW = 1.0 + ((1.0 - b.radiation.emissW) / b.radiation.emissW)*1.0;
		doublereal apE = 1.0 + ((1.0 - b.radiation.emissE) / b.radiation.emissE)*1.0;
		doublereal apS = 1.0 + ((1.0 - b.radiation.emissS) / b.radiation.emissS)*1.0;
		doublereal apN = 1.0 + ((1.0 - b.radiation.emissN) / b.radiation.emissN)*1.0;
		doublereal apB = 1.0 + ((1.0 - b.radiation.emissB) / b.radiation.emissB)*1.0;
		doublereal apT = 1.0 + ((1.0 - b.radiation.emissT) / b.radiation.emissT)*1.0;

		//doublereal residual = 1.0; // невязка.

		doublereal JW_new, JE_new, JS_new, JN_new, JB_new, JT_new;
		for (integer k = 0; k <= maxiter_Radiation_block; ++k)
		{
			// Возможно недоэтерированность это очень хорошо а то иначе мы получим глобальную расходимость.

			// ВНИМАНИЕ !!! b.radiation.Temp в Кельвинах.

			JW_new = alpha_Radiation_block*((STEFAN_BOLCMAN_CONST*b.radiation.TempW*b.radiation.TempW*b.radiation.TempW*b.radiation.TempW + ((1.0 - b.radiation.emissW) / b.radiation.emissW)*
				(b.radiation.FWE*b.radiation.JE + b.radiation.FWS*b.radiation.JS + b.radiation.FWN*b.radiation.JN + b.radiation.FWB*b.radiation.JB +
				b.radiation.FWT*b.radiation.JT)) / apW) + (1.0 - alpha_Radiation_block)*b.radiation.JW;

			JE_new = alpha_Radiation_block*((STEFAN_BOLCMAN_CONST*b.radiation.TempE*b.radiation.TempE*b.radiation.TempE*b.radiation.TempE + ((1.0 - b.radiation.emissE) / b.radiation.emissE)*
				(b.radiation.FEW*b.radiation.JW + b.radiation.FES*b.radiation.JS + b.radiation.FEN*b.radiation.JN + b.radiation.FEB*b.radiation.JB +
				b.radiation.FET*b.radiation.JT)) / apE) + (1.0 - alpha_Radiation_block)*b.radiation.JE;

			JS_new = alpha_Radiation_block*((STEFAN_BOLCMAN_CONST*b.radiation.TempS*b.radiation.TempS*b.radiation.TempS*b.radiation.TempS + ((1.0 - b.radiation.emissS) / b.radiation.emissS)*
				(b.radiation.FSW*b.radiation.JW + b.radiation.FSE*b.radiation.JE + b.radiation.FSN*b.radiation.JN + b.radiation.FSB*b.radiation.JB +
				b.radiation.FST*b.radiation.JT)) / apS) + (1.0 - alpha_Radiation_block)*b.radiation.JS;

			JN_new = alpha_Radiation_block*((STEFAN_BOLCMAN_CONST*b.radiation.TempN*b.radiation.TempN*b.radiation.TempN*b.radiation.TempN + ((1.0 - b.radiation.emissN) / b.radiation.emissN)*
				(b.radiation.FNW*b.radiation.JW + b.radiation.FNE*b.radiation.JE + b.radiation.FNS*b.radiation.JS + b.radiation.FNB*b.radiation.JB +
				b.radiation.FNT*b.radiation.JT)) / apN) + (1.0 - alpha_Radiation_block)*b.radiation.JN;

			JB_new = alpha_Radiation_block*((STEFAN_BOLCMAN_CONST*b.radiation.TempB*b.radiation.TempB*b.radiation.TempB*b.radiation.TempB + ((1.0 - b.radiation.emissB) / b.radiation.emissB)*
				(b.radiation.FBW*b.radiation.JW + b.radiation.FBE*b.radiation.JE + b.radiation.FBS*b.radiation.JS + b.radiation.FBN*b.radiation.JN +
				b.radiation.FBT*b.radiation.JT)) / apB) + (1.0 - alpha_Radiation_block)*b.radiation.JB;

			JT_new = alpha_Radiation_block*((STEFAN_BOLCMAN_CONST*b.radiation.TempT*b.radiation.TempT*b.radiation.TempT*b.radiation.TempT + ((1.0 - b.radiation.emissT) / b.radiation.emissT)*
				(b.radiation.FTW*b.radiation.JW + b.radiation.FTE*b.radiation.JE + b.radiation.FTS*b.radiation.JS + b.radiation.FTN*b.radiation.JN +
				b.radiation.FTB*b.radiation.JB)) / apT) + (1.0 - alpha_Radiation_block)*b.radiation.JT;

			// Вычисление невязки.

			// update:
			b.radiation.JW = JW_new;
			b.radiation.JE = JE_new;
			b.radiation.JS = JS_new;
			b.radiation.JN = JN_new;
			b.radiation.JB = JB_new;
			b.radiation.JT = JT_new;
		}

		if ((steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)&&
			((steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T)))
		{
			printf("tW=%e tE=%e tS=%e tN=%e tB=%e tT=%e\n", b.radiation.TempW-273.15, b.radiation.TempE - 273.15, b.radiation.TempS - 273.15, b.radiation.TempN - 273.15, b.radiation.TempB - 273.15, b.radiation.TempT - 273.15);
			printf("JW=%e JE=%e JS=%e JN=%e JB=%e JT=%e\n", b.radiation.JW, b.radiation.JE, b.radiation.JS, b.radiation.JN, b.radiation.JB, b.radiation.JT);
		}
	    //	system("PAUSE");		
		
	}
} // calculation_density_radiation_heat_flux




// Пересекаются ли два отрезка 02.08.2019
bool b_is_intersect(doublereal x1, doublereal y1,
	doublereal x2, doublereal y2, doublereal x3, doublereal y3, 
	doublereal x4, doublereal y4)
{
	bool bintersect = true;
	doublereal Ua, Ub, numerator_a, numerator_b, denominator;

	denominator = (y4 - y3)*(x1 - x2) - (x4 - x3)*(y1 - y2);
    
	if (fabs(denominator) < 1.0e-36) {
		if (fabs((x1*y2 - x2 * y1)*(x4 - x3) - (x3*y4 - x4 * y3)*(x2 - x1)) <= 1.0e-36
			&& fabs((x1*y2 - x2 * y1)*(y4 - y3) - (x3*y4 - x4 * y3)*(y2 - y1)) <= 1.0e-36)
		{
			
			if (sqrt((x4 - x1)*(x4 - x1) + (y4 - y1)*(y4 - y1)) > 1.0e-36) {
				// Четвертая точка не должна совпадать с первой. Она может совпадать если у нас всего пять точек и
				// первый и последний отрезок лежат на одной прямой и имеют общую точку. Нам нужно корректно 
				// обработать такой случай.

				if ((fabs(x1 - x2) < 1.0e-36) && (fabs(x1 - x3) < 1.0e-36) && (fabs(x1 - x4) < 1.0e-36)) {
					if (((y1 < y3) && (y2 < y3) && (y1 < y4) && (y2 < y4)) ||
						((y1 > y3) && (y2 > y3) && (y1 > y4) && (y2 > y4))) {
						// Лежат на одной прямой но отрезки не пересекаются.
						bintersect = false; // не пересекаются
					}
					else {
						printf("denominator=%e\n", denominator);
						printf("%e %e\n", fabs((x1*y2 - x2 * y1)*(x4 - x3) - (x3*y4 - x4 * y3)*(x2 - x1)),
							fabs((x1*y2 - x2 * y1)*(y4 - y3) - (x3*y4 - x4 * y3)*(y2 - y1)));
						printf("x1=%e y1=%e\n", x1, y1);
						printf("x2=%e y2=%e\n", x2, y2);
						printf("x3=%e y3=%e\n", x3, y3);
						printf("x4=%e y4=%e\n", x4, y4);
						//bintersect = true; // пересекаются
					}
				}
				else if ((fabs(y1 - y2) < 1.0e-36) && (fabs(y1 - y3) < 1.0e-36) && (fabs(y1 - y4) < 1.0e-36)) {
					if (((x1 < x3) && (x2 < x3) && (x1 < x4) && (x2 < x4)) ||
						((x1 > x3) && (x2 > x3) && (x1 > x4) && (x2 > x4))) {
						// Лежат на одной прямой но отрезки не пересекаются.
						bintersect = false; // не пересекаются
					}
					else {
						printf("denominator=%e\n", denominator);
						printf("%e %e\n", fabs((x1*y2 - x2 * y1)*(x4 - x3) - (x3*y4 - x4 * y3)*(x2 - x1)),
							fabs((x1*y2 - x2 * y1)*(y4 - y3) - (x3*y4 - x4 * y3)*(y2 - y1)));
						printf("x1=%e y1=%e\n", x1, y1);
						printf("x2=%e y2=%e\n", x2, y2);
						printf("x3=%e y3=%e\n", x3, y3);
						printf("x4=%e y4=%e\n", x4, y4);
						//bintersect = true; // пересекаются
					}
				}
				else {
					printf("denominator=%e\n", denominator);
					printf("%e %e\n", fabs((x1*y2 - x2 * y1)*(x4 - x3) - (x3*y4 - x4 * y3)*(x2 - x1)),
						fabs((x1*y2 - x2 * y1)*(y4 - y3) - (x3*y4 - x4 * y3)*(y2 - y1)));
					printf("x1=%e y1=%e\n", x1, y1);
					printf("x2=%e y2=%e\n", x2, y2);
					printf("x3=%e y3=%e\n", x3, y3);
					printf("x4=%e y4=%e\n", x4, y4);
					//bintersect = true; // пересекаются
				}
			}
			else {
				bintersect = false; // не пересекаются
			}
		}
		else {
			bintersect = false; // не пересекаются
		}
	}
	else {
		numerator_a = (x4 - x2)*(y4 - y3) - (x4 - x3)*(y4 - y2);
		numerator_b = (x1 - x2)*(y4 - y2) - (x4 - x2)*(y1 - y2);
		Ua = numerator_a / denominator;
		Ub = numerator_b / denominator;
		
		if (Ua > 1.0e-36 && Ua < 1.0-1.0e-36 && Ub > 1.0e-36 && Ub < 1.0-1.0e-36) {
			printf("denominator=%e\n", denominator);
			printf("Ua=%e Ub=%e\n", Ua, Ub);
			//bintersect = true;
		}
		else {
			bintersect = false;
		}
	}

	return bintersect;
}



// Параметры которые используются для 
// настройки модели Смагоринского.
typedef struct TSMAGORINSKYINFO {
	doublereal Cs; // константа Смагоринского.
	bool bDynamic_Stress; // Для определения константы Смагоринского Cs используется динамическая модель Германо.
	bool bLimiters_Cs; // ограничивать ли постоянную Смагоринского ?
	doublereal minCs, maxCs; // минимальное и максимальное значения константы Смагоринского.
	// 2-SIMPSON
	integer itypeFiltrGermano; // тип фильтра который используется для осреднения в модели Германо.
	// 10 micron
	doublereal roughness; // значение шероховатости на твёрдой неподвижной стенке.
						  // показатель степени для учёта шероховатости на стенке.
	integer ipowerroughness; // может принимать значения только 1 или 2.
	bool bfdelta; // использовать ли поправку дающую улучшение на неравномерной сетке ?
	bool bSmagorinsky_Lilly; // использовать ли модель Смагоринского-Лиллу ?
	bool bsurface_roughness; // использовать ли поправку учитывающую шероховатость стенки ?
	bool bSelectiveSmagorinsky; // использовать ли Selective Smagorinsky Model ?
	// 2-SIMPSON
	integer itypeFILTRSelectiveSmagorinsky; // тип фильтра который используется для осреднения в модели Selective Smagorinsky.
	doublereal SSangle; // угол между вихрем и осреднённым вихрем в модели Selective Smagorinsky;
	bool bRichardsonCorrect; // использовать ли поправку связанную с числом Ричардсона для течений с кривизной линий тока.
	doublereal rRichardsonMultiplyer; // коэффициент в поправочной формуле связанной с кривизной линий тока.
	TSMAGORINSKYINFO() {
		 Cs=0.151; // константа Смагоринского.
		bDynamic_Stress=false; // Для определения константы Смагоринского Cs используется динамическая модель Германо.
		bLimiters_Cs=false; // ограничивать ли постоянную Смагоринского ?
		minCs=-1.0e30; maxCs=1.0e30; // минимальное и максимальное значения константы Смагоринского.
		// 2-SIMPSON
		itypeFiltrGermano=2; // тип фильтра который используется для осреднения в модели Германо.
		// 10 micron
		roughness=10.0; // значение шероховатости на твёрдой неподвижной стенке.
						  // показатель степени для учёта шероховатости на стенке.
		ipowerroughness=1; // может принимать значения только 1 или 2.
		bfdelta=true; // использовать ли поправку дающую улучшение на неравномерной сетке ?
		bSmagorinsky_Lilly=false; // использовать ли модель Смагоринского-Лиллу ?
		bsurface_roughness=false; // использовать ли поправку учитывающую шероховатость стенки ?
		bSelectiveSmagorinsky=false; // использовать ли Selective Smagorinsky Model ?
		// 2-SIMPSON
		itypeFILTRSelectiveSmagorinsky=2; // тип фильтра который используется для осреднения в модели Selective Smagorinsky.
		SSangle=15.0; // угол между вихрем и осреднённым вихрем в модели Selective Smagorinsky;
		bRichardsonCorrect=false; // использовать ли поправку связанную с числом Ричардсона для течений с кривизной линий тока.
		rRichardsonMultiplyer=1.0; // коэффициент в поправочной формуле связанной с кривизной линий тока.
	}
} SMAGORINSKYINFO;

// Для полилинейного метода:
typedef struct TNODELR {
	integer id; // идентификатор внутреннего узла
	struct TNODELR *next; // ссылка на следующий узел или nullptr
	TNODELR() {
		id=-1; // идентификатор внутреннего узла
	    next = nullptr; // ссылка на следующий узел или nullptr
	}
} NODELR;

typedef struct TNODELR_BASE {
	integer ilineid; // идентификатор сеточной линии
	integer iN; // количество узловых точек включая граничные
	struct TNODELR *root; // корень новой сеточной линии
	struct TNODELR_BASE *next; // указатель на следующую сеточную линию или nullptr

	// Особый случай:
	// обработка плоского бесконечно 
	// тонкого источника тепла.
	// Переменные равны истине если связь с источником 
	// прервана (в случае если воздух граничит с плоским источником).
	bool bNeimanStart;
	bool bNeimanEnd;
	TNODELR_BASE() {
		ilineid=-1; // идентификатор сеточной линии
		iN=-1; // количество узловых точек включая граничные
		root = nullptr; // корень новой сеточной линии
		next = nullptr; // указатель на следующую сеточную линию или nullptr

		// Особый случай:
		// обработка плоского бесконечно 
		// тонкого источника тепла.
		// Переменные равны истине если связь с источником 
		// прервана (в случае если воздух граничит с плоским источником).
		bNeimanStart=false;
		bNeimanEnd=false;
	}
} NODELR_BASE;


typedef struct TBOUND {
	// граничный узел и узлы вглубь расчётной области
	int iB, iI, iII;
	int iI1, iI2;
	integer Norm; // внутренняя нормаль
				  // marker boundary:
				  /*
				  * О маркере границы:
				  * if (MCB < ls) значит источник с номером MCB;
				  * else if (MCB < ls+lw) значит стенка c номером MCB-ls;
				  * else значение по умолчанию. MCB==ls+lw.
				  */
	integer MCB; // marker boundary

				 // соседи на границе области,
				 // их позиции нужны в матрице:
				 // На этих позициях будут стоять нули.
				 // Если позиция отсутствует, то стоит -1.
				 // Всего возможно не более 4 позиций, но 
				 // здесь сделано больше из-за 
				 // алгоритмического удобства.
				 // Позиции соответствуют сторонам света.
	integer iW[6];

	// Для внутреннего источника указывает
	// на соседний внутренний узел с обратной стороны
	// источника противоположной текущей стороне, или на -1
	// если источник на границе hollow блока.
	// TODO 6 мая 2016.

	// для теплообмена излучением хранит emissivity.
	doublereal emissivity; // излучательная способность границы.

						   // Площадь грани.
						   // Это новое поле введённое лишь 20 сентября 2016.
	doublereal dS;

	// координаты центра грани.
	TOCHKA p_c;

	TBOUND() {
		// граничный узел и узлы вглубь расчётной области
		iB=-1; iI=-1; iII=-1;
		iI1=-1; iI2=-1;
		Norm=-10000; // внутренняя нормаль
		// marker boundary:
		/*
		* О маркере границы:
		* if (MCB < ls) значит источник с номером MCB;
		* else if (MCB < ls+lw) значит стенка c номером MCB-ls;
		* else значение по умолчанию. MCB==ls+lw.
		*/
		MCB=-1; // marker boundary

		 // соседи на границе области,
		 // их позиции нужны в матрице:
		 // На этих позициях будут стоять нули.
		 // Если позиция отсутствует, то стоит -1.
		 // Всего возможно не более 4 позиций, но 
		 // здесь сделано больше из-за 
		 // алгоритмического удобства.
		 // Позиции соответствуют сторонам света.
		iW[0] = -1;
		iW[1] = -1;
		iW[2] = -1;
		iW[3] = -1;
		iW[4] = -1;
		iW[5] = -1;
	

		// Для внутреннего источника указывает
		// на соседний внутренний узел с обратной стороны
		// источника противоположной текущей стороне, или на -1
		// если источник на границе hollow блока.
		// TODO 6 мая 2016.

		// для теплообмена излучением хранит emissivity.
		emissivity=0.8; // излучательная способность границы.

		// Площадь грани.
		// Это новое поле введённое лишь 20 сентября 2016.
		dS=0.0;

		// координаты центра грани.
		//TOCHKA p_c;
	}
} BOUND;

enum class FLOW_REGIME {LAMINAR=0, TURBULENT=1};
//const unsigned char LAMINAR = 0; // ламинарное течение
//const unsigned char ZEROEQMOD = 1; // турбулентное течение - алгебраическая модель турбулентности Zero Equation Model
//const unsigned char SMAGORINSKY = 2; // турбулентное течение - LES моделирование одна из разновидностей модели Смагоринского.
//const unsigned char RNG_LES = 3; // Based on Renormalization Group Theory. (модель соответствует описанию CFD-Wiki).
//const unsigned char RANS_SPALART_ALLMARES = 4; // RANS модель Спаларта - Аллмареса.
//const unsigned char RANS_MENTER_SST = 5; //  Shear Stress Model SST Ментера.
//const unsigned char RANS_STANDART_K_EPS = 6; // Двухслойная k-epsilon модель на основе стандартной k-epsilon модели.
//const unigned char RANS_LANGTRY_MENTOR_SST = 7; // Модель Ментора Лангтрии 2009.
// Динамическая модель Германо 1991 года. (основывается на модели Смагоринского и реализуется в виде её опции - bDynamic_Stress).
enum class VISCOSITY_MODEL { LAMINAR = 0,
	ZEROEQMOD = 1,
	SMAGORINSKY = 2, 
	RNG_LES = 3,
	RANS_SPALART_ALLMARES = 4, 
	RANS_MENTER_SST = 5,
	RANS_STANDART_K_EPS = 6,
    RANS_LANGTRY_MENTOR_SST = 7
};
enum class TURBULENT_MODEL { ZEROEQMOD = 0,
	SMAGORINSKY = 1,
	RNG_LES = 2,
	RANS_SPALART_ALLMARES = 3,
	RANS_MENTER_SST = 4, 
	RANS_STANDART_K_EPS = 5,
	RANS_LANGTRY_MENTOR_SST = 6
};

// информация о жидкой зоне
typedef struct TFLOWINFO {
	// опорная точка
	doublereal xc, yc, zc;
	// нужно ли рассчитывать поле течения
	integer iflow;
	// режим течения: 0 - ламинарный, 1 - турбулентный
	FLOW_REGIME iflowregime;
	// модель турбулентности
	// 0 - Zero Equation Turbulence Model (RANS).
	// 1 - модель Смагоринского (LES). (как опция модели, модель Германо 1991 (LES)).
	// 2 - RNG (LES).
	// 3 - Spalart-Allmares (RANS). 19.04.2019
	// 4 - K - Omega SST (RANS). 06.10.2019
	// 5 - Двухслойная модель на основе стандартной K-Epsilon модели (RANS). 31.10.2019
	// 6 - Модель Ламинарно турбулентного перехода Ментора Лангтрии (RANS). 15.01.2021
	TURBULENT_MODEL iturbmodel;
	// параметры модели Смагоринского:
	doublereal Cs; // постоянная Смагоринского.
	bool bDynamic_Stress ; // если   то включает динамическую модель Германо для определения квадрата постоянной Смагоринского.
	bool bLimiters_Cs; // включает ограничение на постоянную Смагоринского: ограничение на минимальное и максимальное значения.
	doublereal rminCs, rmaxCs; // минимальное и максимальное ограничение для константы Смагоринского.
	// тестовый фильтр который используется для осреднения в модели Германо 1991 года.
	integer itypeFiltrGermano; // SIMPSON filter	
	// поправка на неравномерной сетке,
	// модель Смагоринского-Лиллу,
	// учёт шероховатости стенки.
	integer bfdelta;
	integer bSmagorinskyLilly;
	integer bsurface_roughness;
	doublereal roughness; // шероховатость стенки 10micron.
	integer ipowerroughness; // показатель степени в модели учёта шероховатости стенки.
	doublereal rRi_mult; // корректирующий множитель подправляющий турбулентное число Ричардсона.
	doublereal rSelectiveAngle; // пороговое значение угла в модели Selective Smagorinsky.
	integer itypeSelectiveSmagorinsky_filtr; // SIMPSON тип фильтра в модели Смагоринского.
	// поправка на течения с кривизной линий тока,
	// избирательная модель Смагоринского.
	bool bSwirlAmendment, bSelectiveSmagorinsky;

	// Spalart-Allmares (RANS) [1992]. 19.04.2019
	// Эмпирические константы, определяющие SA модель турбулентности.
	doublereal sigma_nu;
	doublereal c_b1;
	doublereal c_b2;
	doublereal karman;
	doublereal c_w1;
	doublereal c_w2;
	doublereal c_w3;
	doublereal c_nu1;
	doublereal C_t3;
	doublereal C_t4;

	// SST Ментер (RANS) [1993]. 03.10.2019
	// А.В. Гарбарук Моделирование турбулентности в расчётах сложных течений. Политех 2012.
	// Эмпирические константы, определяющие SST Ментера модель турбулентности.
	doublereal sigma_k1;
	doublereal sigma_omega1;
	doublereal beta1;
	doublereal sigma_k2;
	doublereal sigma_omega2;
	doublereal beta2;
	doublereal beta_zvezda;
	doublereal menter_a1;
	doublereal alpha1;
	doublereal alpha2;

	// Модель ламинарно турбулентного перехода Лантгрии Ментера.
	// The calibration constants for the Langtry-Menter model are:
	doublereal Ca1;
	doublereal Ca2;
	doublereal Cepsilon1;
	doublereal Cepsilon2;
	doublereal CthetaT;
	doublereal s1;
	doublereal sigmaf;
	doublereal sigma_theta_T;


	// Стандартная k-epsilon модель не работает вблизи стенок.
	// Здесь используется двухслойная модель на основе стандартной
	// k-epsilon модели которая работает также и в пристеночной области.
	// А.В. Кузьминов, В.Н. Лапин, С.Г. Черный
	// Метод расчёта турбулентных течений 
	// несжимаемой жидкости на основе двухслойной (k-epsilon)-модели.
	// Новосибирск, 2001.
	// Граница между вязким подслоем 
	// и развитым турбулентным течением
	// характеризуется турбулентным числом Рейнольдса:
	// Standart k-epsilon (RANS) [2001] 23.10.2019
	doublereal Rey_zvezda;//100.0;
	doublereal A_switch; // A_switch == 1.0 --- 10.0;
	doublereal C_mu_std_ke;
	doublereal kar_std_ke;
	doublereal Cl_std_ke;
	doublereal Aepsilon_std_ke;
	doublereal CD_std_ke;
	doublereal Anu_std_ke;
	doublereal sigma_epsilon_std_ke;
	doublereal C_epsilon1_std_ke;
	doublereal C_epsilon2_std_ke;

	TFLOWINFO() {
		// опорная точка
		xc=1.0e30; yc=1.0e30; zc=1.0e30;
		// нужно ли рассчитывать поле течения
		iflow=0;
		// режим течения: 0 - ламинарный, 1 - турбулентный
		iflowregime= FLOW_REGIME::LAMINAR;
		// модель турбулентности
		// 0 - Zero Equation Turbulence Model (RANS).
		// 1 - модель Смагоринского (LES). (как опция модели, модель Германо 1991 (LES)).
		// 2 - RNG (LES).
		// 3 - Spalart-Allmares (RANS). 19.04.2019
		// 4 - K - Omega SST (RANS). 06.10.2019
		// 5 - Двухслойная модель на основе стандартной K-Epsilon модели (RANS). 31.10.2019
		iturbmodel= TURBULENT_MODEL::ZEROEQMOD;
		// параметры модели Смагоринского:
		Cs=0.151; // постоянная Смагоринского.
		bDynamic_Stress = false; // если   то включает динамическую модель Германо для определения квадрата постоянной Смагоринского.
		bLimiters_Cs = false; // включает ограничение на постоянную Смагоринского: ограничение на минимальное и максимальное значения.
		rminCs=-1.0e20; rmaxCs=1.0e23; // минимальное и максимальное ограничение для константы Смагоринского.
		// тестовый фильтр который используется для осреднения в модели Германо 1991 года.
		itypeFiltrGermano=2; // SIMPSON filter	
		// поправка на неравномерной сетке,
		// модель Смагоринского-Лиллу,
		// учёт шероховатости стенки.
		bfdelta=0; bSmagorinskyLilly=0; bsurface_roughness=0;
		roughness=10.0e-6; // шероховатость стенки 10micron.
		ipowerroughness=2; // показатель степени в модели учёта шероховатости стенки.
		rRi_mult=1.0; // корректирующий множитель подправляющий турбулентное число Ричардсона.
		rSelectiveAngle = 15.0; // пороговое значение угла в модели Selective Smagorinsky.
		itypeSelectiveSmagorinsky_filtr = 2; // SIMPSON тип фильтра в модели Смагоринского.
		// поправка на течения с кривизной линий тока,
		// избирательная модель Смагоринского.
		bSwirlAmendment = false; bSelectiveSmagorinsky=false;

		// Spalart-Allmares (RANS) [1992]. 19.04.2019
		// Эмпирические константы, определяющие SA модель турбулентности.
		sigma_nu = 2.0 / 3.0;
		c_b1 = 0.1355;
		c_b2 = 0.622;
		karman = 0.41;
		c_w1 = (c_b1/(karman*karman))+((1.0+c_b2)/sigma_nu);
		c_w2 = 0.3;
		c_w3 = 2.0;
		c_nu1 = 7.1;
		C_t3 = 1.2;
		C_t4 = 0.5;

		// SST Ментер (RANS) [1993]. 03.10.2019
		// А.В. Гарбарук Моделирование турбулентности в расчётах сложных течений. Политех 2012.
		// Эмпирические константы, определяющие SST Ментера модель турбулентности.
		sigma_k1 = 0.85;
		sigma_omega1 = 0.5;
		beta1 = 0.075;
		sigma_k2 = 1.0;
		sigma_omega2 = 0.856;
		beta2 = 0.0828;
		beta_zvezda = 0.09;
		menter_a1 = 0.31;
		alpha1 = 0.555;
		alpha2 = 0.44;

		// Модель ламинарно турбулентного перехода Лантгрии Ментера.
	    // The calibration constants for the Langtry-Menter model are:
		Ca1 = 2.0;
		Ca2 = 0.06;
		Cepsilon1 = 1.0;
		Cepsilon2 = 50.0;
		CthetaT = 0.03;
		s1 = 2.0;
		sigmaf = 1.0;
		sigma_theta_T = 2.0;


		// Стандартная k-epsilon модель не работает вблизи стенок.
		// Здесь используется двухслойная модель на основе стандартной
		// k-epsilon модели которая работает также и в пристеночной области.
		// А.В. Кузьминов, В.Н. Лапин, С.Г. Черный
		// Метод расчёта турбулентных течений 
		// несжимаемой жидкости на основе двухслойной (k-epsilon)-модели.
		// Новосибирск, 2001.
		// Граница между вязким подслоем 
		// и развитым турбулентным течением
		// характеризуется турбулентным числом Рейнольдса:
		// Standart k-epsilon (RANS) [2001] 23.10.2019
		Rey_zvezda = 100.0;//100.0;
		A_switch = 5.0; // A_switch == 1.0 --- 10.0;
		C_mu_std_ke = 0.09;
		kar_std_ke = 0.42;
		Cl_std_ke = kar_std_ke*pow(C_mu_std_ke,-3.0/4.0);
		Aepsilon_std_ke = 2.0*Cl_std_ke;
		CD_std_ke = 1.0;
		Anu_std_ke = 70.0;
		sigma_epsilon_std_ke = 1.3;
		C_epsilon1_std_ke = 1.44;
		C_epsilon2_std_ke = 1.92;
	}

	// функция переключатель между пристеночной областью и 
	// областью развитого турбулентного потока. 
	doublereal lambda_switch(doublereal Re_y) {
		/*if (Re_y > 10.0) {
			printf("Re=%e\n", Re_y);
			system("pause");
		}*/
		return 0.5*(1.0+tanh((Re_y - Rey_zvezda)/ A_switch));
	}

} FLOWINFO;

// вся информация о полном наборе решаемых уравнений
typedef struct TEQUATIONINFO {
	// нужно ли решать уравнение теплопроводности
	// 0 - ненужно, 1 - нужно методом контрольного объема, 2  - методом конечных элементов.
	integer itemper;
	// максимальное количество изолированных по жидкости
	// жидких зон (FLUID interior)
	integer imaxflD;
	// база данных с информацией о жидких зонах.
	FLOWINFO* fluidinfo;

	TEQUATIONINFO() {
		// нужно ли решать уравнение теплопроводности
		// 0 - ненужно, 1 - нужно методом контрольного объема, 2  - методом конечных элементов.
		itemper=1;
		// максимальное количество изолированных по жидкости
		// жидких зон (FLUID interior)
		imaxflD=0;
		// база данных с информацией о жидких зонах.
		fluidinfo=nullptr;
	}
} EQUATIONINFO;

EQUATIONINFO eqin; // информация о наборе решаемых уравнений

struct Tdatabase {
	float *x, *y, *z; // координаты узлов.
	integer maxelm;
	int** nvtxcell;
	integer ncell;
	// связь теплопередачи с гидродинамикой.
	int **ptr;

	Tdatabase() {
		x = nullptr; y = nullptr; z = nullptr; // координаты узлов.
		maxelm=0;
		nvtxcell = nullptr;
		ncell=0;
		// связь теплопередачи с гидродинамикой.
		ptr = nullptr;// для тестирования алгебраического многосеточного метода
	}
};

// одна строка матрицы СЛАУ в 3D варианте
// для внутреннего КО.
typedef struct Tequation3D {
	doublereal ap, ae, an, aw, as, at, ab, b;
	int iP, iE, iN, iT, iW, iS, iB;
	// Расширение структуры данных под АЛИС сетку.
	// На АЛИС сетке шаблон имеет переменное число связей и 
	// их число колеблется от семиточечного шаблона до
	// 25 точечного максимально.
	// На АЛИС сетке даже чисто диффузионная матрица не имеет диагонального преобладания и 
	// для решения СЛАУ с такой матрицей нужны специальные "робастые" методы типа BiCGStab + ILU2.
	// Признак существования дополнительных связей:
	// true если коэффициент существует.
	bool bE2, bN2, bT2, bW2, bS2, bB2;
	bool bE3, bN3, bT3, bW3, bS3, bB3;
	bool bE4, bN4, bT4, bW4, bS4, bB4;
	// Значение матричного коэффициента:
	doublereal ae2, an2, aw2, as2, at2, ab2;
	doublereal ae3, an3, aw3, as3, at3, ab3;
	doublereal ae4, an4, aw4, as4, at4, ab4;
	// индексация для дополнительных связей.
	int iE2, iN2, iT2, iW2, iS2, iB2;
	int iE3, iN3, iT3, iW3, iS3, iB3;
	int iE4, iN4, iT4, iW4, iS4, iB4;

	Tequation3D() {
		ap=0.0; ae=0.0; an=0.0; aw=0.0; as=0.0; at=0.0; ab=0.0; b=0.0;
		iP = -1; iE = -1; iN = -1; iT = -1; iW = -1; iS = -1; iB = -1;
		// Расширение структуры данных под АЛИС сетку.
		// На АЛИС сетке шаблон имеет переменное число связей и 
		// их число колеблется от семиточечного шаблона до
		// 25 точечного максимально.
		// На АЛИС сетке даже чисто диффузионная матрица не имеет диагонального преобладания и 
		// для решения СЛАУ с такой матрицей нужны специальные "робастые" методы типа BiCGStab + ILU2.
		// Признак существования дополнительных связей:
		// true если коэффициент существует.
		 bE2 = false; bN2 = false; bT2 = false; bW2 = false; bS2 = false; bB2 = false;
		 bE3 = false; bN3 = false; bT3 = false; bW3 = false; bS3 = false; bB3 = false;
		 bE4 = false; bN4 = false; bT4 = false; bW4 = false; bS4 = false; bB4 = false;
		 // Значение матричного коэффициента:
		 ae2 = 0.0; an2 = 0.0; aw2 = 0.0; as2 = 0.0; at2 = 0.0; ab2 = 0.0;
		 ae3 = 0.0; an3 = 0.0; aw3 = 0.0; as3 = 0.0; at3 = 0.0; ab3 = 0.0;
		 ae4 = 0.0; an4 = 0.0; aw4 = 0.0; as4 = 0.0; at4 = 0.0; ab4 = 0.0;
		 // индексация для дополнительных связей.
		 iE2 = -1; iN2 = -1; iT2 = -1; iW2 = -1; iS2 = -1; iB2 = -1;
		 iE3 = -1; iN3 = -1; iT3 = -1; iW3 = -1; iS3 = -1; iB3 = -1;
		 iE4 = -1; iN4 = -1; iT4 = -1; iW4 = -1; iS4 = -1; iB4 = -1;
	}
	
} equation3D;

// одна строка матрицы СЛАУ в 3D варианте
// для граничного КО.
typedef struct Tequation3D_bon {
	doublereal aw, ai, b; // wall, internal, правая часть
	int iW, iI;

	// соседи на границе области,
	// их позиции нужны в матрице:
	// На этих позициях будут стоять нули.
	integer iW1, iW2, iW3, iW4;

	Tequation3D_bon() {
		aw = 0.0; ai = 0.0; b = 0.0; // wall, internal, правая часть
		iW = -1; iI = -1;

		// соседи на границе области,
		// их позиции нужны в матрице:
		// На этих позициях будут стоять нули.
		iW1 = -1; iW2 = -1; iW3 = -1; iW4 = -1;
	}
} equation3D_bon;

typedef struct TTEMPER {
	// флаг отвечающий за освобождение 
	// памяти первого уровня. Пояснение:
	// После сборки матрицы СЛАУ можно 
	// уничтожить почти все структуры необходимые
	// для сборки матрицы. Это можно сделать 
	// только в том случае если матрица собирается
	// лишь единожды, как например в случае статики
	// для уравнения чистой теплопроводности.
	// true если освобождаем.
	bool free_temper_level1;
	// флаг отвечающий за освобождение памяти
	// второго уровня. Когда матрица СЛАУ перезаписывается
	// в формат SIMPLESPARSE. Исходную матрицу в формате 
	// equation3D можно убрать из оперативной памяти компьютера.
	bool free_temper_level2;

	int maxnod; // максимальный номер узла (размерность массива)
					// pa[0..maxnod-1];
	TOCHKA* pa; // координаты узлов сетки принадлежащие расчётной области

	int maxelm; // число ненулевых контрольных объёмов
					// nvtx[0..7][0..maxelm-1]
	int **nvtx; // список узлов для каждого элемента (ненулевого контрольного объёма)
						   // neighbors_for_the_internal_node[0..11][0..maxelm-1]
						   //integer **neighbors_for_the_internal_node; // соседние контрольные объёмы для каждого внутреннего контрольного объёма
						   // AliceMesh
	int*** neighbors_for_the_internal_node;// соседние контрольные объёмы для каждого внутреннего контрольного объёма
	int maxbound; // число граничных узлов
	int maxp; // maxp == maxelm + maxbound;
				  // border_neighbor[0..maxbound-1];
	BOUND* border_neighbor; // граничные узлы расчётной области 
						  // для всех граничных КО. Равно истине если имеем дело с
						  // границей строго внутри расчётной области причём на ней 
						  // расположен именно плоский бесконечно тонкий источник и по 
						  // одну его сторону расположена жидкость а по другую твёрдое тело.
	bool* binternalsource;

	// какому блоку принадлежит внутренний КО
	int* whot_is_block;

	int **ptr; // Связь с гидродинамикой.

						  // для графической визуализации
	integer ncell; // количество связей для контрольных объёмов.
	int **nvtxcell; // связи для контрольных объёмов.

							   // Для АЛИС сетки хранит номер уровня ячейки, это 
							   // потребуется при сборке матрицы.
							   // ilevel_alice[0..maxelm-1].
	integer *ilevel_alice;


	doublereal *potent; // массив узловых потенциалов (искомых функций)
	doublereal **total_deformation; // Полная деформация.

										   // Свойства материалов разделяются для  
										   // внутренних и для граничных КО.
										   // сначала prop[0..2][0..maxelm-1]
	bool *bActiveShearModule; // Задает ли пользователь модуль сдвига.
	float **prop; // свойства материалов для внутренних КО.
							  // сначала prop_b[0..2][0..maxbound-1]
	float **prop_b; // свойства для граничных КО.


	doublereal *Sc; // объёмное тепловыделение приходящееся на один выбранный контрольный объём.
	POWER_TIME_DEPEND  *ipower_time_depend; // закон зависимости мощности тепловыделения от времени.

	doublereal alpha; // параметр релаксации для температуры
	equation3D *slau; // коэффициенты матрицы СЛАУ для внутренних КО
	equation3D_bon *slau_bon; // коэффициенты матрицы СЛАУ для граничных КО

									 // для полилинейного метода:
									 // полилинейный метод рекомендован проф. Минесотского университета С. Патанкаром.
									 // полилинейный метод обладает фирменной особенностью - за первые несколько итераций 
									 // невязка падает на  несколько порядков. Это говорит о том что полилинейный метод 
									 // может быть использован как предобуславливатель в алгоритме Хенка ван дер Ворста - BiCGStab.
	NODELR_BASE *rootWE;
	NODELR_BASE *rootSN;
	NODELR_BASE *rootBT;

	integer iWE, iSN, iBT; // число сеточных линий вдоль каждого из направлений.

						   // Для полилинейного метода LR1:
						   // память будет выделяться и 
						   // уничтожаться один раз а не 
						   // каждый раз в цикле. 
						   // Это должно ускорить вычисления.
						   // Т.к. вычисления требуется распараллелить то память будет выделяться
						   // и уничтожаться многократно для каждой прогонки - одно выделение и одно
						   // уничтожение памяти. В каждой сеточной линии относительно небольшое число 
						   // узлов поэтому выделение памяти наверно не должно занять много времени по сравнению
						   // со временем вычисления.

						   // выходная невязка для температуры
						   // согласованная с точностью аппроксимации уравнения.
	doublereal resLR1sk; // O(h!3)

						 // Копия сеточных размеров для восстановления данных.
	int inx_copy, iny_copy, inz_copy;
	doublereal operatingtemperature_copy;

	doublereal *xpos_copy, *ypos_copy, *zpos_copy;

	// 9 августа 2015.
	// Для распараллеливания 
	int *ifrontregulationgl;
	int *ibackregulationgl; // обратное преобразование.

	doublereal operatingtemperature;

	Tdatabase database;

	TTEMPER() {
		// флаг отвечающий за освобождение 
		// памяти первого уровня. Пояснение:
		// После сборки матрицы СЛАУ можно 
		// уничтожить почти все структуры необходимые
		// для сборки матрицы. Это можно сделать 
		// только в том случае если матрица собирается
		// лишь единожды, как например в случае статики
		// для уравнения чистой теплопроводности.
		// true если освобождаем.
		free_temper_level1=true;
		// флаг отвечающий за освобождение памяти
		// второго уровня. Когда матрица СЛАУ перезаписывается
		// в формат SIMPLESPARSE. Исходную матрицу в формате 
		// equation3D можно убрать из оперативной памяти компьютера.
		 free_temper_level2=true;

		 maxnod=0; // максимальный номер узла (размерность массива)
					// pa[0..maxnod-1];
		pa = nullptr; // координаты узлов сетки принадлежащие расчётной области

		maxelm=0; // число ненулевых контрольных объёмов
					// nvtx[0..7][0..maxelm-1]
		nvtx = nullptr; // список узлов для каждого элемента (ненулевого контрольного объёма)
						   // neighbors_for_the_internal_node[0..11][0..maxelm-1]
						   //integer **neighbors_for_the_internal_node; // соседние контрольные объёмы для каждого внутреннего контрольного объёма
						   // AliceMesh
		neighbors_for_the_internal_node = nullptr;// соседние контрольные объёмы для каждого внутреннего контрольного объёма
		maxbound=0; // число граничных узлов
		maxp=0; // maxp == maxelm + maxbound;
				  // border_neighbor[0..maxbound-1];
		border_neighbor = nullptr; // граничные узлы расчётной области 
						  // для всех граничных КО. Равно истине если имеем дело с
						  // границей строго внутри расчётной области причём на ней 
						  // расположен именно плоский бесконечно тонкий источник и по 
						  // одну его сторону расположена жидкость а по другую твёрдое тело.
		binternalsource = nullptr;

		// какому блоку принадлежит внутренний КО
		whot_is_block = nullptr;

		ptr = nullptr; // Связь с гидродинамикой.

						  // для графической визуализации
		ncell=0; // количество связей для контрольных объёмов.
		nvtxcell = nullptr; // связи для контрольных объёмов.

							   // Для АЛИС сетки хранит номер уровня ячейки, это 
							   // потребуется при сборке матрицы.
							   // ilevel_alice[0..maxelm-1].
		ilevel_alice = nullptr;


		potent = nullptr; // массив узловых потенциалов (искомых функций)
		total_deformation = nullptr; // Полная деформация.

										   // Свойства материалов разделяются для  
										   // внутренних и для граничных КО.
										   // сначала prop[0..2][0..maxelm-1]
		bActiveShearModule = nullptr; // задаёт ли пользователь ЭВМ значение модуля сдвига.
		prop = nullptr; // свойства материалов для внутренних КО.
							  // сначала prop_b[0..2][0..maxbound-1]
		prop_b = nullptr; // свойства для граничных КО.


		Sc = nullptr; // объёмное тепловыделение приходящееся на один выбранный контрольный объём.
		ipower_time_depend = nullptr; // закон зависимости мощности тепловыделения от времени.

		alpha=1.0; // параметр релаксации для температуры
		slau = nullptr; // коэффициенты матрицы СЛАУ для внутренних КО
		slau_bon = nullptr; // коэффициенты матрицы СЛАУ для граничных КО

									 // для полилинейного метода:
									 // полилинейный метод рекомендован проф. Минесотского университета С. Патанкаром.
									 // полилинейный метод обладает фирменной особенностью - за первые несколько итераций 
									 // невязка падает на  несколько порядков. Это говорит о том что полилинейный метод 
									 // может быть использован как предобуславливатель в алгоритме Хенка ван дер Ворста - BiCGStab.
		rootWE = nullptr;
		rootSN = nullptr;
		rootBT = nullptr;

		iWE=-1000; iSN = -1000; iBT = -1000; // число сеточных линий вдоль каждого из направлений.

						   // Для полилинейного метода LR1:
						   // память будет выделяться и 
						   // уничтожаться один раз а не 
						   // каждый раз в цикле. 
						   // Это должно ускорить вычисления.
						   // Т.к. вычисления требуется распараллелить то память будет выделяться
						   // и уничтожаться многократно для каждой прогонки - одно выделение и одно
						   // уничтожение памяти. В каждой сеточной линии относительно небольшое число 
						   // узлов поэтому выделение памяти наверно не должно занять много времени по сравнению
						   // со временем вычисления.

						   // выходная невязка для температуры
						   // согласованная с точностью аппроксимации уравнения.
		resLR1sk=0.0; // O(h!3)

						 // Копия сеточных размеров для восстановления данных.
		inx_copy = 0; iny_copy = 0; inz_copy=0;
		operatingtemperature_copy=0.0;

		xpos_copy = nullptr; ypos_copy = nullptr; zpos_copy = nullptr;

		// 9 августа 2015.
		// Для распараллеливания 
		ifrontregulationgl = nullptr;
		ibackregulationgl = nullptr; // обратное преобразование.

		operatingtemperature = 20.0;

		//Tdatabase database;
	}

}  TEMPER;




typedef struct TFLOW {
	int maxnod; // максимальный номер узла (размерность массива)
	int maxelm; // число внутренних контрольных объёмов
					// nvtx[0..7][0..maxelm-1]
	int **nvtx; // список узлов для каждого внутреннего элемента (контрольного объёма)

						   // pa[0..maxnod-1]
	TOCHKA* pa; // координаты узлов сетки принадлежащие расчётной области

					   // neighbors_for_the_internal_node[0..11][0..maxelm-1]
					   //integer **neighbors_for_the_internal_node; // соседние контрольные объёмы для каждого внутреннего КО
					   // Для ALICEMESH сетки.
	int ***neighbors_for_the_internal_node;// соседние контрольные объёмы для каждого внутреннего КО
	int maxbound; // число граничных КО
	int maxp; // число уравнений
				  // border_neighbor[0..maxbound-1];
	BOUND* border_neighbor; // граничные узлы расчётной области

	int *ptr; // Связь с теплопроводностью


	// какому блоку принадлежит внутренний КО
	int* whot_is_block;

	TOCHKA* center_coord;
	TOCHKA* volume;

	// potent[iVar][0..maxp-1]
	doublereal **potent; // массив узловых потенциалов (искомых функций)
								// prop[0..2][0..maxelm-1]

	doublereal** turbulent_parameters_old_time_step;
	
	float **prop; // свойства материалов
							  // prop_b[0..2][0..maxbound-1]
	float **prop_b; // свойства материалов для граничных КО

	doublereal *alpha; // параметры нижней релаксации
	equation3D **slau; // коэффициенты матрицы СЛАУ для внутренних КО.
	equation3D_bon **slau_bon; // коэффициенты матрицы СЛАУ для граничных КО
									  // для реализации монотонизатора Рхи-Чоу 1983 требуется хранить диагональные коэффициенты.
	doublereal **diag_coef;
	doublereal OpTemp; // Operating Temperature


	bool bactive; // нужно-ли рассчитывать поле течения
	bool bPressureFix; // нужно ли фиксировать давление в одной точке
	bool bLR1free; // нужно ли применять плавающий полилинейный солвер (он показывает более быструю сходимость).

				   // для полилинейного метода:
				   // полилинейный метод рекомендован проф. Минесотского университета С. Патанкаром.
				   // полилинейный метод обладает фирменной особенностью - за первые несколько итераций 
				   // невязка падает на  несколько порядков. Это говорит о том что полилинейный метод 
				   // может быть использован как предобуславливатель в алгоритме Хенка ван дер Ворста - BiCGStab 1992.

	integer iWE, iSN, iBT; // число сеточных линий вдоль каждого из направлений.
	integer** iN; //iN[3][max3(iWE,iSN,iBT)];
	integer*** id; //id[3][max3(iWE,iSN,iBT)][max(iN)]; 

						  // Для полилинейного метода LR1:
						  // память будет выделяться и 
						  // уничтожаться один раз а не 
						  // каждый раз в цикле. 
						  // Это должно ускорить вычисления.
						  // Т.к. вычисления требуется распараллелить то память будет выделяться
						  // и уничтожаться многократно для каждой прогонки - одно выделение и одно
						  // уничтожение памяти. В каждой сеточной линии относительно небольшое число 
						  // узлов поэтому выделение памяти наверно не должно занять много времени по сравнению
						  // со временем вычисления.

						  // В реальности большинство течений Турбулентны,
						  // здесь содержится некоторая вспомогательная информация для расчёта турбулентных течений.
						  // режим течения для данной зоны FLUID
	VISCOSITY_MODEL iflowregime; // default LAMINAR
						 // Кратчайшее расстояние до ближайшей стенки [0..maxelm-1]
	doublereal* rdistWall; // расстояние до ближайшей твёрдой стенки.
								  // Толщина пограничного слоя в формуле Эскудиера
	doublereal rdistWallmax;
	// S инвариант тензора скоростей-деформаций
	doublereal* SInvariantStrainRateTensor; // [0..maxelm+maxbound-1]; // инициализируется нулём.

												   // массовый поток через грани КО:
	doublereal** mf;

	SMAGORINSKYINFO smaginfo; // параметры модели Смагоринского.

							  // выходная невязка для поправки давления
							  // согласованная с точностью аппроксимации уравнения.
	doublereal resICCG; // O(h!2)
	doublereal resLR1sk; // O(h!3)

						 // Для правильной работы mass balance для естественно конвективных задач
						 // нужен идентификатор разных по связности гидродинамических областей для
						 // внутренних КО.
	integer *icolor_different_fluid_domain; // Разные гидродинамические подобласти имеют разные цвета.

												   // 9 августа 2015.
												   // Для распараллеливания 
	int *ifrontregulationgl;
	int *ibackregulationgl; // обратное преобразование.

	TFLOW() {
		 maxnod=0; // максимальный номер узла (размерность массива)
		maxelm = 0; // число внутренних контрольных объёмов
					// nvtx[0..7][0..maxelm-1]
		nvtx = nullptr; // список узлов для каждого внутреннего элемента (контрольного объёма)

						   // pa[0..maxnod-1]
		pa = nullptr; // координаты узлов сетки принадлежащие расчётной области

					   // neighbors_for_the_internal_node[0..11][0..maxelm-1]
					   //integer **neighbors_for_the_internal_node; // соседние контрольные объёмы для каждого внутреннего КО
					   // Для ALICEMESH сетки.
		neighbors_for_the_internal_node = nullptr;// соседние контрольные объёмы для каждого внутреннего КО
		maxbound = 0; // число граничных КО
		maxp = 0; // число уравнений
				  // border_neighbor[0..maxbound-1];
		border_neighbor = nullptr; // граничные узлы расчётной области

		ptr = nullptr; // Связь с теплопроводностью


						 // какому блоку принадлежит внутренний КО
		whot_is_block = nullptr;

		center_coord = nullptr;
		volume = nullptr;

		// potent[iVar][0..maxp-1]
		potent = nullptr; // массив узловых потенциалов (искомых функций)
								// prop[0..2][0..maxelm-1]

		turbulent_parameters_old_time_step = nullptr;

		
		prop = nullptr; // свойства материалов
							  // prop_b[0..2][0..maxbound-1]
		prop_b = nullptr; // свойства материалов для граничных КО.

		alpha = nullptr; // параметры нижней релаксации.
		slau = nullptr; // коэффициенты матрицы СЛАУ для внутренних КО.
		slau_bon = nullptr; // коэффициенты матрицы СЛАУ для граничных КО
									  // для реализации монотонизатора Рхи-Чоу 1983 требуется хранить диагональные коэффициенты.
		diag_coef = nullptr;
		 OpTemp=0.0; // Operating Temperature


		bactive=true; // нужно-ли рассчитывать поле течения.
		bPressureFix = false; // нужно ли фиксировать давление в одной точке.
		bLR1free=false; // нужно ли применять плавающий полилинейный солвер (он показывает более быструю сходимость).

				   // для полилинейного метода:
				   // полилинейный метод рекомендован проф. Минесотского университета С. Патанкаром.
				   // полилинейный метод обладает фирменной особенностью - за первые несколько итераций 
				   // невязка падает на  несколько порядков. Это говорит о том что полилинейный метод 
				   // может быть использован как предобуславливатель в алгоритме Хенка ван дер Ворста - BiCGStab 1992.

		iWE=-1000; iSN = -1000; iBT = -1000; // число сеточных линий вдоль каждого из направлений.
		iN = nullptr; //iN[3][max3(iWE,iSN,iBT)];
		id = nullptr; //id[3][max3(iWE,iSN,iBT)][max(iN)]; 

						  // Для полилинейного метода LR1:
						  // память будет выделяться и 
						  // уничтожаться один раз а не 
						  // каждый раз в цикле. 
						  // Это должно ускорить вычисления.
						  // Т.к. вычисления требуется распараллелить то память будет выделяться
						  // и уничтожаться многократно для каждой прогонки - одно выделение и одно
						  // уничтожение памяти. В каждой сеточной линии относительно небольшое число 
						  // узлов поэтому выделение памяти наверно не должно занять много времени по сравнению
						  // со временем вычисления.

						  // В реальности большинство течений Турбулентны,
						  // здесь содержится некоторая вспомогательная информация для расчёта турбулентных течений.
						  // режим течения для данной зоны FLUID
		iflowregime= VISCOSITY_MODEL::LAMINAR; // default LAMINAR
						 // Кратчайшее расстояние до ближайшей стенки [0..maxelm-1]
		rdistWall = nullptr; // расстояние до ближайшей твёрдой стенки.
								  // Толщина пограничного слоя в формуле Эскудиера
		 rdistWallmax=0.0;
		// S инвариант тензора скоростей-деформаций
		SInvariantStrainRateTensor = nullptr; // [0..maxelm+maxbound-1]; // инициализируется нулём.

												   // массовый поток через грани КО:
		mf = nullptr;

		//SMAGORINSKYINFO smaginfo; // параметры модели Смагоринского.

							  // выходная невязка для поправки давления
							  // согласованная с точностью аппроксимации уравнения.
		resICCG=0.0; // O(h!2)
		resLR1sk=0.0; // O(h!3)

						 // Для правильной работы mass balance для естественно конвективных задач
						 // нужен идентификатор разных по связности гидродинамических областей для
						 // внутренних КО.
		icolor_different_fluid_domain = nullptr; // Разные гидродинамические подобласти имеют разные цвета.

												   // 9 августа 2015.
												   // Для распараллеливания 
		ifrontregulationgl = nullptr;
		ibackregulationgl = nullptr; // обратное преобразование.
	}

} FLOW;




// Объединение
// Для блочно структурированной расчётной сетки.
// Начало разработки 25.04.2018.
typedef struct TUNION {
	// id передаётся из интерфейса.
	integer id; // Уникальный номер объединения.

	// -1 Cabinet
	integer iunion_parent; // идентификатор родительского объединения.

				// Из кабинета юнион видится как Hollow блок
				// в виде прямой прямоугольной призмы.
				// размеры передаются из интерфейса.
	doublereal xS, xE, yS, yE, zS, zE;

	// Внутренняя сетка union.
	// Union является кабинетом для своих внутренних блоков.
	// для внутреннего пользования.
	doublereal *xpos, *ypos, *zpos;
	doublereal *xposadd, *yposadd, *zposadd;

	// Для внутренней сетки (размерности).
	// Передаётся из интерфейса.
	int inx, iny, inz;
	// для внутреннего пользования.
	int inxadd, inyadd, inzadd;

	//  Тип сеточного генератора
	// передаётся из интерфейса.
	CONFORMAL_MESH_GENERATOR_SELECTOR iswitchMeshGenerator; // 2 - CoarseMeshGen

															// Локальные объявления.
	TEMPER t;
	int flow_interior; // Суммарное число FLUID зон
	FLOW* f;

	bool active;

	char name[80]; // имя асемблеса.

	TUNION() {


		name[0] = '\0'; // Текстовое имя асемблеса.
		iunion_parent = -1; // -1 Cabinet

		// id передаётся из интерфейса.
		id = -1; // Уникальный номер объединения.

				 // Из кабинета юнион видится как Hollow блок
				 // в виде прямой прямоугольной призмы.
				 // размеры передаются из интерфейса.
		xS = 0.0; xE = 0.0; yS = 0.0;
		yE = 0.0; zS = 0.0; zE = 0.0;

		// Внутренняя сетка union.
		// Union является кабинетом для своих внутренних блоков.
		// для внутреннего пользования.
		xpos = nullptr; ypos = nullptr; zpos = nullptr;
		xposadd = nullptr; yposadd = nullptr; zposadd = nullptr;

		// Для внутренней сетки (размерности).
		// Передаётся из интерфейса.
		inx = -1; iny = -1; inz = -1;
		// для внутреннего пользования.
		inxadd = -1; inyadd = -1; inzadd = -1;

		//  Тип сеточного генератора
		// передаётся из интерфейса.
		iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; // 2 - CoarseMeshGen

													 // Локальные объявления.
													 //TEMPER t;
		flow_interior = -1; // Суммарное число FLUID зон
		f = nullptr;

		active = false;
	}

} UNION;

// Проверяет если ли выход за пределы кабинета
// среди блоков, стенок и источников тепла. 02.08.2019.
// Переделал с printf("..."); на std::cout<<"..."<<std::endl; 14.08.2019.
void BODY_CHECK(BLOCK* &b, int lb, WALL* &w, int lw, SOURCE* &s, int ls, UNION* &u, int &lu) {
	bool bOk = true;
	const doublereal dcheck_eps = 1e-10; // 1.0e-10 подобрано для tgf*
	
	// Нельзя распараллеливать т.к. там вложенный параллелизм при вычислении объёма полигона.

	for (integer i = lb - 1; i >= 0; i--) {
		// Проверка типа геометрии.

		// Излучение в вакуумной призме не работает совместно с температурным солвером 
		// на основе графового метода решения уравнения теплопередачи.
		if (0&&(b[i].radiation.binternalRadiation)&&(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T)) {
			std::cout << "Error:  body[" << i << "].name = " << b[i].name << " body[" << i << "] radiation.binternalRadiation is ON." << std::endl;
			std::cout << "steady_or_unsteady_global_determinant == NETWORK_T" << std::endl;
			std::cout << "synopsis: incompatible solver settings..." << std::endl;
			std::cout << "recomended: set the control volume method as a temperature solver and a structured grid." << std::endl;
			system("pause");
			bOk = false;
			exit(1);
		}

		// Излучение в вакуумной призме не работает совместно с температурным солвером
		// на основе метода конечных элементов.
		if ((b[i].radiation.binternalRadiation) && (eqin.itemper == 2)) {
			std::cout << "Error:  body[" << i << "].name = " << b[i].name << " body[" << i << "] radiation.binternalRadiation is ON." << std::endl;
			std::cout << "eqin.itemper == 2 Finite Element Method for Temperature solver." << std::endl;
			std::cout << "synopsis: incompatible solver settings..." << std::endl;
			std::cout << "recomended: set the control volume method as a temperature solver and a structured grid." << std::endl;
			system("pause");
			bOk = false;
			exit(1);
		}

		// Излучение в вакуумной призме не работает совместно с гидродинамическим расчётом.
		if ((b[i].radiation.binternalRadiation) && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY)) {
			std::cout << "Error:  body[" << i << "].name = " << b[i].name << " body[" << i << "] radiation.binternalRadiation is ON." << std::endl;
			std::cout << "steady_or_unsteady_global_determinant == CFD_STEADY." << std::endl;
			std::cout << "synopsis: incompatible solver settings..." << std::endl;
			std::cout << "recomended: set the control volume method as a temperature solver and a structured grid. turn OFF cfd." << std::endl;
			system("pause");
			bOk = false;
			exit(1);
		}
		
		// Излучение в вакуумной призме не работает с Адаптивными Локально Измельченными Сетками.
		if ((b[i].radiation.binternalRadiation) && (b_on_adaptive_local_refinement_mesh)&&
			(!(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY))) {
			std::cout << "Error:  body[" << i << "].name = " << b[i].name << " body[" << i << "] radiation.binternalRadiation is ON." << std::endl;
			std::cout << "b_on_adaptive_local_refinement_mesh   Adaptive Local Refinement Mesh is ON." << std::endl;
			std::cout << "synopsis: incompatible solver settings..." << std::endl;
			std::cout << "recomended: set the control volume method as a temperature solver and a structured grid. turn OFF cfd." << std::endl;
			system("pause");
			bOk = false;
			exit(1);
		}


		CHECK_TYPE_GEOM(b[i].g.itypegeom);
		if (b[i].g.itypegeom == PRISM) {
			// Ошибка не число.
			if (b[i].g.xS != b[i].g.xS) {
				//printf("Error body[%lld] xS NOT VALUE=%e.\n", i, b[i].g.xS);
				std::cout << "Error:  body[" << i << "].name = " << b[i].name << " body[" << i << "] xS NOT VALUE=" << b[i].g.xS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.xE != b[i].g.xE) {
				//printf("Error body[%lld] xE NOT VALUE=%e.\n", i, b[i].g.xE);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] xE NOT VALUE=" << b[i].g.xE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.yS != b[i].g.yS) {
				//printf("Error body[%lld] yS NOT VALUE=%e.\n", i, b[i].g.yS);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] yS NOT VALUE=" << b[i].g.yS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.yE != b[i].g.yE) {
				//printf("Error body[%lld] yE NOT VALUE=%e.\n", i, b[i].g.yE);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] yE NOT VALUE=" << b[i].g.yE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.zS != b[i].g.zS) {
				//printf("Error body[%lld] zS NOT VALUE=%e.\n", i, b[i].g.zS);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] zS NOT VALUE=" << b[i].g.zS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.zE != b[i].g.zE) {
				//printf("Error body[%lld] zE NOT VALUE=%e.\n", i, b[i].g.zE);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] zE NOT VALUE=" << b[i].g.zE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (b[i].g.xS >= b[i].g.xE) {
				//printf("body[%lld].xS=%e >= body[%lld].xE=%e\n", i, b[i].g.xS, i, b[i].g.xE);
				std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "].xS=" << b[i].g.xS << " >= body[" << i << "].xE=" << b[i].g.xE << std::endl;
				system("pause");
				doublereal temp = b[i].g.xE;
				b[i].g.xE = b[i].g.xS;
				b[i].g.xS = temp;
				bOk = false;
			}
			if (b[i].g.yS >= b[i].g.yE) {
				//printf("body[%lld].yS=%e >= body[%lld].yE=%e\n", i, b[i].g.yS, i, b[i].g.yE);
				std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "].yS=" << b[i].g.yS << " >= body[" << i << "].yE=" << b[i].g.yE << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				doublereal temp = b[i].g.yE;
				b[i].g.yE = b[i].g.yS;
				b[i].g.yS = temp;
				bOk = false;
			}
			if (b[i].g.zS >= b[i].g.zE) {
				//printf("body[%lld].zS=%e >= body[%lld].zE=%e\n", i, b[i].g.zS, i, b[i].g.zE);
				std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "].zS=" << b[i].g.zS << " >= body[" << i << "].zE=" << b[i].g.zE << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				doublereal temp = b[i].g.zE;
				b[i].g.zE = b[i].g.zS;
				b[i].g.zS = temp;
				bOk = false;
			}
			if (i != 0) {
				// Не cabinet
				if (b[i].g.xE > b[0].g.xE) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.xE=%e > cabinet.xE=%e.\n", i, b[i].g.xE, b[0].g.xE);
					std::cout << "body[" << i << "].name = " << b[i].name << " b[" << i << "].g.xE=" << b[i].g.xE << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
					system("pause");
					b[i].g.xE = b[0].g.xE; // исправление.
					bOk = false;
				}
				if (b[i].g.yE > b[0].g.yE) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.yE=%e > cabinet.yE=%e.\n", i, b[i].g.yE, b[0].g.yE);
					std::cout << "body[" << i << "].name = " << b[i].name << " b[" << i << "].g.yE=" << b[i].g.yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
					system("pause");
					b[i].g.yE = b[0].g.yE; // исправление.
					bOk = false;
				}
				if (b[i].g.zE > b[0].g.zE) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.zE=%e > cabinet.zE=%e.\n", i, b[i].g.zE, b[0].g.zE);
					std::cout << "body[" << i << "].name = " << b[i].name << " b[" << i << "].g.zE=" << b[i].g.zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
					system("pause");
					b[i].g.zE = b[0].g.zE; // исправление.
					bOk = false;
				}
				if (b[i].g.xS < b[0].g.xS) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.xS=%e < cabinet.xS=%e.\n", i, b[i].g.xS, b[0].g.xS);
					std::cout << "body[" << i << "].name = " << b[i].name << " b[" << i << "].g.xS=" << b[i].g.xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
					system("pause");
					b[i].g.xS = b[0].g.xS; // исправление.
					bOk = false;
				}
				if (b[i].g.yS < b[0].g.yS) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.yS=%e < cabinet.yS=%e.\n", i, b[i].g.yS, b[0].g.yS);
					std::cout << "body[" << i << "].name = " << b[i].name << " b[" << i << "].g.yS=" << b[i].g.yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
					std::cout << "" << std::endl;
					system("pause");
					b[i].g.yS = b[0].g.yS; // исправление.
					bOk = false;
				}
				if (b[i].g.zS < b[0].g.zS) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.zS=%e < cabinet.zS=%e.\n", i, b[i].g.zS, b[0].g.zS);
					std::cout << "body[" << i << "].name = " << b[i].name << " b[" << i << "].g.zS=" << b[i].g.zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
					system("pause");
					b[i].g.zS = b[0].g.zS; // исправление.
					bOk = false;
				}
			}
		}
		if (b[i].g.itypegeom == CYLINDER) {
			if (b[i].g.xC != b[i].g.xC) {
				//printf("Error body[%lld] xC NOT VALUE=%e.\n",i, b[i].g.xC);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] xC NOT VALUE=" << b[i].g.xC << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.yC != b[i].g.yC) {
				//printf("Error body[%lld] yC NOT VALUE=%e.\n", i, b[i].g.yC);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] yC NOT VALUE=" << b[i].g.yC << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.zC != b[i].g.zC) {
				//printf("Error body[%lld] zC NOT VALUE=%e.\n", i, b[i].g.zC);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] zC NOT VALUE=" << b[i].g.zC << "." << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.Hcyl != b[i].g.Hcyl) {
				//printf("Error body[%lld] Hcyl NOT VALUE=%e.\n", i, b[i].g.Hcyl);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] Hcyl NOT VALUE=" << b[i].g.Hcyl << "." << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_out_cyl != b[i].g.R_out_cyl) {
				//printf("Error body[%lld] R_out_cyl NOT VALUE=%e.\n", i, b[i].g.R_out_cyl);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] R_out_cyl NOT VALUE=" << b[i].g.R_out_cyl << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl != b[i].g.R_in_cyl) {
				//printf("Error body[%lld] R_in_cyl NOT VALUE=%e.\n", i, b[i].g.R_in_cyl);
				std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] R_in_cyl NOT VALUE=" << b[i].g.R_in_cyl << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.Hcyl <= 0.0) {
				//printf("Error: negative or zero Hcyl=%e of cylinder body[%lld]\n", b[i].g.Hcyl, i);
				std::cout << "Error: body[" << i << "].name = " << b[i].name << " negative or zero Hcyl=" << b[i].g.Hcyl << " of cylinder body[" << i << "]" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_out_cyl <= 0.0) {
				//printf("Error: negative or zero R_out_cyl=%e of cylinder body[%lld]\n", b[i].g.R_out_cyl, i);
				std::cout << "Error: body[" << i << "].name = " << b[i].name << " negative or zero R_out_cyl=" << b[i].g.R_out_cyl << " of cylinder body[" << i << "]" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl < 0.0) {
				//printf("Error: negative R_in_cyl=%e of cylinder body[%lld]\n", b[i].g.R_in_cyl, i);
				std::cout << "Error: body[" << i << "].name = " << b[i].name << " negative  R_in_cyl=" << b[i].g.R_in_cyl << " of cylinder body[" << i << "]" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl > b[i].g.R_out_cyl) {
				//printf("Error: b[%lld].g.R_in_cyl=%e > b[%lld].g.R_out_cyl=%e\n", i, b[i].g.R_in_cyl, i, b[i].g.R_out_cyl);
				std::cout << "Error: body[" << i << "].name = " << b[i].name << " b[" << i << "].g.R_in_cyl=" << b[i].g.R_in_cyl << " > b[" << i << "].g.R_out_cyl=" << b[i].g.R_out_cyl << std::endl;
				system("pause");
				bOk = false;
			}
			if (i != 0) {
				// Не cabinet

				// Внешние границы цилиндра -
				// окаймляющая прямая призма.
				doublereal xS_box_cyl = b[i].g.xC;
				doublereal xE_box_cyl = b[i].g.xC;
				doublereal yS_box_cyl = b[i].g.yC;
				doublereal yE_box_cyl = b[i].g.yC;
				doublereal zS_box_cyl = b[i].g.zC;
				doublereal zE_box_cyl = b[i].g.zC;
				switch (b[i].g.iPlane) {
				case XY_PLANE:
					xS_box_cyl -= b[i].g.R_out_cyl;
					xE_box_cyl += b[i].g.R_out_cyl;
					yS_box_cyl -= b[i].g.R_out_cyl;
					yE_box_cyl += b[i].g.R_out_cyl;
					zE_box_cyl += b[i].g.Hcyl;
					break;
				case XZ_PLANE:
					xS_box_cyl -= b[i].g.R_out_cyl;
					xE_box_cyl += b[i].g.R_out_cyl;
					zS_box_cyl -= b[i].g.R_out_cyl;
					zE_box_cyl += b[i].g.R_out_cyl;
					yE_box_cyl += b[i].g.Hcyl;
					break;
				case YZ_PLANE:
					zS_box_cyl -= b[i].g.R_out_cyl;
					zE_box_cyl += b[i].g.R_out_cyl;
					yS_box_cyl -= b[i].g.R_out_cyl;
					yE_box_cyl += b[i].g.R_out_cyl;
					xE_box_cyl += b[i].g.Hcyl;
					break;
				}

				if (xE_box_cyl > b[0].g.xE + dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] xE_box_cyl=%e > cabinet.xE=%e.\n", i, xE_box_cyl, b[0].g.xE);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] xE_box_cyl=" << xE_box_cyl << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
					printf("tolerance=%e\n", xE_box_cyl - b[0].g.xE);
					system("pause");

					bOk = false;
				}
				if (yE_box_cyl > b[0].g.yE + dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] yE_box_cyl=%e > cabinet.yE=%e.\n", i, yE_box_cyl, b[0].g.yE);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] yE_box_cyl=" << yE_box_cyl << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
					printf("tolerance=%e\n", yE_box_cyl - b[0].g.yE);
					system("pause");

					bOk = false;
				}
				if (zE_box_cyl > b[0].g.zE + dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] zE_box_cyl=%e > cabinet.zE=%e.\n", i, zE_box_cyl, b[0].g.zE);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] zE_box_cyl=" << zE_box_cyl << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
					printf("tolerance=%e\n", zE_box_cyl - b[0].g.zE);
					system("pause");

					bOk = false;
				}
				if (xS_box_cyl < b[0].g.xS - dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] xS_box_cyl=%e < cabinet.xS=%e.\n", i, xS_box_cyl, b[0].g.xS);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] xS_box_cyl=" << xS_box_cyl << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
					printf("tolerance=%e\n", b[0].g.xS - xS_box_cyl);
					system("pause");

					bOk = false;
				}
				if (yS_box_cyl < b[0].g.yS - dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] yS_box_cyl=%e < cabinet.yS=%e.\n", i, yS_box_cyl, b[0].g.yS);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] yS_box_cyl=" << yS_box_cyl << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
					printf("tolerance=%e\n", b[0].g.yS - yS_box_cyl);
					system("pause");

					bOk = false;
				}
				if (zS_box_cyl < b[0].g.zS - dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] zS_box_cyl=%e < cabinet.zS=%e.\n", i, zS_box_cyl, b[0].g.zS);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] zS_box_cyl=" << zS_box_cyl << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
					printf("tolerance=%e\n", b[0].g.zS - zS_box_cyl);
					system("pause");

					bOk = false;
				}
			}
		}

		if (b[i].g.itypegeom == POLYGON) {
			if (i != 0) {
				// Не cabinet

				// Внешние границы полигона -
				// окаймляющая прямая призма.
				doublereal xS_box_polygon = 1.0e30;
				doublereal xE_box_polygon = -1.0e30;
				doublereal yS_box_polygon = 1.0e30;
				doublereal yE_box_polygon = -1.0e30;
				doublereal zS_box_polygon = 1.0e30;
				doublereal zE_box_polygon = -1.0e30;

				if (b[i].g.nsizei<2) {
					//printf("Error number points on polygon very small: b[%lld].g.nsizei=%lld<2\n",i, b[i].g.nsizei);
					std::cout << "Error number points on polygon very small: body[" << i << "].name = " << b[i].name << " b[" << i << "].g.nsizei=" << b[i].g.nsizei << "<2." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.xi == nullptr) {
					//printf("Error: no memory allocate for b[%lld].g.xi array.\n",i);
					std::cout << "Error: body[" << i << "].name = " << b[i].name << " no memory allocate for b[" << i << "].g.xi array." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.yi == nullptr) {
					//printf("Error: no memory allocate for b[%lld].g.yi array.\n", i);
					std::cout << "Error: body[" << i << "].name = " << b[i].name << " no memory allocate for b[" << i << "].g.yi array." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.zi == nullptr) {
					//printf("Error: no memory allocate for b[%lld].g.zi array.\n", i);
					std::cout << "Error: body[" << i << "].name = " << b[i].name << " no memory allocate for b[" << i << "].g.zi array." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.hi == nullptr) {
					//printf("Error: no memory allocate for b[%lld].g.hi array.\n", i);
					std::cout << "Error: body[" << i << "].name = " << b[i].name << " no memory allocate for b[" << i << "].g.hi array." << std::endl;
					system("pause");
					bOk = false;
				}

				for (integer ip = 0; ip < b[i].g.nsizei; ip++) {
					if (b[i].g.xi[ip] != b[i].g.xi[ip]) {
						//printf("Error body[%lld] xi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.xi[ip]);
						std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] xi[" << ip << "] NOT VALUE=" << b[i].g.xi[ip] << "." << std::endl;
						system("pause");
						bOk = false;
					}
					if (b[i].g.yi[ip] != b[i].g.yi[ip]) {
						//printf("Error body[%lld] yi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.yi[ip]);
						std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] yi[" << ip << "] NOT VALUE=" << b[i].g.yi[ip] << "." << std::endl;
						system("pause");
						bOk = false;
					}
					if (b[i].g.zi[ip] != b[i].g.zi[ip]) {
						//printf("Error body[%lld] zi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.zi[ip]);
						std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] zi[" << ip << "] NOT VALUE=" << b[i].g.zi[ip] << "." << std::endl;
						system("pause");
						bOk = false;
					}
					if (b[i].g.hi[ip] != b[i].g.hi[ip]) {
						//printf("Error body[%lld] hi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.hi[ip]);
						std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "] hi[" << ip << "] NOT VALUE=" << b[i].g.hi[ip] << "." << std::endl;
						system("pause");
						bOk = false;
					}
				}

				integer i_t1, i_t2, i_t3, it_c = 0;
				switch (b[i].g.iPlane_obj2) {
				case XY_PLANE:
					for (integer ip = 0; ip < b[i].g.nsizei; ip++) {

						// Детектирование ОШИБКИ самопересечения отрезков в полигоне.
						// Мы проверяем каждый с каждым
						it_c = 0;
						for (integer ip_t = 2; ip_t < b[i].g.nsizei; ip_t++) {
							i_t1 = ip + 1;
							if (i_t1 >= b[i].g.nsizei) i_t1 -= b[i].g.nsizei;
							i_t2 = ip + it_c + 2;
							if (i_t2 >= b[i].g.nsizei) i_t2 -= b[i].g.nsizei;
							i_t3 = ip + it_c + 3;
							if (i_t3 >= b[i].g.nsizei) i_t3 -= b[i].g.nsizei;

							// Мы попарно проверяем есть ли самопересечения каждый с каждым.
							if (b_is_intersect(b[i].g.xi[ip], b[i].g.yi[ip], b[i].g.xi[i_t1], b[i].g.yi[i_t1],
								b[i].g.xi[i_t2], b[i].g.yi[i_t2], b[i].g.xi[i_t3], b[i].g.yi[i_t3])) {
								// ОШИБКА: обнаружено самопересечение отрезков в полигоне.
								std::cout << "Error in POLYGON: body[" << i << "].name = " << b[i].name << " Self-intersection of segments in body[" << i << "]\n";
								system("pause");
								bOk = false;
							}

							it_c++;
						}

						if (b[i].g.xi[ip] < xS_box_polygon) {
							xS_box_polygon = b[i].g.xi[ip];
						}
						if (b[i].g.xi[ip] > xE_box_polygon) {
							xE_box_polygon = b[i].g.xi[ip];
						}
						if (b[i].g.yi[ip] < yS_box_polygon) {
							yS_box_polygon = b[i].g.yi[ip];
						}
						if (b[i].g.yi[ip] > yE_box_polygon) {
							yE_box_polygon = b[i].g.yi[ip];
						}
						if (ip != 0) {
							if (fabs(b[i].g.zi[ip] - b[i].g.zi[0]) > 1.0e-36) {
								std::cout << "Error: NON planar polygon. body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								std::cout << "Diagnostic: fabs(b[" << i << "].g.zi[" << ip << "]- b[" << i << "].g.zi[0]) > 1.0e-36\n";
								system("pause");
								bOk = false;
							}
							if (fabs(b[i].g.hi[ip] - b[i].g.hi[0]) > 1.0e-36) {
								std::cout << "Error: NON planar polygon. body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								std::cout << "Diagnostic: fabs(b[" << i << "].g.hi[" << ip << "]- b[" << i << "].g.hi[0]) > 1.0e-36\n";
								system("pause");
								bOk = false;
							}
						}
						if (b[i].g.hi[ip] <= 0.0) {
							//printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "].g.hi[" << ip << "]=" << b[i].g.hi[ip] << " is NEGATIV or ZERO." << std::endl;
							system("pause");
							bOk = false;
						}
						zS_box_polygon = b[i].g.zi[ip];
						zE_box_polygon = b[i].g.zi[ip] + b[i].g.hi[ip];
					}
					break;
				case XZ_PLANE:
					for (integer ip = 0; ip < b[i].g.nsizei; ip++) {

						// Детектирование ОШИБКИ самопересечения отрезков в полигоне.
						// Мы проверяем каждый с каждым
						it_c = 0;
						for (integer ip_t = 2; ip_t < b[i].g.nsizei; ip_t++) {
							i_t1 = ip + 1;
							if (i_t1 >= b[i].g.nsizei) i_t1 -= b[i].g.nsizei;
							i_t2 = ip + it_c + 2;
							if (i_t2 >= b[i].g.nsizei) i_t2 -= b[i].g.nsizei;
							i_t3 = ip + it_c + 3;
							if (i_t3 >= b[i].g.nsizei) i_t3 -= b[i].g.nsizei;

							// Мы попарно проверяем есть ли самопересечения каждый с каждым.
							if (b_is_intersect(b[i].g.xi[ip], b[i].g.zi[ip], b[i].g.xi[i_t1], b[i].g.zi[i_t1],
								b[i].g.xi[i_t2], b[i].g.zi[i_t2], b[i].g.xi[i_t3], b[i].g.zi[i_t3])) {
								// ОШИБКА: обнаружено самопересечение отрезков в полигоне.
								std::cout << "Error in POLYGON: Self-intersection of segments in body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								system("pause");
								bOk = false;
							}

							it_c++;
						}


						if (b[i].g.xi[ip] < xS_box_polygon) {
							xS_box_polygon = b[i].g.xi[ip];
						}
						if (b[i].g.xi[ip] > xE_box_polygon) {
							xE_box_polygon = b[i].g.xi[ip];
						}
						if (b[i].g.zi[ip] < zS_box_polygon) {
							zS_box_polygon = b[i].g.zi[ip];
						}
						if (b[i].g.zi[ip] > zE_box_polygon) {
							zE_box_polygon = b[i].g.zi[ip];
						}
						if (ip != 0) {
							if (fabs(b[i].g.yi[ip] - b[i].g.yi[0]) > 1.0e-36) {
								std::cout << "Error: NON planar polygon. body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								std::cout << "Diagnostic: fabs(b[" << i << "].g.yi[" << ip << "]- b[" << i << "].g.yi[0]) > 1.0e-36\n";
								system("pause");
								bOk = false;
							}
							if (fabs(b[i].g.hi[ip] - b[i].g.hi[0]) > 1.0e-36) {
								std::cout << "Error: NON planar polygon. body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								std::cout << "Diagnostic: fabs(b[" << i << "].g.hi[" << ip << "]- b[" << i << "].g.hi[0]) > 1.0e-36\n";
								system("pause");
								bOk = false;
							}
						}
						if (b[i].g.hi[ip] <= 0.0) {
							//printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							std::cout << "Error body[" << i << "].name = " << b[i].name << "  body[" << i << "].g.hi[" << ip << "]=" << b[i].g.hi[ip] << " is NEGATIV or ZERO." << std::endl;
							system("pause");
							bOk = false;
						}
						yS_box_polygon = b[i].g.yi[ip];
						yE_box_polygon = b[i].g.yi[ip] + b[i].g.hi[ip];
					}
					break;
				case YZ_PLANE:
					for (integer ip = 0; ip < b[i].g.nsizei; ip++) {


						// Детектирование ОШИБКИ самопересечения отрезков в полигоне.
						// Мы проверяем каждый с каждым
						it_c = 0;
						for (integer ip_t = 2; ip_t < b[i].g.nsizei; ip_t++) {
							i_t1 = ip + 1;
							if (i_t1 >= b[i].g.nsizei) i_t1 -= b[i].g.nsizei;
							i_t2 = ip + it_c + 2;
							if (i_t2 >= b[i].g.nsizei) i_t2 -= b[i].g.nsizei;
							i_t3 = ip + it_c + 3;
							if (i_t3 >= b[i].g.nsizei) i_t3 -= b[i].g.nsizei;

							// Мы попарно проверяем есть ли самопересечения каждый с каждым.
							if (b_is_intersect(b[i].g.yi[ip], b[i].g.zi[ip], b[i].g.yi[i_t1], b[i].g.zi[i_t1],
								b[i].g.yi[i_t2], b[i].g.zi[i_t2], b[i].g.yi[i_t3], b[i].g.zi[i_t3])) {
								// ОШИБКА: обнаружено самопересечение отрезков в полигоне.
								std::cout << "Error in POLYGON: Self-intersection of segments in body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								system("pause");
								bOk = false;
							}

							it_c++;
						}

						if (b[i].g.zi[ip] < zS_box_polygon) {
							zS_box_polygon = b[i].g.zi[ip];
						}
						if (b[i].g.zi[ip] > zE_box_polygon) {
							zE_box_polygon = b[i].g.zi[ip];
						}
						if (b[i].g.yi[ip] < yS_box_polygon) {
							yS_box_polygon = b[i].g.yi[ip];
						}
						if (b[i].g.yi[ip] > yE_box_polygon) {
							yE_box_polygon = b[i].g.yi[ip];
						}
						if (ip != 0) {
							if (fabs(b[i].g.xi[ip] - b[i].g.xi[0]) > 1.0e-36) {
								std::cout << "Error: NON planar polygon. body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								std::cout << "Diagnostic: fabs(b[" << i << "].g.xi[" << ip << "]- b[" << i << "].g.xi[0]) > 1.0e-36\n";
								system("pause");
								bOk = false;
							}
							if (fabs(b[i].g.hi[ip] - b[i].g.hi[0]) > 1.0e-36) {
								std::cout << "Error: NON planar polygon. body[" << i << "] body[" << i << "].name = " << b[i].name << " \n";
								std::cout << "Diagnostic: fabs(b[" << i << "].g.hi[" << ip << "]- b[" << i << "].g.hi[0]) > 1.0e-36\n";
								system("pause");
								bOk = false;
							}
						}
						if (b[i].g.hi[ip] <= 0.0) {
							//printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							std::cout << "Error body[" << i << "].name = " << b[i].name << " body[" << i << "].g.hi[" << ip << "]=" << b[i].g.hi[ip] << " is NEGATIV or ZERO." << std::endl;
							system("pause");
							bOk = false;
						}
						xS_box_polygon = b[i].g.xi[ip];
						xE_box_polygon = b[i].g.xi[ip] + b[i].g.hi[ip];
					}
					break;
				}

				if (xE_box_polygon > b[0].g.xE + dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] xE_box_polygon=%e > cabinet.xE=%e.\n", i, xE_box_polygon, b[0].g.xE);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] xE_box_polygon=" << xE_box_polygon << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (yE_box_polygon > b[0].g.yE + dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] yE_box_polygon=%e > cabinet.yE=%e.\n", i, yE_box_polygon, b[0].g.yE);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] yE_box_polygon=" << yE_box_polygon << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (zE_box_polygon > b[0].g.zE + dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] zE_box_polygon=%e > cabinet.zE=%e.\n", i, zE_box_polygon, b[0].g.zE);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] zE_box_polygon=" << zE_box_polygon << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (xS_box_polygon < b[0].g.xS - dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] xS_box_polygon=%e < cabinet.xS=%e.\n", i, xS_box_polygon, b[0].g.xS);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] xS_box_polygon=" << xS_box_polygon << " > cabinet.xS=" << b[0].g.xS << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (yS_box_polygon < b[0].g.yS - dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] yS_box_polygon=%e < cabinet.yS=%e.\n", i, yS_box_polygon, b[0].g.yS);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] yS_box_polygon=" << yS_box_polygon << " > cabinet.yS=" << b[0].g.yS << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (zS_box_polygon < b[0].g.zS - dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] zS_box_polygon=%e < cabinet.zS=%e.\n", i, zS_box_polygon, b[0].g.zS);
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] zS_box_polygon=" << zS_box_polygon << " > cabinet.zS=" << b[0].g.zS << "." << std::endl;
					system("pause");

					bOk = false;
				}
				doublereal vol_poly = Volume_polygon(b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2);
				if (vol_poly < 1.0e-30) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					std::cout << "body[" << i << "].name = " << b[i].name << " body[" << i << "] volume polygon=" << vol_poly << " is very small." << std::endl;
					system("pause");

					bOk = false;
				}


			}

		}
		if (b[i].g.itypegeom == CAD_STL) {
			// 15.11.2020

			if (i != 0) {
				if (b[i].g.root_CAD_STL == nullptr) {
					std::cout << "Error tringle faces CAD geom is empty: body[" << i << "].name = " << b[i].name << std::endl;
					system("pause");
					bOk = false;
				}
				else {
					integer icount_triangles = 0;

					gCAD_STL* tmp = b[i].g.root_CAD_STL;

					while (tmp != nullptr) {

						icount_triangles++;

						tmp = tmp->next;
					}

					if (icount_triangles < 4) {
						std::cout << "Error number tringles faces CAD geom  < 4: body[" << i << "].name = " << b[i].name << std::endl;
						system("pause");
						bOk = false;
					}

				}
			}
		}
	}

	for (integer i = 0; i < lw; ++i) {

		// Ошибка не число.
		if (w[i].g.xS != w[i].g.xS) {
			//printf("Error wall[%lld] xS NOT VALUE=%e.\n", i, w[i].g.xS);
			std::cout << "Error wall[" << i << "].name = " << w[i].name << " wall[" << i << "] xS NOT VALUE=" << w[i].g.xS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.xE != w[i].g.xE) {
			//printf("Error wall[%lld] xE NOT VALUE=%e.\n", i, w[i].g.xE);
			std::cout << "Error  wall[" << i << "].name = " << w[i].name << " wall[" << i << "] xE NOT VALUE=" << w[i].g.xE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.yS != w[i].g.yS) {
			//printf("Error wall[%lld] yS NOT VALUE=%e.\n", i, w[i].g.yS);
			std::cout << "Error  wall[" << i << "].name = " << w[i].name << " wall[" << i << "] yS NOT VALUE=" << w[i].g.yS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.yE != w[i].g.yE) {
			//printf("Error wall[%lld] yE NOT VALUE=%e.\n", i, w[i].g.yE);
			std::cout << "Error  wall[" << i << "].name = " << w[i].name << " wall[" << i << "] yE NOT VALUE=" << w[i].g.yE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.zS != w[i].g.zS) {
			//printf("Error wall[%lld] zS NOT VALUE=%e.\n", i, w[i].g.zS);
			std::cout << "Error  wall[" << i << "].name = " << w[i].name << " wall[" << i << "] zS NOT VALUE=" << w[i].g.zS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.zE != w[i].g.zE) {
			//printf("Error wall[%lld] zE NOT VALUE=%e.\n", i, w[i].g.zE);
			std::cout << "Error  wall[" << i << "].name = " << w[i].name << " wall[" << i << "] zE NOT VALUE=" << w[i].g.zE << "." << std::endl;
			system("pause");
			bOk = false;
		}

		switch (w[i].iPlane) {
		case XY_PLANE:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.xS >= w[i].g.xE) {
				//printf("wall[%lld].xS=%e >= wall[%lld].xE=%e\n", i, w[i].g.xS, i, w[i].g.xE);
				std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].xS=" << w[i].g.xS << " >= wall[" << i << "].xE=" << w[i].g.xE << std::endl;
				system("pause");
				doublereal temp = w[i].g.xE;
				w[i].g.xE = w[i].g.xS;
				w[i].g.xS = temp;
				bOk = false;
			}
			if (w[i].g.yS >= w[i].g.yE) {
				//printf("wall[%lld].yS=%e >= wall[%lld].yE=%e\n", i, w[i].g.yS, i, w[i].g.yE);
				std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].yS=" << w[i].g.yS << " >= wall[" << i << "].yE=" << w[i].g.yE << std::endl;
				system("pause");
				doublereal temp = w[i].g.yE;
				w[i].g.yE = w[i].g.yS;
				w[i].g.yS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.zS - w[i].g.zE) > 1.0e-36) {
				//printf("non union iso position on plane wall[%lld]: zS=%e zE=%e\n",i, w[i].g.zS, w[i].g.zE);
				std::cout << "non union iso position on plane   wall[" << i << "].name = " << w[i].name << " wall[" << i << "]: zS=" << w[i].g.zS << " zE=" << w[i].g.zE << std::endl;
				system("pause");
				w[i].g.zS = w[i].g.zE;
				bOk = false;
			}
			break;
		case XZ_PLANE:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.xS >= w[i].g.xE) {
				//printf("wall[%lld].xS=%e >= wall[%lld].xE=%e\n", i, w[i].g.xS, i, w[i].g.xE);
				std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].xS=" << w[i].g.xS << " >= wall[" << i << "].xE=" << w[i].g.xE << std::endl;
				system("pause");
				doublereal temp = w[i].g.xE;
				w[i].g.xE = w[i].g.xS;
				w[i].g.xS = temp;
				bOk = false;
			}
			if (w[i].g.zS >= w[i].g.zE) {
				//printf("wall[%lld].zS=%e >= wall[%lld].zE=%e\n", i, w[i].g.zS, i, w[i].g.zE);
				std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].zS=" << w[i].g.zS << " >= wall[" << i << "].zE=" << w[i].g.zE << std::endl;
				system("pause");
				doublereal temp = w[i].g.zE;
				w[i].g.zE = w[i].g.zS;
				w[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.yS - w[i].g.yE) > 1.0e-36) {
				//printf("non union iso position on plane wall[%lld]: yS=%e yE=%e\n", i, w[i].g.yS, w[i].g.yE);
				std::cout << "non union iso position on plane   wall[" << i << "].name = " << w[i].name << " wall[" << i << "]: yS=" << w[i].g.yS << " yE=" << w[i].g.yE << std::endl;
				system("pause");
				w[i].g.yS = w[i].g.yE;
				bOk = false;
			}
			break;
		case YZ_PLANE:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.yS >= w[i].g.yE) {
				//printf("wall[%lld].yS=%e >= wall[%lld].yE=%e\n", i, w[i].g.yS, i, w[i].g.yE);
				std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].yS=" << w[i].g.yS << " >= wall[" << i << "].yE=" << w[i].g.yE << std::endl;
				system("pause");
				doublereal temp = w[i].g.yE;
				w[i].g.yE = w[i].g.yS;
				w[i].g.yS = temp;
				bOk = false;
			}
			if (w[i].g.zS >= w[i].g.zE) {
				//printf("wall[%lld].zS=%e >= wall[%lld].zE=%e\n", i, w[i].g.zS, i, w[i].g.zE);
				std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].zS=" << w[i].g.zS << " >= wall[" << i << "].zE=" << w[i].g.zE << std::endl;
				system("pause");
				doublereal temp = w[i].g.zE;
				w[i].g.zE = w[i].g.zS;
				w[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.xS - w[i].g.xE) > 1.0e-36) {
				//printf("non union iso position on plane wall[%lld]: xS=%e xE=%e\n", i, w[i].g.xS, w[i].g.xE);
				std::cout << "non union iso position on plane   wall[" << i << "].name = " << w[i].name << " wall[" << i << "]: xS=" << w[i].g.xS << " xE=" << w[i].g.xE << std::endl;
				system("pause");
				w[i].g.xS = w[i].g.xE;
				bOk = false;
			}
			break;
		}
		if (w[i].g.xE > b[0].g.xE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.xE=%e > cabinet.xE=%e.\n", i, w[i].g.xE, b[0].g.xE);
			std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].g.xE=" << w[i].g.xE << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
			system("pause");
			w[i].g.xE = b[0].g.xE; // исправление.
			bOk = false;
		}
		if (w[i].g.yE > b[0].g.yE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.yE=%e > cabinet.yE=%e.\n", i, w[i].g.yE, b[0].g.yE);
			std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].g.yE=" << w[i].g.yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
			system("pause");
			w[i].g.yE = b[0].g.yE; // исправление.
			bOk = false;
		}
		if (w[i].g.zE > b[0].g.zE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.zE=%e > cabinet.zE=%e.\n", i, w[i].g.zE, b[0].g.zE);
			std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].g.zE=" << w[i].g.zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
			system("pause");
			w[i].g.zE = b[0].g.zE; // исправление.
			bOk = false;
		}
		if (w[i].g.xS < b[0].g.xS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.xS=%e < cabinet.xS=%e.\n", i, w[i].g.xS, b[0].g.xS);
			std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].g.xS=" << w[i].g.xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
			system("pause");
			w[i].g.xS = b[0].g.xS; // исправление.
			bOk = false;
		}
		if (w[i].g.yS < b[0].g.yS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.yS=%e < cabinet.yS=%e.\n", i, w[i].g.yS, b[0].g.yS);
			std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].g.yS=" << w[i].g.yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
			system("pause");
			w[i].g.yS = b[0].g.yS; // исправление.
			bOk = false;
		}
		if (w[i].g.zS < b[0].g.zS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.zS=%e < cabinet.zS=%e.\n", i, w[i].g.zS, b[0].g.zS);
			std::cout << "  wall[" << i << "].name = " << w[i].name << " wall[" << i << "].g.zS=" << w[i].g.zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
			system("pause");
			w[i].g.zS = b[0].g.zS; // исправление.
			bOk = false;
		}
	}

	for (integer i = 0; i < ls; ++i) {

		// Ошибка не число.
		if (s[i].g.xS != s[i].g.xS) {
			//printf("Error source[%lld] xS NOT VALUE=%e.\n", i, s[i].g.xS);
			std::cout << "Error   source[" << i << "].name = " << s[i].name << " source[" << i << "] xS NOT VALUE=" << s[i].g.xS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.xE != s[i].g.xE) {
			//printf("Error source[%lld] xE NOT VALUE=%e.\n", i, s[i].g.xE);
			std::cout << "Error   source[" << i << "].name = " << s[i].name << " source[" << i << "] xE NOT VALUE=" << s[i].g.xE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.yS != s[i].g.yS) {
			//printf("Error source[%lld] yS NOT VALUE=%e.\n", i, s[i].g.yS);
			std::cout << "Error   source[" << i << "].name = " << s[i].name << " source[" << i << "] yS NOT VALUE=" << s[i].g.yS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.yE != s[i].g.yE) {
			//printf("Error source[%lld] yE NOT VALUE=%e.\n", i, s[i].g.yE);
			std::cout << "Error   source[" << i << "].name = " << s[i].name << " source[" << i << "] yE NOT VALUE=" << s[i].g.yE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.zS != s[i].g.zS) {
			//printf("Error source[%lld] zS NOT VALUE=%e.\n", i, s[i].g.zS);
			std::cout << "Error   source[" << i << "].name = " << s[i].name << " source[" << i << "] zS NOT VALUE=" << s[i].g.zS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.zE != s[i].g.zE) {
			//printf("Error source[%lld] zE NOT VALUE=%e.\n", i, s[i].g.zE);
			std::cout << "Error   source[" << i << "].name = " << s[i].name << " source[" << i << "] zE NOT VALUE=" << s[i].g.zE << "." << std::endl;
			system("pause");
			bOk = false;
		}

		switch (s[i].iPlane) {
		case XY_PLANE:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.xS >= s[i].g.xE) {
				//printf("source[%lld].xS=%e >= source[%lld].xE=%e\n", i, s[i].g.xS, i, s[i].g.xE);
				std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].xS=" << s[i].g.xS << " >= source[" << i << "].xE=" << s[i].g.xE << std::endl;
				system("pause");
				doublereal temp = s[i].g.xE;
				s[i].g.xE = s[i].g.xS;
				s[i].g.xS = temp;
				bOk = false;
			}
			if (s[i].g.yS >= s[i].g.yE) {
				//printf("source[%lld].yS=%e >= source[%lld].yE=%e\n", i, s[i].g.yS, i, s[i].g.yE);
				std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].yS=" << s[i].g.yS << " >= source[" << i << "].yE=" << s[i].g.yE << std::endl;
				system("pause");
				doublereal temp = s[i].g.yE;
				s[i].g.yE = s[i].g.yS;
				s[i].g.yS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.zS - s[i].g.zE) > 1.0e-36) {
				//printf("non union iso position on plane source[%lld]: zS=%e zE=%e\n", i, s[i].g.zS, s[i].g.zE);
				std::cout << "non union iso position on plane   source[" << i << "].name = " << s[i].name << " source[" << i << "]: zS=" << s[i].g.zS << " zE=" << s[i].g.zE << std::endl;
				system("pause");
				s[i].g.zS = s[i].g.zE;
				bOk = false;
			}
			break;
		case XZ_PLANE:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.xS >= s[i].g.xE) {
				//printf("source[%lld].xS=%e >= source[%lld].xE=%e\n", i, s[i].g.xS, i, s[i].g.xE);
				std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].xS=" << s[i].g.xS << " >= source[" << i << "].xE=" << s[i].g.xE << std::endl;
				system("pause");
				doublereal temp = s[i].g.xE;
				s[i].g.xE = s[i].g.xS;
				s[i].g.xS = temp;
				bOk = false;
			}
			if (s[i].g.zS >= s[i].g.zE) {
				//printf("source[%lld].zS=%e >= source[%lld].zE=%e\n", i, s[i].g.zS, i, s[i].g.zE);
				std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].zS=" << s[i].g.zS << " >= source[" << i << "].zE=" << s[i].g.zE << std::endl;
				system("pause");
				doublereal temp = s[i].g.zE;
				s[i].g.zE = s[i].g.zS;
				s[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.yS - s[i].g.yE) > 1.0e-36) {
				//printf("non union iso position on plane source[%lld]: yS=%e yE=%e\n", i, s[i].g.yS, s[i].g.yE);
				std::cout << "non union iso position on plane   source[" << i << "].name = " << s[i].name << " source[" << i << "]: yS=" << s[i].g.yS << " yE=" << s[i].g.yE << std::endl;
				system("pause");
				s[i].g.yS = s[i].g.yE;
				bOk = false;
			}
			break;
		case YZ_PLANE:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.yS >= s[i].g.yE) {
				//printf("source[%lld].yS=%e >= source[%lld].yE=%e\n", i, s[i].g.yS, i, s[i].g.yE);
				std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].yS=" << s[i].g.yS << " >= source[" << i << "].yE=" << s[i].g.yE << std::endl;
				system("pause");
				doublereal temp = s[i].g.yE;
				s[i].g.yE = s[i].g.yS;
				s[i].g.yS = temp;
				bOk = false;
			}
			if (s[i].g.zS >= s[i].g.zE) {
				//printf("source[%lld].zS=%e >= source[%lld].zE=%e\n", i, s[i].g.zS, i, s[i].g.zE);
				std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].zS=" << s[i].g.zS << " >= source[" << i << "].zE=" << s[i].g.zE << std::endl;
				system("pause");
				doublereal temp = s[i].g.zE;
				s[i].g.zE = s[i].g.zS;
				s[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.xS - s[i].g.xE) > 1.0e-36) {
				//printf("non union iso position on plane source[%lld]: xS=%e xE=%e\n", i, s[i].g.xS, s[i].g.xE);
				std::cout << "non union iso position on plane   source[" << i << "].name = " << s[i].name << " source[" << i << "]: xS=" << s[i].g.xS << " xE=" << s[i].g.xE << std::endl;
				system("pause");
				s[i].g.xS = s[i].g.xE;
				bOk = false;
			}
			break;
		}
		if (s[i].g.xE > b[0].g.xE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.xE=%e > cabinet.xE=%e.\n", i, s[i].g.xE, b[0].g.xE);
			std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].g.xE=" << s[i].g.xE << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
			system("pause");
			s[i].g.xE = b[0].g.xE; // исправление.
			bOk = false;
		}
		if (s[i].g.yE > b[0].g.yE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.yE=%e > cabinet.yE=%e.\n", i, s[i].g.yE, b[0].g.yE);
			std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].g.yE=" << s[i].g.yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
			system("pause");
			s[i].g.yE = b[0].g.yE; // исправление.
			bOk = false;
		}
		if (s[i].g.zE > b[0].g.zE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.zE=%e > cabinet.zE=%e.\n", i, s[i].g.zE, b[0].g.zE);
			std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].g.zE=" << s[i].g.zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
			system("pause");
			s[i].g.zE = b[0].g.zE; // исправление.
			bOk = false;
		}
		if (s[i].g.xS < b[0].g.xS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.xS=%e < cabinet.xS=%e.\n", i, s[i].g.xS, b[0].g.xS);
			std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].g.xS=" << s[i].g.xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
			system("pause");
			s[i].g.xS = b[0].g.xS; // исправление.
			bOk = false;
		}
		if (s[i].g.yS < b[0].g.yS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.yS=%e < cabinet.yS=%e.\n", i, s[i].g.yS, b[0].g.yS);
			std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].g.yS=" << s[i].g.yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
			system("pause");
			s[i].g.yS = b[0].g.yS; // исправление.
			bOk = false;
		}
		if (s[i].g.zS < b[0].g.zS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.zS=%e < cabinet.zS=%e.\n", i, s[i].g.zS, b[0].g.zS);
			std::cout << "  source[" << i << "].name = " << s[i].name << " source[" << i << "].g.zS=" << s[i].g.zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
			system("pause");
			s[i].g.zS = b[0].g.zS; // исправление.
			bOk = false;
		}
	}


	if ((b_on_adaptive_local_refinement_mesh) && (lu > 0)) {
		// Не допустимо совместное использование
		// Адаптивных локально измельченных сеток и
		// асемблесов 24.07.2020.
		std::cout << "ERROR!!! ((b_on_adaptive_local_refinement_mesh) && (lu > 0))==true." << std::endl;
		std::cout << "dignostics: not compatible with user settings." << std::endl;
		system("PAUSE");
		exit(1);
	}

	if ((lu==0)||((lu > 0) && ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT) ||
		(((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE)||(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE)) && (eqin.itemper == 2))||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL)||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))))
	{
		// Всё впорядке настройки совместимы.
		// Асемблес может существовать при построении расчётной сетки,
		// а также может существовать для температурного солвера на основе метода конечных элементов МКЭ.
		// Либо асемблес может быть неактивен.
	}
	else {
		std::wcout << "not if ((lu > 0) && ((steady_or_unsteady_global_determinant == MESHER_ONLY) ||" << std::endl;
		std::wcout << "	((steady_or_unsteady_global_determinant == STEADY_TEMPERATURE) && (eqin.itemper == 2))))" << std::endl;
		std::cout << "dignostics: not compatible with user settings." << std::endl;
		system("PAUSE");
		exit(1);
	}

	for (integer i = 0; i < lu; ++i) {
		if (u[i].active) {

			if (u[i].xS != u[i].xS) {
				//printf("Error my_union[%lld] xS NOT VALUE=%e.\n", i, u[i].xS);
				std::cout << "Error:  my_union[" << i << "].name = " << "union" << i + 1 << "  "<< u[i].name  << "  my_union[" << i << "] xS NOT VALUE=" << u[i].xS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (u[i].xE != u[i].xE) {
				//printf("Error my_union[%lld] xE NOT VALUE=%e.\n", i, u[i].xE);
				std::cout << "Error:  my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "] xE NOT VALUE=" << u[i].xE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (u[i].yS != u[i].yS) {
				//printf("Error my_union[%lld] yS NOT VALUE=%e.\n", i, u[i].yS);
				std::cout << "Error:  my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "] yS NOT VALUE=" << u[i].yS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (u[i].yE != u[i].yE) {
				//printf("Error my_union[%lld] yE NOT VALUE=%e.\n", i, u[i].yE);
				std::cout << "Error:  my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "] yE NOT VALUE=" << u[i].yE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (u[i].zS != u[i].zS) {
				//printf("Error my_union[%lld] zS NOT VALUE=%e.\n", i, u[i].zS);
				std::cout << "Error:  my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "] zS NOT VALUE=" << u[i].zS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (u[i].zE != u[i].zE) {
				//printf("Error my_union[%lld] zE NOT VALUE=%e.\n", i, u[i].zE);
				std::cout << "Error:  my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "] zE NOT VALUE=" << u[i].zE << "." << std::endl;
				system("pause");
				bOk = false;
			}

			// Ошибка порядка следования Start должен быть строго меньше End.
			if (u[i].xS >= u[i].xE) {
				//printf("my_union[%lld].xS=%e >= my_union[%lld].xE=%e\n", i, u[i].xS, i, u[i].xE);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].xS=" << u[i].xS << " >= my_union[" << i << "].xE=" << u[i].xE << std::endl;
				system("pause");
				doublereal temp = u[i].xE;
				u[i].xE = u[i].xS;
				u[i].xS = temp;
				bOk = false;
			}

			if (u[i].yS >= u[i].yE) {
				//printf("my_union[%lld].yS=%e >= my_union[%lld].yE=%e\n", i, u[i].yS, i, u[i].yE);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].yS=" << u[i].yS << " >= my_union[" << i << "].yE=" << u[i].yE << std::endl;
				system("pause");
				doublereal temp = u[i].yE;
				u[i].yE = u[i].yS;
				u[i].yS = temp;
				bOk = false;
			}

			if (u[i].zS >= u[i].zE) {
				//printf("my_union[%lld].zS=%e >= my_union[%lld].zE=%e\n", i, u[i].zS, i, u[i].zE);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].zS=" << u[i].zS << " >= my_union[" << i << "].zE=" << u[i].zE << std::endl;
				system("pause");
				doublereal temp = u[i].zE;
				u[i].zE = u[i].zS;
				u[i].zS = temp;
				bOk = false;
			}

			// Асемблес не вылазит за границы кабинета
			if (u[i].xE > b[0].g.xE) {
				printf("ERROR. Your model is incorrect.\n");
				//printf("my_union[%lld].xE=%e > cabinet.xE=%e.\n", i, u[i].xE, b[0].g.xE);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].xE=" << u[i].xE << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
				system("pause");
				u[i].xE = b[0].g.xE; // исправление.
				bOk = false;
			}

			if (u[i].yE > b[0].g.yE) {
				printf("ERROR. Your model is incorrect.\n");
				//printf("my_union[%lld].yE=%e > cabinet.yE=%e.\n", i, u[i].yE, b[0].g.yE);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].yE=" << u[i].yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
				system("pause");
				u[i].yE = b[0].g.yE; // исправление.
				bOk = false;
			}

			if (u[i].zE > b[0].g.zE) {
				printf("ERROR. Your model is incorrect.\n");
				//printf("my_union[%lld].zE=%e > cabinet.zE=%e.\n", i, u[i].zE, b[0].g.zE);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].zE=" << u[i].zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
				system("pause");
				u[i].zE = b[0].g.zE; // исправление.
				bOk = false;
			}

			if (u[i].xS < b[0].g.xS) {
				printf("ERROR. Your model is incorrect.\n");
				//printf("my_union[%lld].xS=%e < cabinet.xS=%e.\n", i, u[i].xS, b[0].g.xS);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].xS=" << u[i].xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
				system("pause");
				u[i].xS = b[0].g.xS; // исправление.
				bOk = false;
			}

			if (u[i].yS < b[0].g.yS) {
				printf("ERROR. Your model is incorrect.\n");
				//printf("my_union[%lld].yS=%e < cabinet.yS=%e.\n", i, u[i].yS, b[0].g.yS);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].yS=" << u[i].yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
				system("pause");
				u[i].yS = b[0].g.yS; // исправление.
				bOk = false;
			}

			if (u[i].zS < b[0].g.zS) {
				printf("ERROR. Your model is incorrect.\n");
				//printf("my_union[%lld].zS=%e < cabinet.zS=%e.\n", i, u[i].zS, b[0].g.zS);
				std::cout << "my_union[" << i << "].name = " << "union" << i + 1 << "  " << u[i].name << " my_union[" << i << "].zS=" << u[i].zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
				system("pause");
				u[i].zS = b[0].g.zS; // исправление.
				bOk = false;
			}
		}
	}


	if (bOk) {
		printf("GEOMETRY CHECKED. GEOMETRY IS CORRECT.\n");
	}
	else {
		printf("GEOMETRY IS INCORRECT. YOUR MODEL IS INCORRECT.\n");
		system("pause");
		exit(1);
	}
} // body_check()


// считывание параметров из 
// входного файла premeshin.txt
void premeshin_old(const char *fname, integer &lmatmax, int &lb, int &ls, int &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w, 
	           doublereal &dgx, doublereal &dgy, doublereal &dgz, int &inx, int &iny, int &inz, doublereal &operatingtemperature, 
			   integer &ltdp, TEMP_DEP_POWER* &gtdps, int &lu, UNION* &my_union) {


	


	// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
	// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
	// при повторном обращении идет обращение к ячейки хеш-таблицы.
	// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
	// Хеш-таблицы для automatic
	// Инициализация хеш-таблицы.
	for (integer i_1 = 0; i_1 < isize_shorter_hash; ++i_1) {
		bshorter_hash_X[i_1] = false;
		bshorter_hash_Y[i_1] = false;
		bshorter_hash_Z[i_1] = false;
		shorter_hash_X[i_1] = 1.0e-10;
		shorter_hash_Y[i_1] = 1.0e-10;
		shorter_hash_Z[i_1] = 1.0e-10;
	}

	doublereal dmult = 1.0 / rdivision_interval;

#ifdef MINGW_COMPILLER

	

	// eqin - информация о наборе решаемых уравнений.

	// dgx, dgy, dgz - вектор силы тяжести.
	// inx, iny, inz - количество точек по каждой из осей.	

	FILE *fp;
	int err1=0;
	fp=fopen64(fname, "r");

	if (fp != nullptr) {
		err1 = 0;
	}
	else {
		err1 = 1; // ошибка открытия.
	}
	if (err1 != 0) {
		printf("Error!!! No input File premeshin.txt.\n");
		printf("You must use the graphical user interface\n");
		printf("AliceMesh_v0_45.exe in Delphi, which will \n");
		printf("prepare the model and write the premeshin.txt\n");
		printf("file in the required format.\n");
		//system("PAUSE");
		system("pause");
		// Если файла premeshin.txt нет то мы не можем ничего обрабатывать
		// т.к. элементарно отсутствуют входные данные. Поэтому мы выходим из 
		// приложения.
		exit(1);
	}
	else
	{
		if (fp != nullptr) {
			float fin = 0.0;
			int din = 0;
			doublereal scale = 1.0;
			doublereal scale_all = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
			case 1: ionly_solid_visible= WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE; break;
			default: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
			}
			//ionly_solid_visible=din;
			switch (ionly_solid_visible) {
			case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE:
				std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
				break;
			case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE:
				std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
				break;
			}
			fscanf(fp, "%f", &fin);
			scale = fin;
			scale_all = 1.0f / scale;
			fscanf(fp, "%d", &din);
			lmatmax = din;
			fscanf(fp, "%d", &din);
			lb = din;
			fscanf(fp, "%d", &din);
			ls = din;
			fscanf(fp, "%d", &din);
			lw = din;
			fscanf(fp, "%d", &din);
			ltdp = din; // количество уникальных данных с табличными данными по зависимости расеиваемой мощности от температуры.
			
						// Считываем значение вектора силы тяжести:
			fscanf(fp, "%f", &fin);
			dgx = fin;
			fscanf(fp, "%f", &fin);
			dgy = fin;
			fscanf(fp, "%f", &fin);
			dgz = fin;

			// считываем количество точек на каждой координатной оси
			fscanf(fp, "%d", &din);
			inx = din;
			fscanf(fp, "%d", &din);
			iny = din;
			fscanf(fp, "%d", &din);
			inz = din;

			fscanf(fp, "%f", &fin);
			operatingtemperature = fin; // Operating Temperature
			operating_temperature_for_film_coeff=fin;


			// инициализация компонент скорости константой.
			// Единые значения для всей расчётной области.
			// initialization value.
			fscanf(fp, "%f", &fin);
			starting_speed_Vx = fin;
			fscanf(fp, "%f", &fin);
			starting_speed_Vy = fin;
			fscanf(fp, "%f", &fin);
			starting_speed_Vz = fin;

			// Считываем координаты опорной точки через которую проходит пользовательская линия Variation Plot.
			fscanf(fp, "%f", &fin);
			Tochka_position_X0_for_XY_Plot = scale*fin;
			fscanf(fp, "%f", &fin);
			Tochka_position_Y0_for_XY_Plot = scale*fin;
			fscanf(fp, "%f", &fin);
			Tochka_position_Z0_for_XY_Plot = scale*fin;
			// Направление линии совпадает с направлением 
			// одной из осей декартовой прямоугольной системы координат:
			// 0 - Ox, 1 - Oy, 2 - Oz.
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
			case 1: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
			case 2: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
			default: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
			}
			//idirectional_for_XY_Plot = din;

			fscanf(fp, "%f", &fin);
			etalon_max_size_ratio = fin; // подробность расчётной сетки.

			fscanf(fp, "%f", &fin);
			etalon_max_size_ratio2 = fin; // Критерий качества расчётной сетки на основе FlowVision.

			// 0.	none
            // 1.	Snap to grid
            // 2.	Snap to grid ALICE
            // 3.	Snap to grid ++
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: bsnap_TO_global = 0;  break;
			case 1: bsnap_TO_global = 1;  break;
			case 2: bsnap_TO_global = 2;  break;
			case 3: bsnap_TO_global = 3;  break;
			default: bsnap_TO_global = 1;  break;
			}
			

			fscanf(fp, "%d", &din);
			iswitchsolveramg_vs_BiCGstab_plus_ILU2 = din; // Выбор решающего устройства: либо amg1r5 либо BiCGStab+ILU2.

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
			case 1: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER; break;
			case 2: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::CAMG_RUMBA_v0_14_SECOND_T_SOLVER; break;
			case 3: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMG1R5_SECOND_T_SOLVER; break;
			case 4: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMGCL_SECONT_T_SOLVER; break;
			default: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
			}
			//iswitchsolveramg_vs_BiCGstab_plus_ILU6 = din; // Выбор решающего устройства: либо РУМБА0.14 либо BiCGStab+ILU6.

			fscanf(fp, "%d", &din);
			if (din == 1) {
				// SIMPLEC algorithm.
				iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby;
			}
			else {
				// SIMPLE algorithm 1972.
				iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;
			}

			fscanf(fp, "%d", &din);
			
				switch (din) {
				case 0:  iprefix_Scheme_Flow = CR;
					break;
				case 1:  iprefix_Scheme_Flow = UDS;
					break;
				case 2:  iprefix_Scheme_Flow = COMB;
					break;
				case 3:  iprefix_Scheme_Flow = POLY;
					break;
				case 4:  iprefix_Scheme_Flow = EXP;
					break;
				case 5:  iprefix_Scheme_Flow = BULG;
					break;
				case 6:  iprefix_Scheme_Flow = POW;
					break;
				default:
					iprefix_Scheme_Flow = UDS;
					break;
				}
			

			fscanf(fp, "%d", &din);
			// выбор схемы для потока жидкости.
			// Внимание эти определения должны полностью соответствовать 
			// определениям в файле my_approx_convective2.c
			switch (din) {
			case 1: iFLOWScheme = UNEVEN_MUSCL; break; // MUSCL 2
			case 2: iFLOWScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
			case 3: iFLOWScheme = UNEVEN_HLPA; break; // HLPA 2
			case 4: iFLOWScheme = UNEVEN_SMART; break; // SMART 3
			case 5: iFLOWScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
			case 6: iFLOWScheme = UNEVEN_SMARTER; break; // SMARTER 3
			case 7: iFLOWScheme = UNEVEN_LPPA; break; // LPPA 3
			case 8: iFLOWScheme = UNEVEN_VONOS; break; // VONOS 3
			case 9: iFLOWScheme = UNEVEN_STOIC; break; // STOIC
			case 10: iFLOWScheme = UNEVEN_CLAM; break; // CLAM
			case 11: iFLOWScheme = UNEVEN_OSHER; break; // OSHER
			case 12: iFLOWScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
			case 13: iFLOWScheme = UNEVEN_SUPER_C; break; // SUPER_C
			case 14: iFLOWScheme = UNEVEN_ISNAS; break; // ISNAS
			case 15: iFLOWScheme = UNEVEN_CUBISTA; break; // CUBISTA
			case 16: iFLOWScheme = UNEVEN_GAMMA; break; // GAMMA
			case 17: iFLOWScheme = UNEVEN_COPLA; break; // COPLA
			case 18: iFLOWScheme = UNEVEN_SECBC; break; // SECBC
			case 19: iFLOWScheme = UNEVEN_SGSD; break; // SGSD
			case 20: iFLOWScheme = UNEVEN_WENO5; break; // WENO5
			default: iFLOWScheme = UDS; break; // UDS самая стабильная схема.
			}

			fscanf(fp, "%d", &din);
			// выбор схемы для температуры в потоке жидкости.
			// Внимание эти определения должны полностью соответствовать 
			// определениям в файле my_approx_convective2.c
			switch (din) {
			case 1: iTEMPScheme = UNEVEN_MUSCL; break; // MUSCL 2
			case 2: iTEMPScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
			case 3: iTEMPScheme = UNEVEN_HLPA; break; // HLPA 2
			case 4: iTEMPScheme = UNEVEN_SMART; break; // SMART 3
			case 5: iTEMPScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
			case 6: iTEMPScheme = UNEVEN_SMARTER; break; // SMARTER 3
			case 7: iTEMPScheme = UNEVEN_LPPA; break; // LPPA 3
			case 8: iTEMPScheme = UNEVEN_VONOS; break; // VONOS 3
			case 9: iTEMPScheme = UNEVEN_STOIC; break; // STOIC
			case 10: iTEMPScheme = UNEVEN_CLAM; break; // CLAM
			case 11: iTEMPScheme = UNEVEN_OSHER; break; // OSHER
			case 12: iTEMPScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
			case 13: iTEMPScheme = UNEVEN_SUPER_C; break; // SUPER_C
			case 14: iTEMPScheme = UNEVEN_ISNAS; break; // ISNAS
			case 15: iTEMPScheme = UNEVEN_CUBISTA; break; // CUBISTA
			case 16: iTEMPScheme = UNEVEN_GAMMA; break; // GAMMA
			case 17: iTEMPScheme = UNEVEN_COPLA; break; // COPLA
			case 18: iTEMPScheme = UNEVEN_SECBC; break; // SECBC
			case 19: iTEMPScheme = UNEVEN_SGSD; break; // SGSD
			case 20: iTEMPScheme = UNEVEN_WENO5; break; // WENO5
			default: iTEMPScheme = UDS; break; // UDS самая стабильная схема.
			}


			// Выбор сеточного генератора.
			fscanf(fp, "%d", &din);
			switch (din) {
				case 0: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; break;
				case 1: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER; break;
				case 2: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
				default: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
			}
			//iswitchMeshGenerator = din;


			fscanf(fp, "%d", &din);
			steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY;
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 5) || (din == 6)|| 
				(din == 7) || (din == 8) || (din == 9)||(din == 10)||(din == 11)||(din==12)||(din==13)) {
				// 0 - thermal only steady state calculation,
				// 1 - thermal only unsteady calculation,
				// 2 - mesh generator only.
				// 3 - fluid dynamic steady state.
				// 5 - Static Structural (Thermal solver #2)
				// 6 - Thermal Stress
				// 7 - Unsteady thermal solver #2
				// 8 - Visualisation only
				// 9 - cfd unsteady fluid dynamic.
				// 10 - NETWORK_T Графовый метод решения уравнения теплопроводности.
				// 11 - UNSTEADY NETWORK_T Нестационарный графовый метод решения уравнения теплопроводности.
				// 12 - Нестационарная механика,
				// 13 - Нестационарная механика совместно с нестационарной теплопередачей.
				switch(din) {
				case 0: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE; break;
				case 1: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE; break;
				case 2: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
				case 3: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_STEADY; break;
				case 5: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL; break;
				case 6: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
				case 7: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER; break;
				case 8: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT; break;
				case 9: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY; break;
				case 10: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T;  break;
				case 11: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY; break;
				case 12: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL; break;
				case 13: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
				default: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
				}
				//steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady - 1 calculation.
			}
			else {
				printf("error input parametr steady or unsteady calculation\n");
				system("PAUSE");
				exit(1);
			}

			fscanf(fp, "%d", &din);
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3)||(din == 4)) {
				switch (din) {
				case 0: glTSL.id_law = TIME_STEP_lAW_SELECTOR::LINEAR;
					break;
				case 1: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE;
					break;
				case 2: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE2;
					break;
				case 3: glTSL.id_law = TIME_STEP_lAW_SELECTOR::HOT_COLD;
					break;
				case 4: glTSL.id_law = TIME_STEP_lAW_SELECTOR::PIECEWISE_CONSTANT;
					break;
				}
			}
			else {
				printf("error input parametr timestep law\n");
				system("PAUSE");
				exit(1);
			}
			fscanf(fp, "%f", &fin);
			glTSL.Factor_a_for_Linear = fin; // Factor_a
			if ((fin < 0.0)||(fin>=1.0)) {
				printf("error input parametr timestep law Factor a\n");
				system("PAUSE");
				exit(1);
            }
			fscanf(fp, "%f", &fin);
			if (fin< 0.0) {
				printf("error input parametr timestep law tau must be strongly positive\n");
				system("PAUSE");
				exit(1);
            }
			glTSL.tau = fin; // длительность импульса.
			fscanf(fp, "%f", &fin);
			glTSL.Q = fin;  // Скважность.
			// Параметры импульсного режима для SquareWave 2 режима.
			fscanf(fp, "%f", &fin);
			if ((fin< 0.0)||(fin>=1.0)) {
				printf("error input parametr timestep law SquareWave2 multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf(fp, "%f", &fin);
			if (fin< 0.0) {
				printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;			
			fscanf(fp, "%f", &fin);
			if (fin< 0.0) {
				printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau2 = fin;
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau_pause = fin;

			fscanf(fp, "%f", &fin);
			if ((fin < 0.0)||(fin>1.0)) {
				printf("error input parametr timestep law SquareWave2 off_multiplyer must be [0..1]\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.off_multiplyer = fin;

			fscanf(fp, "%d", &din);
			glTSL.n_cycle = din;
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law SquareWave2 Period must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.T_all = fin;
			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl <= 0.0) {
				printf("error in parameters Square Wave SquareWave2 time step law.\n");
				system("PAUSE");
				exit(1);
			}

			fscanf(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr on_time_double_linear law hot cold reshime must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.on_time_double_linear = fin;


			// Время окончания нестационарного моделирования при расчёте теплопередачи в твёрдом теле.
			fscanf(fp, "%f", &fin);
			globalEndTimeUnsteadyTemperatureCalculation = fin;

			// Newton-Richman condition.
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
            case 1: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC; break;
	        case 2: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC; break;
		    case 3: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC; break;
			default :adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
			}
			//adiabatic_vs_heat_transfer_coeff = din;  // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
			fscanf(fp, "%f", &fin);
			film_coefficient = fin;
			// AЛИС сетка
			fscanf(fp, "%d", &din);
			if (din==0) {
				// Обычная структурированная сетка.
				b_on_adaptive_local_refinement_mesh = false;
			}
			else {
				// АЛИС
				b_on_adaptive_local_refinement_mesh = true;
			}
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
			case 1: itype_ALICE_Mesh= TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH; break;
			default: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
			}
			//itype_ALICE_Mesh=din;
			fscanf(fp, "%d", &din);
			my_amg_manager.m_restart = static_cast<int>(din);
			// classical algebraic multigrid parameters:
			// only for my_agregat_amg.cu.
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
				break;
			case 1: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::QUICK_SORT;
				break;
			case 2: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::HEAP_SORT;
				break;
			case 3: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::TIM_SORT;
				break;
			default:
				my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
				break;
			}
			fscanf(fp, "%d", &din);
			//my_amg_manager.maximum_levels = din;
			my_amg_manager.maximum_delete_levels_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Stress = din;

			// type interpolation procedure:
			//fscanf(fp, "%d", &din);
			//my_amg_manager.number_interpolation_procedure = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Stress = din;

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}


			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}


			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}


			fscanf(fp, "%d", &din);
			//my_amg_manager.itypemodifyinterpol = din;
			//my_amg_manager.baglomeration_with_consistency_scaling = din;
			my_amg_manager.bdiagonal_dominant = din;
			fscanf(fp, "%d", &din);
			//my_amg_manager.inumberadaptpass = din;

			// 23.02.2018
			// print matrix portrait
			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bTemperatureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bTemperatureMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}
			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bSpeedMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bSpeedMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}
			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bPressureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bPressureMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}
			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bStressMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bStressMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}

			// 01.05.2017
			// truncation of interpolation:
			fscanf(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Stress = din;
			fscanf(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Temperature = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Speed = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Pressure = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Stress = fin;

			// number nFinnest sweeps:
			fscanf(fp, "%d", &din);
			my_amg_manager.nFinnest_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nFinnest_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nFinnest_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nFinnest_Stress = din;

			// number presweeps:
			fscanf(fp, "%d", &din);
			my_amg_manager.nu1_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu1_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu1_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu1_Stress = din;

			// number postsweeps:
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Stress = din;

			// memory size:
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Stress = din;


			// Параметр верхней релаксации в сглаживателе.
			//fscanf(fp, "%f", &fin);
			//my_amg_manager.gold_const_Temperature = fin;
			//fscanf(fp, "%f", &fin);
			//my_amg_manager.gold_const_Speed = fin;
			//fscanf(fp, "%f", &fin);
			//my_amg_manager.gold_const_Pressure = fin;
			//fscanf(fp, "%f", &fin);
			//my_amg_manager.gold_const_Stress = fin;

			// использовать ли ilu2 smoother.
			fscanf(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Temperature = din;
			
			/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Temperature = din;
				my_amg_manager.b_gmresTemp = true;
			}
			*/
			
			fscanf(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Speed = din;
			/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Speed = din;
				my_amg_manager.b_gmresSpeed = true;
			}
			*/
			

			fscanf(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Pressure = din;
/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Pressure = din;
				my_amg_manager.b_gmresPressure = true;
			}
			*/
					
			fscanf(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Stress = din;
/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Stress = din;
				my_amg_manager.b_gmresStress = true;
			}
	*/		
			my_amg_manager.gold_const_Temperature = return_gold_const(my_amg_manager.ilu2_smoother_Temperature);
			my_amg_manager.gold_const_Speed = return_gold_const(my_amg_manager.ilu2_smoother_Speed);
			my_amg_manager.gold_const_Pressure = return_gold_const(my_amg_manager.ilu2_smoother_Pressure);
			my_amg_manager.gold_const_Stress = return_gold_const(my_amg_manager.ilu2_smoother_Stress);


			fscanf(fp, "%d", &din);
			my_amg_manager.Chebyshev_degree = din;

			if (!((din == 5) || (din == 8) || (din == 16) || (din == 25) || (din == 32) || (din == 64) || (din == 128) || (din == 256) || (din == 512))) {
				std::cout << "Chebyshev degree undefined. Chebyshev degree mast be equal: 5, 8, 16, 25, 32, 64, 128, 256, 512 in this programm. Chebyshev degree == " << din << std::endl;
				system("pause");
				my_amg_manager.Chebyshev_degree = 5;
			}

			// strength threshold:
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Temperature = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Speed = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Pressure = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Stress = fin;

			// magic threshold:
			//fscanf(fp, "%f", &fin);
			//my_amg_manager.magic = fin;
			// magic <=> F_to_F
			fscanf(fp, "%f", &fin);
			my_amg_manager.F_to_F_Temperature = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.F_to_F_Speed = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.F_to_F_Pressure = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.F_to_F_Stress = fin;

			// AMG Splitting (coarsening)
			// Способ построения C-F разбиения: 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
			// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}
			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}

			// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
			// Если din==1, Stabilization BiCGStab.
			// 8.01.2017 Метод Хенка ван дер Ворста BiCGStab 1992. 
			// предобусловленный алгебраическим многосеточным методом.
			// 9.01.2018 Если din==2, FGMRes предобусловленный алгебраическим многосеточным методом.
			fscanf(fp, "%d", &din);
			my_amg_manager.istabilizationTemp = din; // 0 - none
			fscanf(fp, "%d", &din);
			my_amg_manager.istabilizationSpeed = din; // 0 - none
			fscanf(fp, "%d", &din);
			my_amg_manager.istabilizationPressure = din; // 0 - none
			fscanf(fp, "%d", &din);
			my_amg_manager.istabilizationStress = din; // 0 - none
			fscanf(fp, "%d", &din);
			my_amg_manager.ipatch_number = din; // 0 - патч не применяется.

			// Печать лога на консоль.
			fscanf(fp, "%d", &din);
			my_amg_manager.iprint_log_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.iprint_log_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.iprint_log_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.iprint_log_Stress = din;

			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Temperature_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Temperature_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}

			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Speed_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Speed_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}

			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Pressure_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Pressure_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}

			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Stress_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Stress_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}


			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Temperature = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Temperature = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}

			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Speed = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Speed = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}


			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Pressure = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Pressure = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}

			fscanf(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Stress = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Stress = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_smoother = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_selector = din;

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
				break;
			case 1: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::FGMRes;
				break;
			default: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
				break;
			}

			fscanf(fp, "%d", &din);
			my_amg_manager.lfil = din;

			fscanf(fp, "%d", &din);
			ireconstruction_free_construct_alloc = din;

			fscanf(fp, "%d", &din);
			ianimation_write_on = din;

			// выделение оперативной памяти.
			gtdps = new TEMP_DEP_POWER[ltdp];
			matlist = new TPROP[lmatmax];
			b = new BLOCK[lb];
			s = new SOURCE[ls];
			w = new WALL[lw+1];// +1 это гран условие по умолчанию, заглушка, чтобы избежать неопределенности.
			integer i = 0; // счётчик цикла for

			for (i = 0; i < ltdp; ++i) {
				// считывание имён файлов.
				gtdps[i].sname = new char[100]; // выделение памяти
				fscanf(fp, "%s", gtdps[i].sname, 100);
				//printf("%s",gtdps[i].sname);
				//system("PAUSE");
				// построение таблицы в памяти.
				my_read_power_table(gtdps[i].sname, gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
				printf("printeger table\n"); // debug
				my_print_table(gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
				printf("Please, press any key to continue...\n");
				//system("PAUSE");
				system("pause");

			}

		
			// считывание базы материалов
			for (i = 0; i < lmatmax; ++i) {
				// свойства материалов:
				// плотность
				fscanf(fp, "%f", &fin);
				matlist[i].rho = fin;
				// теплоёмкость при постоянном давлении
				//fscanf(fp, "%f", &fin);
				//matlist[i].cp = fin;
				fscanf(fp, "%d", &din);
				matlist[i].n_cp = din;
				matlist[i].arr_cp = nullptr;
				matlist[i].temp_cp = nullptr;
				matlist[i].arr_cp = new float[matlist[i].n_cp];
				matlist[i].temp_cp= new float[matlist[i].n_cp];
				if (matlist[i].temp_cp == nullptr) {
					printf("problem memory allocation for temp_cp\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_cp == nullptr) {
					printf("problem memory allocation for arr_cp\n");
					system("pause");
					exit(1);
				}				
				for (integer i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					matlist[i].temp_cp[i_4] = fin;
					fscanf(fp, "%f", &fin);
					matlist[i].arr_cp[i_4] = fin;
				}
				// теплопроводность
				//fscanf(fp, "%f", &fin);
				//matlist[i].lam = fin;
				fscanf(fp, "%d", &din);
				matlist[i].n_lam = din;
				matlist[i].arr_lam = nullptr;
				matlist[i].temp_lam = nullptr;
				matlist[i].arr_lam = new float[matlist[i].n_lam];
				matlist[i].temp_lam = new float[matlist[i].n_lam];
				if (matlist[i].temp_lam == nullptr) {
					printf("problem memory allocation for temp_lam\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_lam == nullptr) {
					printf("problem memory allocation for arr_lam\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					matlist[i].temp_lam[i_4] = fin;
					fscanf(fp, "%f", &fin);
					matlist[i].arr_lam[i_4] = fin;
				}
				// ортотропность теплопроводности:
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x=fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y=fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z=fin;


				// 28.08.2020
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x_beta_t_solid = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y_beta_t_solid = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z_beta_t_solid = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x_Young_Module = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y_Young_Module = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z_Young_Module = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_xy = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_xz = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_yz = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_yx = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_zx = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_zy = fin;
				fscanf(fp, "%d", &din);
				if (din == 1) {
					matlist[i].bActive_ShearModule = true;
				}
				else {
					matlist[i].bActive_ShearModule = false;
				}
				fscanf(fp, "%f", &fin);
				matlist[i].ShearModule_xy = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].ShearModule_yz = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].ShearModule_xz = fin;



				// 5.08.2017.
				// Коэффициенты для задачи упругости.
				// Модуль Юнга и коэффициент Пуассона.
				//doublereal Poissonratio = 0.154;
				//doublereal Youngmodule = 217.5e9;
				//fscanf(fp, "%f", &fin);
				//Poissonratio = fin;
				//if (fabs(fin) < 1.0e-30) {
					//printf("ERROR !!! Zero Poisson Ratio in model. \n");
					//printf("ERROR !!! Model is incorrect...\n");
					//system("PAUSE");
					//exit(1);
				//}
				//matlist[i].Poisson_ratio = fin; 

				fscanf(fp, "%d", &din);
				matlist[i].n_Poisson_ratio = din;
				matlist[i].arr_Poisson_ratio = nullptr;
				matlist[i].temp_Poisson_ratio = nullptr;
				matlist[i].arr_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
				matlist[i].temp_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
				if (matlist[i].temp_Poisson_ratio == nullptr) {
					printf("problem memory allocation for temp_Poisson_ratio\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_Poisson_ratio == nullptr) {
					printf("problem memory allocation for arr_Poisson_ratio\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_Poisson_ratio; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					matlist[i].temp_Poisson_ratio[i_4] = fin;
					fscanf(fp, "%f", &fin);
					matlist[i].arr_Poisson_ratio[i_4] = fin;
				}


				//fscanf(fp, "%f", &fin);
				//Youngmodule = fin*1e9;
				//if (fabs(fin) < 1.0e-30) {
					//printf("ERROR !!! Zero Young Module in model. \n");
					//printf("ERROR !!! Model is incorrect...\n");
					//system("PAUSE");
					//exit(1);
				//}
				//matlist[i].Young_Module= fin * 1e9; // Модуль Юнга

				fscanf(fp, "%d", &din);
				matlist[i].n_YoungModule = din;
				matlist[i].arr_Young_Module = nullptr;
				matlist[i].temp_Young_Module = nullptr;
				matlist[i].arr_Young_Module = new float[matlist[i].n_YoungModule];
				matlist[i].temp_Young_Module = new float[matlist[i].n_YoungModule];
				if (matlist[i].temp_Young_Module == nullptr) {
					printf("problem memory allocation for temp_Young_Module\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_Young_Module == nullptr) {
					printf("problem memory allocation for arr_Young_Module\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_YoungModule; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					matlist[i].temp_Young_Module[i_4] = fin;
					fscanf(fp, "%f", &fin);
					matlist[i].arr_Young_Module[i_4] = fin * 1.0e+9;
				}


				//fscanf(fp, "%f", &fin);
				//matlist[i].n_beta_t_solid = 1;
				//matlist[i].temp_beta_t_solid = new doublereal[1];
				//matlist[i].temp_beta_t_solid[0] = 25.0;
				//matlist[i].arr_beta_t_solid = new doublereal[1];
				//matlist[i].arr_beta_t_solid[0]= fin * 1E-6;

				fscanf(fp, "%d", &din);
				matlist[i].n_beta_t_solid = din;
				matlist[i].arr_beta_t_solid = nullptr;
				matlist[i].temp_beta_t_solid = nullptr;
				matlist[i].arr_beta_t_solid = new float[matlist[i].n_beta_t_solid];
				matlist[i].temp_beta_t_solid = new float[matlist[i].n_beta_t_solid];
				if (matlist[i].temp_beta_t_solid == nullptr) {
					printf("problem memory allocation for temp_beta_t_solid\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_beta_t_solid == nullptr) {
					printf("problem memory allocation for arr_beta_t_solid\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_beta_t_solid; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					matlist[i].temp_beta_t_solid[i_4] = fin;
					fscanf(fp, "%f", &fin);
					matlist[i].arr_beta_t_solid[i_4] = fin * 1.0e-6;
				}
			
				// Коэффициенты Ламе не используются в данном коде.
				
				// коэффициент динамической вязкости
				fscanf(fp, "%f", &fin);
				matlist[i].mu = fin;
				// коэффициент линейного температурного расширения
				fscanf(fp, "%f", &fin);
				matlist[i].beta_t = fin;
				// признак библиотечности материала
				fscanf(fp, "%d", &din);
				matlist[i].blibmat = din;
				// номер материала в библиотеке
				fscanf(fp, "%d", &din);
				matlist[i].ilibident = din;
				
				// для каждого КО данного блока,
				// если только он не перекрывается другими блоками
				// может быть использовано приближение 
				// Обербека-Буссинеска с соответствующей опорной температурой Tref.
				fscanf(fp, "%d", &din);
				switch (din) {
				case 0: matlist[i].bBussineskApproach = false; break;
				case 1: matlist[i].bBussineskApproach = true; break;
				default: matlist[i].bBussineskApproach = false; break;
				}
				// номер закона для зависимости динамической вязкости от напряжения сдвига
				fscanf(fp, "%d", &din);
				matlist[i].ilawmu = din;
				// минимальное значение динамической вязкости
				fscanf(fp, "%f", &fin);
				matlist[i].mumin = fin;
				// максимальное значение динамической вязкости
				fscanf(fp, "%f", &fin);
				matlist[i].mumax = fin;
				// параметры модельных законов для зависимости вязкости от напряжения сдвига
				fscanf(fp, "%f", &fin);
				matlist[i].Amu = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].Bmu = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].Cmu = fin;
				// показатель степени
				fscanf(fp, "%f", &fin);
				matlist[i].degreennmu = fin;

				// печать считанных значений на консоль
				//printf("%e %e %e %e %e\n", matlist[i].rho, matlist[i].arr_cp[0], matlist[i].arr_lam[0], matlist[i].mu, matlist[i].beta_t);
				//std::cout << matlist[i].rho << " " << matlist[i].arr_cp[0] << " " << matlist[i].arr_lam[0] << " " << matlist[i].mu << " " << matlist[i].beta_t << std::endl;
				//printf("%d %d %d\n", matlist[i].blibmat, matlist[i].ilibident, matlist[i].ilawmu); // bBoussinesq не печатается
				//printf("%e %e %e %e %e %e\n", matlist[i].mumin, matlist[i].mumax, matlist[i].Amu, matlist[i].Bmu, matlist[i].Cmu, matlist[i].degreennmu);
				//std::cout << matlist[i].mumin << " " << matlist[i].mumax << " " << matlist[i].Amu << " " << matlist[i].Bmu << " " << matlist[i].Cmu << " " << matlist[i].degreennmu << std::endl;

			}

			// считывание блоков
			for (i = 0; i<lb; ++i) {

				fscanf(fp, "%d", &din);
				b[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				// геометрия


				fscanf(fp, "%d", &din);
				CHECK_TYPE_GEOM(din);
				b[i].g.itypegeom = din;
				fscanf(fp, "%d", &din);
				if (din == 1) {
					b[i].bvisible = true;
				}
				else {
					b[i].bvisible = false;
				}


				fscanf(fp, "%f", &fin);
				b[i].g.xS = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.yS = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.zS = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.xE = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.yE = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.zE = scale*fin;
				// swap
				if (b[i].g.xS>b[i].g.xE) {
					dbuf = b[i].g.xS;
					b[i].g.xS = b[i].g.xE;
					b[i].g.xE = dbuf;
				}
				if (b[i].g.yS > b[i].g.yE) {
					dbuf = b[i].g.yS;
					b[i].g.yS = b[i].g.yE;
					b[i].g.yE = dbuf;
				}
				if (b[i].g.zS > b[i].g.zE) {
					dbuf = b[i].g.zS;
					b[i].g.zS = b[i].g.zE;
					b[i].g.zE = dbuf;
				}


				// Cylinder

				fscanf(fp, "%d", &din);
				b[i].g.iPlane = din;

				fscanf(fp, "%f", &fin);
				b[i].g.xC= scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.yC = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.zC = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.Hcyl = scale*fin;
				if (b[i].g.Hcyl < 0.0) {
					// ликвидируем отрицательную высоту цилиндра.
					switch (b[i].g.iPlane) {
					case XY_PLANE:
						b[i].g.zC += b[i].g.Hcyl;
						break;
					case XZ_PLANE:
						b[i].g.yC += b[i].g.Hcyl;
						break;
					case YZ_PLANE:
						b[i].g.xC += b[i].g.Hcyl;
						break;
					}
					b[i].g.Hcyl = fabs(b[i].g.Hcyl);
				}


				fscanf(fp, "%f", &fin);
				b[i].g.R_out_cyl = scale*fin;
				fscanf(fp, "%f", &fin);
				b[i].g.R_in_cyl = scale*fin;

				// Polygon
				fscanf(fp, "%d", &din);
				b[i].g.iPlane_obj2 = din;
				fscanf(fp, "%d", &din);
				b[i].g.nsizei = din;
				b[i].g.hi = new doublereal[b[i].g.nsizei];
				b[i].g.xi = new doublereal[b[i].g.nsizei];
				b[i].g.yi = new doublereal[b[i].g.nsizei];
				b[i].g.zi = new doublereal[b[i].g.nsizei];
				for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
					fscanf(fp, "%f", &fin);
					b[i].g.hi[i73] = scale*fin;
					fscanf(fp, "%f", &fin);
					b[i].g.xi[i73] = scale*fin;
					fscanf(fp, "%f", &fin);
					b[i].g.yi[i73] = scale*fin;
					fscanf(fp, "%f", &fin);
					b[i].g.zi[i73] = scale*fin;
					if (b[i].g.hi[i73] < 0.0) {
						// ликвидируем отрицательную высоту цилиндра.
						switch (b[i].g.iPlane_obj2) {
						case XY_PLANE:
							b[i].g.zi[i73] += b[i].g.hi[i73];
							break;
						case XZ_PLANE:
							b[i].g.yi[i73] += b[i].g.hi[i73];
							break;
						case YZ_PLANE:
							b[i].g.xi[i73] += b[i].g.hi[i73];
							break;
						}
						b[i].g.hi[i73] = fabs(b[i].g.hi[i73]);
					}
				}
				if ((b[i].g.itypegeom == POLYGON) && (b[i].g.nsizei > 0)) {
					// Мы заносим в окаймляющую призму границы полигона чтобы ускорить 
					// поиски в inmodel_temp быстро отсекая ячейки вне границ полигона.
					doublereal xmin53 = 1.0e30;
					doublereal ymin53 = 1.0e30;
					doublereal zmin53 = 1.0e30;
					doublereal xmax53 = -1.0e30;
					doublereal ymax53 = -1.0e30;
					doublereal zmax53 = -1.0e30;
					for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
						if (b[i].g.xi[i73] > xmax53) xmax53 = b[i].g.xi[i73];
						if (b[i].g.yi[i73] > ymax53) ymax53 = b[i].g.yi[i73];
						if (b[i].g.zi[i73] > zmax53) zmax53 = b[i].g.zi[i73];
						if (b[i].g.xi[i73] < xmin53) xmin53 = b[i].g.xi[i73];
						if (b[i].g.yi[i73] < ymin53) ymin53 = b[i].g.yi[i73];
						if (b[i].g.zi[i73] < zmin53) zmin53 = b[i].g.zi[i73];
					}
					switch (b[i].g.iPlane_obj2) {
					case XY_PLANE:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmin53 + b[i].g.hi[0];
						break;
					case XZ_PLANE:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymin53 + b[i].g.hi[0];
						break;
					case YZ_PLANE:
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.xS = xmin53;
						b[i].g.xE = xmin53 + b[i].g.hi[0];
						break;
					}
				}

				if (b[i].g.itypegeom == CAD_STL) {
					// CAD_STL

					// Считывание бинарного STL файла.
					b[i].g.ReadSTL_binary(lb, b[0].g);

				}

				// Ввод предполагается корректным.
				// emissivity
				fscanf(fp, "%f", &fin);
				b[i].radiation.emissW = fin;
				fscanf(fp, "%f", &fin);
				b[i].radiation.emissE = fin;
				fscanf(fp, "%f", &fin);
				b[i].radiation.emissS = fin;
				fscanf(fp, "%f", &fin);
				b[i].radiation.emissN = fin;
				fscanf(fp, "%f", &fin);
				b[i].radiation.emissB = fin;
				fscanf(fp, "%f", &fin);
				b[i].radiation.emissT = fin;
				fscanf(fp, "%d", &din);
				if (din==0) {
					// Блок не является вакуумным промежутком.
					b[i].radiation.binternalRadiation = false;
				}
				else {
					// блок является вакуумным промежутком.
					b[i].radiation.binternalRadiation = true;
					if (bvacuumPrism) {
						bdouble_vacuum_PRISM = true;
					}
					bvacuumPrism=true;
					// Вычисляем View Factors.
					calculate_view_factors(b[i]);
				}
				b[i].radiation.nodelistW = nullptr;
				b[i].radiation.nodelistE = nullptr;
				b[i].radiation.nodelistS = nullptr;
				b[i].radiation.nodelistN = nullptr;
				b[i].radiation.nodelistB = nullptr;
				b[i].radiation.nodelistT = nullptr;
				b[i].radiation.nodelistWsize = 0;
				b[i].radiation.nodelistEsize = 0;
				b[i].radiation.nodelistSsize = 0;
				b[i].radiation.nodelistNsize = 0;
				b[i].radiation.nodelistBsize = 0;
				b[i].radiation.nodelistTsize = 0;

				// идентификатор материала в базе материалов
				fscanf(fp, "%d", &din);
				b[i].imatid = din;

				fscanf(fp, "%d", &din);

				// bCylinderFixed
				if (din == 1) {
					b[i].CylinderFixed = true;
				}
				else {
					b[i].CylinderFixed = false;
				}

				// мощность тепловыделения
				//fscanf(fp, "%f", &fin);
				//b[i].Sc = fin;
				// 19 november 2016 температурно зависимая мощность тепловыделения.
				fscanf(fp, "%d", &din);
				b[i].n_Sc = din;
				b[i].arr_Sc = nullptr;
				b[i].temp_Sc = nullptr;
				b[i].arr_Sc = new doublereal[b[i].n_Sc];
				b[i].temp_Sc = new doublereal[b[i].n_Sc];
				if (b[i].temp_Sc == nullptr) {
					printf("problem memory allocation for temp_Sc\n");
					system("pause");
					exit(1);
				}
				if (b[i].arr_Sc == nullptr) {
					printf("problem memory allocation for arr_Sc\n");
					system("pause");
					exit(1);
				}
				// Объём полигона.
				doublereal vol_poly = Volume_polygon(b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2);
				for (integer i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					b[i].temp_Sc[i_4] = fin;
					fscanf(fp, "%f", &fin);
					if (fin != fin) {
						b[i].arr_Sc[i_4] = 0.0;
					}
					else {

						if (b[i].g.itypegeom == POLYGON) {

							// Для полигона передается из интерфейса просто мощность, а не удельная мощность.
							// Т.к. интерфейс не содержит функцию расчёта объёма полигона.
							// Для единообразия здесь мощность преобразуется в удельную мощность.
							if (vol_poly > 1.0e-30) {
								b[i].arr_Sc[i_4] = fin / vol_poly;
							}
							else {
								printf("error zero volume in polygon %s...\n", b[i].name);
								system("PAUSE");
								exit(1);
							}
						}
						else {
							b[i].arr_Sc[i_4] = fin;
						}
					}
				}

                // стиль зависимости мощности тепловыделения в блоке от времени.
				fscanf(fp, "%d", &din);
				switch (din) {
				case 0: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
					break;
				case 1: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE;
					break;
				case 2: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE2;
					break;
				case 3: b[i].ipower_time_depend = POWER_TIME_DEPEND::HOT_COLD;
					break;
				case 4: b[i].ipower_time_depend = POWER_TIME_DEPEND::PIECEWISE_CONST;
					break;
				default: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
					break;
				}

				// тип блока
				fscanf(fp, "%d", &din);
				switch (din) {
				case 1: b[i].itype = PHYSICS_TYPE_IN_BODY::SOLID;
					break;
				case 2: b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
					break;
				case 3: b[i].itype = PHYSICS_TYPE_IN_BODY::FLUID;
					break;
				default:
					b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
					break;
				}

				// печать считанных значений на консоль
				//printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				//std::cout << b[i].g.xS << " " << b[i].g.yS << " " << b[i].g.zS << " " << b[i].g.xE << " " << b[i].g.yE << " " << b[i].g.zE << std::endl;
				//printf("%d %d\n", b[i].imatid,  b[i].itype);
				//printf("temperature depend power\n");
				//printf("t_C power_W\n");
				for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					//printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
					//std::cout << b[i].temp_Sc[i_54] << " " << b[i].arr_Sc[i_54] << std::endl;
				}
			}

			// считывание источников тепла
			for (i = 0; i < ls; ++i) {

				fscanf(fp, "%d", &din);
				s[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				fscanf(fp, "%f", &fin);
				s[i].power = fin;
				fscanf(fp, "%d", &din);
				if (din == 0) {
					// const
					// задаётся постоянное значение рассеиваемой в тепло мощности.
					s[i].power_multiplyer = 1.0;
					s[i].bgarber_depend = false;
				}
				else if (din == 1) {
					// рассеиваемая в тепло мощность задаётся таблично
					// в зависимости от максимальной температуры и рабочего 
					// значения смещения стока.
					s[i].bgarber_depend = true;
					s[i].power_multiplyer = s[i].power;
					// мощность будет вычислена при ambient Temperature несколькими строчками позже.
				}
				fscanf(fp, "%d", &din);
				s[i].igarber_depend = din; // уникальный номер таблицы.
				fscanf(fp, "%f", &fin);
				s[i].roperation_offset_drain = fin; // рабочее значение смещения стока.
													//printf("offset drain is %e\n",s[i].roperation_offset_drain);
													//system("PAUSE");
				bool bsplinereadOk = true;
				if (s[i].bgarber_depend) {
					s[i].power = my_splain_interpol_power_table(gtdps[s[i].igarber_depend].intemp,
						gtdps[s[i].igarber_depend].inoffset_drain,
						gtdps[s[i].igarber_depend].rtemp,
						gtdps[s[i].igarber_depend].roffset_drain,
						gtdps[s[i].igarber_depend].rpower_table,
						operatingtemperature,
						s[i].roperation_offset_drain);
					if (bsplinereadOk) {
						// одиночная тестовая проверка сплайновой аппроксимации.
						printf("single test validation spline approximation...\n");
						//printf("calculate initial power=%e\n", s[i].power);
						std::cout << "calculate initial power=" << s[i].power << std::endl;
						printf("please, press any key to continue...");
						// system("PAUSE");
						system("pause");

						bsplinereadOk = false;
					}
					s[i].power *= s[i].power_multiplyer; // домножение на корректирующий множитель.
				}



				fscanf(fp, "%d", &din);
				s[i].iPlane = din;
				// геометрия
				fscanf(fp, "%f", &fin);
				s[i].g.xS = scale*fin;
				fscanf(fp, "%f", &fin);
				s[i].g.yS = scale*fin;
				fscanf(fp, "%f", &fin);
				s[i].g.zS = scale*fin;
				fscanf(fp, "%f", &fin);
				s[i].g.xE = scale*fin;
				fscanf(fp, "%f", &fin);
				s[i].g.yE = scale*fin;
				fscanf(fp, "%f", &fin);
				s[i].g.zE = scale*fin;



				// swap
				if (s[i].g.xS > s[i].g.xE) {
					dbuf = s[i].g.xS;
					s[i].g.xS = s[i].g.xE;
					s[i].g.xE = dbuf;
				}
				if (s[i].g.yS > s[i].g.yE) {
					dbuf = s[i].g.yS;
					s[i].g.yS = s[i].g.yE;
					s[i].g.yE = dbuf;
				}
				if (s[i].g.zS > s[i].g.zE) {
					dbuf = s[i].g.zS;
					s[i].g.zS = s[i].g.zE;
					s[i].g.zE = dbuf;
				}
				switch (s[i].iPlane) {
				case XY_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ_PLANE: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				//printf("%e %d %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
				//std::cout << s[i].power << " " << s[i].iPlane << " " << s[i].g.xS << " " << s[i].g.yS << " " << s[i].g.zS << " " << s[i].g.xE << " " << s[i].g.yE << " " << s[i].g.zE << " " << s[i].square << std::endl;
			}

			// считывание твёрдых стенок

			// 19.01.2021
		    // Стенка заглушка с условием прилипания по скорости(нулевой вектор скорости), 
		    // нулевым тепловым потоком по температуре, 
		    // свободной free границе в механике без приложенной силы.
			w[lw].iunion_id = 0; // cabinet
			w[lw].ifamily = WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY;
			w[lw].Tamb = 0.0;
			w[lw].emissivity = 0.0;
			w[lw].film_coefficient = 0.0;
			w[lw].ViewFactor = 1.0;
			w[lw].hf = 0.0;
			w[lw].bsymmetry = false;
			w[lw].bpressure = false;
			w[lw].bopening = false;
			w[lw].Vx = 0.0;
			w[lw].Vy = 0.0;
			w[lw].Vz = 0.0;
			w[lw].P = 0.0;
			w[lw].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
			w[lw].xForce = 0.0;
			w[lw].yForce = 0.0;
			w[lw].zForce = 0.0;
			w[lw].iPlane = XY_PLANE;
			w[lw].g.xS = 0.0;
			w[lw].g.yS = 0.0;
			w[lw].g.zS = 0.0;
			w[lw].g.xE = 0.0;
			w[lw].g.yE = 0.0;
			w[lw].g.zE = 0.0;

			for (i = 0; i < lw; ++i) {

				fscanf(fp, "%d", &din);
				w[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				fscanf(fp, "%d", &din);
				switch (din) {
				case 1 : w[i].ifamily= WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY; break;
				case 2: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
				case 3: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY; break;
				case 4: w[i].ifamily= WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY; break;
				default: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
				}
				//w[i].ifamily = din;
				switch (din) {
				case 1:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// termostability wall
					w[i].emissivity=0.0;
					w[i].film_coefficient=0.0;
					fscanf(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = 1.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // первого рода
				case 2:  fscanf(fp, "%f", &fin);
					w[i].Tamb = 0.0;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// adiabatic wall
					w[i].emissivity=0.0;
					w[i].film_coefficient=0.0;
					fscanf(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = 1.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // однородное условие Неймана
				case 3:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// Newton-Richman condition, film coefficient.
					w[i].emissivity=0.0;
					w[i].film_coefficient=fin;
					fscanf(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = 1.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Ньютон-Рихман.
				case 4:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// Stefan - Bolcman condition
					w[i].emissivity=fin;
					w[i].film_coefficient=0.0;
					fscanf(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = fmax(0.0,fmin(fin,1.0));
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Стефан-Больцман.
				default: break;
				}
				fscanf(fp, "%d", &din);
				if (din == 1) w[i].bsymmetry = true; else w[i].bsymmetry = false;
				fscanf(fp, "%d", &din);
				if (din == 1) w[i].bpressure = true; else w[i].bpressure = false;
				fscanf(fp, "%d", &din);
				if (din == 1) w[i].bopening = true; else w[i].bopening = false;

				fscanf(fp, "%f", &fin);
				w[i].Vx = fin;
				fscanf(fp, "%f", &fin);
				w[i].Vy = fin;
				fscanf(fp, "%f", &fin);
				w[i].Vz = fin;
				fscanf(fp, "%f", &fin);
				w[i].P = fin;
				fscanf(fp, "%d", &din);
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
				switch (din) {
				case 0: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
					break;
				case 1: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FIXIT;
					break;
				case 2: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FIXIT;
					break;
				case 3: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FIXIT;
					break;
				case 4: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XY_FIXIT;
					break;
				case 5: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XZ_FIXIT;
					break;
				case 6: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::YZ_FIXIT;
					break;
				case 7: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::ALL_FIXIT;
					break;
				case 8: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FORCE;
					break;
				case 9: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FORCE;
					break;
				case 10: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FORCE;
					break;
				default:
					w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // Free all
					break;
				}
				fscanf(fp, "%f", &fin);
				w[i].xForce = fin;
				fscanf(fp, "%f", &fin);
				w[i].yForce = fin;
				fscanf(fp, "%f", &fin);
				w[i].zForce = fin;
				fscanf(fp, "%d", &din);
				w[i].iPlane = din;
				// геометрия
				fscanf(fp, "%f", &fin);
				w[i].g.xS = scale*fin;
				fscanf(fp, "%f", &fin);
				w[i].g.yS = scale*fin;
				fscanf(fp, "%f", &fin);
				w[i].g.zS = scale*fin;
				fscanf(fp, "%f", &fin);
				w[i].g.xE = scale*fin;
				fscanf(fp, "%f", &fin);
				w[i].g.yE = scale*fin;
				fscanf(fp, "%f", &fin);
				w[i].g.zE = scale*fin;
				// swap
				if (w[i].g.xS > w[i].g.xE) {
					dbuf = w[i].g.xS;
					w[i].g.xS = w[i].g.xE;
					w[i].g.xE = dbuf;
				}
				if (w[i].g.yS > w[i].g.yE) {
					dbuf = w[i].g.yS;
					w[i].g.yS = w[i].g.yE;
					w[i].g.yE = dbuf;
				}
				if (w[i].g.zS > w[i].g.zE) {
					dbuf = w[i].g.zS;
					w[i].g.zS = w[i].g.zE;
					w[i].g.zE = dbuf;
				}
				//w[i].bfixboundary = false;// Свободная граница.
				//printf("%d %e %e %d %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
				//std::cout << w[i].ifamily << " " << w[i].Tamb << " " << w[i].hf << " " << w[i].iPlane << " " << w[i].g.xS << " " << w[i].g.yS << " " << w[i].g.zS << " " << w[i].g.xE << " " << w[i].g.yE << " " << w[i].g.zE << std::endl;
            }


			// АСЕМБЛЕСЫ.
			// Новое считывание работает только с новым форматом файла premeshin.txt 
			// с текстовыми метками, данный код здесь устарел. 31.10.2021.
			fscanf(fp, "%d", &din);
			lu = din;
			if (lu == 0) {
				my_union = nullptr;
			}
			else {
				my_union = new UNION[lu];
				// инициализация.
				for (i = 0; i < lu; ++i) {
					my_union[i].f = nullptr;
					my_union[i].xpos = nullptr;
					my_union[i].ypos = nullptr;
					my_union[i].zpos = nullptr;
					my_union[i].xposadd = nullptr;
					my_union[i].yposadd = nullptr;
					my_union[i].zposadd = nullptr;
					my_union[i].iunion_parent = -1;
					my_union[i].iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; // 2 - CoarseMeshGen
					my_union[i].inxadd = -1;
					my_union[i].inyadd = -1;
					my_union[i].inzadd = -1;
					my_union[i].flow_interior = 0;
					my_union[i].active = false;
				}
			}
			for (i = 0; i < lu; ++i) {

				fscanf(fp, "%d", &din);
				my_union[i].iunion_parent = din;

				fscanf(fp, "%f", &fin);
				my_union[i].xS = scale*fin;
				fscanf(fp, "%f", &fin);
				my_union[i].xE = scale*fin;
				fscanf(fp, "%f", &fin);
				my_union[i].yS = scale*fin;
				fscanf(fp, "%f", &fin);
				my_union[i].yE = scale*fin;
				fscanf(fp, "%f", &fin);
				my_union[i].zS = scale*fin;
				fscanf(fp, "%f", &fin);
				my_union[i].zE = scale*fin;

				fscanf(fp, "%d", &din);
				my_union[i].id = din; // Уникальный идентификатор АССЕМБЛЕСА.
				fscanf(fp, "%d", &din);
				my_union[i].inx = din;
				fscanf(fp, "%d", &din);
				my_union[i].iny = din;
				fscanf(fp, "%d", &din);
				my_union[i].inz = din;

			}

			// считывание информации о наборе решаемых уравнений
			fscanf(fp, "%d", &din);
			eqin.itemper = din;
			fscanf(fp, "%d", &din);
			eqin.imaxflD = din;
			if (eqin.imaxflD == 0) {
				eqin.fluidinfo = nullptr;
			}
			else
			{
				// выделение оперативной памяти
				if (eqin.fluidinfo != nullptr) {
					delete eqin.fluidinfo;
					eqin.fluidinfo = nullptr;
				}
				eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
				for (i = 0; i < eqin.imaxflD; ++i) {
					// Считывание координат опорной точки
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].xc = scale*fin;
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].yc = scale*fin;
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].zc = scale*fin;
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].iflow = din;
					fscanf(fp, "%d", &din);
					switch (din) {
					case 0: eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
						break;
					case 1:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::TURBULENT;
						break;
					default:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
						break;
					}
					fscanf(fp, "%d", &din);
					switch (din) {
					case 0: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
						break;
					case 1: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::SMAGORINSKY;
						break;
					case 2: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RNG_LES;
						break;
					case 3: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_SPALART_ALLMARES;
						break;
					case 4: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_MENTER_SST;
						break;
					case 5: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_STANDART_K_EPS;
						break;
					case 6: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_LANGTRY_MENTOR_SST;
						break;
					default: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
						break;
					}
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].Cs = fin; // постоянная Смагоринского.
					fscanf(fp, "%d", &din);
					// учёт динамического определения квадрата постоянной Смагоринского.
					// включает Dynamic Subgrid Scale Model Германо 1991 года.
					if (din == 1) {
						eqin.fluidinfo[i].bDynamic_Stress = true;
					}
					else {
						eqin.fluidinfo[i].bDynamic_Stress = false;
					}
					fscanf(fp, "%d", &din);
					// включает ограничение сверху и снизу на возможные значения постоянной Смагоринского.
					if (din == 1) {
						eqin.fluidinfo[i].bLimiters_Cs = true;
					}
					else {
						eqin.fluidinfo[i].bLimiters_Cs = false;
					}
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].rminCs = fin; // минимальное возможное значение постоянной Смагоринского.
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].rmaxCs = fin; // максимальное возможное значение постоянной Смагоринского.
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].itypeFiltrGermano = din; // тип фильтра в модели Германо 1991 года.
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].roughness = 1.0e-6*fin; // шероховатость стенки в м.
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].rRi_mult = fin; // множитель корректирующий турбулентное число Ричардсона.
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].rSelectiveAngle = fin; // пороговое значение угла в модели Selective Smagorinsky.
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].ipowerroughness = din; // показатель степени в модели учёта шероховатости.
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = din; // тип фильтра в модели Selective Smagorinsky.
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].bfdelta = din; // учёт неравномерности сетки.
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].bSmagorinskyLilly = din; // модель Смагоринского-Лиллу.
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].bsurface_roughness = din; // учёт шероховатости стенки.
					fscanf(fp, "%d", &din);
					// учёт течений с кривизной линий тока.
					if (din == 1) {
						eqin.fluidinfo[i].bSwirlAmendment = true;
					}
					else {
						eqin.fluidinfo[i].bSwirlAmendment = false;
					}
					fscanf(fp, "%d", &din);
					// учёт избирательности в модели Смагоринского
					if (din == 1) {
						eqin.fluidinfo[i].bSelectiveSmagorinsky = true;
					}
					else {
						eqin.fluidinfo[i].bSelectiveSmagorinsky = false;
					}

					// Параметры преобразователя картинок для отчетов.
					// 5.01.2018
					fscanf(fp, "%f", &fin);
					pfpir.fminimum = scale * fin;
					fscanf(fp, "%f", &fin);
					pfpir.fmaximum = scale * fin;
#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					switch (din) {
					case 0: pfpir.idir = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
					case 1: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					case 2: pfpir.idir = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
					default: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					}
					//pfpir.idir = din;

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					AMG1R6_LABEL = din;

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					if (din > 9999) {
						b_iluk_amg1r5_LABEL_D = true;
						nrd_LABEL = din - 10000;
					}
					else {
						nrd_LABEL = din;
					}

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					if (din > 9999) {
						b_iluk_amg1r5_LABEL_U = true;
						nru_LABEL = din - 10000;
					}
					else {
						nru_LABEL = din;
					}

#if doubleintprecision == 1
					fscanf(fp, "%f", &fin);
#else
					fscanf(fp, "%f", &fin);
#endif
					ecg2_LABEL = din;

#if doubleintprecision == 1
					fscanf(fp, "%f", &fin);
#else
					fscanf(fp, "%f", &fin);
#endif
					ewt2_LABEL = din;

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					number_processors_global_var = (int)(din);


#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					idevice_Tesla = (int)(din);

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					if (din>=0) {
						number_iteration_SIMPLE_algorithm = (unsigned int)(din);
					}
					else {
						number_iteration_SIMPLE_algorithm =0;
						std::cout<<"WARNING : number_iteration_SIMPLE_algorithm =0;"<<std::endl;
				    }

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					switch (din) {
					case 0: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::NONE_only_amg1r5;
						break;
					case 1: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
						break;
					case 2: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::FGMRes_plus_amg1r5;
						break;
					case 3: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5;
						break;
					default:
						stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
						break;
					}


#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					BonLevelDrobim = din;

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					isizearrAdaptRegion = din;


					if (isizearrAdaptRegion > 0) {

						AdaptRegion = new ADAPTREGION[isizearrAdaptRegion];

						AdaptRegion[0].ilevelMeshAdaptRegion = -1;
						AdaptRegion[1].ilevelMeshAdaptRegion = -1;
						AdaptRegion[2].ilevelMeshAdaptRegion = -1;

						for (int i_63 = 0; i_63 < isizearrAdaptRegion; ++i_63) {

							// Адаптация сетки региона номер 1.
							if (i_63 == 0) {

#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[0].xS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[0].xE = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[0].yS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[0].yE = scale * fin;

#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[0].zS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[0].zE = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%lld", &din);
#else
								fscanf(fp, "%d", &din);
#endif
								AdaptRegion[0].ilevelMeshAdaptRegion = din;

							}


							if (i_63 == 1) {


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[1].xS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[1].xE = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[1].yS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[1].yE = scale * fin;

#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[1].zS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[1].zE = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%lld", &din);
#else
								fscanf(fp, "%d", &din);
#endif
								AdaptRegion[1].ilevelMeshAdaptRegion = din;

							}


							if (i_63 == 2) {

#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[2].xS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[2].xE = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[2].yS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[2].yE = scale * fin;

#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[2].zS = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%f", &fin);
#else
								fscanf(fp, "%f", &fin);
#endif
								AdaptRegion[2].zE = scale * fin;


#if doubleintprecision == 1
								fscanf(fp, "%lld", &din);
#else
								fscanf(fp, "%d", &din);
#endif
								AdaptRegion[2].ilevelMeshAdaptRegion = din;


		}

	}

}




#if doubleintprecision == 1
					fscanf(fp, "%f", &fin);
#else
					fscanf(fp, "%f", &fin);
#endif
					 
					// Найдено успешно.
					// Свободный параметр передаваемый из интерфейса и используемый для отладки.
					free_debug_parametr1 = static_cast<doublereal>(fin);
					
					 
				}
			}

			if (fp != NULL) {
				fclose(fp); // закрытие файла
				fp = NULL;
			}
		}
	
	}


	integer ilb_p = 0;// Количество блоков внутри которых задана тепловая мощность.
	doublereal dpower = 0.0; // Суммарная тепловая мощность в блоках.
	integer ipoly = 0, icyl = 0, iprism = 0;
	integer ihol = 0, isol = 0, iflui = 0;
	for (integer i_1 = 0; i_1 < lb; ++i_1) {
		if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			ihol++;
		}
		if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
			iflui++;
		}
		if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
			isol++;
		}
		if (b[i_1].g.itypegeom == PRISM) {
			// 0 - PRISM object

			// 13.08.2019 Автомат настройки допусков сетки.
			// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
			// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
			// Она не может быть больше стороны кабинета деленной на 15.
			if (0.1*fabs(b[i_1].g.xE - b[i_1].g.xS) < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*fabs(b[i_1].g.xE - b[i_1].g.xS);
			if (0.1*fabs(b[i_1].g.yE - b[i_1].g.yS) < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*fabs(b[i_1].g.yE - b[i_1].g.yS);
			if (0.1*fabs(b[i_1].g.zE - b[i_1].g.zS) < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*fabs(b[i_1].g.zE - b[i_1].g.zS);

			if (lb == 1) {// Течение в каверне или тест Дэвиса.
				shorter_length_for_simplificationX_BASIC = 1.0e-10;
				shorter_length_for_simplificationY_BASIC = 1.0e-10;
				shorter_length_for_simplificationZ_BASIC = 1.0e-10;
			}
			else {
				if (shorter_length_for_simplificationX_BASIC > 0.067*(fabs(b[0].g.xE - b[0].g.xS))) shorter_length_for_simplificationX_BASIC = 0.067*(fabs(b[0].g.xE - b[0].g.xS));
				if (shorter_length_for_simplificationY_BASIC > 0.067*(fabs(b[0].g.yE - b[0].g.yS))) shorter_length_for_simplificationY_BASIC = 0.067*(fabs(b[0].g.yE - b[0].g.yS));
				if (shorter_length_for_simplificationZ_BASIC > 0.067*(fabs(b[0].g.zE - b[0].g.zS))) shorter_length_for_simplificationZ_BASIC = 0.067*(fabs(b[0].g.zE - b[0].g.zS));
			}

			iprism++;
			if (b[i_1].n_Sc > 0) {
				doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
				doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS)*fabs(b[i_1].g.yE - b[i_1].g.yS)*fabs(b[i_1].g.zE - b[i_1].g.zS);
				if (vol < 1.0e-36) {
					printf("ERROR: zero volume in PRISM block number %lld\n", i_1);
					system("PAUSE");
					exit(1);
				}
				//if (fabs(b[i_1].arr_Sc[0]) > 0.0) {
				if (pdiss > 0.0) {
					ilb_p++;
					//dpower += b[i_1].arr_Sc[0];
					dpower += pdiss*vol;
				}
			}
		}
		if (b[i_1].g.itypegeom == CYLINDER) {
			// Cylinder

			// На тот случай если геометрия состоит только из цилиндров.
			// 13.08.2019 Автомат настройки допусков сетки. 
			// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
			// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
			switch (b[i_1].g.iPlane) {
			case XY_PLANE:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-36) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
				}
				break;
			case XZ_PLANE:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-36) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_in_cyl;
				}
				break;
			case YZ_PLANE:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-36) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_in_cyl;
				}
				break;
			}


			icyl++;
			if (b[i_1].n_Sc > 0) {
				doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
				doublereal vol = 0.0;
				
				const doublereal MPI0 = 3.1415926535;

				vol = b[i_1].g.Hcyl*MPI0*(b[i_1].g.R_out_cyl*b[i_1].g.R_out_cyl - b[i_1].g.R_in_cyl*b[i_1].g.R_in_cyl);
				if (vol < 1.0e-36) {
					printf("ERROR: zero volume in CYLINDER block number %lld\n", i_1);
					system("PAUSE");
					exit(1);
				}
				if (pdiss > 0.0) {
					ilb_p++;

					dpower += pdiss*vol;
					//printf("ERROR: non zero power in cylinder object.\n");
					//system("PAUSE");
					//exit(1);
				}
			}
		}
		if (b[i_1].g.itypegeom == POLYGON) {
			// Polygon
			ipoly++;
			if (b[i_1].n_Sc > 0) {
				doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
				// Объём полигона.
				doublereal vol_poly = Volume_polygon(b[i_1].g.nsizei, b[i_1].g.xi, b[i_1].g.yi, b[i_1].g.zi, b[i_1].g.hi, b[i_1].g.iPlane_obj2);
				if (vol_poly < 1.0e-36) {
					printf("ERROR: zero volume in POLYGON block number %lld\n", i_1);
					system("PAUSE");
					exit(1);
				}
				if (pdiss > 0.0) {
					ilb_p++;

					dpower += pdiss*vol_poly;
					//printf("ERROR: non zero power in polygon object.\n");
					//system("PAUSE");
					//exit(1);
				}
			}
		}
	}

	doublereal dsoupow = 0.0; // интегральная тепловая мощность плоских источников тепла.
	for (integer i_1 = 0; i_1 < ls; ++i_1) {
		dsoupow += s[i_1].power;
	}

	printf("Apriory quick model statistics:\n");
	printf("number of thermal power blocks lb_p=%lld\n", ilb_p);
	//printf("Blocks integral power =%e W\n", dpower);
	std::cout << "Blocks integral power =" << dpower << " W" << std::endl;
	printf("number of sources ls=%lld\n", ls);
	//printf("Sources integral power = %e W\n", dsoupow);
	std::cout <<"Sources integral power ="<<dsoupow<<" W"<< std::endl;
	//printf("Full total power = %e W\n", dpower + dsoupow);
	std::cout <<"Full total power ="<<(dpower + dsoupow)<<" W"<< std::endl;
	// Запоминаем полное тепловыделение в модели.
	d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = dpower + dsoupow;
	printf("number of blocks lb=%lld\n", lb);
	printf("PRISMS = %lld, CYLINDERS = %lld, POLYGONS = %lld\n", iprism, icyl, ipoly);
	printf("SOLID: %lld\n", isol);
	printf("HOLLOW: %lld\n", ihol);
	printf("FLUID: %lld\n", iflui);
	printf("number of walls lw=%lld\n", lw);
	printf("number of units lu=%lld\n", lu);



#endif

#ifndef MINGW_COMPILLER
#if doubleintprecision == 0

	// eqin - информация о наборе решаемых уравнений.

	// dgx, dgy, dgz - вектор силы тяжести.
	// inx, iny, inz - количество точек по каждой из осей.

	FILE* fp;
	errno_t err1;
	err1 = fopen_s(&fp, fname, "r");

	if ((err1) != 0) {
		printf("Error!!! No input File premeshin.txt.\n");
		printf("You must use the graphical user interface\n");
		printf("AliceMesh_v0_45.exe in Delphi, which will \n");
		printf("prepare the model and write the premeshin.txt\n");
		printf("file in the required format.\n");
		//system("PAUSE");
		system("pause");
		// Если файла premeshin.txt нет то мы не можем ничего обрабатывать
		// т.к. элементарно отсутствуют входные данные. Поэтому мы выходим из 
		// приложения.
		exit(1);
	}
	else
	{
		if (fp != nullptr) {
			float fin = 0.0;
			integer din = 0;
			doublereal scale = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания



			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
			case 1: ionly_solid_visible= WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE; break;
			default: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
			}
			//ionly_solid_visible = din;
			switch (ionly_solid_visible) {
			case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE:
				std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
				break;
			case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE:
				std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
				break;
			}

			fscanf_s(fp, "%f", &fin);
			scale = fin;
#ifndef NO_OPENGL_GLFW
			scale_all = 1.0f / scale;
#endif
			fscanf_s(fp, "%d", &din);
			lmatmax = din;
			fscanf_s(fp, "%d", &din);
			lb = din;
			fscanf_s(fp, "%d", &din);
			ls = din;
			fscanf_s(fp, "%d", &din);
			lw = din;
			fscanf_s(fp, "%d", &din);
			ltdp = din; // количество уникальных данных с табличными данными по зависимости расеиваемой мощности от температуры.


			// Считываем значение вектора силы тяжести:
			fscanf_s(fp, "%f", &fin);
			dgx = fin;
			fscanf_s(fp, "%f", &fin);
			dgy = fin;
			fscanf_s(fp, "%f", &fin);
			dgz = fin;

			// считываем количество точек на каждой координатной оси
			fscanf_s(fp, "%d", &din);
			inx = din;
			fscanf_s(fp, "%d", &din);
			iny = din;
			fscanf_s(fp, "%d", &din);
			inz = din;

			fscanf_s(fp, "%f", &fin);
			operatingtemperature = fin; // Operating Temperature
			operating_temperature_for_film_coeff = fin;

			// инициализация компонент скорости константой.
			// Единые значения для всей расчётной области.
			// initialization value.
			fscanf_s(fp, "%f", &fin);
			starting_speed_Vx = fin;
			fscanf_s(fp, "%f", &fin);
			starting_speed_Vy = fin;
			fscanf_s(fp, "%f", &fin);
			starting_speed_Vz = fin;


			// Считываем координаты опорной точки через которую проходит пользовательская линия Variation Plot.
			fscanf_s(fp, "%f", &fin);
			Tochka_position_X0_for_XY_Plot = scale * fin;
			fscanf_s(fp, "%f", &fin);
			Tochka_position_Y0_for_XY_Plot = scale * fin;
			fscanf_s(fp, "%f", &fin);
			Tochka_position_Z0_for_XY_Plot = scale * fin;
			// Направление линии совпадает с направлением 
			// одной из осей декартовой прямоугольной системы координат:
			// 0 - Ox, 1 - Oy, 2 - Oz.
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
			case 1: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
			case 2: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
			default: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL;	break;
			}
			//idirectional_for_XY_Plot = din;

			fscanf_s(fp, "%f", &fin);
			etalon_max_size_ratio = fin; // подробность расчётной сетки.

			fscanf_s(fp, "%f", &fin);
			etalon_max_size_ratio2 = fin; // Критерий качества расчётной сетки на основе FlowVision.

			//printf("etalon_max_size_ratio=%e etalon_max_size_ratio2=%e\n", etalon_max_size_ratio, etalon_max_size_ratio2);
			//system("PAUSE");

			
			// 0.	none
			// 1.	Snap to grid
			// 2.	Snap to grid ALICE
			// 3.	Snap to grid ++
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: bsnap_TO_global = 0;  break;
			case 1: bsnap_TO_global = 1;  break;
			case 2: bsnap_TO_global = 2;  break;
			case 3: bsnap_TO_global = 3;  break;
			default: bsnap_TO_global = 1;  break;
			}


			fscanf_s(fp, "%d", &din);
			iswitchsolveramg_vs_BiCGstab_plus_ILU2 = din; // Выбор решающего устройства: либо amg1r5 либо BiCGStab+ILU2.

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
			case 1: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER; break;
			case 2: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::CAMG_RUMBA_v0_14_SECOND_T_SOLVER; break;
			case 3: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMG1R5_SECOND_T_SOLVER; break;
			case 4: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMGCL_SECONT_T_SOLVER; break;
			default: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
			}
			//iswitchsolveramg_vs_BiCGstab_plus_ILU6 = din; // Выбор решающего устройства: либо РУМБА0.14 либо BiCGStab+ILU6.

			fscanf_s(fp, "%d", &din);
			if (din == 1) {
				// SIMPLEC algorithm.
				iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby;
			}
			else {
				// SIMPLE algorithm 1972.
				iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;
			}


			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0:  iprefix_Scheme_Flow = CR;
				break;
			case 1:  iprefix_Scheme_Flow = UDS;
				break;
			case 2:  iprefix_Scheme_Flow = COMB;
				break;
			case 3:  iprefix_Scheme_Flow = POLY;
				break;
			case 4:  iprefix_Scheme_Flow = EXP;
				break;
			case 5:  iprefix_Scheme_Flow = BULG;
				break;
			case 6:  iprefix_Scheme_Flow = POW;
				break;
			default:
				iprefix_Scheme_Flow = UDS;
				break;
			}

			fscanf_s(fp, "%d", &din);
			// выбор схемы для потока жидкости.
			// Внимание эти определения должны полностью соответствовать 
			// определениям в файле my_approx_convective2.c
			switch (din) {
			case 1: iFLOWScheme = UNEVEN_MUSCL; break; // MUSCL 2
			case 2: iFLOWScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
			case 3: iFLOWScheme = UNEVEN_HLPA; break; // HLPA 2
			case 4: iFLOWScheme = UNEVEN_SMART; break; // SMART 3
			case 5: iFLOWScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
			case 6: iFLOWScheme = UNEVEN_SMARTER; break; // SMARTER 3
			case 7: iFLOWScheme = UNEVEN_LPPA; break; // LPPA 3
			case 8: iFLOWScheme = UNEVEN_VONOS; break; // VONOS 3
			case 9: iFLOWScheme = UNEVEN_STOIC; break; // STOIC
			case 10: iFLOWScheme = UNEVEN_CLAM; break; // CLAM
			case 11: iFLOWScheme = UNEVEN_OSHER; break; // OSHER
			case 12: iFLOWScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
			case 13: iFLOWScheme = UNEVEN_SUPER_C; break; // SUPER_C
			case 14: iFLOWScheme = UNEVEN_ISNAS; break; // ISNAS
			case 15: iFLOWScheme = UNEVEN_CUBISTA; break; // CUBISTA
			case 16: iFLOWScheme = UNEVEN_GAMMA; break; // GAMMA
			case 17: iFLOWScheme = UNEVEN_COPLA; break; // COPLA
			case 18: iFLOWScheme = UNEVEN_SECBC; break; // SECBC
			case 19: iFLOWScheme = UNEVEN_SGSD; break; // SGSD
			case 20: iFLOWScheme = UNEVEN_WENO5; break; // WENO5
			default: iFLOWScheme = UDS; break; // UDS самая стабильная схема.
			}

			fscanf_s(fp, "%d", &din);
			// выбор схемы для температуры в потоке жидкости.
			// Внимание эти определения должны полностью соответствовать 
			// определениям в файле my_approx_convective2.c
			switch (din) {
			case 1: iTEMPScheme = UNEVEN_MUSCL; break; // MUSCL 2
			case 2: iTEMPScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
			case 3: iTEMPScheme = UNEVEN_HLPA; break; // HLPA 2
			case 4: iTEMPScheme = UNEVEN_SMART; break; // SMART 3
			case 5: iTEMPScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
			case 6: iTEMPScheme = UNEVEN_SMARTER; break; // SMARTER 3
			case 7: iTEMPScheme = UNEVEN_LPPA; break; // LPPA 3
			case 8: iTEMPScheme = UNEVEN_VONOS; break; // VONOS 3
			case 9: iTEMPScheme = UNEVEN_STOIC; break; // STOIC
			case 10: iTEMPScheme = UNEVEN_CLAM; break; // CLAM
			case 11: iTEMPScheme = UNEVEN_OSHER; break; // OSHER
			case 12: iTEMPScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
			case 13: iTEMPScheme = UNEVEN_SUPER_C; break; // SUPER_C
			case 14: iTEMPScheme = UNEVEN_ISNAS; break; // ISNAS
			case 15: iTEMPScheme = UNEVEN_CUBISTA; break; // CUBISTA
			case 16: iTEMPScheme = UNEVEN_GAMMA; break; // GAMMA
			case 17: iTEMPScheme = UNEVEN_COPLA; break; // COPLA
			case 18: iTEMPScheme = UNEVEN_SECBC; break; // SECBC
			case 19: iTEMPScheme = UNEVEN_SGSD; break; // SGSD
			case 20: iTEMPScheme = UNEVEN_WENO5; break; // WENO5
			default: iTEMPScheme = UDS; break; // UDS самая стабильная схема.
			}



			// Выбор сеточного генератора.
			fscanf_s(fp, "%d", &din);
			switch (din) {
				case 0: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; break;
				case 1: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER; break;
				case 2: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
				default: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
			}
			//iswitchMeshGenerator = din;

			fscanf_s(fp, "%d", &din);
			steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY;
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3)  || (din == 5) || (din == 6) || 
				(din == 7) || (din == 8) || (din == 9) || (din == 10) || (din == 11)||(din==12)||(din==13)) {
				// 0 - thermal only steady state calculation,
				// 1 - thermal only unsteady calculation,
				// 2 - mesh generator only.
				// 3 - fluid dynamic steady state.
				// 5 - Static Structural (Thermal solver #2)
				// 6 - Thermal Stress
				// 7 - Unsteady thermal solver #2
				// 8 - Visualisation only
				// 9 - cfd unsteady fluid dynamic.
				// 10 - NETWORK_T Графовый метод решения уравнения теплопроводности.
				// 11 - UNSTEADY NETWORK_T Нестационврный графовый метод решения уравнения теплопроводности.
				// 12 - Нестационарная механика,
				// 13 - Нестационарная механика совместно с нестационарной теплопередачей.
				switch (din) {
				case 0: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE; break;
				case 1: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE; break;
				case 2: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
				case 3: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_STEADY; break;
				case 5: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL; break;
				case 6: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
				case 7: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER; break;
				case 8: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT; break;
				case 9: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY; break;
				case 10: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T;  break;
				case 11: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY; break;
				case 12: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL;  break;
				case 13: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
				default: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
				}
				//steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady -1 calculation.
			}
			else {
				printf("error input parametr steady or unsteady calculation\n");
				system("PAUSE");
				exit(1);
			}

			fscanf_s(fp, "%d", &din);
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 4)) {
				switch (din) {
				case 0: glTSL.id_law = TIME_STEP_lAW_SELECTOR::LINEAR;
					break;
				case 1: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE;
					break;
				case 2: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE2;
					break;
				case 3: glTSL.id_law = TIME_STEP_lAW_SELECTOR::HOT_COLD;
					break;
				case 4: glTSL.id_law = TIME_STEP_lAW_SELECTOR::PIECEWISE_CONSTANT;
					break;
				}
			}
			else {
				printf("error input parametr timestep law\n");
				system("PAUSE");
				exit(1);
			}
			fscanf_s(fp, "%f", &fin);
			if ((fin < 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law Factor a\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.Factor_a_for_Linear = fin; // Factor_a
			fscanf_s(fp, "%f", &fin);
			if (fin < 0.0) {
				printf("error input parametr timestep law tau must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau = fin; // длительность импульса.
			fscanf_s(fp, "%f", &fin);
			glTSL.Q = fin;  // Скважность.
			// Параметры импульсного режима для SquareWave 2 режима.
			fscanf_s(fp, "%f", &fin);
			if ((fin < 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law SquareWave2 multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin < 0.0) {
				printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin < 0.0) {
				printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau2 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin < 0.0) {
				printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau_pause = fin;

			fscanf_s(fp, "%f", &fin);
			if ((fin < 0.0) || (fin > 1.0)) {
				printf("error input parametr timestep law SquareWave2 off_multiplyer must be [0..1]\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.off_multiplyer = fin;

			fscanf_s(fp, "%d", &din);
			glTSL.n_cycle = din;
			fscanf_s(fp, "%f", &fin);
			if (fin < 0.0) {
				printf("error input parametr timestep law SquareWave2 Period must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.T_all = fin;
			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle * (2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl < 0.0) {
				printf("error in parameters Square Wave SquareWave2 time step law.\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}

			fscanf_s(fp, "%f", &fin);
			if (fin < 0.0) {
				printf("error input parametr on_time_double_linear law hot cold reshime must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.on_time_double_linear = fin;

			// Время окончания нестационарного моделирования при расчёте теплопередачи в твёрдом теле.
			fscanf_s(fp, "%f", &fin);
			globalEndTimeUnsteadyTemperatureCalculation = fin;




			// Newton-Richman condition.
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
            case 1: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC; break;
	        case 2: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC; break;
		    case 3: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC; break;
			default :adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
			}
			//adiabatic_vs_heat_transfer_coeff = din; // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
			fscanf_s(fp, "%f", &fin);
			film_coefficient = fin;
			// AЛИС сетка
			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				// Обычная структурированная сетка.
				b_on_adaptive_local_refinement_mesh = false;
			}
			else {
				// АЛИС
				b_on_adaptive_local_refinement_mesh = true;
			}
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
			case 1: itype_ALICE_Mesh= TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH; break;
			default: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
			}
			//itype_ALICE_Mesh = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.m_restart = static_cast<int>(din);
			// classical algebraic multigrid parameters:
			// only for my_agregat_amg.cu.
			// only for my_agregat_amg.cu.
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
				break;
			case 1: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::QUICK_SORT;
				break;
			case 2: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::HEAP_SORT;
				break;
			case 3: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::TIM_SORT;
				break;
			default:
				my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
				break;
			}
			fscanf_s(fp, "%d", &din);
			//my_amg_manager.maximum_levels = din;
			my_amg_manager.maximum_delete_levels_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Stress = din;

			// type interpolation procedure:
			//fscanf_s(fp, "%d", &din);
			//my_amg_manager.number_interpolation_procedure = din;

			fscanf_s(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.number_interpolation_procedure_Stress = din;

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
				break;
			case 1: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			case 2: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
				break;
			case 3: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
				break;
			case 4: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
				break;
			case 5: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
				break;
			case 6: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
				break;
			default:
				my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
				break;
			}


			fscanf_s(fp, "%d", &din);
			//my_amg_manager.itypemodifyinterpol = din;
			//my_amg_manager.baglomeration_with_consistency_scaling = din;
			my_amg_manager.bdiagonal_dominant = din;
			fscanf_s(fp, "%d", &din);
			//my_amg_manager.inumberadaptpass = din;


			// 23.02.2018
			// print matrix portrait
			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bTemperatureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bTemperatureMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}
			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bSpeedMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bSpeedMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}
			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bPressureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bPressureMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}
			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bStressMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
			}
			else {
				my_amg_manager.bStressMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
			}


			// 01.05.2017
			// truncation of interpolation:
			fscanf_s(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.itruncation_interpolation_Stress = din;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Temperature = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Speed = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Pressure = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.truncation_interpolation_Stress = fin;

			// number nFinnest sweeps:
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nFinnest_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nFinnest_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nFinnest_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nFinnest_Stress = din;

			// number presweeps:
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu1_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu1_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu1_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu1_Stress = din;

			// number postsweeps:
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Stress = din;

			// memory size:
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Stress = din;


			// Параметр верхней релаксации в сглаживателе.
			//fscanf_s(fp, "%f", &fin);
			//my_amg_manager.gold_const_Temperature = fin;
			//fscanf_s(fp, "%f", &fin);
			//my_amg_manager.gold_const_Speed = fin;
			//fscanf_s(fp, "%f", &fin);
			//my_amg_manager.gold_const_Pressure = fin;
			//fscanf_s(fp, "%f", &fin);
			//my_amg_manager.gold_const_Stress = fin;

			// использовать ли ilu2 smoother.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Temperature = din;
			/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Temperature = din;
				my_amg_manager.b_gmresTemp = true;
			}
			*/

			fscanf_s(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Speed = din;
			/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Speed = din;
				my_amg_manager.b_gmresSpeed = true;
			}
			*/

			fscanf_s(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Pressure = din;
			/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Pressure = din;
				my_amg_manager.b_gmresPressure = true;
			}
			*/

			fscanf_s(fp, "%d", &din);
			my_amg_manager.ilu2_smoother_Stress = din;
			/*
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Stress = din;
				my_amg_manager.b_gmresStress = true;
			}
			*/
			my_amg_manager.gold_const_Temperature = return_gold_const(my_amg_manager.ilu2_smoother_Temperature);
			my_amg_manager.gold_const_Speed = return_gold_const(my_amg_manager.ilu2_smoother_Speed);
			my_amg_manager.gold_const_Pressure = return_gold_const(my_amg_manager.ilu2_smoother_Pressure);
			my_amg_manager.gold_const_Stress = return_gold_const(my_amg_manager.ilu2_smoother_Stress);


			fscanf_s(fp, "%d", &din);
			my_amg_manager.Chebyshev_degree = din;

			if (!((din == 5) || (din == 8) || (din == 16) || (din == 25) || (din == 32) || (din == 64) || (din == 128) || (din == 256) || (din == 512))) {
				std::cout << "Chebyshev degree undefined. Chebyshev degree mast be equal: 5, 8, 16, 25, 32, 64, 128, 256, 512 in this programm. Chebyshev degree == " << din << std::endl;
				system("pause");
				my_amg_manager.Chebyshev_degree = 5;
			}

			// strength threshold:
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Temperature = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Speed = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Pressure = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Stress = fin;

			// magic threshold:
			//fscanf_s(fp, "%f", &fin);
			//my_amg_manager.magic = fin;

			// magic <=> F_to_F
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.F_to_F_Temperature = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.F_to_F_Speed = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.F_to_F_Pressure = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.F_to_F_Stress = fin;

			// AMG Splitting (coarsening)
			// Способ построения C-F разбиения: 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
			// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}
			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
				break;
			case 1: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
				break;
			case 2: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
				break;
			case 3: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
				break;
			case 4: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
				break;
			case 5: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
				break;
			case 6: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
				break;
			case 7: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
				break;
			case 8: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			case 9: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
				break;
			case 10: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
				break;
			case 11: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
				break;
			default:
				my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
				break;
			}

			// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
			// Если din==1, Stabilization BiCGStab.
			// 8.01.2017 Метод Хенка ван дер Ворста BiCGStab 1992.
			// предобусловленный алгебраическим многосеточным методом.
			// 9.01.2018 Если din==2, FGMRes предобусловленный алгебраическим многосеточным методом.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.istabilizationTemp = din; // 0 - none
			fscanf_s(fp, "%d", &din);
			my_amg_manager.istabilizationSpeed = din; // 0 - none
			fscanf_s(fp, "%d", &din);
			my_amg_manager.istabilizationPressure = din; // 0 - none
			fscanf_s(fp, "%d", &din);
			my_amg_manager.istabilizationStress = din; // 0 - none
			fscanf_s(fp, "%d", &din);
			my_amg_manager.ipatch_number = din; // 0 - патч не применяется.

			// Печать лога на консоль.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iprint_log_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iprint_log_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iprint_log_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iprint_log_Stress = din;


			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Temperature_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Temperature_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}

			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Speed_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Speed_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}

			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Pressure_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Pressure_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}

			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bthreshold_Stress_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}
			else {
				my_amg_manager.bthreshold_Stress_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
			}


			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Temperature = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Temperature = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}

			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Speed = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Speed = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}


			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Pressure = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Pressure = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}

			fscanf_s(fp, "%d", &din);
			if (din == 0) {
				my_amg_manager.bcf_reorder_Stress = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}
			else {
				my_amg_manager.bcf_reorder_Stress = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
			}



			fscanf_s(fp, "%d", &din);
			my_amg_manager.amgcl_smoother = din;

			fscanf_s(fp, "%d", &din);
			my_amg_manager.amgcl_selector = din;

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
				break;
			case 1: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::FGMRes;
				break;
			default: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
				break;
			}

			fscanf_s(fp, "%d", &din);
			my_amg_manager.lfil = din;

			fscanf_s(fp, "%d", &din);
			ireconstruction_free_construct_alloc = din;

			fscanf_s(fp, "%d", &din);
			ianimation_write_on = din;

			// выделение оперативной памяти.
			gtdps = new TEMP_DEP_POWER[ltdp];
			matlist = new TPROP[lmatmax];
			b = new BLOCK[lb];
			s = new SOURCE[ls];
			w = new WALL[lw+1];// +1 это гран условие по умолчанию, заглушка, чтобы избежать неопределенности.
			integer i = 0; // счётчик цикла for

			for (i = 0; i < ltdp; ++i) {
				// считывание имён файлов.
				gtdps[i].sname = new char[100]; // выделение памяти
				fscanf_s(fp, "%s", gtdps[i].sname, 100);
				//printf("%s",gtdps[i].sname);
				//system("PAUSE");
				// построение таблицы в памяти.
				my_read_power_table(gtdps[i].sname, gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
				printf("printeger table\n"); // debug
				my_print_table(gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
				printf("Please, press any key to continue...\n");
				//system("PAUSE");
				system("pause");

			}



			// считывание базы материалов
			for (i = 0; i < lmatmax; ++i) {
				// свойства материалов:
				// плотность
				fscanf_s(fp, "%f", &fin);
				matlist[i].rho = fin;
				// теплоёмкость при постоянном давлении
				//fscanf(fp, "%f", &fin);
				//matlist[i].cp = fin;
				fscanf_s(fp, "%d", &din);
				matlist[i].n_cp = din;
				matlist[i].arr_cp = nullptr;
				matlist[i].temp_cp = nullptr;
				matlist[i].arr_cp = new float[matlist[i].n_cp];
				matlist[i].temp_cp = new float[matlist[i].n_cp];
				if (matlist[i].temp_cp == nullptr) {
					printf("problem memory allocation for temp_cp\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_cp == nullptr) {
					printf("problem memory allocation for arr_cp\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					matlist[i].temp_cp[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					matlist[i].arr_cp[i_4] = fin;
				}
				// теплопроводность
				//fscanf(fp, "%f", &fin);
				//matlist[i].lam = fin;
				fscanf_s(fp, "%d", &din);
				matlist[i].n_lam = din;
				matlist[i].arr_lam = nullptr;
				matlist[i].temp_lam = nullptr;
				matlist[i].arr_lam = new float[matlist[i].n_lam];
				matlist[i].temp_lam = new float[matlist[i].n_lam];
				if (matlist[i].temp_lam == nullptr) {
					printf("problem memory allocation for temp_lam\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_lam == nullptr) {
					printf("problem memory allocation for arr_lam\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					matlist[i].temp_lam[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					matlist[i].arr_lam[i_4] = fin;
				}
				// ортотропность теплопроводности:
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z = fin;

				// 28.08.2020
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x_beta_t_solid = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y_beta_t_solid = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z_beta_t_solid = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x_Young_Module = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y_Young_Module = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z_Young_Module = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_xy = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_xz = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_yz = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_yx = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_zx = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_Poisson_ratio_zy = fin;
				fscanf_s(fp, "%d", &din);
				if (din == 1) {
					matlist[i].bActive_ShearModule = true;
				}
				else {
					matlist[i].bActive_ShearModule = false;
				}
				fscanf_s(fp, "%f", &fin);
				matlist[i].ShearModule_xy = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].ShearModule_yz = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].ShearModule_xz = fin;


				// 5.08.2017.
				// Коэффициенты для задачи упругости.
				// Модуль Юнга и коэффициент Пуассона.
				//doublereal Poissonratio = 0.154;
				//doublereal Youngmodule = 217.5e9;
				//fscanf_s(fp, "%f", &fin);
				//Poissonratio = fin;
				//if (fabs(fin) < 1.0e-30) {
					//printf("ERROR !!! Zero Poisson Ratio in model. \n");
					//printf("ERROR !!! Model is incorrect...\n");
					//system("PAUSE");
					//exit(1);
				//}
				//matlist[i].Poisson_ratio = fin;

				fscanf_s(fp, "%d", &din);
				matlist[i].n_Poisson_ratio = din;
				matlist[i].arr_Poisson_ratio = nullptr;
				matlist[i].temp_Poisson_ratio = nullptr;
				matlist[i].arr_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
				matlist[i].temp_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
				if (matlist[i].temp_Poisson_ratio == nullptr) {
					printf("problem memory allocation for temp_Poisson_ratio\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_Poisson_ratio == nullptr) {
					printf("problem memory allocation for arr_Poisson_ratio\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_Poisson_ratio; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					matlist[i].temp_Poisson_ratio[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					matlist[i].arr_Poisson_ratio[i_4] = fin;
				}


				//fscanf_s(fp, "%f", &fin);
				//Youngmodule = fin * 1e9;
				//if (fabs(fin) < 1.0e-30) {
					//printf("ERROR !!! Zero Young Module in model. \n");
					//printf("ERROR !!! Model is incorrect...\n");
					//system("PAUSE");
					//exit(1);
				//}
				//matlist[i].Young_Module = fin * 1e9;

				fscanf_s(fp, "%d", &din);
				matlist[i].n_YoungModule = din;
				matlist[i].arr_Young_Module = nullptr;
				matlist[i].temp_Young_Module = nullptr;
				matlist[i].arr_Young_Module = new float[matlist[i].n_YoungModule];
				matlist[i].temp_Young_Module = new float[matlist[i].n_YoungModule];
				if (matlist[i].temp_Young_Module == nullptr) {
					printf("problem memory allocation for temp_Young_Module\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_Young_Module == nullptr) {
					printf("problem memory allocation for arr_Young_Module\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_YoungModule; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					matlist[i].temp_Young_Module[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					matlist[i].arr_Young_Module[i_4] = fin * 1.0e+9;
				}

				//fscanf_s(fp, "%f", &fin);
				// beta_t_solid*1E-6
				//matlist[i].n_beta_t_solid = 1;
				//matlist[i].temp_beta_t_solid = new doublereal[1];
				//matlist[i].temp_beta_t_solid[0] = 25.0;
				//matlist[i].arr_beta_t_solid = new doublereal[1];
				//matlist[i].arr_beta_t_solid[0] = fin * 1E-6;

				fscanf_s(fp, "%d", &din);
				matlist[i].n_beta_t_solid = din;
				matlist[i].arr_beta_t_solid = nullptr;
				matlist[i].temp_beta_t_solid = nullptr;
				matlist[i].arr_beta_t_solid = new float[matlist[i].n_beta_t_solid];
				matlist[i].temp_beta_t_solid = new float[matlist[i].n_beta_t_solid];
				if (matlist[i].temp_beta_t_solid == nullptr) {
					printf("problem memory allocation for temp_beta_t_solid\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_beta_t_solid == nullptr) {
					printf("problem memory allocation for arr_beta_t_solid\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < matlist[i].n_beta_t_solid; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					matlist[i].temp_beta_t_solid[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					matlist[i].arr_beta_t_solid[i_4] = fin * 1.0e-6;
				}

				// Коэффициенты Ламе не используются в данном коде.
				
				// коэффициент динамической вязкости
				fscanf_s(fp, "%f", &fin);
				matlist[i].mu = fin;
				// коэффициент линейного температурного расширения
				fscanf_s(fp, "%f", &fin);
				matlist[i].beta_t = fin;
				// признак библиотечности материала
				fscanf_s(fp, "%d", &din);
				matlist[i].blibmat = din;
				// номер материала в библиотеке
				fscanf_s(fp, "%d", &din);
				matlist[i].ilibident = din;
				//printf("blibmat=%d ilibident=%d\n", matlist[i].blibmat, matlist[i].ilibident);
				//system("pause");

				// для каждого КО данного блока,
				// если только он не перекрывается другими блоками
				// может быть использовано приближение 
				// Обербека-Буссинеска с соответствующей опорной температурой Tref.
				fscanf_s(fp, "%d", &din);
				switch (din) {
				case 0: matlist[i].bBussineskApproach = false; break;
				case 1: matlist[i].bBussineskApproach = true; break;
				default: matlist[i].bBussineskApproach = false; break;
				}
				// номер закона для зависимости динамической вязкости от напряжения сдвига
				fscanf_s(fp, "%d", &din);
				matlist[i].ilawmu = din;
				// минимальное значение динамической вязкости
				fscanf_s(fp, "%f", &fin);
				matlist[i].mumin = fin;
				// максимальное значение динамической вязкости
				fscanf_s(fp, "%f", &fin);
				matlist[i].mumax = fin;
				// параметры модельных законов для зависимости вязкости от напряжения сдвига
				fscanf_s(fp, "%f", &fin);
				matlist[i].Amu = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].Bmu = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].Cmu = fin;
				// показатель степени
				fscanf_s(fp, "%f", &fin);
				matlist[i].degreennmu = fin;

				// печать считанных значений на консоль
				//printf("%e %e %e %e %e\n", matlist[i].rho, matlist[i].cp, matlist[i].lam, matlist[i].mu, matlist[i].beta_t);
				printf("HEAT_CAPACITY\n");
				printf("t_C HEAT_CAPACITY\n");
				for (integer i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
					printf("%e %e\n", matlist[i].temp_cp[i_4], matlist[i].arr_cp[i_4]);
				}
				printf("lam\n");
				printf("t_C lam\n");
				for (integer i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
					printf("%e %e\n", matlist[i].temp_lam[i_4], matlist[i].arr_lam[i_4]);
				}
				printf("%e %e %e\n", matlist[i].rho, matlist[i].mu, matlist[i].beta_t);
				printf("%d %d %d\n", matlist[i].blibmat, matlist[i].ilibident, matlist[i].ilawmu); // bBoussinesq не печатается
				printf("%e %e %e %e %e %e\n", matlist[i].mumin, matlist[i].mumax, matlist[i].Amu, matlist[i].Bmu, matlist[i].Cmu, matlist[i].degreennmu);
			}

			// считывание блоков
			for (i = 0; i < lb; ++i) {

				fscanf_s(fp, "%d", &din);
				b[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				fscanf_s(fp, "%d", &din);
				CHECK_TYPE_GEOM(din);
				b[i].g.itypegeom = din; // 0 - Prism, 1 - Cylinder, 2 - Polygon
				fscanf_s(fp, "%d", &din);
				if (din == 1) {
					b[i].bvisible = true;
				}
				else {
					b[i].bvisible = false;
				}


				// геометрия
				fscanf_s(fp, "%f", &fin);
				b[i].g.xS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.yS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.zS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.xE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.yE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.zE = scale * fin;
				// swap
				if (b[i].g.xS > b[i].g.xE) {
					dbuf = b[i].g.xS;
					b[i].g.xS = b[i].g.xE;
					b[i].g.xE = dbuf;
				}
				if (b[i].g.yS > b[i].g.yE) {
					dbuf = b[i].g.yS;
					b[i].g.yS = b[i].g.yE;
					b[i].g.yE = dbuf;
				}
				if (b[i].g.zS > b[i].g.zE) {
					dbuf = b[i].g.zS;
					b[i].g.zS = b[i].g.zE;
					b[i].g.zE = dbuf;
				}

				// Cylinder

				fscanf_s(fp, "%d", &din);
				b[i].g.iPlane = din;

				fscanf_s(fp, "%f", &fin);
				b[i].g.xC = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.yC = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.zC = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.Hcyl = scale * fin;
				if (b[i].g.Hcyl < 0.0) {
					// ликвидируем отрицательную высоту цилиндра.
					switch (b[i].g.iPlane) {
					case XY_PLANE:
						b[i].g.zC += b[i].g.Hcyl;
						break;
					case XZ_PLANE:
						b[i].g.yC += b[i].g.Hcyl;
						break;
					case YZ_PLANE:
						b[i].g.xC += b[i].g.Hcyl;
						break;
					}
					b[i].g.Hcyl = fabs(b[i].g.Hcyl);
				}

				fscanf_s(fp, "%f", &fin);
				b[i].g.R_out_cyl = scale * fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.R_in_cyl = scale * fin;

				// Polygon
				fscanf_s(fp, "%d", &din);
				b[i].g.iPlane_obj2 = din;
				fscanf_s(fp, "%d", &din);
				b[i].g.nsizei = din;
				b[i].g.hi = new doublereal[b[i].g.nsizei];
				b[i].g.xi = new doublereal[b[i].g.nsizei];
				b[i].g.yi = new doublereal[b[i].g.nsizei];
				b[i].g.zi = new doublereal[b[i].g.nsizei];
				for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
					fscanf_s(fp, "%f", &fin);
					b[i].g.hi[i73] = scale * fin;
					fscanf_s(fp, "%f", &fin);
					b[i].g.xi[i73] = scale * fin;
					fscanf_s(fp, "%f", &fin);
					b[i].g.yi[i73] = scale * fin;
					fscanf_s(fp, "%f", &fin);
					b[i].g.zi[i73] = scale * fin;
					if (b[i].g.hi[i73] < 0.0) {
						// ликвидируем отрицательную высоту цилиндра.
						switch (b[i].g.iPlane_obj2) {
						case XY_PLANE:
							b[i].g.zi[i73] += b[i].g.hi[i73];
							break;
						case XZ_PLANE:
							b[i].g.yi[i73] += b[i].g.hi[i73];
							break;
						case YZ_PLANE:
							b[i].g.xi[i73] += b[i].g.hi[i73];
							break;
						}
						b[i].g.hi[i73] = fabs(b[i].g.hi[i73]);
					}
				}

				if ((b[i].g.itypegeom == POLYGON) && (b[i].g.nsizei > 0)) {
					// Мы заносим в окаймляющую призму границы полигона чтобы ускорить 
					// поиски в inmodel_temp быстро отсекая ячейки вне границ полигона.
					doublereal xmin53 = 1.0e30;
					doublereal ymin53 = 1.0e30;
					doublereal zmin53 = 1.0e30;
					doublereal xmax53 = -1.0e30;
					doublereal ymax53 = -1.0e30;
					doublereal zmax53 = -1.0e30;
					for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
						if (b[i].g.xi[i73] > xmax53) xmax53 = b[i].g.xi[i73];
						if (b[i].g.yi[i73] > ymax53) ymax53 = b[i].g.yi[i73];
						if (b[i].g.zi[i73] > zmax53) zmax53 = b[i].g.zi[i73];
						if (b[i].g.xi[i73] < xmin53) xmin53 = b[i].g.xi[i73];
						if (b[i].g.yi[i73] < ymin53) ymin53 = b[i].g.yi[i73];
						if (b[i].g.zi[i73] < zmin53) zmin53 = b[i].g.zi[i73];
					}
					switch (b[i].g.iPlane_obj2) {
					case XY_PLANE:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmin53 + b[i].g.hi[0];
						break;
					case XZ_PLANE:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymin53 + b[i].g.hi[0];
						break;
					case YZ_PLANE:
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.xS = xmin53;
						b[i].g.xE = xmin53 + b[i].g.hi[0];
						break;
					}
				}

				if (b[i].g.itypegeom == CYLINDER) {
					// Cylinder
					//printf("%e %e %e %e %e %e %e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE, b[i].g.xC, b[i].g.yC, b[i].g.zC, b[i].g.Hcyl, b[i].g.R_out_cyl, b[i].g.R_in_cyl);
					//system("PAUSE");
				}


				// Ввод предполагается корректным.
				// emissivity
				fscanf_s(fp, "%f", &fin);
				b[i].radiation.emissW = fin;
				fscanf_s(fp, "%f", &fin);
				b[i].radiation.emissE = fin;
				fscanf_s(fp, "%f", &fin);
				b[i].radiation.emissS = fin;
				fscanf_s(fp, "%f", &fin);
				b[i].radiation.emissN = fin;
				fscanf_s(fp, "%f", &fin);
				b[i].radiation.emissB = fin;
				fscanf_s(fp, "%f", &fin);
				b[i].radiation.emissT = fin;
				fscanf_s(fp, "%d", &din);
				if (din == 0) {
					// Блок не является вакуумным промежутком.
					b[i].radiation.binternalRadiation = false;
				}
				else {
					// блок является вакуумным промежутком.
					b[i].radiation.binternalRadiation = true;
					if (bvacuumPrism) {
						bdouble_vacuum_PRISM = true;
					}
					bvacuumPrism = true;
					// Вычисляем View Factors.
					calculate_view_factors(b[i]);
				}
				b[i].radiation.nodelistW = nullptr;
				b[i].radiation.nodelistE = nullptr;
				b[i].radiation.nodelistS = nullptr;
				b[i].radiation.nodelistN = nullptr;
				b[i].radiation.nodelistB = nullptr;
				b[i].radiation.nodelistT = nullptr;
				b[i].radiation.nodelistWsize = 0;
				b[i].radiation.nodelistEsize = 0;
				b[i].radiation.nodelistSsize = 0;
				b[i].radiation.nodelistNsize = 0;
				b[i].radiation.nodelistBsize = 0;
				b[i].radiation.nodelistTsize = 0;

				// идентификатор материала в базе материалов
				fscanf_s(fp, "%d", &din);
				b[i].imatid = din;

				fscanf_s(fp, "%d", &din);

				// bCylinderFixed
				if (din == 1) {
					b[i].CylinderFixed = true;
				}
				else {
					b[i].CylinderFixed = false;
				}

				// мощность тепловыделения
				//fscanf_s(fp, "%f", &fin);
				//b[i].Sc = fin;
				// 19 november 2016 температурно зависимая мощность тепловыделения.
				fscanf_s(fp, "%d", &din);
				b[i].n_Sc = din;
				b[i].arr_Sc = nullptr;
				b[i].temp_Sc = nullptr;
				b[i].arr_Sc = new doublereal[b[i].n_Sc];
				b[i].temp_Sc = new doublereal[b[i].n_Sc];
				if (b[i].temp_Sc == nullptr) {
					printf("problem memory allocation for temp_Sc\n");
					system("pause");
					exit(1);
				}
				if (b[i].arr_Sc == nullptr) {
					printf("problem memory allocation for arr_Sc\n");
					system("pause");
					exit(1);
				}
				// Объём полигона.
				doublereal vol_poly = Volume_polygon(b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2);

				for (integer i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					b[i].temp_Sc[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					if (fin != fin) {
						b[i].arr_Sc[i_4] = 0.0;
					}
					else {

						if (b[i].g.itypegeom == POLYGON) {

							// Для полигона передается из интерфейса просто мощность, а не удельная мощность.
							// Т.к. интерфейс не содержит функцию расчёта объёма полигона.
							// Для единообразия здесь мощность преобразуется в удельную мощность.
							if (vol_poly > 1.0e-30) {
								b[i].arr_Sc[i_4] = fin / vol_poly;
							}
							else {
								printf("error zero volume in polygon %s...\n", b[i].name);
								system("PAUSE");
								exit(1);
							}
						}
						else {
							b[i].arr_Sc[i_4] = fin;
						}
					}
				}

				// debug
				//if (fabs(b[i].Sc)>1.0e-30) {
					//printf("%e\n", b[i].Sc);
					//system("PAUSE");
				//}
				// стиль зависимости мощности тепловыделения в блоке от времени.
				// 0 - не зависит от времени, 1 - square wave зависимость, 
				// 2 - square wave 2 зависимость, 3 - hot cold режим,
				// 4 - piecewise const.
				fscanf_s(fp, "%d", &din);
				switch (din) {
				case 0: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
					break;
				case 1: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE;
					break;
				case 2: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE2;
					break;
				case 3: b[i].ipower_time_depend = POWER_TIME_DEPEND::HOT_COLD;
					break;
				case 4: b[i].ipower_time_depend = POWER_TIME_DEPEND::PIECEWISE_CONST;
					break;
				default: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
					break;
				}
				// тип блока
				fscanf_s(fp, "%d", &din);
				switch (din) {
				case 1: b[i].itype = PHYSICS_TYPE_IN_BODY::SOLID;
					break;
				case 2: b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
					break;
				case 3: b[i].itype = PHYSICS_TYPE_IN_BODY::FLUID;
					break;
				default:
					b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
					break;
				}

				// печать считанных значений на консоль
				printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				printf("%d %d\n", b[i].imatid, b[i].itype);
				printf("temperature depend power\n");
				printf("t_C power_W\n");
				for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
				}
			}

			// считывание источников тепла
			for (i = 0; i < ls; ++i) {

				fscanf_s(fp, "%d", &din);
				s[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				fscanf_s(fp, "%f", &fin);
				s[i].power = fin;
				fscanf_s(fp, "%d", &din);
				if (din == 0) {
					// const
					// задаётся постоянное значение рассеиваемой в тепло мощности.
					s[i].power_multiplyer = 1.0;
					s[i].bgarber_depend = false;
				}
				else if (din == 1) {
					// рассеиваемая в тепло мощность задаётся таблично
					// в зависимости от максимальной температуры и рабочего 
					// значения смещения стока.
					s[i].bgarber_depend = true;
					s[i].power_multiplyer = s[i].power;
					// мощность будет вычислена при ambient Temperature несколькими строчками позже.
				}
				fscanf_s(fp, "%d", &din);
				s[i].igarber_depend = din; // уникальный номер таблицы.
				fscanf_s(fp, "%f", &fin);
				s[i].roperation_offset_drain = fin; // рабочее значение смещения стока.
				//printf("offset drain is %e\n",s[i].roperation_offset_drain);
				//system("PAUSE");
				bool bsplinereadOk = true;
				if (s[i].bgarber_depend) {
					s[i].power = my_splain_interpol_power_table(gtdps[s[i].igarber_depend].intemp,
						gtdps[s[i].igarber_depend].inoffset_drain,
						gtdps[s[i].igarber_depend].rtemp,
						gtdps[s[i].igarber_depend].roffset_drain,
						gtdps[s[i].igarber_depend].rpower_table,
						operatingtemperature,
						s[i].roperation_offset_drain);
					if (bsplinereadOk) {
						// одиночная тестовая проверка сплайновой аппроксимации.
						printf("single test validation spline approximation...\n");
						printf("calculate initial power=%e\n", s[i].power);
						printf("please, press any key to continue...");
						// system("PAUSE");
						system("pause");

						bsplinereadOk = false;
					}
					s[i].power *= s[i].power_multiplyer; // домножение на корректирующий множитель.
				}



				fscanf_s(fp, "%d", &din);
				s[i].iPlane = din;
				// геометрия
				fscanf_s(fp, "%f", &fin);
				s[i].g.xS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				s[i].g.yS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				s[i].g.zS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				s[i].g.xE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				s[i].g.yE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				s[i].g.zE = scale * fin;



				// swap
				if (s[i].g.xS > s[i].g.xE) {
					dbuf = s[i].g.xS;
					s[i].g.xS = s[i].g.xE;
					s[i].g.xE = dbuf;
				}
				if (s[i].g.yS > s[i].g.yE) {
					dbuf = s[i].g.yS;
					s[i].g.yS = s[i].g.yE;
					s[i].g.yE = dbuf;
				}
				if (s[i].g.zS > s[i].g.zE) {
					dbuf = s[i].g.zS;
					s[i].g.zS = s[i].g.zE;
					s[i].g.zE = dbuf;
				}
				switch (s[i].iPlane) {
				case XY_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS) * fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS) * fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ_PLANE: s[i].square = fabs(s[i].g.yE - s[i].g.yS) * fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				printf("%e %d %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
			}

			// считывание твёрдых стенок

			// 19.01.2021
		    // Стенка заглушка с условием прилипания по скорости(нулевой вектор скорости), 
		    // нулевым тепловым потоком по температуре, 
		    // свободной free границе в механике без приложенной силы.
			w[lw].iunion_id = 0; // cabinet
			w[lw].ifamily = WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY;
			w[lw].Tamb = 0.0;
			w[lw].emissivity = 0.0;
			w[lw].film_coefficient = 0.0;
			w[lw].ViewFactor = 1.0;
			w[lw].hf = 0.0;
			w[lw].bsymmetry = false;
			w[lw].bpressure = false;
			w[lw].bopening = false;
			w[lw].Vx = 0.0;
			w[lw].Vy = 0.0;
			w[lw].Vz = 0.0;
			w[lw].P = 0.0;
			w[lw].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
			w[lw].xForce = 0.0;
			w[lw].yForce = 0.0;
			w[lw].zForce = 0.0;
			w[lw].iPlane = XY_PLANE;
			w[lw].g.xS = 0.0;
			w[lw].g.yS = 0.0;
			w[lw].g.zS = 0.0;
			w[lw].g.xE = 0.0;
			w[lw].g.yE = 0.0;
			w[lw].g.zE = 0.0;

			for (i = 0; i < lw; ++i) {

				fscanf_s(fp, "%d", &din);
				w[i].iunion_id = din; // 0 == Кабинет или >0 - значит номер АССЕМБЛЕСА которому принадлежит.

				fscanf_s(fp, "%d", &din);
				switch (din) {
				case 1 : w[i].ifamily= WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY; break;
				case 2: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
				case 3: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY; break;
				case 4: w[i].ifamily= WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY; break;
				default: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
				}
				//w[i].ifamily = din;
				switch (din) {
				case 1:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf_s(fp, "%f", &fin); // Stefan Bolcman
					// termostability wall
					w[i].emissivity = 0.0;
					w[i].film_coefficient = 0.0;
					fscanf_s(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = 1.0;
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // первого рода
				case 2:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = 0.0;
					fscanf_s(fp, "%f", &fin);  // Stefan Bolcman
					// adiabatic wall
					w[i].emissivity = 0.0;
					w[i].film_coefficient = 0.0;
					fscanf_s(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = 1.0;
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // однородное условие Неймана
				case 3:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf_s(fp, "%f", &fin); // Stefan Bolcman
					// Newton-Richman condition, film coefficient.
					w[i].emissivity = 0.0;
					w[i].film_coefficient = fin;
					fscanf_s(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = 1.0;
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Ньютон-Рихман.
				case 4:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf_s(fp, "%f", &fin); // Stefan Bolcman
					// Stefan - Bolcman condition
					w[i].emissivity = fin;
					w[i].film_coefficient = 0.0;
					fscanf_s(fp, "%f", &fin); // ViewFactor
					w[i].ViewFactor = fmax(0.0,fmin(fin,1.0));
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Стефан-Больцман.
				default: break;
				}
				fscanf_s(fp, "%d", &din);
				if (din == 1) w[i].bsymmetry = true; else w[i].bsymmetry = false;
				fscanf_s(fp, "%d", &din);
				if (din == 1) w[i].bpressure = true; else w[i].bpressure = false;
				fscanf_s(fp, "%d", &din);
				if (din == 1) w[i].bopening = true; else w[i].bopening = false;
				fscanf_s(fp, "%f", &fin);
				w[i].Vx = fin;
				fscanf_s(fp, "%f", &fin);
				w[i].Vy = fin;
				fscanf_s(fp, "%f", &fin);
				w[i].Vz = fin;
				fscanf_s(fp, "%f", &fin);
				w[i].P = fin;
				fscanf_s(fp, "%d", &din);
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
				switch (din) {
				case 0: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
					break;
				case 1: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FIXIT;
					break;
				case 2: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FIXIT;
					break;
				case 3: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FIXIT;
					break;
				case 4: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XY_FIXIT;
					break;
				case 5: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XZ_FIXIT;
					break;
				case 6: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::YZ_FIXIT;
					break;
				case 7: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::ALL_FIXIT;
					break;
				case 8: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FORCE;
					break;
				case 9: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FORCE;
					break;
				case 10: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FORCE;
					break;
				default:
					w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // Free all
					break;
				}
				fscanf_s(fp, "%f", &fin);
				w[i].xForce = fin;
				fscanf_s(fp, "%f", &fin);
				w[i].yForce = fin;
				fscanf_s(fp, "%f", &fin);
				w[i].zForce = fin;
				fscanf_s(fp, "%d", &din);
				w[i].iPlane = din;
				// геометрия
				fscanf_s(fp, "%f", &fin);
				w[i].g.xS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				w[i].g.yS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				w[i].g.zS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				w[i].g.xE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				w[i].g.yE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				w[i].g.zE = scale * fin;
				// swap
				if (w[i].g.xS > w[i].g.xE) {
					dbuf = w[i].g.xS;
					w[i].g.xS = w[i].g.xE;
					w[i].g.xE = dbuf;
				}
				if (w[i].g.yS > w[i].g.yE) {
					dbuf = w[i].g.yS;
					w[i].g.yS = w[i].g.yE;
					w[i].g.yE = dbuf;
				}
				if (w[i].g.zS > w[i].g.zE) {
					dbuf = w[i].g.zS;
					w[i].g.zS = w[i].g.zE;
					w[i].g.zE = dbuf;
				}

				printf("%d %e %e %d %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
			}


			// АСЕМБЛЕСЫ.
			// Новое считывание работает только с новым форматом файла premeshin.txt 
			// с текстовыми метками, данный код здесь устарел. 31.10.2021.
			fscanf_s(fp, "%d", &din);
			lu = din;
			if (lu == 0) {
				my_union = nullptr;
			}
			else {
				my_union = new UNION[lu];
				// инициализация.
				for (i = 0; i < lu; ++i) {
					my_union[i].f = nullptr;
					my_union[i].xpos = nullptr;
					my_union[i].ypos = nullptr;
					my_union[i].zpos = nullptr;
					my_union[i].xposadd = nullptr;
					my_union[i].yposadd = nullptr;
					my_union[i].zposadd = nullptr;
					my_union[i].iunion_parent = -1; // Cabinet.
					my_union[i].iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; // 2 - CoarseMeshGen
					my_union[i].inxadd = -1;
					my_union[i].inyadd = -1;
					my_union[i].inzadd = -1;
					my_union[i].flow_interior = 0;
					my_union[i].active = false;
				}
			}
			for (i = 0; i < lu; ++i) {

				fscanf(fp, "%d", &din);
				my_union[i].iunion_parent = din;

				fscanf_s(fp, "%f", &fin);
				my_union[i].xS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				my_union[i].xE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				my_union[i].yS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				my_union[i].yE = scale * fin;
				fscanf_s(fp, "%f", &fin);
				my_union[i].zS = scale * fin;
				fscanf_s(fp, "%f", &fin);
				my_union[i].zE = scale * fin;

				fscanf_s(fp, "%d", &din);
				my_union[i].id = din; // Уникальный идентификатор АССЕМБЛЕСА.
				fscanf_s(fp, "%d", &din);
				my_union[i].inx = din;
				fscanf_s(fp, "%d", &din);
				my_union[i].iny = din;
				fscanf_s(fp, "%d", &din);
				my_union[i].inz = din;

			}


			// считывание информации о наборе решаемых уравнений
			fscanf_s(fp, "%d", &din);
			eqin.itemper = din;
			fscanf_s(fp, "%d", &din);
			eqin.imaxflD = din;
			//printf("eqin.imaxflD= %d\n", eqin.imaxflD); system("pause");
			if (eqin.imaxflD == 0) {
				eqin.fluidinfo = nullptr;
			}
			else
			{
				// выделение оперативной памяти
				if (eqin.fluidinfo != nullptr) {
					delete eqin.fluidinfo;
					eqin.fluidinfo = nullptr;
				}
				eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
				for (i = 0; i < eqin.imaxflD; ++i) {
					// Считывание координат опорной точки
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].xc = scale * fin;
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].yc = scale * fin;
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].zc = scale * fin;
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].iflow = din;
					fscanf_s(fp, "%d", &din);
					switch (din) {
					case 0: eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
						break;
					case 1:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::TURBULENT;
						break;
					default:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
						break;
					}
					fscanf_s(fp, "%d", &din);
					switch (din) {
					case 0: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
						break;
					case 1: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::SMAGORINSKY;
						break;
					case 2: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RNG_LES;
						break;
					case 3: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_SPALART_ALLMARES;
						break;
					case 4: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_MENTER_SST;
						break;
					case 5: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_STANDART_K_EPS;
						break;
					case 6: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_LANGTRY_MENTOR_SST;
						break;
					default: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
						break;
					}
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].Cs = fin; // постоянная Смагоринского.
					fscanf_s(fp, "%d", &din);
					// учёт динамического определения квадрата постоянной Смагоринского.
					// включает Dynamic Subgrid Scale Model Германо 1991 года.
					if (din == 1) {
						eqin.fluidinfo[i].bDynamic_Stress = true;
					}
					else {
						eqin.fluidinfo[i].bDynamic_Stress = false;
					}
					fscanf_s(fp, "%d", &din);
					// включает ограничение сверху и снизу на возможные значения постоянной Смагоринского.
					if (din == 1) {
						eqin.fluidinfo[i].bLimiters_Cs = true;
					}
					else {
						eqin.fluidinfo[i].bLimiters_Cs = false;
					}
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].rminCs = fin; // минимальное возможное значение постоянной Смагоринского.
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].rmaxCs = fin; // максимальное возможное значение постоянной Смагоринского.
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].itypeFiltrGermano = din; // тип фильтра в модели Германо 1991 года.
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].roughness = 1.0e-6 * fin; // шероховатость стенки в м.
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].rRi_mult = fin; // множитель корректирующий турбулентное число Ричардсона.
					fscanf_s(fp, "%f", &fin);
					eqin.fluidinfo[i].rSelectiveAngle = fin; // пороговое значение угла в модели Selective Smagorinsky.
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].ipowerroughness = din; // показатель степени в модели учёта шероховатости.
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = din; // тип фильтра в модели Selective Smagorinsky.
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].bfdelta = din; // учёт неравномерности сетки.
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].bSmagorinskyLilly = din; // модель Смагоринского-Лиллу.
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].bsurface_roughness = din; // учёт шероховатости стенки.
					fscanf_s(fp, "%d", &din);
					// учёт течений с кривизной линий тока.
					if (din == 1) {
						eqin.fluidinfo[i].bSwirlAmendment = true;
					}
					else {
						eqin.fluidinfo[i].bSwirlAmendment = false;
					}
					fscanf_s(fp, "%d", &din);
					// учёт избирательности в модели Смагоринского
					if (din == 1) {
						eqin.fluidinfo[i].bSelectiveSmagorinsky = true;
					}
					else {
						eqin.fluidinfo[i].bSelectiveSmagorinsky = false;
					}

					// Параметры преобразователя картинок для отчетов.
					// 5.01.2018
					fscanf_s(fp, "%f", &fin);
					pfpir.fminimum = scale * fin;
					fscanf_s(fp, "%f", &fin);
					pfpir.fmaximum = scale * fin;
#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					switch (din) {
					case 0: pfpir.idir = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
					case 1: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					case 2: pfpir.idir = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
					default: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					}
					//pfpir.idir = din;

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					AMG1R6_LABEL = din;

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					if (din > 9999) {
						b_iluk_amg1r5_LABEL_D = true;
						nrd_LABEL = din-10000;
					}
					else {
						nrd_LABEL = din;
					}

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					if (din > 9999) {
						b_iluk_amg1r5_LABEL_U = true;
						nru_LABEL = din-10000;
					}
					else {
						nru_LABEL = din;
					}

#if doubleintprecision == 1
					fscanf_s(fp, "%f", &fin);
#else
					fscanf_s(fp, "%f", &fin);
#endif
					ecg2_LABEL = din;

#if doubleintprecision == 1
					fscanf_s(fp, "%f", &fin);
#else
					fscanf_s(fp, "%f", &fin);
#endif
					ewt2_LABEL = din;

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					number_processors_global_var = (int)(din);

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					idevice_Tesla = (int)(din);

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					if (din>=0) {
						number_iteration_SIMPLE_algorithm = (unsigned int)(din);
					}
					else {
						number_iteration_SIMPLE_algorithm =0;
						std::cout<<"WARNING : number_iteration_SIMPLE_algorithm =0;"<<std::endl;
				    }

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					switch (din) {
					case 0: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::NONE_only_amg1r5;
						break;
					case 1: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
						break;
					case 2: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::FGMRes_plus_amg1r5;
						break;
					case 3: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5;
						break;
					default:
						stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
						break;
					}


#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					BonLevelDrobim = din;

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					isizearrAdaptRegion = din;


					if (isizearrAdaptRegion > 0) {

						AdaptRegion = new ADAPTREGION[isizearrAdaptRegion];

						AdaptRegion[0].ilevelMeshAdaptRegion = -1;
						AdaptRegion[1].ilevelMeshAdaptRegion = -1;
						AdaptRegion[2].ilevelMeshAdaptRegion = -1;

						for (int i_63 = 0; i_63 < isizearrAdaptRegion; ++i_63) {

							// Адаптация сетки региона номер 1.
							if (i_63 == 0) {

#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[0].xS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[0].xE = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[0].yS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[0].yE = scale * fin;

#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[0].zS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[0].zE = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%lld", &din);
#else
								fscanf_s(fp, "%d", &din);
#endif
								AdaptRegion[0].ilevelMeshAdaptRegion = din;

							}


							if (i_63 == 1) {


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[1].xS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[1].xE = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[1].yS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[1].yE = scale * fin;

#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[1].zS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[1].zE = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%lld", &din);
#else
								fscanf_s(fp, "%d", &din);
#endif
								AdaptRegion[1].ilevelMeshAdaptRegion = din;

							}


							if (i_63 == 2) {

#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[2].xS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[2].xE = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[2].yS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[2].yE = scale * fin;

#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[2].zS = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%f", &fin);
#else
								fscanf_s(fp, "%f", &fin);
#endif
								AdaptRegion[2].zE = scale * fin;


#if doubleintprecision == 1
								fscanf_s(fp, "%lld", &din);
#else
								fscanf_s(fp, "%d", &din);
#endif
								AdaptRegion[2].ilevelMeshAdaptRegion = din;


							}

	}

}



#if doubleintprecision == 1
					fscanf_s(fp, "%f", &fin);
#else
					fscanf_s(fp, "%f", &fin);
#endif

					// Найдено успешно.
					// Свободный параметр передаваемый из интерфейса и используемый для отладки.
					free_debug_parametr1 = static_cast<doublereal>(fin);
					

				}
			}

			if (fp != NULL) {
				fclose(fp); // закрытие файла
				fp = NULL;
			}
		}
	}


#else
// eqin - информация о наборе решаемых уравнений.

// dgx, dgy, dgz - вектор силы тяжести.
// inx, iny, inz - количество точек по каждой из осей.

FILE* fp;
//errno_t err1;
//if ((err1 = fopen_s(&fp, fname, "r")) != 0) 
if ((fopen_s(&fp, fname, "r")) != 0)
{
	printf("Error!!! No input File premeshin.txt.\n");
	printf("You must use the graphical user interface\n");
	printf("AliceMesh_v0_45.exe in Delphi, which will \n");
	printf("prepare the model and write the premeshin.txt\n");
	printf("file in the required format.\n");
	//system("PAUSE");
	system("pause");
	// Если файла premeshin.txt нет то мы не можем ничего обрабатывать
	// т.к. элементарно отсутствуют входные данные. Поэтому мы выходим из 
	// приложения.
	exit(1);
}
else
{
	if (fp != nullptr) {
		float fin = 0.0;
		integer din = 0;
		//doublereal scale = 1.0;
		doublereal dbuf; // для упорядочивания в порядке возрастания



		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
		case 1: ionly_solid_visible= WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE; break;
		default: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
		}
		//ionly_solid_visible = din;
		switch (ionly_solid_visible) {
		case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE:
			std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
			break;
		case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE:
			std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
			break;
		}
		fscanf_s(fp, "%f", &fin);
		const doublereal scale = fin;
#ifndef NO_OPENGL_GLFW
		scale_all = (GLfloat) (1.0f / scale);
#endif
		fscanf_s(fp, "%lld", &din);
		lmatmax = din;
		fscanf_s(fp, "%lld", &din);
		lb = (int)(din);
		fscanf_s(fp, "%lld", &din);
		ls = (int)(din);
		fscanf_s(fp, "%lld", &din);
		lw = (int)(din);
		fscanf_s(fp, "%lld", &din);
		ltdp = din; // количество уникальных данных с табличными данными по зависимости рассеиваемой мощности от температуры.
		

					// Считываем значение вектора силы тяжести:
		fscanf_s(fp, "%f", &fin);
		dgx = fin;
		fscanf_s(fp, "%f", &fin);
		dgy = fin;
		fscanf_s(fp, "%f", &fin);
		dgz = fin;

		// считываем количество точек на каждой координатной оси
		fscanf_s(fp, "%lld", &din);
		inx = (int)(din);
		fscanf_s(fp, "%lld", &din);
		iny = (int)(din);
		fscanf_s(fp, "%lld", &din);
		inz = (int)(din);

		fscanf_s(fp, "%f", &fin);
		operatingtemperature = fin; // Operating Temperature
		operating_temperature_for_film_coeff = fin;

		// инициализация компонент скорости константой.
		// Единые значения для всей расчётной области.
		// initialization value.
		fscanf_s(fp, "%f", &fin);
		starting_speed_Vx = fin;
		fscanf_s(fp, "%f", &fin);
		starting_speed_Vy = fin;
		fscanf_s(fp, "%f", &fin);
		starting_speed_Vz = fin;

		// Считываем координаты опорной точки через которую проходит пользовательская линия Variation Plot.
		fscanf_s(fp, "%f", &fin);
		Tochka_position_X0_for_XY_Plot = scale*fin;
		fscanf_s(fp, "%f", &fin);
		Tochka_position_Y0_for_XY_Plot = scale*fin;
		fscanf_s(fp, "%f", &fin);
		Tochka_position_Z0_for_XY_Plot = scale*fin;
		// Направление линии совпадает с направлением 
		// одной из осей декартовой прямоугольной системы координат:
		// 0 - Ox, 1 - Oy, 2 - Oz.
		fscanf_s(fp, "%lld", &din);
		switch (din) {
			case 0: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
			case 1: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
			case 2: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
			default: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL;	break;
		}
		//idirectional_for_XY_Plot = din;

		fscanf_s(fp, "%f", &fin);
		etalon_max_size_ratio = fin; // подробность расчётной сетки.

		fscanf_s(fp, "%f", &fin);
		etalon_max_size_ratio2 = fin; // Критерий качества расчётной сетки на основе FlowVision.

		
		// 0.	none
		// 1.	Snap to grid
		// 2.	Snap to grid ALICE
		// 3.	Snap to grid ++
		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: bsnap_TO_global = 0;  break;
		case 1: bsnap_TO_global = 1;  break;
		case 2: bsnap_TO_global = 2;  break;
		case 3: bsnap_TO_global = 3;  break;
		default: bsnap_TO_global = 1;  break;
		}

		fscanf_s(fp, "%lld", &din);
		iswitchsolveramg_vs_BiCGstab_plus_ILU2 = din; // Выбор решающего устройства: либо amg1r5 либо BiCGStab+ILU2.

		fscanf_s(fp, "%lld", &din);
		switch (din) {
			case 0: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
			case 1: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER; break;
			case 2: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::CAMG_RUMBA_v0_14_SECOND_T_SOLVER; break;
			case 3: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMG1R5_SECOND_T_SOLVER; break;
			case 4: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMGCL_SECONT_T_SOLVER; break;
			default: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
		}
		//iswitchsolveramg_vs_BiCGstab_plus_ILU6 = din; // Выбор решающего устройства: либо РУМБА0.14 либо BiCGStab+ILU6.

		fscanf_s(fp, "%lld", &din);
		if (din == 1) {
			// SIMPLEC algorithm.
			iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby;
		}
		else {
			// SIMPLE algorithm 1972.
			iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;
		}


		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0:  iprefix_Scheme_Flow = CR;
			break;
		case 1:  iprefix_Scheme_Flow = UDS;
			break;
		case 2:  iprefix_Scheme_Flow = COMB;
			break;
		case 3:  iprefix_Scheme_Flow = POLY;
			break;
		case 4:  iprefix_Scheme_Flow = EXP;
			break;
		case 5:  iprefix_Scheme_Flow = BULG;
			break;
		case 6:  iprefix_Scheme_Flow = POW;
			break;
		default:
			iprefix_Scheme_Flow = UDS;
			break;
		}

		fscanf_s(fp, "%lld", &din);
		// выбор схемы для потока жидкости.
		// Внимание эти определения должны полностью соответствовать 
		// определениям в файле my_approx_convective2.c
		switch (din) {
		case 1: iFLOWScheme = UNEVEN_MUSCL; break; // MUSCL 2
		case 2: iFLOWScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
		case 3: iFLOWScheme = UNEVEN_HLPA; break; // HLPA 2
		case 4: iFLOWScheme = UNEVEN_SMART; break; // SMART 3
		case 5: iFLOWScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
		case 6: iFLOWScheme = UNEVEN_SMARTER; break; // SMARTER 3
		case 7: iFLOWScheme = UNEVEN_LPPA; break; // LPPA 3
		case 8: iFLOWScheme = UNEVEN_VONOS; break; // VONOS 3
		case 9: iFLOWScheme = UNEVEN_STOIC; break; // STOIC
		case 10: iFLOWScheme = UNEVEN_CLAM; break; // CLAM
		case 11: iFLOWScheme = UNEVEN_OSHER; break; // OSHER
		case 12: iFLOWScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
		case 13: iFLOWScheme = UNEVEN_SUPER_C; break; // SUPER_C
		case 14: iFLOWScheme = UNEVEN_ISNAS; break; // ISNAS
		case 15: iFLOWScheme = UNEVEN_CUBISTA; break; // CUBISTA
		case 16: iFLOWScheme = UNEVEN_GAMMA; break; // GAMMA
		case 17: iFLOWScheme = UNEVEN_COPLA; break; // COPLA
		case 18: iFLOWScheme = UNEVEN_SECBC; break; // SECBC
		case 19: iFLOWScheme = UNEVEN_SGSD; break; // SGSD
		case 20: iFLOWScheme = UNEVEN_WENO5; break; // WENO5
		default: iFLOWScheme = UDS; break; // UDS самая стабильная схема.
		}

		fscanf_s(fp, "%lld", &din);
		// выбор схемы для температуры в потоке жидкости.
		// Внимание эти определения должны полностью соответствовать 
		// определениям в файле my_approx_convective2.c
		switch (din) {
		case 1: iTEMPScheme = UNEVEN_MUSCL; break; // MUSCL 2
		case 2: iTEMPScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
		case 3: iTEMPScheme = UNEVEN_HLPA; break; // HLPA 2
		case 4: iTEMPScheme = UNEVEN_SMART; break; // SMART 3
		case 5: iTEMPScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
		case 6: iTEMPScheme = UNEVEN_SMARTER; break; // SMARTER 3
		case 7: iTEMPScheme = UNEVEN_LPPA; break; // LPPA 3
		case 8: iTEMPScheme = UNEVEN_VONOS; break; // VONOS 3
		case 9: iTEMPScheme = UNEVEN_STOIC; break; // STOIC
		case 10: iTEMPScheme = UNEVEN_CLAM; break; // CLAM
		case 11: iTEMPScheme = UNEVEN_OSHER; break; // OSHER
		case 12: iTEMPScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
		case 13: iTEMPScheme = UNEVEN_SUPER_C; break; // SUPER_C
		case 14: iTEMPScheme = UNEVEN_ISNAS; break; // ISNAS
		case 15: iTEMPScheme = UNEVEN_CUBISTA; break; // CUBISTA
		case 16: iTEMPScheme = UNEVEN_GAMMA; break; // GAMMA
		case 17: iTEMPScheme = UNEVEN_COPLA; break; // COPLA
		case 18: iTEMPScheme = UNEVEN_SECBC; break; // SECBC
		case 19: iTEMPScheme = UNEVEN_SGSD; break; // SGSD
		case 20: iTEMPScheme = UNEVEN_WENO5; break; // WENO5	
		default: iTEMPScheme = UDS; break; // UDS самая стабильная схема.
		}



		// Выбор сеточного генератора.
		fscanf_s(fp, "%lld", &din);
		switch (din) {
				case 0: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; break;
				case 1: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER; break;
				case 2: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
				default: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
		}
		//iswitchMeshGenerator = din;

		fscanf_s(fp, "%lld", &din);
		steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY;
		if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 5) || (din == 6) ||
			(din == 7) || (din == 8) || (din == 9) || (din == 10) || (din == 11) || (din == 12) || (din == 13)) {
			// 0 - thermal only steady state calculation,
			// 1 - thermal only unsteady calculation,
			// 2 - mesh generator only.
			// 3 - fluid dynamic steady state.
			// 5 - Static Structural (Thermal solver #2)
			// 6 - Thermal Stress
			// 7 - Unsteady thermal solver #2
			// 8 - Visualisation only
			// 9 - cfd unsteady fluid dynamic.
			// 10 - NETWORK_T Графовый метод решения уравнения теплопроводности.
			// 11 - UNSTEADY NETWORK_T Нестационарный графовый метод решения уравнения теплопроводности.
			// 12 - Нестационарная механика,
			// 13 - Нестационарная механика совместно с нестационарной теплопередачей.
			switch(din) {
				case 0: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE; break;
				case 1: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE; break;
				case 2: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
				case 3: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_STEADY; break;
				case 5: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL; break;
				case 6: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
				case 7: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER; break;
				case 8: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT; break;
				case 9: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY; break;
				case 10: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T;  break;
				case 11: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY; break;
				case 12: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL;  break;
				case 13: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
				default: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
			}
			//steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady - 1 calculation.
			//printf("steady_or_unsteady_global_determinant =%lld\n",din);
			//system("PAUSE");
		}
		else {
			printf("error input parametr steady or unsteady calculation\n");
			system("PAUSE");
			exit(1);
		}

		fscanf_s(fp, "%lld", &din);
		if ((din == 0) || (din == 1) || (din == 2) || (din == 3)|| (din == 4)) {
			switch (din) {
			case 0: glTSL.id_law = TIME_STEP_lAW_SELECTOR::LINEAR;
				break;
			case 1: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE;
				break;
			case 2: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE2;
				break;
			case 3: glTSL.id_law = TIME_STEP_lAW_SELECTOR::HOT_COLD;
				break;
			case 4: glTSL.id_law = TIME_STEP_lAW_SELECTOR::PIECEWISE_CONSTANT;
				break;
			}
			
		}
		else {
			printf("error input parametr timestep law\n");
			system("PAUSE");
			exit(1);
		}
		fscanf_s(fp, "%f", &fin);
		if ((fin < 0.0) || (fin >= 1.0)) {
			printf("error input parametr timestep law Factor a\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.Factor_a_for_Linear = fin; // Factor_a
		fscanf_s(fp, "%f", &fin);
		if (fin < 0.0) {
			printf("error input parametr timestep law tau must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau = fin; // длительность импульса.
		fscanf_s(fp, "%f", &fin);
		glTSL.Q = fin;  // Скважность.
						// Параметры импульсного режима для SquareWave 2 режима.
		fscanf_s(fp, "%f", &fin);
		if ((fin < 0.0) || (fin >= 1.0)) {
			printf("error input parametr timestep law SquareWave2 multiplyer\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.m1 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin < 0.0) {
			printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau1 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin < 0.0) {
			printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau2 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin < 0.0) {
			printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau_pause = fin;

		fscanf_s(fp, "%f", &fin);
		if ((fin < 0.0) || (fin > 1.0)) {
			printf("error input parametr timestep law SquareWave2 off_mulptiplyer must be [0..1]\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.off_multiplyer = fin;

		fscanf_s(fp, "%lld", &din);
		glTSL.n_cycle = din;
		fscanf_s(fp, "%f", &fin);
		if (fin < 0.0) {
			printf("error input parametr timestep law SquareWave2 Period must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.T_all = fin;
		doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
		if (t_pause_gl < 0.0) {
			printf("error in parameters Square Wave 2 time step law.\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}

		fscanf_s(fp, "%f", &fin);
		if (fin < 0.0) {
			printf("error input parametr on_time_double_linear law hot cold reshime must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.on_time_double_linear = fin;

		// Время окончания нестационарного моделирования при расчёте теплопередачи в твёрдом теле.
		fscanf_s(fp, "%f", &fin);
		globalEndTimeUnsteadyTemperatureCalculation = fin;


		// Newton-Richman condition.
		fscanf_s(fp, "%lld", &din);
		switch (din) {
			case 0: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
            case 1: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC; break;
	        case 2: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC; break;
		    case 3: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC; break;
			default :adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
		}
		//adiabatic_vs_heat_transfer_coeff = din; // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
		fscanf_s(fp, "%f", &fin);
		film_coefficient = fin;
		// AЛИС сетка
		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			// Обычная структурированная сетка.
			b_on_adaptive_local_refinement_mesh = false;
		}
		else {
			// АЛИС
			b_on_adaptive_local_refinement_mesh = true;
		}
		fscanf_s(fp, "%lld", &din);
		switch (din) {
			case 0: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
			case 1: itype_ALICE_Mesh= TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH; break;
			default: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
		}
		//itype_ALICE_Mesh = din;

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.m_restart = static_cast<int>(din);
		// classical algebraic multigrid parameters:
		// only for my_agregat_amg.cu.
		// only for my_agregat_amg.cu.
		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
			break;
		case 1: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::QUICK_SORT;
			break;
		case 2: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::HEAP_SORT;
			break;
		case 3: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::TIM_SORT;
			break;
		default:
			my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
			break;
		}
		 
		fscanf_s(fp, "%lld", &din);
		//my_amg_manager.maximum_levels = din;
		my_amg_manager.maximum_delete_levels_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.maximum_delete_levels_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.maximum_delete_levels_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.maximum_delete_levels_Stress = din;

		// type interpolation procedure:
		//fscanf_s(fp, "%d", &din);
		//my_amg_manager.number_interpolation_procedure = din;

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.number_interpolation_procedure_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.number_interpolation_procedure_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.number_interpolation_procedure_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.number_interpolation_procedure_Stress = din;

		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
			break;
		case 1: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		case 2: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
			break;
		case 3: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
			break;
		case 4: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
			break;
		case 5: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
			break;
		case 6: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
			break;
		default:
			my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		}
		
		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
			break;
		case 1: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		case 2: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
			break;
		case 3: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
			break;
		case 4: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
			break;
		case 5: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
			break;
		case 6: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
			break;
		default:
			my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		}

		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
			break;
		case 1: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		case 2: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
			break;
		case 3: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
			break;
		case 4: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
			break;
		case 5: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
			break;
		case 6: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
			break;
		default:
			my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		}

		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
			break;
		case 1: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		case 2: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
			break;
		case 3: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
			break;
		case 4: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
			break;
		case 5: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
			break;
		case 6: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
			break;
		default:
			my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
			break;
		}

		fscanf_s(fp, "%lld", &din);
		//my_amg_manager.itypemodifyinterpol = din;
		//my_amg_manager.baglomeration_with_consistency_scaling = din;
		my_amg_manager.bdiagonal_dominant = din;
		fscanf_s(fp, "%lld", &din);
		//my_amg_manager.inumberadaptpass = din;


		// 23.02.2018
		// print matrix portrait
		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bTemperatureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
		}
		else {
			my_amg_manager.bTemperatureMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
		}
		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bSpeedMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
		}
		else {
			my_amg_manager.bSpeedMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
		}		
		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bPressureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
		}
		else {
			my_amg_manager.bPressureMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
		}
		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bStressMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
		}
		else {
			my_amg_manager.bStressMatrixPortrait = true; // false - NO_PRINT, true - PRINT.
		}

		// truncation of interpolation:
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.itruncation_interpolation_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.itruncation_interpolation_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.itruncation_interpolation_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.itruncation_interpolation_Stress = din;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.truncation_interpolation_Temperature = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.truncation_interpolation_Speed = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.truncation_interpolation_Pressure = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.truncation_interpolation_Stress = fin;

		// number nFinnest sweeps:
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nFinnest_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nFinnest_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nFinnest_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nFinnest_Stress = din;

		// number presweeps:
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu1_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu1_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu1_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu1_Stress = din;

		// number postsweeps:
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Stress = din;

		// memory size:
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Stress = din;


		// Параметр верхней релаксации в сглаживателе.
		//fscanf_s(fp, "%f", &fin);
		//my_amg_manager.gold_const_Temperature = fin;
		//fscanf_s(fp, "%f", &fin);
		//my_amg_manager.gold_const_Speed = fin;
		//fscanf_s(fp, "%f", &fin);
		//my_amg_manager.gold_const_Pressure = fin;
		//fscanf_s(fp, "%f", &fin);
		//my_amg_manager.gold_const_Stress = fin;

		// использовать ли ilu2 smoother.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.ilu2_smoother_Temperature = din;
		/*
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Temperature = din;
			my_amg_manager.b_gmresTemp = true;
		}
		*/

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.ilu2_smoother_Speed = din;
		/*
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Speed = din;
			my_amg_manager.b_gmresSpeed = true;
		}
		*/

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.ilu2_smoother_Pressure = din;
		/*
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Pressure = din;
			my_amg_manager.b_gmresPressure = true;
		}
		*/

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.ilu2_smoother_Stress = din;
		/*
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Stress = din;
			my_amg_manager.b_gmresStress = true;
		}
		*/
		my_amg_manager.gold_const_Temperature = return_gold_const(my_amg_manager.ilu2_smoother_Temperature);
		my_amg_manager.gold_const_Speed = return_gold_const(my_amg_manager.ilu2_smoother_Speed);
		my_amg_manager.gold_const_Pressure = return_gold_const(my_amg_manager.ilu2_smoother_Pressure);
		my_amg_manager.gold_const_Stress = return_gold_const(my_amg_manager.ilu2_smoother_Stress);


		fscanf_s(fp, "%lld", &din);
		my_amg_manager.Chebyshev_degree = din;

		if (!((din == 5) || (din == 8) || (din == 16) || (din == 25) || (din == 32) || (din == 64) || (din == 128) || (din == 256) || (din == 512))) {
			std::cout << "Chebyshev degree undefined. Chebyshev degree mast be equal: 5, 8, 16, 25, 32, 64, 128, 256, 512 in this programm. Chebyshev degree == " << din << std::endl;
			system("pause");
			my_amg_manager.Chebyshev_degree = 5;
		}

		// strength threshold:
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Temperature = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Speed = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Pressure = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Stress = fin;

		// magic threshold:
		//fscanf_s(fp, "%f", &fin);
		//my_amg_manager.magic = fin;

		// magic <=> F_to_F
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.F_to_F_Temperature = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.F_to_F_Speed = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.F_to_F_Pressure = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.F_to_F_Stress = fin;

		// AMG Splitting (coarsening)
		// Способ построения C-F разбиения: 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
		// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
			break;
		case 1: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
			break;
		case 2: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
			break;
		case 3: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
			break;
		case 4: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
			break;
		case 5: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
			break;
		case 6: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
			break;
		case 7: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
			break;
		case 8: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		case 9: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
			break;
		case 10: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
			break;
		case 11: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
			break;
		default:
			my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		}
		
		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
			break;
		case 1: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
			break;
		case 2: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
			break;
		case 3: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
			break;
		case 4: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
			break;
		case 5: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
			break;
		case 6: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
			break;
		case 7: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
			break;
		case 8: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		case 9: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
			break;
		case 10: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
			break;
		case 11: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
			break;
		default:
			my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		}

		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
			break;
		case 1: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
			break;
		case 2: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
			break;
		case 3: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
			break;
		case 4: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
			break;
		case 5: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
			break;
		case 6: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
			break;
		case 7: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
			break;
		case 8: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		case 9: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
			break;
		case 10: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
			break;
		case 11: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
			break;
		default:
			my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		}


		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
			break;
		case 1: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
			break;
		case 2: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
			break;
		case 3: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
			break;
		case 4: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
			break;
		case 5: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
			break;
		case 6: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
			break;
		case 7: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
			break;
		case 8: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		case 9: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
			break;
		case 10: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
			break;
		case 11: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
			break;
		default:
			my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
			break;
		}

		// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
		// Если din==1, Stabilization BiCGStab.
		// 8.01.2017 Метод Хенка ван дер Ворста BiCGStab 1992.
		// предобусловленный алгебраическим многосеточным методом.
		// 9.01.2018 Если din==2, FGMRes предобусловленный алгебраическим многосеточным методом.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.istabilizationTemp = din; // 0 - none
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.istabilizationSpeed = din; // 0 - none
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.istabilizationPressure = din; // 0 - none
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.istabilizationStress = din; // 0 - none
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.ipatch_number = din; // 0 - патч не применяется.

		// Печать лога на консоль.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iprint_log_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iprint_log_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iprint_log_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iprint_log_Stress = din;


		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bthreshold_Temperature_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}
		else {
			my_amg_manager.bthreshold_Temperature_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}

		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bthreshold_Speed_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}
		else {
			my_amg_manager.bthreshold_Speed_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}

		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bthreshold_Pressure_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}
		else {
			my_amg_manager.bthreshold_Pressure_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}

		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bthreshold_Stress_auto = false; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}
		else {
			my_amg_manager.bthreshold_Stress_auto = true; // false - используется константа пользователя заданная вручную, true - автоматическая настройка на каждом уровне.
		}


		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bcf_reorder_Temperature = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}
		else {
			my_amg_manager.bcf_reorder_Temperature = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}

		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bcf_reorder_Speed = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}
		else {
			my_amg_manager.bcf_reorder_Speed = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}


		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bcf_reorder_Pressure = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}
		else {
			my_amg_manager.bcf_reorder_Pressure = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}

		fscanf_s(fp, "%lld", &din);
		if (din == 0) {
			my_amg_manager.bcf_reorder_Stress = false; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}
		else {
			my_amg_manager.bcf_reorder_Stress = true; // false - CF упорядочивание узлов при итерированиии не используется, true - Используется CF упорядочивание узлов при итерировании.
		}


		fscanf_s(fp, "%lld", &din);
		my_amg_manager.amgcl_smoother = din;

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.amgcl_selector = din;

		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
			break;
		case 1: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::FGMRes;
			break;
		default: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
			break;
		}
		 

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.lfil = din;

		fscanf_s(fp, "%lld", &din);
		ireconstruction_free_construct_alloc = din;

		fscanf_s(fp, "%lld", &din);
		ianimation_write_on = din;

		// выделение оперативной памяти.
		gtdps = new TEMP_DEP_POWER[ltdp];
		matlist = new TPROP[lmatmax];
		b = new BLOCK[lb];
		s = new SOURCE[ls];
		w = new WALL[lw+1];// +1 это гран условие по умолчанию, заглушка, чтобы избежать неопределенности.
		integer i = 0; // счётчик цикла for

		for (i = 0; i < ltdp; ++i) {
			// считывание имён файлов.
			gtdps[i].sname = new char[100]; // выделение памяти
			fscanf_s(fp, "%s", gtdps[i].sname, 100);
			//printf("%s",gtdps[i].sname);
			//system("PAUSE");
			// построение таблицы в памяти.
			my_read_power_table(gtdps[i].sname, gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
			printf("printeger table\n"); // debug
			my_print_table(gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
			printf("Please, press any key to continue...\n");
			//system("PAUSE");
			system("pause");

		}



		// считывание базы материалов
		for (i = 0; i < lmatmax; ++i) {
			// свойства материалов:
			// плотность
			fscanf_s(fp, "%f", &fin);
			matlist[i].rho = fin;
			// теплоёмкость при постоянном давлении
			//fscanf(fp, "%f", &fin);
			//matlist[i].cp = fin;
			fscanf_s(fp, "%lld", &din);
			matlist[i].n_cp = din;
			matlist[i].arr_cp = nullptr;
			matlist[i].temp_cp = nullptr;
			matlist[i].arr_cp = new float[matlist[i].n_cp];
			matlist[i].temp_cp = new float[matlist[i].n_cp];
			if (matlist[i].temp_cp == nullptr) {
				printf("problem memory allocation for temp_cp\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_cp == nullptr) {
				printf("problem memory allocation for arr_cp\n");
				system("pause");
				exit(1);
			}
			for (integer i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				matlist[i].temp_cp[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].arr_cp[i_4] = fin;
			}
			// теплопроводность
			//fscanf(fp, "%f", &fin);
			//matlist[i].lam = fin;
			fscanf_s(fp, "%lld", &din);
			matlist[i].n_lam = din;
			matlist[i].arr_lam = nullptr;
			matlist[i].temp_lam = nullptr;
			matlist[i].arr_lam = new float[matlist[i].n_lam];
			matlist[i].temp_lam = new float[matlist[i].n_lam];
			if (matlist[i].temp_lam == nullptr) {
				printf("problem memory allocation for temp_lam\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_lam == nullptr) {
				printf("problem memory allocation for arr_lam\n");
				system("pause");
				exit(1);
			}
			for (integer i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				matlist[i].temp_lam[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].arr_lam[i_4] = fin;
			}
			// ортотропность теплопроводности:
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_x = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_y = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_z = fin;

			// 28.08.2020
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_x_beta_t_solid = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_y_beta_t_solid = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_z_beta_t_solid = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_x_Young_Module = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_y_Young_Module = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_z_Young_Module = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_Poisson_ratio_xy = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_Poisson_ratio_xz = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_Poisson_ratio_yz = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_Poisson_ratio_yx = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_Poisson_ratio_zx = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_Poisson_ratio_zy = fin;
			fscanf_s(fp, "%lld", &din);
			if (din == 1) {
				matlist[i].bActive_ShearModule = true;
			}
			else {
				matlist[i].bActive_ShearModule = false;
			}
			fscanf_s(fp, "%f", &fin);
			matlist[i].ShearModule_xy = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].ShearModule_yz = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].ShearModule_xz = fin;

			// 5.08.2017.
			// Коэффициенты для задачи упругости.
			// Модуль Юнга и коэффициент Пуассона.
			// Алюминий стр. 232 В.Н.Сидоров, В.В. Вершинин. Метод конечных элементов в расчёте сооружений.
			//doublereal Poissonratio = 0.33;
			//doublereal Youngmodule = 71.7e9;
			// steel
			// Tensile Yield Strength 4.12E+8
			// Compressive Yield Strength 4.12E+8
			// Tensile Ultimate Strength 5.7E+8
			//doublereal Poissonratio = 0.154;
			//doublereal Youngmodule = 217.5E9;
			//doublereal Poissonratio = 0.3;
			//doublereal Youngmodule = 200.0E9;
			//fscanf_s(fp, "%f", &fin);
			//doublereal Poissonratio = fin;
			//if (fabs(fin) < 1.0e-30) {
				//printf("ERROR !!! Zero Poisson Ratio in model. \n");
				//printf("ERROR !!! Model is incorrect...\n");
				//system("PAUSE");
				//exit(1);
			//}
			//matlist[i].Poisson_ratio = fin;

			fscanf_s(fp, "%lld", &din);
			matlist[i].n_Poisson_ratio = din;
			matlist[i].arr_Poisson_ratio = nullptr;
			matlist[i].temp_Poisson_ratio = nullptr;
			matlist[i].arr_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
			matlist[i].temp_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
			if (matlist[i].temp_Poisson_ratio == nullptr) {
				printf("problem memory allocation for temp_Poisson_ratio\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_Poisson_ratio == nullptr) {
				printf("problem memory allocation for arr_Poisson_ratio\n");
				system("pause");
				exit(1);
			}
			for (integer i_4 = 0; i_4 < matlist[i].n_Poisson_ratio; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				matlist[i].temp_Poisson_ratio[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].arr_Poisson_ratio[i_4] = fin;
			}


			//fscanf_s(fp, "%f", &fin);
			//doublereal Youngmodule = fin*1e9;
			//if (fabs(fin) < 1.0e-30) {
				//printf("ERROR !!! Zero Young Module in model. \n");
				//printf("ERROR !!! Model is incorrect...\n");
				//system("PAUSE");
				//exit(1);
			//}
			//matlist[i].Young_Module = fin * 1e9;

			fscanf_s(fp, "%lld", &din);
			matlist[i].n_YoungModule = din;
			matlist[i].arr_Young_Module = nullptr;
			matlist[i].temp_Young_Module = nullptr;
			matlist[i].arr_Young_Module = new float[matlist[i].n_YoungModule];
			matlist[i].temp_Young_Module = new float[matlist[i].n_YoungModule];
			if (matlist[i].temp_Young_Module == nullptr) {
				printf("problem memory allocation for temp_Young_Module\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_Young_Module == nullptr) {
				printf("problem memory allocation for arr_Young_Module\n");
				system("pause");
				exit(1);
			}
			for (integer i_4 = 0; i_4 < matlist[i].n_YoungModule; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				matlist[i].temp_Young_Module[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].arr_Young_Module[i_4] = static_cast <float>(fin * 1.0e+9);
			}

			//fscanf_s(fp, "%f", &fin);
			// beta_t_solid*1E-6
			//matlist[i].n_beta_t_solid = 1;
			//matlist[i].temp_beta_t_solid = new doublereal[1];
			//matlist[i].temp_beta_t_solid[0] = 25.0;
			//matlist[i].arr_beta_t_solid = new doublereal[1];
			//matlist[i].arr_beta_t_solid[0] = fin * 1E-6;

			fscanf_s(fp, "%lld", &din);
			matlist[i].n_beta_t_solid = din;
			matlist[i].arr_beta_t_solid = nullptr;
			matlist[i].temp_beta_t_solid = nullptr;
			matlist[i].arr_beta_t_solid = new float[matlist[i].n_beta_t_solid];
			matlist[i].temp_beta_t_solid = new float[matlist[i].n_beta_t_solid];
			if (matlist[i].temp_beta_t_solid == nullptr) {
				printf("problem memory allocation for temp_beta_t_solid\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_beta_t_solid == nullptr) {
				printf("problem memory allocation for arr_beta_t_solid\n");
				system("pause");
				exit(1);
			}
			for (integer i_4 = 0; i_4 < matlist[i].n_beta_t_solid; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				matlist[i].temp_beta_t_solid[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].arr_beta_t_solid[i_4] = static_cast <float>(fin * 1.0e-6);
			}


			// Коэффициенты Ламе не используются в данном коде.
			
			// коэффициент динамической вязкости
			fscanf_s(fp, "%f", &fin);
			matlist[i].mu = fin;
			// коэффициент линейного температурного расширения
			fscanf_s(fp, "%f", &fin);
			matlist[i].beta_t = fin;
			// признак библиотечности материала
			fscanf_s(fp, "%lld", &din);
			matlist[i].blibmat = din;
			// номер материала в библиотеке
			fscanf_s(fp, "%lld", &din);
			matlist[i].ilibident = din;
			//printf("blibmat=%d ilibident=%d\n", matlist[i].blibmat, matlist[i].ilibident);
			//system("pause");

			// для каждого КО данного блока,
			// если только он не перекрывается другими блоками
			// может быть использовано приближение 
			// Обербека-Буссинеска с соответствующей опорной температурой Tref.
			fscanf_s(fp, "%lld", &din);
			switch (din) {
			case 0: matlist[i].bBussineskApproach = false; break;
			case 1: matlist[i].bBussineskApproach = true; break;
			default: matlist[i].bBussineskApproach = false; break;
			}
			// номер закона для зависимости динамической вязкости от напряжения сдвига
			fscanf_s(fp, "%lld", &din);
			matlist[i].ilawmu = din;
			// минимальное значение динамической вязкости
			fscanf_s(fp, "%f", &fin);
			matlist[i].mumin = fin;
			// максимальное значение динамической вязкости
			fscanf_s(fp, "%f", &fin);
			matlist[i].mumax = fin;
			// параметры модельных законов для зависимости вязкости от напряжения сдвига
			fscanf_s(fp, "%f", &fin);
			matlist[i].Amu = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].Bmu = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].Cmu = fin;
			// показатель степени
			fscanf_s(fp, "%f", &fin);
			matlist[i].degreennmu = fin;

			// печать считанных значений на консоль
			//printf("%e %e %e %e %e\n", matlist[i].rho, matlist[i].cp, matlist[i].lam, matlist[i].mu, matlist[i].beta_t);
			if (0) {
				printf("HEAT_CAPACITY\n");
				printf("t_C HEAT_CAPACITY\n");
				for (integer i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
					printf("%e %e\n", matlist[i].temp_cp[i_4], matlist[i].arr_cp[i_4]);
				}
				printf("lam\n");
				printf("t_C lam\n");
				for (integer i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
					printf("%e %e\n", matlist[i].temp_lam[i_4], matlist[i].arr_lam[i_4]);
				}
				printf("%e %e %e\n", matlist[i].rho, matlist[i].mu, matlist[i].beta_t);
				printf("%lld %lld %lld\n", matlist[i].blibmat, matlist[i].ilibident, matlist[i].ilawmu); // bBoussinesq не печатается
				printf("%e %e %e %e %e %e\n", matlist[i].mumin, matlist[i].mumax, matlist[i].Amu, matlist[i].Bmu, matlist[i].Cmu, matlist[i].degreennmu);
			}
		}

		// считывание блоков
		for (i = 0; i<lb; ++i) {

			fscanf_s(fp, "%lld", &din);
			b[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

			fscanf_s(fp, "%lld", &din);
			CHECK_TYPE_GEOM(din);
			b[i].g.itypegeom = din; // 0 - Prism, 1 - Cylinder, 2 - Polygon
			fscanf_s(fp, "%lld", &din);
			if (din == 1) {
				b[i].bvisible = true;
			}
			else {
				b[i].bvisible = false;
			}


			// геометрия
			fscanf_s(fp, "%f", &fin);
			b[i].g.xS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.yS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.zS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.xE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.yE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.zE = scale*fin;
			// swap
			if (b[i].g.xS>b[i].g.xE) {
				dbuf = b[i].g.xS;
				b[i].g.xS = b[i].g.xE;
				b[i].g.xE = dbuf;
			}
			if (b[i].g.yS > b[i].g.yE) {
				dbuf = b[i].g.yS;
				b[i].g.yS = b[i].g.yE;
				b[i].g.yE = dbuf;
			}
			if (b[i].g.zS > b[i].g.zE) {
				dbuf = b[i].g.zS;
				b[i].g.zS = b[i].g.zE;
				b[i].g.zE = dbuf;
			}

			// Cylinder

			fscanf_s(fp, "%lld", &din);
			b[i].g.iPlane = din;

			fscanf_s(fp, "%f", &fin);
			b[i].g.xC = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.yC = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.zC = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.Hcyl = scale*fin;
			if (b[i].g.Hcyl < 0.0) {
				// ликвидируем отрицательную высоту цилиндра.
				switch (b[i].g.iPlane) {
				case XY_PLANE:
					b[i].g.zC += b[i].g.Hcyl;
					break;
				case XZ_PLANE:
					b[i].g.yC += b[i].g.Hcyl;
					break;
				case YZ_PLANE:
					b[i].g.xC += b[i].g.Hcyl;
					break;
				}
				b[i].g.Hcyl = fabs(b[i].g.Hcyl);
			}

			fscanf_s(fp, "%f", &fin);
			b[i].g.R_out_cyl = scale*fin;
			fscanf_s(fp, "%f", &fin);
			b[i].g.R_in_cyl = scale*fin;

			// Polygon
			//printf("Polygon\n");
			fscanf_s(fp, "%lld", &din);
			b[i].g.iPlane_obj2 = din;
#if doubleintprecision == 1
			//printf("iPlane_obj2=%lld\n", b[i].g.iPlane_obj2);
#else
			//printf("iPlane_obj2=%d\n", b[i].g.iPlane_obj2);
#endif			
			fscanf_s(fp, "%lld", &din);
			b[i].g.nsizei = din;
#if doubleintprecision == 1
			//printf("nsizei=%lld\n", b[i].g.nsizei);
#else
			//printf("nsizei=%d\n", b[i].g.nsizei);
#endif			
			b[i].g.hi = new doublereal[b[i].g.nsizei];
			b[i].g.xi = new doublereal[b[i].g.nsizei];
			b[i].g.yi = new doublereal[b[i].g.nsizei];
			b[i].g.zi = new doublereal[b[i].g.nsizei];
			for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
				fscanf_s(fp, "%f", &fin);
				b[i].g.hi[i73] = scale*fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.xi[i73] = scale*fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.yi[i73] = scale*fin;
				fscanf_s(fp, "%f", &fin);
				b[i].g.zi[i73] = scale*fin;
				if (b[i].g.hi[i73] < 0.0) {
					// ликвидируем отрицательную высоту цилиндра.
					switch (b[i].g.iPlane_obj2) {
					case XY_PLANE:
						b[i].g.zi[i73] += b[i].g.hi[i73];
						break;
					case XZ_PLANE:
						b[i].g.yi[i73] += b[i].g.hi[i73];
						break;
					case YZ_PLANE:
						b[i].g.xi[i73] += b[i].g.hi[i73];
						break;
					}
					b[i].g.hi[i73] = fabs(b[i].g.hi[i73]);
				}

#if doubleintprecision == 1
				//printf("%lld h=%e x=%e y=%e z=%e",i73, b[i].g.hi[i73], b[i].g.xi[i73], b[i].g.yi[i73], b[i].g.zi[i73]);
#else
				//printf("%d h=%e x=%e y=%e z=%e", i73, b[i].g.hi[i73], b[i].g.xi[i73], b[i].g.yi[i73], b[i].g.zi[i73]);
#endif
				//system("PAUSE");
			}

			if ((b[i].g.itypegeom== POLYGON)&&(b[i].g.nsizei > 0)) {
				// Мы заносим в окаймляющую призму границы полигона чтобы ускорить 
				// поиски в inmodel_temp быстро отсекая ячейки вне границ полигона.
				doublereal xmin53 = 1.0e30;
				doublereal ymin53 = 1.0e30;
				doublereal zmin53 = 1.0e30;
				doublereal xmax53 = -1.0e30;
				doublereal ymax53 = -1.0e30;
				doublereal zmax53 = -1.0e30;
				for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
					if (b[i].g.xi[i73] > xmax53) xmax53 = b[i].g.xi[i73];
					if (b[i].g.yi[i73] > ymax53) ymax53 = b[i].g.yi[i73];
					if (b[i].g.zi[i73] > zmax53) zmax53 = b[i].g.zi[i73];
					if (b[i].g.xi[i73] < xmin53) xmin53 = b[i].g.xi[i73];
					if (b[i].g.yi[i73] < ymin53) ymin53 = b[i].g.yi[i73];
					if (b[i].g.zi[i73] < zmin53) zmin53 = b[i].g.zi[i73];
				}
				switch (b[i].g.iPlane_obj2) {
				case XY_PLANE:
					b[i].g.xS = xmin53;
					b[i].g.xE = xmax53;
					b[i].g.yS = ymin53;
					b[i].g.yE = ymax53;
					b[i].g.zS = zmin53;
					b[i].g.zE = zmin53 + b[i].g.hi[0];
					break;
				case XZ_PLANE:
					b[i].g.xS = xmin53;
					b[i].g.xE = xmax53;
					b[i].g.zS = zmin53;
					b[i].g.zE = zmax53;
					b[i].g.yS = ymin53;
					b[i].g.yE = ymin53 + b[i].g.hi[0];
					break;
				case YZ_PLANE:
					b[i].g.yS = ymin53;
					b[i].g.yE = ymax53;
					b[i].g.zS = zmin53;
					b[i].g.zE = zmax53;
					b[i].g.xS = xmin53;
					b[i].g.xE = xmin53 + b[i].g.hi[0];
					break;
				}
			}

			if (b[i].g.itypegeom == CYLINDER) {
				// Cylinder
				
				//	printf("%e %e %e %e %e %e %e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE, b[i].g.xC, b[i].g.yC, b[i].g.zC, b[i].g.Hcyl, b[i].g.R_out_cyl, b[i].g.R_in_cyl);
				//system("PAUSE");
				
			}


			// Ввод предполагается корректным.
			// emissivity
			fscanf_s(fp, "%f", &fin);
			b[i].radiation.emissW = fin;
			fscanf_s(fp, "%f", &fin);
			b[i].radiation.emissE = fin;
			fscanf_s(fp, "%f", &fin);
			b[i].radiation.emissS = fin;
			fscanf_s(fp, "%f", &fin);
			b[i].radiation.emissN = fin;
			fscanf_s(fp, "%f", &fin);
			b[i].radiation.emissB = fin;
			fscanf_s(fp, "%f", &fin);
			b[i].radiation.emissT = fin;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			if (din == 0) {
				// Блок не является вакуумным промежутком.
				b[i].radiation.binternalRadiation = false;
			}
			else {
				// блок является вакуумным промежутком.
				b[i].radiation.binternalRadiation = true;
				if (bvacuumPrism) {
					bdouble_vacuum_PRISM = true;
				}
				bvacuumPrism = true;
				// Вычисляем View Factors.
				calculate_view_factors(b[i]);
			}
			b[i].radiation.nodelistW = nullptr;
			b[i].radiation.nodelistE = nullptr;
			b[i].radiation.nodelistS = nullptr;
			b[i].radiation.nodelistN = nullptr;
			b[i].radiation.nodelistB = nullptr;
			b[i].radiation.nodelistT = nullptr;
			b[i].radiation.nodelistWsize = 0;
			b[i].radiation.nodelistEsize = 0;
			b[i].radiation.nodelistSsize = 0;
			b[i].radiation.nodelistNsize = 0;
			b[i].radiation.nodelistBsize = 0;
			b[i].radiation.nodelistTsize = 0;

			// идентификатор материала в базе материалов
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			b[i].imatid = din;

#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			// bCylinderFixed
			if (din == 1) {
				b[i].CylinderFixed = true;
			}
			else {
				b[i].CylinderFixed = false;
			}

			// мощность тепловыделения
			//fscanf_s(fp, "%f", &fin);
			//b[i].Sc = fin;
			// 19 november 2016 температурно зависимая мощность тепловыделения.
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			b[i].n_Sc = din;
			b[i].arr_Sc = nullptr;
			b[i].temp_Sc = nullptr;
			b[i].arr_Sc = new doublereal[b[i].n_Sc];
			b[i].temp_Sc = new doublereal[b[i].n_Sc];
			if (b[i].temp_Sc == nullptr) {
				printf("problem memory allocation for temp_Sc\n");
				system("pause");
				exit(1);
			}
			if (b[i].arr_Sc == nullptr) {
				printf("problem memory allocation for arr_Sc\n");
				system("pause");
				exit(1);
			}
			// Объём полигона.
			doublereal vol_poly = Volume_polygon(b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2);

			for (integer i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				b[i].temp_Sc[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				if (fin != fin) {
					b[i].arr_Sc[i_4] = 0.0;
				}
				else {

					if (b[i].g.itypegeom == POLYGON) {

						// Для полигона передается из интерфейса просто мощность, а не удельная мощность.
						// Т.к. интерфейс не содержит функцию расчёта объёма полигона.
						// Для единообразия здесь мощность преобразуется в удельную мощность.
						if (vol_poly > 1.0e-30) {
							b[i].arr_Sc[i_4] = fin / vol_poly;
						}
						else {
							printf("error zero volume in polygon %s...\n", b[i].name);
							system("PAUSE");
							exit(1);
						}

					}
					else {
						b[i].arr_Sc[i_4] = fin;
					}
				}
			}

			// debug
			//if (fabs(b[i].Sc)>1.0e-30) {
			//printf("%e\n", b[i].Sc);
			//system("PAUSE");
			//}
			// стиль зависимости мощности тепловыделения в блоке от времени.
			// 0 - не зависит от времени, 1 - square wave зависимость, 
			// 2 - square wave 2 зависимость, 3 - hot cold режим,
			// 4 - piecewise const режим.
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			switch (din) {
			case 0: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
				break;
			case 1: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE;
				break;
			case 2: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE2;
				break;
			case 3: b[i].ipower_time_depend = POWER_TIME_DEPEND::HOT_COLD;
				break;
			case 4: b[i].ipower_time_depend = POWER_TIME_DEPEND::PIECEWISE_CONST;
				break;
			default: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
				break;
			}
			
			// тип блока
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			switch (din) {
			case 1: b[i].itype = PHYSICS_TYPE_IN_BODY::SOLID;
				break;
			case 2: b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
				break;
			case 3: b[i].itype = PHYSICS_TYPE_IN_BODY::FLUID;
				break;
			default:
				b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
				break;
			}
			

			// печать считанных значений на консоль
			//printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
			//printf("%lld %lld\n", b[i].imatid, b[i].itype);
			//printf("temperature depend power\n");
			//printf("t_C power_W\n");
			//for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
				//printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
			//}
		}

		// считывание источников тепла
		for (i = 0; i < ls; ++i) {

			fscanf_s(fp, "%lld", &din);
			s[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

			fscanf_s(fp, "%f", &fin);
			s[i].power = fin;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			if (din == 0) {
				// const
				// задаётся постоянное значение рассеиваемой в тепло мощности.
				s[i].power_multiplyer = 1.0;
				s[i].bgarber_depend = false;
			}
			else if (din == 1) {
				// рассеиваемая в тепло мощность задаётся таблично
				// в зависимости от максимальной температуры и рабочего 
				// значения смещения стока.
				s[i].bgarber_depend = true;
				s[i].power_multiplyer = s[i].power;
				// мощность будет вычислена при ambient Temperature несколькими строчками позже.
			}
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			s[i].igarber_depend = din; // уникальный номер таблицы.
			fscanf_s(fp, "%f", &fin);
			s[i].roperation_offset_drain = fin; // рабочее значение смещения стока.
												//printf("offset drain is %e\n",s[i].roperation_offset_drain);
												//system("PAUSE");
			
			if (s[i].bgarber_depend) {
				s[i].power = my_splain_interpol_power_table(gtdps[s[i].igarber_depend].intemp,
					gtdps[s[i].igarber_depend].inoffset_drain,
					gtdps[s[i].igarber_depend].rtemp,
					gtdps[s[i].igarber_depend].roffset_drain,
					gtdps[s[i].igarber_depend].rpower_table,
					operatingtemperature,
					s[i].roperation_offset_drain);
				{
					// одиночная тестовая проверка сплайновой аппроксимации.
					printf("single test validation spline approximation...\n");
					printf("calculate initial power=%e\n", s[i].power);
					printf("please, press any key to continue...");
					// system("PAUSE");
					system("pause");
				}
				s[i].power *= s[i].power_multiplyer; // домножение на корректирующий множитель.
			}



#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			s[i].iPlane = din;
			// геометрия
			fscanf_s(fp, "%f", &fin);
			s[i].g.xS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			s[i].g.yS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			s[i].g.zS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			s[i].g.xE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			s[i].g.yE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			s[i].g.zE = scale*fin;



			// swap
			if (s[i].g.xS > s[i].g.xE) {
				dbuf = s[i].g.xS;
				s[i].g.xS = s[i].g.xE;
				s[i].g.xE = dbuf;
			}
			if (s[i].g.yS > s[i].g.yE) {
				dbuf = s[i].g.yS;
				s[i].g.yS = s[i].g.yE;
				s[i].g.yE = dbuf;
			}
			if (s[i].g.zS > s[i].g.zE) {
				dbuf = s[i].g.zS;
				s[i].g.zS = s[i].g.zE;
				s[i].g.zE = dbuf;
			}
			switch (s[i].iPlane) {
			case XY_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
			case XZ_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
			case YZ_PLANE: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
			default: break;
			}
			//printf("source %e %lld %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
		}

		// считывание твёрдых стенок

		// 19.01.2021
		// Стенка заглушка с условием прилипания по скорости(нулевой вектор скорости), 
		// нулевым тепловым потоком по температуре, 
		// свободной free границе в механике без приложенной силы.
		w[lw].iunion_id = 0; // cabinet
		w[lw].ifamily = WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY;
		w[lw].Tamb = 0.0;
		w[lw].emissivity = 0.0;
		w[lw].film_coefficient = 0.0;
		w[lw].ViewFactor = 1.0;
		w[lw].hf = 0.0;
		w[lw].bsymmetry = false;
		w[lw].bpressure = false;
		w[lw].bopening = false;
		w[lw].Vx = 0.0;
		w[lw].Vy = 0.0;
		w[lw].Vz = 0.0;
		w[lw].P = 0.0;
		w[lw].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
		w[lw].xForce = 0.0;
		w[lw].yForce = 0.0;
		w[lw].zForce = 0.0;
		w[lw].iPlane = XY_PLANE;
		w[lw].g.xS = 0.0;
		w[lw].g.yS = 0.0;
		w[lw].g.zS = 0.0;
		w[lw].g.xE = 0.0;
		w[lw].g.yE = 0.0;
		w[lw].g.zE = 0.0;

		for (i = 0; i < lw; ++i) {

			fscanf_s(fp, "%lld", &din);
			w[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			switch (din) {
			case 1 : w[i].ifamily= WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY; break;
				case 2: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
				case 3: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY; break;
				case 4: w[i].ifamily= WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY; break;
				default: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
			}
			//w[i].ifamily = din;
			switch (din) {
			case 1:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = fin;
				fscanf_s(fp, "%f", &fin); // Stefan Bolcman
										  // termostability wall
				w[i].emissivity = 0.0;
				w[i].film_coefficient = 0.0;
				fscanf_s(fp, "%f", &fin); // ViewFactor
				w[i].ViewFactor = 1.0;
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // первого рода
			case 2:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = 0.0;
				fscanf_s(fp, "%f", &fin);  // Stefan Bolcman
										   // adiabatic wall
				w[i].emissivity = 0.0;
				w[i].film_coefficient = 0.0;
				fscanf_s(fp, "%f", &fin); // ViewFactor
				w[i].ViewFactor = 1.0;
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // однородное условие Неймана
			case 3:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = fin;
				fscanf_s(fp, "%f", &fin); // Stefan Bolcman
				// Newton-Richman condition, film coefficient.
				w[i].emissivity = 0.0;
				w[i].film_coefficient = fin;
				fscanf_s(fp, "%f", &fin); // ViewFactor
				w[i].ViewFactor = 1.0;
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // Ньютон-Рихман.
			case 4:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = fin;
				fscanf_s(fp, "%f", &fin); // Stefan Bolcman
										  // Stefan - Bolcman condition
				w[i].emissivity = fin;
				w[i].film_coefficient = 0.0;
				fscanf_s(fp, "%f", &fin); // ViewFactor
				w[i].ViewFactor = fmax(0.0,fmin(fin,1.0));
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // Стефан-Больцман.
			default: 
				printf("error: wall unlnown boundary condition type.\n");
				system("PAUSE");
				break;
			}
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			if (din == 1) w[i].bsymmetry = true; else w[i].bsymmetry = false;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			if (din == 1) w[i].bpressure = true; else w[i].bpressure = false;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			if (din == 1) w[i].bopening = true; else w[i].bopening = false;
			fscanf_s(fp, "%f", &fin);
			w[i].Vx = fin;
			fscanf_s(fp, "%f", &fin);
			w[i].Vy = fin;
			fscanf_s(fp, "%f", &fin);
			w[i].Vz = fin;
			fscanf_s(fp, "%f", &fin);
			w[i].P = fin;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif		
			if (din >= 0 && din < 11) {
				// 0- FREE

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
				switch (din) {
				case 0: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
					break;
				case 1: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FIXIT;
					break;
				case 2: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FIXIT;
					break;
				case 3: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FIXIT;
					break;
				case 4: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XY_FIXIT;
					break;
				case 5: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XZ_FIXIT;
					break;
				case 6: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::YZ_FIXIT;
					break;
				case 7: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::ALL_FIXIT;
					break;
				case 8: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FORCE;
					break;
				case 9: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FORCE;
					break;
				case 10: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FORCE;
					break;
				default:
					w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // Free all
					break;
				}
			}
			else {
				printf("error: unknown ithermal_Stress_boundary_condition\n");
				printf("ithermal_Stress_boundary_condition=%lld\n", din);
				system("PAUSE");
				w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // Free all
			}
			fscanf_s(fp, "%f", &fin);
			w[i].xForce = fin;
			fscanf_s(fp, "%f", &fin);
			w[i].yForce = fin;
			fscanf_s(fp, "%f", &fin);
			w[i].zForce = fin;
			//printf("Force Fx=%e Fy=%e Fz=%e\n", w[i].xForce, w[i].yForce, w[i].zForce);
			//system("PAUSE");
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			w[i].iPlane = din;
			// геометрия
			fscanf_s(fp, "%f", &fin);
			w[i].g.xS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			w[i].g.yS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			w[i].g.zS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			w[i].g.xE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			w[i].g.yE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			w[i].g.zE = scale*fin;
			// swap
			if (w[i].g.xS > w[i].g.xE) {
				dbuf = w[i].g.xS;
				w[i].g.xS = w[i].g.xE;
				w[i].g.xE = dbuf;
			}
			if (w[i].g.yS > w[i].g.yE) {
				dbuf = w[i].g.yS;
				w[i].g.yS = w[i].g.yE;
				w[i].g.yE = dbuf;
			}
			if (w[i].g.zS > w[i].g.zE) {
				dbuf = w[i].g.zS;
				w[i].g.zS = w[i].g.zE;
				w[i].g.zE = dbuf;
			}
			//printf("wall %lld %e %e %lld %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
		}


		// АСЕМБЛЕСЫ.
		// Новое считывание работает только с новым форматом файла premeshin.txt 
		// с текстовыми метками, данный код здесь устарел. 31.10.2021.
#if doubleintprecision == 1
		fscanf_s(fp, "%lld", &din);
#else
		fscanf_s(fp, "%d", &din);
#endif
		lu = (int)(din);
		if (lu == 0) {
			my_union = nullptr;
		}
		else {
			my_union = new UNION[lu];
			// инициализация.
			for (i = 0; i < lu; ++i) {
				my_union[i].f = nullptr;
				my_union[i].xpos = nullptr;
				my_union[i].ypos = nullptr;
				my_union[i].zpos = nullptr;
				my_union[i].xposadd = nullptr;
				my_union[i].yposadd = nullptr;
				my_union[i].zposadd = nullptr;
				my_union[i].iunion_parent = -1; // Cabinet.
				my_union[i].iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; // 2 - CoarseMeshGen
				my_union[i].inxadd = -1;
				my_union[i].inyadd = -1;
				my_union[i].inzadd = -1;
				my_union[i].flow_interior = 0;
				my_union[i].active = false;
			}
		}
		for (i = 0; i < lu; ++i) {

#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].iunion_parent = din;

			fscanf_s(fp, "%f", &fin);
			my_union[i].xS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			my_union[i].xE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			my_union[i].yS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			my_union[i].yE = scale*fin;
			fscanf_s(fp, "%f", &fin);
			my_union[i].zS = scale*fin;
			fscanf_s(fp, "%f", &fin);
			my_union[i].zE = scale*fin;

#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].id = din; // Уникальный идентификатор АССЕМБЛЕСА.
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].inx = (int)(din);
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].iny = (int)(din);
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].inz = (int)(din);

		}

		// считывание информации о наборе решаемых уравнений
#if doubleintprecision == 1
		fscanf_s(fp, "%lld", &din);
#else
		fscanf_s(fp, "%d", &din);
#endif
		eqin.itemper = din;
#if doubleintprecision == 1
		fscanf_s(fp, "%lld", &din);
#else
		fscanf_s(fp, "%d", &din);
#endif
		eqin.imaxflD = din;
		//printf("itemper=%lld eqin.imaxflD=%lld\n", eqin.itemper, eqin.imaxflD);
		//system("PAUSE");
		if (eqin.imaxflD == 0) {
			eqin.fluidinfo = nullptr;
		}
		else
		{
			// выделение оперативной памяти
			if (eqin.fluidinfo != nullptr) {
				delete eqin.fluidinfo;
				eqin.fluidinfo = nullptr;
			}
			eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
			for (i = 0; i < eqin.imaxflD; ++i) {
				// Считывание координат опорной точки
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].xc = scale * fin;
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].yc = scale * fin;
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].zc = scale * fin;
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].iflow = din;
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				switch (din) {
				case 0: eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
					break;
				case 1:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::TURBULENT;
					break;
				default:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
					break;
				}
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				// Режим течения: ламинарный или конкретная модель турбулентности.
				switch (din) {
				case 0: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
					break;
				case 1: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::SMAGORINSKY;
					break;
				case 2: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RNG_LES;
					break;
				case 3: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_SPALART_ALLMARES;
					break;
				case 4: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_MENTER_SST;
					break;
				case 5: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_STANDART_K_EPS;
					break;
				case 6: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_LANGTRY_MENTOR_SST;
					break;
				default: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
					break;
				}
				
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].Cs = fin; // постоянная Смагоринского.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				// учёт динамического определения квадрата постоянной Смагоринского.
				// включает Dynamic Subgrid Scale Model Германо 1991 года.
				if (din == 1) {
					eqin.fluidinfo[i].bDynamic_Stress = true;
				}
				else {
					eqin.fluidinfo[i].bDynamic_Stress = false;
				}
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				// включает ограничение сверху и снизу на возможные значения постоянной Смагоринского.
				if (din == 1) {
					eqin.fluidinfo[i].bLimiters_Cs = true;
				}
				else {
					eqin.fluidinfo[i].bLimiters_Cs = false;
				}
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].rminCs = fin; // минимальное возможное значение постоянной Смагоринского.
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].rmaxCs = fin; // максимальное возможное значение постоянной Смагоринского.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].itypeFiltrGermano = din; // тип фильтра в модели Германо 1991 года.
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].roughness = 1.0e-6 * fin; // шероховатость стенки в м.
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].rRi_mult = fin; // множитель корректирующий турбулентное число Ричардсона.
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].rSelectiveAngle = fin; // пороговое значение угла в модели Selective Smagorinsky.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].ipowerroughness = din; // показатель степени в модели учёта шероховатости.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = din; // тип фильтра в модели Selective Smagorinsky.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].bfdelta = din; // учёт неравномерности сетки.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].bSmagorinskyLilly = din; // модель Смагоринского-Лиллу.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].bsurface_roughness = din; // учёт шероховатости стенки.
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				// учёт течений с кривизной линий тока.
				if (din == 1) {
					eqin.fluidinfo[i].bSwirlAmendment = true;
				}
				else {
					eqin.fluidinfo[i].bSwirlAmendment = false;
				}
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				// учёт избирательности в модели Смагоринского
				if (din == 1) {
					eqin.fluidinfo[i].bSelectiveSmagorinsky = true;
				}
				else {
					eqin.fluidinfo[i].bSelectiveSmagorinsky = false;
				}

				// Параметры преобразователя картинок для отчетов.
				// 5.01.2018
				fscanf_s(fp, "%f", &fin);
				pfpir.fminimum = static_cast<doublereal>(scale * fin);
				fscanf_s(fp, "%f", &fin);
				pfpir.fmaximum = static_cast<doublereal>(scale * fin);
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				switch (din) {
					case 0: pfpir.idir = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
					case 1: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					case 2: pfpir.idir = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
					default: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
				}
				//pfpir.idir = din;

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				AMG1R6_LABEL = din;

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				if (din > 9999) {
					b_iluk_amg1r5_LABEL_D = true;
					nrd_LABEL = din-10000;
				}
				else {
					nrd_LABEL = din;
			    }

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				if (din >9999) {
					b_iluk_amg1r5_LABEL_U = true;
					nru_LABEL = din-10000;
				}
				else {
					nru_LABEL = din;
			    }

#if doubleintprecision == 1
				fscanf_s(fp, "%f", &fin);
#else
				fscanf_s(fp, "%f", &fin);
#endif
				ecg2_LABEL = fin;

#if doubleintprecision == 1
				fscanf_s(fp, "%f", &fin);
#else
				fscanf_s(fp, "%f", &fin);
#endif
				ewt2_LABEL = fin;

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				number_processors_global_var = (int)(din);

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				idevice_Tesla = (int)(din);

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				if (din>=0) {
					number_iteration_SIMPLE_algorithm = (unsigned int)(din);
				}
				else {
					number_iteration_SIMPLE_algorithm =0;
					std::cout<<"WARNING : number_iteration_SIMPLE_algorithm =0;"<<std::endl;
				}

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				switch (din) {
				case 0: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::NONE_only_amg1r5;
					break;
				case 1: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
					break;
				case 2: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::FGMRes_plus_amg1r5;
					break;
				case 3: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5;
					break;
				default:
					stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
					break;
				}


#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				BonLevelDrobim = din;

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				isizearrAdaptRegion = din;


				if (isizearrAdaptRegion > 0) {

					AdaptRegion = new ADAPTREGION[isizearrAdaptRegion];

					AdaptRegion[0].ilevelMeshAdaptRegion = -1;
					AdaptRegion[1].ilevelMeshAdaptRegion = -1;
					AdaptRegion[2].ilevelMeshAdaptRegion = -1;

					for (int i_63 = 0; i_63 < isizearrAdaptRegion; ++i_63) {

						// Адаптация сетки региона номер 1.
						if (i_63 == 0) {

#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[0].xS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[0].xE = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[0].yS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[0].yE = scale * fin;

#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[0].zS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[0].zE = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%lld", &din);
#else
							fscanf_s(fp, "%d", &din);
#endif
							AdaptRegion[0].ilevelMeshAdaptRegion = din;

						}


						if (i_63 == 1) {


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[1].xS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[1].xE = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[1].yS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[1].yE = scale * fin;

#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[1].zS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[1].zE = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%lld", &din);
#else
							fscanf_s(fp, "%d", &din);
#endif
							AdaptRegion[1].ilevelMeshAdaptRegion = din;

						}


						if (i_63 == 2) {

#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[2].xS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[2].xE = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[2].yS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[2].yE = scale * fin;

#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[2].zS = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%f", &fin);
#else
							fscanf_s(fp, "%f", &fin);
#endif
							AdaptRegion[2].zE = scale * fin;


#if doubleintprecision == 1
							fscanf_s(fp, "%lld", &din);
#else
							fscanf_s(fp, "%d", &din);
#endif
							AdaptRegion[2].ilevelMeshAdaptRegion = din;


						}

					}

				}

#if doubleintprecision == 1
				fscanf_s(fp, "%f", &fin);
#else
				fscanf_s(fp, "%f", &fin);
#endif

				// Найдено успешно.
				// Свободный параметр передаваемый из интерфейса и используемый для отладки.
				free_debug_parametr1 = static_cast<doublereal>(fin);

			}
		}


		integer ilb_p = 0;// Количество блоков внутри которых задана тепловая мощность.
		doublereal dpower = 0.0; // Суммарная тепловая мощность в блоках.
		integer ipoly = 0, icyl = 0, iprism = 0, icad_stl=0;
		integer ihol = 0, isol = 0, iflui = 0;
		for (integer i_1 = 0; i_1 < lb; ++i_1) {
			if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
				ihol++;
			}
			if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				iflui++;
			}
			if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				isol++;
			}
			if (b[i_1].g.itypegeom == CAD_STL)
			{
				icad_stl++;
			}
			if (b[i_1].g.itypegeom == PRISM) {
				// 0 - PRISM object

				// 13.08.2019 Автомат настройки допусков сетки. 
				// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
				// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
				// Она не может быть больше стороны кабинета деленной на 15.
				if (0.1*fabs(b[i_1].g.xE - b[i_1].g.xS) < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*fabs(b[i_1].g.xE - b[i_1].g.xS);
				if (0.1*fabs(b[i_1].g.yE - b[i_1].g.yS) < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*fabs(b[i_1].g.yE - b[i_1].g.yS);
				if (0.1*fabs(b[i_1].g.zE - b[i_1].g.zS) < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*fabs(b[i_1].g.zE - b[i_1].g.zS);

				if (lb == 1) {// Течение в каверне или тест Дэвиса.
					shorter_length_for_simplificationX_BASIC = 1.0e-10;
					shorter_length_for_simplificationY_BASIC = 1.0e-10;
					shorter_length_for_simplificationZ_BASIC = 1.0e-10;
				}
				else {
					if (shorter_length_for_simplificationX_BASIC > 0.067*(fabs(b[0].g.xE - b[0].g.xS))) shorter_length_for_simplificationX_BASIC = 0.067*(fabs(b[0].g.xE - b[0].g.xS));
					if (shorter_length_for_simplificationY_BASIC > 0.067*(fabs(b[0].g.yE - b[0].g.yS))) shorter_length_for_simplificationY_BASIC = 0.067*(fabs(b[0].g.yE - b[0].g.yS));
					if (shorter_length_for_simplificationZ_BASIC > 0.067*(fabs(b[0].g.zE - b[0].g.zS))) shorter_length_for_simplificationZ_BASIC = 0.067*(fabs(b[0].g.zE - b[0].g.zS));
				}

				iprism++;
				if (b[i_1].n_Sc > 0) {
					doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
					doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS)*fabs(b[i_1].g.yE - b[i_1].g.yS)*fabs(b[i_1].g.zE - b[i_1].g.zS);
					if (vol < 1.0e-36) {
						printf("ERROR: zero volume in PRISM block number %lld\n",i_1);
						system("PAUSE");
						exit(1);
					}
					//if (fabs(b[i_1].arr_Sc[0]) > 0.0) {
					if (pdiss > 0.0) {
						ilb_p++;
						//dpower += b[i_1].arr_Sc[0];
						dpower += pdiss*vol;
					}
				}
		    }
			if (b[i_1].g.itypegeom == CYLINDER) {
				// Cylinder

				// На тот случай если геометрия состоит только из цилиндров.
				// 13.08.2019 Автомат настройки допусков сетки. 
				// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
				// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
				switch (b[i_1].g.iPlane) {
				case XY_PLANE:
					if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.Hcyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-36) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_in_cyl;
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
					}
					break;
				case XZ_PLANE:
					if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.Hcyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-36) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_in_cyl;
					}
					break;
				case YZ_PLANE:
					if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.Hcyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-36) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_in_cyl;
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_in_cyl;
					}
					break;
				}
				

				icyl++;
				if (b[i_1].n_Sc > 0) {
					doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
					doublereal vol = 0.0;
					//const doublereal M_PI = 3.1415926;// math.h
					vol = b[i_1].g.Hcyl*M_PI*(b[i_1].g.R_out_cyl*b[i_1].g.R_out_cyl - b[i_1].g.R_in_cyl*b[i_1].g.R_in_cyl);				
					if (vol < 1.0e-36) {
						printf("ERROR: zero volume in CYLINDER block number %lld\n", i_1);
						system("PAUSE");
						exit(1);
					}
					if (pdiss > 0.0) {
						ilb_p++;

						dpower += pdiss*vol;
						//printf("pdiss=%e vol=%e dpower=%e\n",pdiss, vol, dpower);
						//system("pause");
						//printf("ERROR: non zero power in cylinder object.\n");
						//system("PAUSE");
						//exit(1);
					}
				}
			}
			if (b[i_1].g.itypegeom == POLYGON) {
				// Polygon
				ipoly++;
				if (b[i_1].n_Sc > 0) {
					doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);

					// Объём полигона.
					doublereal vol_poly = Volume_polygon(b[i_1].g.nsizei, b[i_1].g.xi, b[i_1].g.yi, b[i_1].g.zi, b[i_1].g.hi, b[i_1].g.iPlane_obj2);
					if (vol_poly < 1.0e-36) {
						printf("ERROR: zero volume in POLYGON block number %lld\n", i_1);
						system("PAUSE");
						exit(1);
					}

					if (pdiss > 0.0) {
						ilb_p++;
						
						dpower += pdiss*vol_poly;
						
						//printf("ERROR: non zero power in polygon object.\n");
						//system("PAUSE");
						//exit(1);
					}
				}
			}
		}

		doublereal dsoupow = 0.0; // интегральная тепловая мощность плоских источников тепла.
		for (integer i_1 = 0; i_1 < ls; ++i_1) {
			dsoupow += s[i_1].power;
		}


		printf("Apriory quick model statistics:\n");
		printf("number of thermal power blocks lb_p=%lld\n", ilb_p);
		printf("Blocks integral power =%e W\n", dpower);
		printf("number of sources ls=%d\n",ls);
		printf("Sources integral power = %e W\n", dsoupow);
		printf("Full total power = %e W\n", dpower + dsoupow);
		// Запоминаем полное тепловыделение в модели.
		d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = dpower + dsoupow;
		printf("number of blocks lb=%d\n",lb);
		printf("PRISMS = %lld, CYLINDERS = %lld, POLYGONS = %lld, CAD_STL = %lld\n", iprism, icyl, ipoly, icad_stl);
		printf("SOLID: %lld\n", isol);
		printf("HOLLOW: %lld\n", ihol);
		printf("FLUID: %lld\n", iflui);
		printf("number of walls lw=%d\n",lw);
		printf("number of units lu=%d\n", lu);

		if (fp != NULL) {
			fclose(fp); // закрытие файла
			fp = NULL;
		}
	}
}

#endif
#endif
	printf("OK. \n");
} // premeshin_old

// Код 2009 года. Из моей первой программы. Функции необхдимые для работы парсера.


// проверяет совпадают ли две строки:
// имя name и образец pattern
bool isname(char * name, char * pattern) {
	bool bret = false;
	if (name != nullptr) {
		int i; // счётчик
		int il1, il2; // длины строк

		i = 0;
		while ((name[i] != '\0')&&((name[i] != '\n'))) i++;
		il1 = i;
		i = 0;
		while ((pattern[i] != '\0')&&(pattern[i] != '\n')) i++;
		il2 = i;

		if (il1 == il2) {
			//printf("%s \n", pattern);
			//system("pause");
			bret = true;
			for (i = 0; i < il1; ++i) {
				if (name[i] != pattern[i]) {
					bret = false;
					return bret;
				}
			}
		}
		//printf("%d \n", (int)(bret));
	}
	return bret;
} // isname 2009





// Массив строк для ускорения обработки файлов.
const integer ilimit_StringList = 10000000;
integer icurrentSize_StringList = 0;
integer icurrent_position__StringList = 0;
char ** StringList = new char*[ilimit_StringList];

// Загружает данные из входного файла в массив StringList.
// premeshin.txt
// реализовано 25 08 2019.
void loadFromFile()
{
	FILE *fp = NULL;
	

#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("premeshin.txt", "r");
	if (fp != nullptr) {
		err = 0;
	}
	else {
		err = 1; // ошибка открытия.
	}
#else
	errno_t err;
	err = fopen_s(&fp, "premeshin.txt", "r");
#endif

	if (fp != NULL) {
		// создание файла для чтения.
		if (err != 0) {
			printf("Open File Error\n");
			exit(0);
			// return bfound;
		}
		else
		{
			char c;
			c = ' ';
			bool bbeginstring = true; // строка ещё не кончилась
			bool b1 = true;
			int k = 0;
			char * buf = (char *)malloc(1024); // строка буфер
			if (buf != nullptr) {
				c = fgetc(fp); // читает один символ
				buf[k++] = c;
				
				
				while (!feof(fp)) {
					b1 = true;
					c = fgetc(fp); // читает один символ
					if (c == '#') bbeginstring = false;
					if (c == '\n') {
						bbeginstring = true;
						buf[k++] = '\n';
						buf[k] = '\0';

						StringList[icurrentSize_StringList] = new char[k + 1];
						for (integer i_11 = 0; i_11 <= k; i_11++) {
							StringList[icurrentSize_StringList][i_11] = buf[i_11];
						}
						icurrentSize_StringList++;

						if (icurrentSize_StringList > ilimit_StringList) {
							printf("Error!!! buffer StringList is overflow...\n");
							printf("icurrentSize_StringList > ilimit_StringList==4000K\n");
							system("pause");
							exit(1);
						}

						//printf("%s \n",buf);
						//system("pause");
						b1 = false;
						k = 0;
					}
					if ((b1) && (bbeginstring  )) {
						buf[k++] = c;
						// printf("%c \n",c);
					}

				}// while

				buf[k++] = '\n';
				buf[k] = '\0';

				StringList[icurrentSize_StringList] = new char[k + 1];
				for (integer i_11 = 0; i_11 <= k; i_11++) {
					StringList[icurrentSize_StringList][i_11] = buf[i_11];
				}
				icurrentSize_StringList++;

				//printf("%s \n",buf);
				//system("pause");	

				free(buf);
			}
			if (fp != NULL) {
				fclose(fp);
				fp = NULL;
			}

		}
	}
} // loadFromFile

// Освобождает память из под массива StringList.
void freeStringList() {
	
	for (integer i_11 = 0; i_11 < icurrentSize_StringList; i_11++) {
		delete[] StringList[i_11];
	}
	icurrentSize_StringList=0;
	icurrent_position__StringList = 0;

} // freeStringList


// анализирует строку символов buf длиной ilen.
bool sanalizestring(char* name0, char* buf, int ilen, char* strret)
{
	// состав анализируемой строки:
	// ключевое слово, символ равенства, число без знака.
	bool bfound = false;

	bool bstring = false;
	int inum = -1;
	int i; // счётчик цикла for
	for (i = 0; i < ilen; ++i) {
		if (buf[i] == '=') {
			bstring = true;
			inum = i;
		}
	}
	if (bstring) {
		//
		//printf("%d\n",inum);

		// если строка нам подходит по структуре
		char* name = (char*)malloc(1024); // строка имя переменной
		if (name != nullptr) {
			for (i = 0; i < inum; ++i) name[i] = buf[i];
			name[i] = '\n'; // конец имени
			i++;
			name[i] = '\0'; // конец строки
			//printf("%s \n", name); // debug
		}
		char* value = (char*)malloc(1024); // строка значение переменной с имененм name
		i = inum + 1;
		int j = 0;
		if (value != nullptr) {
			while (buf[i] != '\n') {
				value[j] = buf[i];
				j++;
				i++;
			}
			value[j] = '\n';
			j++;
			value[j] = '\0'; // символ конца строки
			//printf("value=%s \n",value); // контроль считанного значения
			//printf("name=%s \n", name);
			//system("pause");

			if (isname(name, name0)) {
				for (integer j45 = 0; j45 < j; j45++) {
					strret[j45] = value[j45];
				}
				strret[j-1]= '\0'; // символ конца строки

				bfound = true;
			}
		}

		if (name != nullptr) {
			free(name);
		}
		if (value != nullptr) {
			free(value);
		}
		return bfound;

	}

	return bfound;
} // sanalizestring по мотивам 2009



// анализирует строку символов buf длиной ilen.
bool ianalizestring(char* name0, char* buf, int ilen, int& iret)
{
	// состав анализируемой строки:
	// ключевое слово, символ равенства, число без знака.
	bool bfound = false;

	bool bstring = false;
	int inum = -1;
	int i; // счётчик цикла for
	for (i = 0; i < ilen; ++i) {
		if (buf[i] == '=') {
			bstring = true;
			inum = i;
		}
	}
	if (bstring) {
		//
		//printf("%d\n",inum);

		// если строка нам подходит по структуре
		char* name = (char*)malloc(1024); // строка имя переменной
		if (name != nullptr) {
			for (i = 0; i < inum; ++i) name[i] = buf[i];
			name[i] = '\n'; // конец имени
			i++;
			name[i] = '\0'; // конец строки
			//printf("%s \n", name); // debug
		}
		char* value = (char*)malloc(1024); // строка значение переменной с имененм name
		i = inum + 1;
		int j = 0;
		if (value != nullptr) {
			while (buf[i] != '\n') {
				value[j] = buf[i];
				j++;
				i++;
			}
			value[j] = '\n';
			j++;
			value[j] = '\0'; // символ конца строки
			//printf("value=%s \n",value); // контроль считанного значения
			//printf("name=%s \n", name);
			//system("pause");

			if (isname(name, name0)) {
				iret = atoi(value);
				bfound = true;
			}
		}

		if (name != nullptr) {
			free(name);
		}
		if (value != nullptr) {
			free(value);
		}
		return bfound;

	}

	return bfound;
} // ianalizestring по мотивам 2009

// анализирует строку символов buf длиной ilen.
bool fanalizestring(char* name0, char* buf, int ilen, double& fret)
{
	// состав анализируемой строки:
	// ключевое слово, символ равенства, число без знака.
	bool bfound = false;

	bool bstring = false;
	int inum = -1;
	int i; // счётчик цикла for
	for (i = 0; i < ilen; ++i) {
		if (buf[i] == '=') {
			bstring = true;
			inum = i;
		}
	}
	if (bstring) {
		//
		//printf("%d\n",inum);

		// если строка нам подходит по структуре
		char* name = (char*)malloc(1024); // строка имя переменной
		if (name != nullptr) {
			for (i = 0; i < inum; ++i) name[i] = buf[i];
			name[i] = '\n'; // конец имени
			i++;
			name[i] = '\0'; // конец строки
			//printf("%s \n", name); // debug
			char* value = (char*)malloc(1024); // строка значение переменной с именем name
			if (value != nullptr) {
				i = inum + 1;
				int j = 0;
				while (buf[i] != '\n') {
					value[j] = buf[i];
					j++;
					i++;
				}
				value[j] = '\n';
				j++;
				value[j] = '\0'; // символ конца строки
				//printf("value=%s \n",value); // контроль считанного значения
				//printf("name=%s \n", name);
				//system("pause");

				if (isname(name, name0)) {
					struct lconv* loc = localeconv();
					char* separator = loc->decimal_point; //этот самый separator и будет искать atof в строке
					if (separator[0] == ',') {
						for (int j29 = 0; j29 < j - 1; j29++) {
							if (value[j29] == '.') {
								value[j29] = ',';
							}
						}
					}
					if (separator[0] == '.') {
						for (int j29 = 0; j29 < j - 1; j29++) {
							if (value[j29] == ',') {
								value[j29] = '.';
							}
						}
					}

					fret = atof(value);
					bfound = true;
				}

				free(value);
			}
			free(name);
		}
		return bfound;

	}

	return bfound;
} // fanalizestring по мотивам 2009

// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool imakesource(const char* name0_const, int& iret)
{
	bool bfound = false;
	char* name0 = new char[strlen(name0_const)+1];
#ifdef MINGW_COMPILLER
	strcpy(name0, name0_const);
#else
	strcpy_s(name0, (strlen(name0_const) + 1)*sizeof(char), name0_const);
#endif

	char* buf = (char*)malloc(1024); // строка буфер
	if (buf != nullptr) {

		bool bfound_loc = false;
		int iret_loc = -1;
		for (integer i_1 = icurrent_position__StringList; i_1 < icurrentSize_StringList; ++i_1) {
			if (!bfound) {
				// Если не найдена.

				bool bbeginstring = true; // строка ещё не кончилась
				int k = 0;
				char c;
				c = ' ';
				while (StringList[i_1][k] != '\0') {
					c = StringList[i_1][k];
					if (c == '#') {
						bbeginstring = false;
						buf[k++] = '\n';
						buf[k] = '\0';
						//bfound_loc = false;
						bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
						if (bfound_loc) {
							bfound = true;
							iret = iret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						break;
						//printf("%s \n",buf);
						//system("pause");	
					}
					if (c == '\n') {
						bbeginstring = false;
						buf[k++] = '\n';
						buf[k] = '\0';

						//bfound_loc = false;
						bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
						if (bfound_loc) {
							bfound = true;
							iret = iret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						//printf("%s \n",buf);
						//system("pause");
						break;
					}
					if (bbeginstring) {
						buf[k++] = c;
						// printf("%c \n",c);
					}
				} // while
			}
			else {
				break;
			}
		}


		if (!bfound) {

			// Сканируем всё сначала.
			for (integer i_1 = 0; i_1 < icurrent_position__StringList; ++i_1) {
				if (!bfound) {
					// Если не найдена.

					bool bbeginstring = true; // строка ещё не кончилась
					int k = 0;
					char c;
					c = ' ';
					while (StringList[i_1][k] != '\0') {
						c = StringList[i_1][k];
						if (c == '#') {
							bbeginstring = false;
							buf[k++] = '\n';
							buf[k] = '\0';
							//bfound_loc = false;
							bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
							if (bfound_loc) {
								bfound = true;
								iret = iret_loc;
								icurrent_position__StringList = i_1 + 1;
							}
							break;
							//printf("%s \n",buf);
							//system("pause");	
						}
						if (c == '\n') {
							bbeginstring = false;
							buf[k++] = '\n';
							buf[k] = '\0';

							//bfound_loc = false;
							bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
							if (bfound_loc) {
								bfound = true;
								iret = iret_loc;
								icurrent_position__StringList = i_1 + 1;
							}
							//printf("%s \n",buf);
							//system("pause");
							break;
						}
						if (bbeginstring) {
							buf[k++] = c;
							// printf("%c \n",c);
						}
					} // while
				}
				else {
					break;
				}
			}
		}


		free(buf);
	}

	delete[] name0;
	return bfound;


} // imakesource по мотивам 2010

// считывает данные из входного файла
// premeshin.txt
// реализовано 23 апреля 2010.
// revised 24 мая 2020.
bool smakesource(const char* name0_const, char* name)
{
	bool bfound = false;
	char* name0 = new char[strlen(name0_const) + 1];
#ifdef MINGW_COMPILLER
	strcpy(name0, name0_const);
#else
	strcpy_s(name0, (strlen(name0_const) + 1) * sizeof(char), name0_const);
#endif

	char* buf = (char*)malloc(1024); // строка буфер
	if (buf != nullptr) {

		bool bfound_loc = false;
		
		for (integer i_1 = icurrent_position__StringList; i_1 < icurrentSize_StringList; ++i_1) {
			if (!bfound) {
				// Если не найдена.

				bool bbeginstring = true; // строка ещё не кончилась
				int k = 0;
				char c;
				c = ' ';
				while (StringList[i_1][k] != '\0') {
					c = StringList[i_1][k];
					if (c == '#') {
						bbeginstring = false;
						buf[k++] = '\n';
						buf[k] = '\0';

						bfound_loc = sanalizestring(name0, buf, k - 1, name);
						if (bfound_loc) {
							bfound = true;
							icurrent_position__StringList = i_1 + 1;
						}
												
						//printf("%s # \n",buf);
						//system("pause");	
						break;
					}
					if (c == '\n') {
						bbeginstring = false;
						buf[k++] = '\n';
						buf[k] = '\0';

						bfound_loc = sanalizestring(name0, buf, k - 1, name);
						if (bfound_loc) {
							bfound = true;
							icurrent_position__StringList = i_1 + 1;
						}

						//printf("%s k=%d n \n",buf,k);
						//system("pause");
						break;
					}
					if (bbeginstring) {
						buf[k++] = c;
						// printf("%c \n",c);
					}
				} // while
			}
			else {
				break;
			}
		}


		if (!bfound) {

			// Сканируем всё сначала.
			for (integer i_1 = 0; i_1 < icurrent_position__StringList; ++i_1) {
				if (!bfound) {
					// Если не найдена.

					bool bbeginstring = true; // строка ещё не кончилась
					int k = 0;
					char c;
					c = ' ';
					while (StringList[i_1][k] != '\0') {
						c = StringList[i_1][k];
						if (c == '#') {
							bbeginstring = false;
							buf[k++] = '\n';
							buf[k] = '\0';

							bfound_loc = sanalizestring(name0, buf, k - 1, name);
							if (bfound_loc) {
								bfound = true;
								icurrent_position__StringList = i_1 + 1;
							}

							//printf("%s \n",buf);
							//system("pause");	
							break;
						}
						if (c == '\n') {
							bbeginstring = false;
							buf[k++] = '\n';
							buf[k] = '\0';

							bfound_loc = sanalizestring(name0, buf, k - 1, name);
							if (bfound_loc) {
								bfound = true;
								icurrent_position__StringList = i_1 + 1;
							}
							
							//printf("%s \n",buf);
							//system("pause");
							break;
						}
						if (bbeginstring) {
							buf[k++] = c;
							// printf("%c \n",c);
						}
					} // while
				}
				else {
					break;
				}
			}
		}


		free(buf);
	}

	delete[] name0;
	return bfound;


} // smakesource по мотивам 2010

// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool fmakesource(const char* name0_const, double& fret)
{
	bool bfound = false;
	char* name0 = new char[strlen(name0_const) + 1];
#ifdef MINGW_COMPILLER
	strcpy(name0, name0_const);
#else
	strcpy_s(name0, (strlen(name0_const) + 1)*sizeof(char), name0_const);
#endif

	const integer sizeN = 1024;
	char* buf = (char*)malloc(sizeN); // строка буфер

	bool bfound_loc = false;
	double fret_loc = 0.0;
	for (integer i_1 = icurrent_position__StringList; i_1 < icurrentSize_StringList; ++i_1) {
		if (!bfound) {
			// Если не найдена.

			bool bbeginstring = true; // строка ещё не кончилась
			int k = 0;
			char c;
			c = ' ';
			while (StringList[i_1][k] != '\0') {
				c = StringList[i_1][k];
				if (c == '#') {
					bbeginstring = false;
					if (buf != nullptr) {
						buf[k++] = '\n';
						buf[k] = '\0';
						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = fret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
					}
					break;
					//printf("%s \n",buf);
					//system("pause");	
				}
				if (c == '\n') {
					bbeginstring = false;
					if (buf != nullptr) {
						if (k < sizeN - 2) {
					        buf[k++] = '\n';
					        buf[k] = '\0';
						}
						else {
							printf("size vector buf is greater %lld.\n", sizeN);
							system("pause");
							exit(1);
						}
					}
					else {
						printf("buf is nullptr pointer.\n");
						system("pause");
						exit(1);
					}

					//bfound_loc = false;
					bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
					if (bfound_loc) {
						bfound = true;
						fret = fret_loc;
						icurrent_position__StringList = i_1 + 1;
					}
					//printf("%s \n",buf);
					//system("pause");
					break;
				}
				if (bbeginstring) {
					if (buf != nullptr) {
						if (k < sizeN - 1) {
					        buf[k++] = c;
						}
						else {
							printf("size vector buf is greater %lld.\n", sizeN);
							system("pause");
							exit(1);
						}
					}
					else {
						printf("buf is nullptr pointer.\n");
						system("pause");
						exit(1);
					}
					// printf("%c \n",c);
				}
			} // while
		}
		else {
			break;
		}
	}

	if (!bfound) {

		// Повторяем всё снова.
		for (integer i_1 = 0; i_1 < icurrent_position__StringList; ++i_1) {
			if (!bfound) {
				// Если не найдена.

				bool bbeginstring = true; // строка ещё не кончилась
				int k = 0;
				char c;
				c = ' ';
				while (StringList[i_1][k] != '\0') {
					c = StringList[i_1][k];
					if (c == '#') {
						bbeginstring = false;
						if (buf != nullptr) {
							if (k < sizeN - 2) {
								buf[k++] = '\n';
								buf[k] = '\0';
							}
							else {
								printf("size vector buf is greater %lld.\n",sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							printf("buf is nullptr pointer.\n");
							system("pause");
							exit(1);
						}
						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = fret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						break;
						//printf("%s \n",buf);
						//system("pause");	
					}
					if (c == '\n') {
						bbeginstring = false;
						if (buf != nullptr) {
							if (k < sizeN - 2) {
						       buf[k++] = '\n';
					           buf[k] = '\0';
							}
							else {
								printf("error!!! string length > %lld\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							printf("error!!! string is nullptr \n");
							system("pause");
							exit(1);
						}


						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = fret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						//printf("%s \n",buf);
						//system("pause");
						break;
					}
					if (bbeginstring) {
						if (buf != nullptr) {
							if (k < sizeN - 1) {
								buf[k++] = c;
							}
							else {
								printf("error!!! string length > %lld\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							k++;
						}
						// printf("%c \n",c);
					}
				} // while
			}
			else {
				break;
			}
		}
	}

	delete[] name0;
	free(buf);

	return bfound;


} // fmakesource по мотивам 2010

// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool fmakesource_float_version(const char* name0_const, float& fret)
{
	bool bfound = false;
	char* name0 = new char[strlen(name0_const) + 1];
#ifdef MINGW_COMPILLER
	strcpy(name0, name0_const);
#else
	strcpy_s(name0, (strlen(name0_const) + 1) * sizeof(char), name0_const);
#endif

	const integer sizeN = 1024;
	char* buf = (char*)malloc(sizeN); // строка буфер

	bool bfound_loc = false;
	double fret_loc = 0.0;
	for (integer i_1 = icurrent_position__StringList; i_1 < icurrentSize_StringList; ++i_1) {
		if (!bfound) {
			// Если не найдена.

			bool bbeginstring = true; // строка ещё не кончилась
			int k = 0;
			char c;
			c = ' ';
			while (StringList[i_1][k] != '\0') {
				c = StringList[i_1][k];
				if (c == '#') {
					bbeginstring = false;
					if (buf != nullptr) {
						buf[k++] = '\n';
						buf[k] = '\0';
						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = static_cast <float>(fret_loc);
							icurrent_position__StringList = i_1 + 1;
						}
					}
					break;
					//printf("%s \n",buf);
					//system("pause");	
				}
				if (c == '\n') {
					bbeginstring = false;
					if (buf != nullptr) {
						if (k < sizeN - 2) {
							buf[k++] = '\n';
							buf[k] = '\0';
						}
						else {
							printf("size vector buf is greater %lld.\n", sizeN);
							system("pause");
							exit(1);
						}
					}
					else {
						printf("buf is nullptr pointer.\n");
						system("pause");
						exit(1);
					}

					//bfound_loc = false;
					bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
					if (bfound_loc) {
						bfound = true;
						fret = static_cast <float>(fret_loc);
						icurrent_position__StringList = i_1 + 1;
					}
					//printf("%s \n",buf);
					//system("pause");
					break;
				}
				if (bbeginstring) {
					if (buf != nullptr) {
						if (k < sizeN - 1) {
							buf[k++] = c;
						}
						else {
							printf("size vector buf is greater %lld.\n", sizeN);
							system("pause");
							exit(1);
						}
					}
					else {
						printf("buf is nullptr pointer.\n");
						system("pause");
						exit(1);
					}
					// printf("%c \n",c);
				}
			} // while
		}
		else {
			break;
		}
	}

	if (!bfound) {

		// Повторяем всё снова.
		for (integer i_1 = 0; i_1 < icurrent_position__StringList; ++i_1) {
			if (!bfound) {
				// Если не найдена.

				bool bbeginstring = true; // строка ещё не кончилась
				int k = 0;
				char c;
				c = ' ';
				while (StringList[i_1][k] != '\0') {
					c = StringList[i_1][k];
					if (c == '#') {
						bbeginstring = false;
						if (buf != nullptr) {
							if (k < sizeN - 2) {
								buf[k++] = '\n';
								buf[k] = '\0';
							}
							else {
								printf("size vector buf is greater %lld.\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							printf("buf is nullptr pointer.\n");
							system("pause");
							exit(1);
						}
						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = static_cast <float>(fret_loc);
							icurrent_position__StringList = i_1 + 1;
						}
						break;
						//printf("%s \n",buf);
						//system("pause");	
					}
					if (c == '\n') {
						bbeginstring = false;
						if (buf != nullptr) {
							if (k < sizeN - 2) {
								buf[k++] = '\n';
								buf[k] = '\0';
							}
							else {
								printf("error!!! string length > %lld\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							printf("error!!! string is nullptr \n");
							system("pause");
							exit(1);
						}


						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = static_cast <float>(fret_loc);
							icurrent_position__StringList = i_1 + 1;
						}
						//printf("%s \n",buf);
						//system("pause");
						break;
					}
					if (bbeginstring) {
						if (buf != nullptr) {
							if (k < sizeN - 1) {
								buf[k++] = c;
							}
							else {
								printf("error!!! string length > %lld\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							k++;
						}
						// printf("%c \n",c);
					}
				} // while
			}
			else {
				break;
			}
		}
	}

	delete[] name0;
	free(buf);

	return bfound;


} // fmakesource_float_version по мотивам 2010


// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool imakesource_old(char *name0, int &iret)
{
	bool bfound = false;
	FILE *fp=NULL;
	

#ifdef MINGW_COMPILLER
	fp = fopen64("premeshin.txt", "r");
	int err = 0;
	if (fp != NULL) {
		err = 0;
	}
	else {
		err = 1; // ошибка открытия.
	}
#else
	errno_t err;
	err = fopen_s(&fp, "premeshin.txt", "r");
#endif

	if (fp != NULL) {
		// создание файла для чтения.
		if (err != 0) {
			printf("Open File Error\n");
			exit(0);
			// return bfound;
		}
		else
		{
			char c;
			c = ' ';
			bool bbeginstring = true; // строка ещё не кончилась
			bool b1 = true;
			int k = 0;
			char* buf = (char*)malloc(1024); // строка буфер
			if (buf != nullptr) {
				c = fgetc(fp); // читает один символ
				buf[k++] = c;
				bool bfound_loc = false;
				int iret_loc = -1;
				while (!feof(fp)) {
					b1 = true;
					c = fgetc(fp); // читает один символ
					if (c == '#') bbeginstring = false;
					if (c == '\n') {
						bbeginstring = true;
						buf[k++] = '\n';
						buf[k] = '\0';
						//bfound_loc = false;
						bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
						if (bfound_loc) {
							bfound = true;
							iret = iret_loc;
							break;
						}
						//printf("%s \n",buf);
						//system("pause");
						b1 = false;
						k = 0;
					}
					if ((b1) && (bbeginstring  )) {
						buf[k++] = c;
						// printf("%c \n",c);
					}

				}// while
				if (!bfound) {
					buf[k++] = '\n';
					buf[k] = '\0';
					//bfound_loc = false;
					bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
					if (bfound_loc) {
						bfound = true;
						iret = iret_loc;
					}
					//printf("%s \n",buf);
					//system("pause");	
				}
				free(buf);
			}
			if (fp != NULL) {
				fclose(fp);
				fp = NULL;
			}

			return bfound;
		}
	}

	return bfound;
} // imakesource по мотивам 2010


// считывает данные из входного файла
// premeshin.txt
// реализовано 23 апреля 2010.
// revised 21 августа 2019.
bool fmakesource_old(char *name0, double &fret)
{
	bool bfound = false;
	FILE *fp = NULL;
	
	// создание файла для чтения.
#ifdef MINGW_COMPILLER
	fp = fopen64("premeshin.txt", "r");
	int err = 0;
	if (fp != NULL) {
		err = 0;
	}
	else {
		err = 1; // ошибка открытия.
	}
#else
	errno_t err;
	err = fopen_s(&fp, "premeshin.txt", "r");
#endif

	if (fp != NULL) {
		if (err != 0) {
			printf("Open File Error\n");
			exit(0);
			// return bfound;
		}
		else
		{
			char c;
			c = ' ';
			bool bbeginstring = true; // строка ещё не кончилась
			bool b1 = true;
			int k = 0;
			const integer sizeN = 1024;
			char* buf = (char*)malloc(sizeN); // строка буфер
			if (buf != nullptr) {
				c = fgetc(fp); // читает один символ
				buf[k++] = c;
				bool bfound_loc = false;
				double fret_loc = 0.0;
				while (!feof(fp)) {
					b1 = true;
					c = fgetc(fp); // читает один символ
					if (c == '#') bbeginstring = false;
					if (c == '\n') {
						bbeginstring = true;
						if (buf != nullptr) {
							if (k < sizeN - 2) {
								buf[k++] = '\n';
								buf[k] = '\0';
							}
							else {
								printf("size vector buf is greater %lld.\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							printf("buf is nullptr pointer.\n");
							system("pause");
							exit(1);
						}
						//bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = fret_loc;
							break;
						}
						//printf("%s \n",buf);
						//system("pause");
						b1 = false;
						k = 0;
					}
					if ((b1) && (bbeginstring  )) {
						if (buf != nullptr) {
							if (k < sizeN - 1) {
								buf[k++] = c;
							}
							else {
								printf("size vector buf is greater %lld.\n", sizeN);
								system("pause");
								exit(1);
							}
						}
						else {
							printf("buf is nullptr pointer.\n");
							system("pause");
							exit(1);
						}
						// printf("%c \n",c);
					}

				}// while
				if (!bfound) {
					if (buf != nullptr) {
						if (k < sizeN - 2) {
							buf[k++] = '\n';
							buf[k] = '\0';
						}
						else {
							printf("size vector buf is greater %lld.\n", sizeN);
							system("pause");
							exit(1);
						}
					}
					else {
						printf("buf is nullptr pointer.\n");
						system("pause");
						exit(1);
					}
					//bfound_loc = false;
					bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
					if (bfound_loc) {
						bfound = true;
						fret = fret_loc;
					}
					//printf("%s \n",buf);
					//system("pause");	
				}
				free(buf);
			}

			if (fp != NULL) {
				fclose(fp);
				fp = NULL;
			}

			return bfound;
		}
	}

	return bfound;

} // fmakesource по мотивам 2010


#include "Mingw_input_laplas.cpp"

// считывание параметров из 
// входного файла premeshin.txt
void premeshin_new(const char *fname, integer &lmatmax, int &lb, int &ls, int &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w,
	doublereal &dgx, doublereal &dgy, doublereal &dgz, int &inx, int &iny, int &inz, doublereal &operatingtemperature,
	integer &ltdp, TEMP_DEP_POWER* &gtdps, int &lu, UNION* &my_union, bool bSTOP_Reading)
{

	

	doublereal dmult = 1.0 / rdivision_interval;


	// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
	// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
	// при повторном обращении идет обращение к ячейки хеш-таблицы.
	// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
	// Хеш-таблицы для automatic
	// Инициализация хеш-таблицы.
	for (integer i_1 = 0; i_1 < isize_shorter_hash; ++i_1) {
		bshorter_hash_X[i_1] = false;
		bshorter_hash_Y[i_1] = false;
		bshorter_hash_Z[i_1] = false;
		shorter_hash_X[i_1] = 1.0e-10;
		shorter_hash_Y[i_1] = 1.0e-10;
		shorter_hash_Z[i_1] = 1.0e-10;
	}

#ifdef MINGW_COMPILLER
	
	mingw_input_new(fname, lmatmax, lb, ls, lw, matlist, b, s, w,
		dgx, dgy, dgz, inx, iny, inz,operatingtemperature,  
		ltdp, gtdps, lu, my_union);

	integer ilb_p = 0;// Количество блоков внутри которых задана тепловая мощность.

	int lb_ic = 0; // Количество пропущенных объектов которые полностью перекрыты
	// другими телами.
	// Список объектов перекрытых другими телами.
	int* b_list_delete_in_scan = new int[lb];
	for (integer i_1 = 1; i_1 < lb; ++i_1) {
		for (integer i_2 = 0; i_2 < i_1; ++i_2) {
			if (b[i_1].g.itypegeom == PRISM && b[i_2].g.itypegeom == PRISM &&
				(b[i_1].g.xS <= b[i_2].g.xS) && (b[i_1].g.xE >= b[i_2].g.xE) &&
				(b[i_1].g.yS <= b[i_2].g.yS) && (b[i_1].g.yE >= b[i_2].g.yE) &&
				(b[i_1].g.zS <= b[i_2].g.zS) && (b[i_1].g.zE >= b[i_2].g.zE))
			{
				if (lb_ic < lb) {
					b_list_delete_in_scan[lb_ic] = i_2;
					lb_ic++;
				}
			}
		}
	}

	doublereal dpower = 0.0; // Суммарная тепловая мощность в блоках.
	integer ipoly = 0, icyl = 0, iprism = 0, icad_stl=0;
	integer ihol = 0, isol = 0, iflui = 0;
	for (integer i_1 = 0; i_1 < lb; ++i_1) {
		bool bscalok = true;
		for (integer i_4 = 0; i_4 < lb_ic; i_4++) {
			if (i_1 == b_list_delete_in_scan[i_4]) {
				bscalok = false; // Тело перекрыто другим телом.
			}
		}
		if (bscalok) {
			if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
				ihol++;
			}
			if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				iflui++;
			}
			if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				isol++;
			}
			if (b[i_1].g.itypegeom == CAD_STL)
			{
				icad_stl++;
			}
			if (b[i_1].g.itypegeom == PRISM) {
				// 0 - PRISM object

				// 13.08.2019 Автомат настройки допусков сетки. 
				// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
				// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
				// Она не может быть больше стороны кабинета деленной на 15.
				if (0.1 * fabs(b[i_1].g.xE - b[i_1].g.xS) < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * fabs(b[i_1].g.xE - b[i_1].g.xS);
				if (0.1 * fabs(b[i_1].g.yE - b[i_1].g.yS) < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * fabs(b[i_1].g.yE - b[i_1].g.yS);
				if (0.1 * fabs(b[i_1].g.zE - b[i_1].g.zS) < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * fabs(b[i_1].g.zE - b[i_1].g.zS);

				if (lb == 1) {// Течение в каверне или тест Дэвиса.
					shorter_length_for_simplificationX_BASIC = 1.0e-10;
					shorter_length_for_simplificationY_BASIC = 1.0e-10;
					shorter_length_for_simplificationZ_BASIC = 1.0e-10;
				}
				else {
					if (shorter_length_for_simplificationX_BASIC > 0.067 * (fabs(b[0].g.xE - b[0].g.xS))) shorter_length_for_simplificationX_BASIC = 0.067 * (fabs(b[0].g.xE - b[0].g.xS));
					if (shorter_length_for_simplificationY_BASIC > 0.067 * (fabs(b[0].g.yE - b[0].g.yS))) shorter_length_for_simplificationY_BASIC = 0.067 * (fabs(b[0].g.yE - b[0].g.yS));
					if (shorter_length_for_simplificationZ_BASIC > 0.067 * (fabs(b[0].g.zE - b[0].g.zS))) shorter_length_for_simplificationZ_BASIC = 0.067 * (fabs(b[0].g.zE - b[0].g.zS));
				}

				iprism++;
				if (b[i_1].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) {
					if (b[i_1].n_Sc > 0) {
						doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
						doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS) * fabs(b[i_1].g.yE - b[i_1].g.yS) * fabs(b[i_1].g.zE - b[i_1].g.zS);
						if (vol < 1.0e-36) {
							printf("ERROR: zero volume in PRISM block number %lld\n", i_1);
							system("PAUSE");
							exit(1);
						}
						//if (fabs(b[i_1].arr_Sc[0]) > 0.0) {
						if (pdiss > 0.0) {
							ilb_p++;
							//dpower += b[i_1].arr_Sc[0];
							dpower += pdiss * vol;
						}
					}
				}
			}
			if (b[i_1].g.itypegeom == CYLINDER) {
				// Cylinder

				// На тот случай если геометрия состоит только из цилиндров.
				// 13.08.2019 Автомат настройки допусков сетки. 
				// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
				// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
				switch (b[i_1].g.iPlane) {
				case XY_PLANE:
					if (0.1 * b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.Hcyl;
					if (0.1 * b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
					if (0.1 * b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-36) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1 * b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
						if (0.1 * b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
					}
					break;
				case XZ_PLANE:
					if (0.1 * b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.Hcyl;
					if (0.1 * b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
					if (0.1 * b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-36) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1 * b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
						if (0.1 * b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
					}
					break;
				case YZ_PLANE:
					if (0.1 * b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.Hcyl;
					if (0.1 * b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
					if (0.1 * b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-36) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1 * b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
						if (0.1 * b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
					}
					break;
				}


				icyl++;
				if (b[i_1].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) {
					if (b[i_1].n_Sc > 0) {
						doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
						doublereal vol = 0.0;
						vol = b[i_1].g.Hcyl * M_PI * (b[i_1].g.R_out_cyl * b[i_1].g.R_out_cyl - b[i_1].g.R_in_cyl * b[i_1].g.R_in_cyl);
						if (vol < 1.0e-36) {
							printf("ERROR: zero volume in CYLINDER block number %lld\n", i_1);
							system("PAUSE");
							exit(1);
						}
						if (pdiss > 0.0) {
							ilb_p++;

							dpower += pdiss * vol;
							//printf("ERROR: non zero power in cylinder object.\n");
							//system("PAUSE");
							//exit(1);
						}
					}
				}
			}
			if (b[i_1].g.itypegeom == POLYGON) {
				// Polygon
				ipoly++;
				if (b[i_1].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) {
					if (b[i_1].n_Sc > 0) {
						doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);

						// Объём полигона.
						doublereal vol_poly = Volume_polygon(b[i_1].g.nsizei, b[i_1].g.xi, b[i_1].g.yi, b[i_1].g.zi, b[i_1].g.hi, b[i_1].g.iPlane_obj2);
						if (vol_poly < 1.0e-36) {
							printf("ERROR: zero volume in POLYGON block number %lld\n", i_1);
							system("PAUSE");
							exit(1);
						}

						if (pdiss > 0.0) {
							ilb_p++;

							dpower += pdiss*vol_poly;

							//printf("ERROR: non zero power in polygon object.\n");
							//system("PAUSE");
							//exit(1);
						}
					}
				}
			}
		}
	}

	delete[] b_list_delete_in_scan;
	b_list_delete_in_scan = nullptr;

	doublereal dsoupow = 0.0; // интегральная тепловая мощность плоских источников тепла.
	for (integer i_1 = 0; i_1 < ls; ++i_1) {
		dsoupow += s[i_1].power;
	}


	printf("Apriory quick model statistics:\n");
	printf("number of thermal power blocks lb_p=%lld\n", ilb_p);
	printf("Blocks integral power =%e W\n", dpower);
	printf("number of sources ls=%lld\n", ls);
	printf("Sources integral power = %e W\n", dsoupow);
	printf("Full total power = %e W\n", dpower + dsoupow);
	// Запоминаем полное тепловыделение в модели.
	d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = dpower + dsoupow;
	printf("number of blocks lb=%lld\n", lb);
	printf("PRISMS = %lld, CYLINDERS = %lld, POLYGONS = %lld, CAD_STL = %lld\n", iprism, icyl, ipoly, icad_stl);
	printf("SOLID: %lld\n", isol);
	printf("HOLLOW: %lld\n", ihol);
	printf("FLUID: %lld\n", iflui);
	printf("number of walls lw=%lld\n", lw);
	printf("number of units lu=%lld\n", lu);


#endif

#ifndef MINGW_COMPILLER

	// eqin - информация о наборе решаемых уравнений.

	// dgx, dgy, dgz - вектор силы тяжести.
	// inx, iny, inz - количество точек по каждой из осей.

	FILE* fp=NULL;
	//errno_t err1;
	//if ((err1 = fopen_s(&fp, fname, "r")) != 0) 
	if ((fopen_s(&fp, fname, "r")) != 0)
	{
		printf("Error!!! No input File premeshin.txt.\n");
		printf("You must use the graphical user interface\n");
		printf("AliceMesh_v0_45.exe in Delphi, which will \n");
		printf("prepare the model and write the premeshin.txt\n");
		printf("file in the required format.\n");
		//system("PAUSE");
		system("pause");
		// Если файла premeshin.txt нет то мы не можем ничего обрабатывать
		// т.к. элементарно отсутствуют входные данные. Поэтому мы выходим из 
		// приложения.
		exit(1);
	}
	else
	{
		if (fp != NULL) {

			
				fclose(fp);
				fp = NULL;
			

			double fin = 0.0;
			integer din = 0;
			int idin = 0;
			doublereal scale = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания

			// Приостанавливает считывание если переменная не найдена.
			//bool bSTOP_Reading = true; 


			if (imakesource("only_solid_visible", idin)) {
			    // Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					switch (idin) {
						case 0: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
						case 1: ionly_solid_visible= WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE; break;
						default: ionly_solid_visible= WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; break;
					}
					//ionly_solid_visible = static_cast<integer>(idin);
					switch (ionly_solid_visible) {
					case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE:
						std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
						break;
					case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE:
						std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
						break;
					}
				}
				else {
					std::cout << "ionly_solid_visible must be equal 0 or 1. now value=" << idin << std::endl;
					system("pause");
					ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;
				}
			}
			else {
				printf("WARNING!!! only_solid_visible not found in file premeshin.txt\n");
				ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE; // Показываем всё
				switch (ionly_solid_visible) {
				case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE:
					std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
					break;
				case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE:
					std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
					break;
				}
				if (bSTOP_Reading) system("pause");
			}
		

			if (fmakesource("mlength", fin)) {
				// Найдено успешно.
				if (scale <= 0.0) {
					printf("ERROR!!! scale must be >0.0\n");
					system("pause");
					exit(1);
				}
				scale = static_cast<doublereal>(fin);
#ifndef NO_OPENGL_GLFW
				scale_all = (GLfloat)( 1.0f/scale);
#endif
			}
			else {
				printf("WARNING!!! mlength not found in file premeshin.txt\n");
				scale = 1.0e-3; // mm
#ifndef NO_OPENGL_GLFW
				scale_all = 1000.0f;
#endif
				printf("scale =%e\n", scale);
				if (bSTOP_Reading) system("pause");
			}

			{// Считываем свободный параметр в самом начале, т.к. он используется сразу в коде.

				char name0[1000] = "egddata";
				name0[0] = '\0'; strcat_s(name0, "free_debug_parametr1");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					// Свободный параметр для выполнения отладки передаваемый из интерфейса внутрь солвера.
					free_debug_parametr1 = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					free_debug_parametr1 = 0.0; // Свободный параметр для выполнения отладки передаваемый из интерфейса внутрь солвера.
					printf(" free_debug_parametr1=%e\n", free_debug_parametr1);
					if (bSTOP_Reading) system("pause");
				}
			}
			
			if (imakesource("lmatmax", idin)) {
				// Найдено успешно.
				if (idin < 1) {
					printf("Error: lmatmax mast be > 0. now value=%d\n",idin);
					system("pause");
					exit(1);
				}
				else {
					lmatmax = static_cast<integer>(idin);
					//printf("lmatmax =%lld\n", lmatmax);
				}
			}
			else {
				printf("WARNING!!! lmatmax not found in file premeshin.txt\n");
				lmatmax = 0; // нет материалов
				printf("lmatmax =%lld\n", lmatmax);
				if (bSTOP_Reading) system("pause");
			}

			
			if (imakesource("lb", idin)) {
				// Найдено успешно.
				if (idin <= 0) {
					printf("Error: lb mast be > 0. now value=%d\n", idin);
					system("pause");
					exit(1);
				}
				else {
					lb = static_cast<integer>(idin);
					//printf("lb =%lld\n", lb);
				}
			}
			else {
				printf("WARNING!!! lb not found in file premeshin.txt\n");
				lb = 0; // нет блоков
				printf("lb =%d\n", lb);
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("ls", idin)) {
				// Найдено успешно.
				if (idin >= 0) {
					ls = static_cast<integer>(idin);
					//printf("ls =%lld\n", ls);
				}
				else {
					printf("Error: ls = %d must be >=0\n",idin);
					system("pause");
					exit(1);
				}
			}
			else {
				printf("WARNING!!! ls not found in file premeshin.txt\n");
				ls = 0; // нет плоских бесконечно тонких источников тепла.
				printf("ls =%d\n", ls);
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("lw", idin)) {
				// Найдено успешно.
				if (idin >= 0) {
					lw = static_cast<integer>(idin);
					//printf("lw =%lld\n", lw);
				}
				else {
					printf("Error: lw = %d must be >=0\n", idin);
					system("pause");
					exit(1);
				}
			}
			else {
				printf("WARNING!!! lw not found in file premeshin.txt\n");
				lw = 0; // нет плоских бесконечно тонких источников тепла.
				printf("lw =%d\n", lw);
				if (bSTOP_Reading) system("pause");
			}


			// количество уникальных данных с табличными данными 
			// по зависимости расеиваемой мощности от температуры.
			if (imakesource("iltdp", idin)) {
				// Найдено успешно.
				if (idin >= 0) {
				    ltdp = static_cast<integer>(idin);
				    //printf("ltdp =%lld\n", ltdp);
				}
				else {
					printf("Error: iltdp = %d must be >=0\n", idin);
					system("pause");
					exit(1);
				}
			}
			else {
				printf("WARNING!!! iltdp not found in file premeshin.txt\n");
				ltdp = 0; // нет плоских бесконечно тонких источников тепла.
				printf("ltdp =%lld\n", ltdp);
				if (bSTOP_Reading) system("pause");
			}

			// Считываем значение вектора силы тяжести:
			if (fmakesource("gx", fin)) {
				// Найдено успешно.
				dgx = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! gx not found in file premeshin.txt\n");
				dgx = 0.0; // Сила тяжести по направлению Ох отсутствует.
				printf("gx =%e\n", dgx);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("gy", fin)) {
				// Найдено успешно.
				dgy = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! gy not found in file premeshin.txt\n");
				dgy = 0.0; // Сила тяжести по направлению Оy отсутствует.
				printf("gy =%e\n", dgy);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("gz", fin)) {
				// Найдено успешно.
				dgz = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! gz not found in file premeshin.txt\n");
				dgz = 0.0; // Сила тяжести по направлению Оz отсутствует.
				printf("gz =%e\n", dgz);
				if (bSTOP_Reading) system("pause");
			}


			// считываем количество точек на каждой координатной оси
			if (imakesource("inx", idin)) {
				// Найдено успешно.
				if (inx > 0) {
					inx = static_cast<integer>(idin);
					//printf("inx =%d\n", inx);
				}
				else {
					printf("ERROR: inx=%d\n",inx);
					system("pause");
				}
			}
			else {
				printf("WARNING!!! inx not found in file premeshin.txt\n");
				inx = 0; // нет информации о количестве узлов сетки в направлении оси Ох.
				printf("inx =%d\n", inx);
				if (bSTOP_Reading) system("pause");
			}

		
			if (imakesource("iny", idin)) {
				// Найдено успешно.
				if (iny > 0) {
					iny = static_cast<integer>(idin);
					//printf("iny =%d\n", iny);
				}
			    else {
				   printf("ERROR: inx=%d\n", inx);
				   system("pause");
		     	}
			}
			else {
				printf("WARNING!!! iny not found in file premeshin.txt\n");
				iny = 0; // нет информации о количестве узлов сетки в направлении оси Оy.
				printf("iny =%d\n", iny);
				if (bSTOP_Reading) system("pause");
			}


			if (imakesource("inz", idin)) {
				// Найдено успешно.
				if (inz>0) {
				   inz = static_cast<integer>(idin);
				   //printf("inz =%d\n", inz);
				}
				else {
					printf("ERROR: inx=%d\n", inx);
					system("pause");
				}
			}
			else {
				printf("WARNING!!! inz not found in file premeshin.txt\n");
				inz = 0; // нет информации о количестве узлов сетки в направлении оси Оz.
				printf("inz =%d\n", inz);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("operatingtemperature", fin)) {
				// Найдено успешно.
				if (fin < -273.15) {
					fin = -273.15;
					printf("Error!!! operatingtemperature < -273.15 C\n");
					system("pause");
				}
				operatingtemperature = static_cast<doublereal>(fin);  // Operating Temperature
				operating_temperature_for_film_coeff = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! operatingtemperature not found in file premeshin.txt\n");
				operatingtemperature = 30.0; // Не задана температура окружающей среды.
				operating_temperature_for_film_coeff = 30.0; 
				printf("operating_temperature_for_film_coeff = operatingtemperature =%e\n", operatingtemperature);
				if (bSTOP_Reading) system("pause");
			}

			// инициализация компонент скорости константой.
			// Единые значения для всей расчётной области.
			// initialization value.
			
			if (fmakesource("SpeedInitializationVx", fin)) {
				// Найдено успешно.
				starting_speed_Vx = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! SpeedInitializationVx not found in file premeshin.txt\n");
				starting_speed_Vx = 0.0; // инициализируем нулевой скоростью в направлении оси Ох.
				printf("starting_speed_Vx =%e\n", starting_speed_Vx);
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("SpeedInitializationVy", fin)) {
				// Найдено успешно.
				starting_speed_Vy = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! SpeedInitializationVy not found in file premeshin.txt\n");
				starting_speed_Vy = 0.0; // инициализируем нулевой скоростью в направлении оси Оy.
				printf("starting_speed_Vy =%e\n", starting_speed_Vy);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("SpeedInitializationVz", fin)) {
				// Найдено успешно.
				starting_speed_Vz = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! SpeedInitializationVz not found in file premeshin.txt\n");
				starting_speed_Vz = 0.0; // инициализируем нулевой скоростью в направлении оси Оz.
				printf("starting_speed_Vz =%e\n", starting_speed_Vz);
				if (bSTOP_Reading) system("pause");
			}


			// Считываем координаты опорной точки через которую проходит пользовательская линия Variation Plot.
			if (fmakesource("XYPlotXo", fin)) {
				// Найдено успешно.
				Tochka_position_X0_for_XY_Plot = scale * (static_cast<doublereal>(fin));
			}
			else {
				printf("WARNING!!! XYPlotXo not found in file premeshin.txt\n");
				Tochka_position_X0_for_XY_Plot = 0.0; 
				printf("Tochka_position_X0_for_XY_Plot =%e\n", Tochka_position_X0_for_XY_Plot);
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("XYPlotYo", fin)) {
				// Найдено успешно.
				Tochka_position_Y0_for_XY_Plot = scale * (static_cast<doublereal>(fin));
			}
			else {
				printf("WARNING!!! XYPlotYo not found in file premeshin.txt\n");
				Tochka_position_Y0_for_XY_Plot = 0.0;
				printf("Tochka_position_Y0_for_XY_Plot =%e\n", Tochka_position_Y0_for_XY_Plot);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("XYPlotZo", fin)) {
				// Найдено успешно.
				Tochka_position_Z0_for_XY_Plot = scale * (static_cast<doublereal>(fin));
			}
			else {
				printf("WARNING!!! XYPlotZo not found in file premeshin.txt\n");
				Tochka_position_Z0_for_XY_Plot = 0.0;
				printf("Tochka_position_Z0_for_XY_Plot =%e\n", Tochka_position_Z0_for_XY_Plot);
				if (bSTOP_Reading) system("pause");
			}

			// Направление линии совпадает с направлением 
			// одной из осей декартовой прямоугольной системы координат:
			// 0 - Ox, 1 - Oy, 2 - Oz.
			if (imakesource("line_directional", idin)) {
				// Найдено успешно.
				switch (idin) {
					case 0: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
					case 1: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					case 2: idirectional_for_XY_Plot = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
					default: idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL;	break;
				}
				//idirectional_for_XY_Plot  = static_cast<integer>(idin);
				//printf("iny =%d\n", iny);
			}
			else {
				printf("WARNING!!! line_directional not found in file premeshin.txt\n");
				idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; // нет информации о количестве узлов сетки в направлении оси Оy.
				printf("idirectional_for_XY_Plot =%d\n", idirectional_for_XY_Plot);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("etalon_max_size_ratio", fin)) {
				// Найдено успешно.
				if (fin < 1.05) {					
					printf("Error!!! etalon_max_size_ratio<1.05. Value = 2 - 10 is recomended.\n");
					printf("current value etalon_max_size_ratio=%e\n",fin);
					fin = 1.05;
					system("pause");
					exit(1);
				}
				etalon_max_size_ratio = static_cast<doublereal>(fin); // подробность расчётной сетки.
			}
			else {
				printf("WARNING!!! etalon_max_size_ratio not found in file premeshin.txt\n");
				etalon_max_size_ratio = 10.0;
				printf("etalon_max_size_ratio =%e\n", etalon_max_size_ratio);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("etalon_max_size_ratio2", fin)) {
				// Найдено успешно.
				if (fin < 10) {
					printf("Error!!! etalon_max_size_ratio2<10. Value = 30 - infinity is recomended.\n");
					printf("current value etalon_max_size_ratio2=%e\n", fin);
					fin = 10;
					system("pause");
					exit(1);
				}
				// Критерий качества расчётной сетки на основе FlowVision.
				etalon_max_size_ratio2 = static_cast<doublereal>(fin); 
			}
			else {
				printf("WARNING!!! etalon_max_size_ratio2 not found in file premeshin.txt\n");
				etalon_max_size_ratio2 = 1.0e9; // ignoring !!!
				printf("etalon_max_size_ratio2 =%e\n", etalon_max_size_ratio2);
				if (bSTOP_Reading) system("pause");
			}

			

			if (imakesource("snap_to_grid", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)|| (idin == 2) || (idin == 3)) {

					// 0.	none
			        // 1.	Snap to grid
			        // 2.	Snap to grid ALICE
			        // 3.	Snap to grid ++					
					switch (idin) {
					case 0: bsnap_TO_global = 0;  break;
					case 1: bsnap_TO_global = 1;  break;
					case 2: bsnap_TO_global = 2;  break;
					case 3: bsnap_TO_global = 3;  break;
					default: bsnap_TO_global = 1;  break;
					}

					//printf("snap_to_grid =%d\n", idin);
				}
				else {
					printf("Error!!! snap_to_grid must be equal 0, 1, 2 or 3. Current value = %d\n",idin);
					system("pause");
					bsnap_TO_global = 0;
				}
			}
			else {
				printf("WARNING!!! snap_to_grid not found in file premeshin.txt\n");
				bsnap_TO_global = 0; // выключен.
				printf("bsnap_TO_global =0\n");
				if (bSTOP_Reading) system("pause");
			}

			// Выбор решающего устройства: либо amg1r5 либо BiCGStab+ILU2.
			if (imakesource("SolverSetting", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin < 14)) {
					iswitchsolveramg_vs_BiCGstab_plus_ILU2 = static_cast<integer>(idin);
					//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);
				}
				else {
					printf("Error !!! must be 0<= iswitchsolveramg_vs_BiCGstab_plus_ILU2 <14. current value iswitchsolveramg_vs_BiCGstab_plus_ILU2=%d\n", idin);
					system("pause");
					idin = 0; // BiCGStab+ilu2 solver default.
					iswitchsolveramg_vs_BiCGstab_plus_ILU2 = static_cast<integer>(idin);
					//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);

				}
			}
			else {
				printf("WARNING!!! SolverSetting not found in file premeshin.txt\n");
				iswitchsolveramg_vs_BiCGstab_plus_ILU2 = 0; //BiCGStab +ILU* solver.
				printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);
				if (bSTOP_Reading) system("pause");
			}

			 
			// Выбор решающего устройства: либо РУМБА0.14 либо BiCGStab+ILU6.
			if (imakesource("StaticStructuralSolverSetting", idin)) {
				// Найдено успешно.
				switch (static_cast<integer>(idin)) {
					case 0: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
					case 1: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER; break;
					case 2: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::CAMG_RUMBA_v0_14_SECOND_T_SOLVER; break;
					case 3: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMG1R5_SECOND_T_SOLVER; break;
					case 4: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::AMGCL_SECONT_T_SOLVER; break;
					case 5: iswitchsolveramg_vs_BiCGstab_plus_ILU6 = SECOND_T_SOLVER_ID_SWITCH::ICCG_SECOND_T_SOLVER; break;
					default: iswitchsolveramg_vs_BiCGstab_plus_ILU6= SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; break;
				}
				//iswitchsolveramg_vs_BiCGstab_plus_ILU6 = static_cast<integer>(idin);
				//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU6 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU6);
			}
			else {
				printf("WARNING!!! StaticStructuralSolverSetting not found in file premeshin.txt\n");
				iswitchsolveramg_vs_BiCGstab_plus_ILU6 = SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; //BiCGStab +ILU* solver.
				printf("iswitchsolveramg_vs_BiCGstab_plus_ILU6 =BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER\n");
				if (bSTOP_Reading) system("pause");
			}
			

			if (imakesource("PressureVelocityCoupling", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 1) {
						// SIMPLEC algorithm.
						iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby;
					}
					else {
						// SIMPLE algorithm 1972.
						iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;
					}
					//printf("iSIMPLE_alg =%d\n", idin);
				}
				else {
					printf("ERROR!!! iSIMPLE_alg must be equal 0 or 1. Current value =%d\n",idin);
					system("pause");
					// SIMPLE algorithm 1972.
					iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;
				}
			}
			else {
				printf("WARNING!!! PressureVelocityCoupling not found in file premeshin.txt\n");
				// SIMPLE algorithm 1972.
				iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;
				printf("iSIMPLE_alg = SIMPLE_Carretto\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("FlowSchemePrefix", idin)) {
				switch (idin) {
				case 0:  iprefix_Scheme_Flow = CR;
					break;
				case 1:  iprefix_Scheme_Flow = UDS;
					break;
				case 2:  iprefix_Scheme_Flow = COMB;
					break;
				case 3:  iprefix_Scheme_Flow = POLY;
					break;
				case 4:  iprefix_Scheme_Flow = EXP;
					break;
				case 5:  iprefix_Scheme_Flow = BULG;
					break;
				case 6:  iprefix_Scheme_Flow = POW;
					break;
				default:
					iprefix_Scheme_Flow = UDS;
					break;
				}
			}
			else {
				//printf("WARNING!!! FlowSchemePrefix not found in file premeshin.txt\n");
				iprefix_Scheme_Flow = UDS; // UDS самая стабильная схема.
				//printf("iprefix_Scheme_Flow = Upwind\n");
				//if (bSTOP_Reading) system("pause");
			}


			if (imakesource("FlowScheme", idin)) {
				// Найдено успешно.
				// выбор схемы для потока жидкости.
			    // Внимание эти определения должны полностью соответствовать 
			    // определениям в файле my_approx_convective2.c
				switch (idin) {
				case 1: iFLOWScheme = UNEVEN_MUSCL; break; // MUSCL 2
				case 2: iFLOWScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
				case 3: iFLOWScheme = UNEVEN_HLPA; break; // HLPA 2
				case 4: iFLOWScheme = UNEVEN_SMART; break; // SMART 3
				case 5: iFLOWScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
				case 6: iFLOWScheme = UNEVEN_SMARTER; break; // SMARTER 3
				case 7: iFLOWScheme = UNEVEN_LPPA; break; // LPPA 3
				case 8: iFLOWScheme = UNEVEN_VONOS; break; // VONOS 3
				case 9: iFLOWScheme = UNEVEN_STOIC; break; // STOIC
				case 10: iFLOWScheme = UNEVEN_CLAM; break; // CLAM
				case 11: iFLOWScheme = UNEVEN_OSHER; break; // OSHER
				case 12: iFLOWScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
				case 13: iFLOWScheme = UNEVEN_SUPER_C; break; // SUPER_C
				case 14: iFLOWScheme = UNEVEN_ISNAS; break; // ISNAS
				case 15: iFLOWScheme = UNEVEN_CUBISTA; break; // CUBISTA
				case 16: iFLOWScheme = UNEVEN_GAMMA; break; // GAMMA
				case 17: iFLOWScheme = UNEVEN_COPLA; break; // COPLA
				case 18: iFLOWScheme = UNEVEN_SECBC; break; // SECBC
				case 19: iFLOWScheme = UNEVEN_SGSD; break; // SGSD
				case 20: iFLOWScheme = UNEVEN_WENO5; break; // WENO5	
				default: iFLOWScheme = UDS; break; // UDS самая стабильная схема.
				}
				//printf("iFLOWScheme =%lld\n", iFLOWScheme);
			}
			else {
				printf("WARNING!!! FlowScheme not found in file premeshin.txt\n");
				iFLOWScheme = UDS; // UDS самая стабильная схема.
				printf("iFLOWScheme = Upwind\n");
				if (bSTOP_Reading) system("pause");
			}
			

			if (imakesource("SchemeTemperature", idin)) {
				// Найдено успешно.
				// выбор схемы для температуры в потоке жидкости.
			    // Внимание эти определения должны полностью соответствовать 
			    // определениям в файле my_approx_convective2.c
				switch (idin) {
				case 1: iTEMPScheme = UNEVEN_MUSCL; break; // MUSCL 2
				case 2: iTEMPScheme = UNEVEN_SOUCUP; break; // SOUCUP [MINMOD] 2
				case 3: iTEMPScheme = UNEVEN_HLPA; break; // HLPA 2
				case 4: iTEMPScheme = UNEVEN_SMART; break; // SMART 3
				case 5: iTEMPScheme = UNEVEN_WACEB; break; // WACEB 3 TVD
				case 6: iTEMPScheme = UNEVEN_SMARTER; break; // SMARTER 3
				case 7: iTEMPScheme = UNEVEN_LPPA; break; // LPPA 3
				case 8: iTEMPScheme = UNEVEN_VONOS; break; // VONOS 3
				case 9: iTEMPScheme = UNEVEN_STOIC; break; // STOIC
				case 10: iTEMPScheme = UNEVEN_CLAM; break; // CLAM
				case 11: iTEMPScheme = UNEVEN_OSHER; break; // OSHER
				case 12: iTEMPScheme = UNEVEN_EXPONENTIAL; break; // EXPONENTIAL
				case 13: iTEMPScheme = UNEVEN_SUPER_C; break; // SUPER_C
				case 14: iTEMPScheme = UNEVEN_ISNAS; break; // ISNAS
				case 15: iTEMPScheme = UNEVEN_CUBISTA; break; // CUBISTA
				case 16: iTEMPScheme = UNEVEN_GAMMA; break; // GAMMA
				case 17: iTEMPScheme = UNEVEN_COPLA; break; // COPLA
				case 18: iTEMPScheme = UNEVEN_SECBC; break; // SECBC
				case 19: iTEMPScheme = UNEVEN_SGSD; break; // SGSD
				case 20: iTEMPScheme = UNEVEN_WENO5; break; // WENO5
				default: iTEMPScheme = UDS; break; // UDS самая стабильная схема.
				}
				//printf("iTEMPScheme  =%lld\n", iTEMPScheme);
			}
			else {
				printf("WARNING!!! SchemeTemperature not found in file premeshin.txt\n");
				iTEMPScheme = UDS; // UDS самая стабильная схема.
				printf("iTEMPScheme = Upwind \n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("mesh_generator_algorithm", idin)) {
				// Найдено успешно.
				// Выбор сеточного генератора.
				switch (idin) {
					case 0: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; break;
					case 1: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER; break;
					case 2: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
					default: iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; break;
				}
				//iswitchMeshGenerator = static_cast<integer>(idin);
				//printf("iswitchMeshGenerator =%lld\n", iswitchMeshGenerator);
			}
			else {
				printf("WARNING!!! mesh_generator_algorithm not found in file premeshin.txt\n");
				iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; // Coarse Mesh.
				printf("iswitchMeshGenerator = COARSEMESHGEN_MESHER\n");
				if (bSTOP_Reading) system("pause");
			}
			

			if (imakesource("Schedule", idin)) {
				// Найдено успешно.
				steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY;
				din = static_cast<integer>(idin);
				if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 5) || (din == 6) || 
					(din == 7) || (din == 8) || (din == 9) || (din == 10) || (din == 11) || (din == 12) || (din == 13)) {
					// 0 - thermal only steady state calculation,
					// 1 - thermal only unsteady calculation,
					// 2 - mesh generator only.
					// 3 - fluid dynamic steady state.
					// 5 - Static Structural (Thermal solver #2)
					// 6 - Thermal Stress
					// 7 - Unsteady thermal solver #2
					// 8 - Visualisation only
					// 9 - cfd unsteady fluid dynamic.
					// 10 - NETWORK_T Графовый метод решения уравнения теплопроводности.
			        // 11 - UNSTEADY NETWORK_T Нестационарный графовый метод решения уравнения теплопроводности.
					// 12 - Нестационарная механика,
					// 13 - Нестационарная механика совместно с нестационарной теплопередачей.
					switch(din) {
						case 0: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE; break;
						case 1: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE; break;
						case 2: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
						case 3: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_STEADY; break;
						case 5: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL; break;
						case 6: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE; break;
						case 7: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER; break;
						case 8: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT; break;
						case 9: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY; break;
						case 10: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T;  break;
						case 11: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY; break;
						case 12: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL;  break;
						case 13: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE;  break;
						default: steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; break;
					}
					//steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady -1 calculation.
					//printf("steady_or_unsteady_global_determinant =%lld\n",din);
					//system("PAUSE");
				}
				else {
					printf("error input parametr steady or unsteady calculation\n");
					system("PAUSE");
					exit(1);
				}
				//printf("steady_or_unsteady_global_determinant  =%lld\n", steady_or_unsteady_global_determinant );
			}
			else {
				printf("WARNING!!! Schedule not found in file premeshin.txt\n");
				steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::MESHER_ONLY; // Mesh Generator only
				printf("steady_or_unsteady_global_determinant =Mesh Generator only\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("TimeStep", idin)) {
				// Найдено успешно.
				// Закон изменения шага по времени.
				din = static_cast<integer>(idin);
				if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 4)) {
					switch (din) {
					case 0: glTSL.id_law = TIME_STEP_lAW_SELECTOR::LINEAR;
						break;
					case 1: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE;
						break;
					case 2: glTSL.id_law = TIME_STEP_lAW_SELECTOR::SQUARE_WAVE2;
						break;
					case 3: glTSL.id_law = TIME_STEP_lAW_SELECTOR::HOT_COLD;
						break;
					case 4: glTSL.id_law = TIME_STEP_lAW_SELECTOR::PIECEWISE_CONSTANT;
						break;
					}
				}
				else {
					printf("error input parametr timestep law\n");
					system("PAUSE");
					exit(1);
				}
				//printf("glTSL.id_law =%lld\n", glTSL.id_law);
			}
			else {
				printf("WARNING!!! TimeStep not found in file premeshin.txt\n");
				glTSL.id_law = TIME_STEP_lAW_SELECTOR::LINEAR; // default Linear.
				printf("glTSL.id_law =Linear\n");
				if (bSTOP_Reading) system("pause");
			}
			
			if (fmakesource("Factor_a_for_Linear", fin)) {
				// Найдено успешно.
				if ((fin < 0.0) || (fin >= 1.0)) {
					printf("error input parametr timestep law Factor a=%e\n",fin);
					system("PAUSE");
					exit(1);
				}
				glTSL.Factor_a_for_Linear = static_cast<doublereal>(fin); // Factor_a
			}
			else {
				printf("WARNING!!! Factor_a_for_Linear not found in file premeshin.txt\n");
				glTSL.Factor_a_for_Linear = 0.3; // Знаменатель геометрической прогрессии.
				printf("glTSL.Factor_a_for_Linear = 0.3\n");
				if (bSTOP_Reading) system("pause");
			}


			
			
			if (fmakesource("tau", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr timestep law tau must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau = static_cast<doublereal>(fin); // длительность импульса.
				
			}
			else {
				printf("WARNING!!! tau not found in file premeshin.txt\n");
				glTSL.tau = 2.1e-3; // Длительность радиоимпульса.
				printf("glTSL.tau = 2.1e-3\n");
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("DutyCycle", fin)) {
				// Найдено успешно.
				if (fin <= 1.0) {
					printf("error input parametr glTSL.Q must be strongly > 1.0\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.Q = static_cast<doublereal>(fin); // Скважность.
			}
			else {
				printf("WARNING!!! DutyCycle not found in file premeshin.txt\n");
				glTSL.Q = 10.0; // Скважность.
				printf("glTSL.Q =%e\n", glTSL.Q);
				if (bSTOP_Reading) system("pause");
			}
			
			// Параметры импульсного режима SquareWave2.
			if (fmakesource("m1", fin)) {
				// Найдено успешно.
				if ((fin < 0.0) || (fin >= 1.0)) {
					printf("error input parametr timestep law SquareWave2 multiplyer\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.m1 = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! m1 not found in file premeshin.txt\n");
				glTSL.m1 = 0.33333; // SquareWave2 parameters.
				printf("glTSL.m1 =%e\n", glTSL.m1);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("tau1", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau1 = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! tau1 not found in file premeshin.txt\n");
				glTSL.tau1 = 180.0; // SquareWave2 parameters.
				printf("glTSL.tau1 =%e\n", glTSL.tau1);
				if (bSTOP_Reading) system("pause");
			}
					   			
			if (fmakesource("tau2", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau2 = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! tau2 not found in file premeshin.txt\n");
				glTSL.tau2 = 240.0; // SquareWave2 parameters.
				printf("glTSL.tau2 =%e\n", glTSL.tau2);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("tau_pause", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau_pause = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! tau_pause not found in file premeshin.txt\n");
				glTSL.tau_pause = 5040.0; // SquareWave2 parameters.
				printf("glTSL.tau_pause =%e\n", glTSL.tau_pause);
				if (bSTOP_Reading) system("pause");
			}

			// Параметры импульсного режима SquareWave2.
			if (fmakesource("off_multiplyer", fin)) {
				// Найдено успешно.
				if ((fin < 0.0) || (fin >= 1.0)) {
					printf("error input parametr timestep law SquareWave2 off_multiplyer\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.off_multiplyer = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!!  off_multiplyer not found in file premeshin.txt\n");
				glTSL.off_multiplyer = 0.0; // SquareWave2 parameters.
				printf("glTSL.off_multiplyer =%e\n", glTSL.off_multiplyer);
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("n", idin)) {
				// Найдено успешно.
				glTSL.n_cycle = static_cast<integer>(idin);
				//printf("glTSL.n_cycle =%lld\n", glTSL.n_cycle);
			}
			else {
				printf("WARNING!!! n not found in file premeshin.txt\n");
				glTSL.n_cycle = 6; // SquareWave2 parameters.
				printf("glTSL.n_cycle =%lld\n", glTSL.n_cycle);
				if (bSTOP_Reading) system("pause");
			}
			
			if (fmakesource("T", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr timestep law SquareWave2 T: Period must be strongly positive=%e\n",fin);
					system("PAUSE");
					exit(1);
				}
				glTSL.T_all = static_cast<doublereal>(fin);				
			}
			else {
				printf("WARNING!!! T not found in file premeshin.txt\n");
				glTSL.T_all = 86400.0; // SquareWave2 parameters.
				printf("glTSL.T_all =%e\n", glTSL.T_all);
				if (bSTOP_Reading) system("pause");
			}

			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl < 0.0) {
				printf("error in parameters Square Wave SquareWave2 time step law.\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}

			
			if (fmakesource("on_time_double_linear", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr on_time_double_linear law hot cold reshime must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.on_time_double_linear = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! on_time_double_linear not found in file premeshin.txt\n");
				glTSL.on_time_double_linear = 4.0; // Hot Cold reshime parameters.
				printf("glTSL.on_time_double_linear =%e\n", glTSL.on_time_double_linear);
				if (bSTOP_Reading) system("pause");
			}

			// Время окончания нестационарного моделирования при расчёте теплопередачи в твёрдом теле.
			if (fmakesource("EndTime", fin)) {
				// Найдено успешно.
				if (fin < 0.0) {
					printf("error input parametr globalEndTimeUnsteadyTemperatureCalculation must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				globalEndTimeUnsteadyTemperatureCalculation = static_cast<doublereal>(fin);
			}
			else {
				printf("WARNING!!! EndTime not found in file premeshin.txt\n");
				globalEndTimeUnsteadyTemperatureCalculation = 2.0; // Продолжительность моделирования. 
				printf("globalEndTimeUnsteadyTemperatureCalculation =%e\n", globalEndTimeUnsteadyTemperatureCalculation);
				if (bSTOP_Reading) system("pause");
			}

			// Piecewise Constant 20.12.2019

			if (imakesource("n_string_PiecewiseConst", idin)) {
				// Найдено успешно.
				// Закон изменения шага по времени.
				din = static_cast<integer>(idin);
				if (din >= 0) {
					glTSL.n_string_PiecewiseConst = din;

					if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY)) {
						glTSL.n_string_PiecewiseConst = 0;// Нету никакого PiecewiseConst, мы просто запускаем мешер.
					}
				}
				else {
					printf("error input parametr timestep law PiecewiseConst\n");
					system("PAUSE");
					exit(1);
				}
			}
			else {
				printf("WARNING!!! timestep law PiecewiseConst not found in file premeshin.txt\n");
				glTSL.n_string_PiecewiseConst = 0; // default Linear.
				printf("glTSL.n_string_PiecewiseConst =0\n");
				if (bSTOP_Reading) system("pause");
			}

			if (glTSL.n_string_PiecewiseConst > 0) {
				if (glTSL.table_law_piecewise_constant != nullptr) {
					delete[] glTSL.table_law_piecewise_constant;
					glTSL.table_law_piecewise_constant = nullptr;
				}
				glTSL.table_law_piecewise_constant = new PiecewiseConstantTimeStepLawTimeStepLaw[glTSL.n_string_PiecewiseConst];
				FILE* fp_piecewise_const = NULL;
				errno_t err_piecewise_const;

#ifdef MINGW_COMPILLER
				fp_piecewise_const = fopen64("premeshin_piecewise_const.txt", "r");
				if (fp_piecewise_const != nullptr) {
					err_piecewise_const = 0;
				}
				else {
					err_piecewise_const = 1; // ошибка открытия.
				}
#else
				err_piecewise_const = fopen_s(&fp_piecewise_const, "premeshin_piecewise_const.txt", "r");
#endif

				if (fp_piecewise_const != NULL) {
					// создание файла для чтения.
					if (err_piecewise_const != 0) {
						printf("Open File Error _piecewise_const\n");
						exit(0);
						// return bfound;
					}
					else
					{

						doublereal time_current_piecewise_const = 0.0;
						// 0 - time (s); 1 - time duration(s);
						int din_piecewise_const = 0;
#ifdef MINGW_COMPILLER
						fscanf(fp_piecewise_const, "%d", &din_piecewise_const);
#else
						fscanf_s(fp_piecewise_const, "%d", &din_piecewise_const);
#endif
						

				        for (integer i_35 = 0; i_35 < glTSL.n_string_PiecewiseConst; i_35++) {
					
							// всего glTSL.n_string_PiecewiseConst троек значений.
							// time, timestep, m
							float fin_piecewise_const = 0.0;
							// time
#ifdef MINGW_COMPILLER
							fscanf(fp_piecewise_const, "%f", &fin_piecewise_const);
#else
							fscanf_s(fp_piecewise_const, "%f", &fin_piecewise_const);
#endif
							
							if (din_piecewise_const == 1) {
								// 1 - time duration(s);
								if (fin_piecewise_const <= 0.0) {
									printf("Zero increment current time. %e %lld\n", fin_piecewise_const, i_35);
									system("pause");
								}
								//printf("increment time=%e\n", fin_piecewise_const);
								time_current_piecewise_const += fin_piecewise_const;
							}
							else {
								// 0 - time(s);
								time_current_piecewise_const=fin_piecewise_const;
							}
							glTSL.table_law_piecewise_constant[i_35].time = time_current_piecewise_const;
								
							if ((i_35>0)&&(glTSL.table_law_piecewise_constant[i_35].time <= glTSL.table_law_piecewise_constant[i_35-1].time)) {
								// Внимание: при возникновении данной ошибки скорее всего произошла следующая ошибка интерфейса AliceMesh_v0_45.
								// Скорее всего забыл после временного интервала указать шаг по времени которым разбивается данный интервал. 
								// Третьим столбцом должен быть написан домножающий мощность множитель. Т.е. столбцов должно быть три а не два.
								printf("ERROR!!! Piecewise constant timestep law time non monotone. i==%lld n_i=%lld %e <=%e\n", i_35, glTSL.n_string_PiecewiseConst, glTSL.table_law_piecewise_constant[i_35].time, glTSL.table_law_piecewise_constant[i_35-1].time);
								system("pause");
								exit(1);
							}
							// timestep
#ifdef MINGW_COMPILLER
							fscanf(fp_piecewise_const, "%f", &fin_piecewise_const);
#else
							fscanf_s(fp_piecewise_const, "%f", &fin_piecewise_const);
#endif
							
							glTSL.table_law_piecewise_constant[i_35].timestep = fin_piecewise_const;
							if (glTSL.table_law_piecewise_constant[i_35].timestep <= 0.0) {
								printf("ERROR!!! Piecewise constant timestep law timestep is zero or negative.\n");
								system("pause");
								exit(1);
							}
							// m - множитель на который домножается значение тепловой мощности на данном шаге по времени. 
#ifdef MINGW_COMPILLER
							fscanf(fp_piecewise_const, "%f", &fin_piecewise_const);
#else
							fscanf_s(fp_piecewise_const, "%f", &fin_piecewise_const);
#endif
							
							glTSL.table_law_piecewise_constant[i_35].m = fin_piecewise_const;
						}
					}				
				}
			}

			// Newton-Richman condition.
			if (imakesource("adiabatic_vs_heat_transfer_coeff", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) 
					|| (idin == 2) || (idin == 3)) {
					// 0 - adiabatic wall, 1 - Newton Richman condition,
					// 2 - Stefan Bolcman condition, 3 - mix condition.
					switch (idin) {
						case 0: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
						case 1: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC; break;
						case 2: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC; break;
						case 3: adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC; break;
						default :adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; break;
					}
					//adiabatic_vs_heat_transfer_coeff = static_cast<integer>(idin);
					//printf("adiabatic_vs_heat_transfer_coeff =%lld\n", idin);
					//system("pause");
				}
				else {
					printf("Error!!! adiabatic_vs_heat_transfer_coeff must be equal 0==ADIABATIC_WALL_BC or 1==NEWTON_RICHMAN_BC or 2==STEFAN_BOLCMAN_BC or 3==MIX_CONDITION_BC. Current value %d\n",idin);
					system("pause");
					adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; // adiabatic wall
				}
			}
			else {
				printf("WARNING!!! adiabatic_vs_heat_transfer_coeff not found in file premeshin.txt\n");
				adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC; //  0 - adiabatic wall.
				printf("adiabatic_vs_heat_transfer_coeff =%d\n", adiabatic_vs_heat_transfer_coeff);
				if (bSTOP_Reading) system("pause");
			}			
			
			if (fmakesource("filmcoefficient", fin)) {
				// Найдено успешно.
				film_coefficient = static_cast<doublereal>(fin);// Коэффициент теплоотдачи.
			}
			else {
				printf("WARNING!!! filmcoefficient not found in file premeshin.txt\n");
				film_coefficient = 0.0; // Коэффициент теплоотдачи.
				printf("film_coefficient =%e\n", film_coefficient);
				if (bSTOP_Reading) system("pause");
			}

			// AЛИС сетка
			if (imakesource("itypeMesh", idin)) {
				// Найдено успешно.
				if (idin == 0) {
					// Обычная структурированная сетка.
					b_on_adaptive_local_refinement_mesh = false;
				}
				else {
					// АЛИС
					b_on_adaptive_local_refinement_mesh = true;
				}
				//printf("b_on_adaptive_local_refinement_mesh =%d\n", idin);
			}
			else {
				printf("WARNING!!! itypeMesh not found in file premeshin.txt\n");
				// Обычная структурированная сетка.
				b_on_adaptive_local_refinement_mesh = false; 
				printf("b_on_adaptive_local_refinement_mesh = Conformal\n");
				if (bSTOP_Reading) system("pause");
			}


			if (imakesource("version_ALICE_Mesh", idin)) {
				// Найдено успешно.
				switch (idin) {
					case 0: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
					case 1: itype_ALICE_Mesh= TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH; break;
					default: itype_ALICE_Mesh= TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; break;
				}
				//itype_ALICE_Mesh = static_cast<integer>(idin);
				//printf("itype_ALICE_Mesh =%lld\n", itype_ALICE_Mesh);
			}
			else {
				printf("WARNING!!! version_ALICE_Mesh not found in file premeshin.txt\n");
				itype_ALICE_Mesh = TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH; // Coarse Mesh.
				printf("itype_ALICE_Mesh =ONE_PASS_COARSE_ALICE_MESH\n");
				if (bSTOP_Reading) system("pause");
			}

			
			if (imakesource("m_restart_gmres", idin)) {
				// Найдено успешно.
				my_amg_manager.m_restart = static_cast<int>(idin);
				//printf("my_amg_manager.m_restart =%d\n", my_amg_manager.m_restart);
			}
			else {
				printf("WARNING!!! m_restart_gmres not found in file premeshin.txt\n");
				my_amg_manager.m_restart = 20; // Количество рестартов агоритма Flexible GMRES.
				printf("my_amg_manager.m_restart =%d\n", my_amg_manager.m_restart);
				if (bSTOP_Reading) system("pause");
			}


			// classical algebraic multigrid parameters:
			// only for my_agregat_amg.cu.
			if (imakesource("amg_manager_sorting_alg", idin)) {
				// Найдено успешно.
				switch (idin) {
				case 0: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
					break;
				case 1: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::QUICK_SORT;
					break;
				case 2: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::HEAP_SORT;
					break;
				case 3: my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::TIM_SORT;
					break;
				default:
					printf("ERROR!!! amg_manager_sorting_alg out of diapazon <0 && >3\n");
					system("PAUSE");
					my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT;
					break;
				}
				//printf("my_amg_manager.imySortAlgorithm =%lld\n", my_amg_manager.imySortAlgorithm);
			}
			else {
				printf("WARNING!!! amg_manager_sorting_alg not found in file premeshin.txt\n");
				my_amg_manager.imySortAlgorithm = MY_SORT_ALGORITHM::COUNTING_SORT; // Counting Sort Х.Г. Сьювард 1954г.
				printf("my_amg_manager.imySortAlgorithm =MY_SORT_ALGORITHM::COUNTING_SORT\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("inum_reduced_levels", idin)) {
				// Найдено успешно.
				my_amg_manager.maximum_delete_levels_Temperature = static_cast<integer>(idin);
				//printf("my_amg_manager.maximum_delete_levels_Temperature =%lld\n", my_amg_manager.maximum_delete_levels_Temperature);
			}
			else {
				printf("WARNING!!! inum_reduced_levels not found in file premeshin.txt\n");
				my_amg_manager.maximum_delete_levels_Temperature = 0; // нет удаляемых уровней.
				printf("my_amg_manager.maximum_delete_levels_Temperature =%lld\n", my_amg_manager.maximum_delete_levels_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("inum_reduced_levels_Speed", idin)) {
				// Найдено успешно.
				my_amg_manager.maximum_delete_levels_Speed  = static_cast<integer>(idin);
				//printf("my_amg_manager.maximum_delete_levels_Speed =%lld\n", my_amg_manager.maximum_delete_levels_Speed);
			}
			else {
				printf("WARNING!!! inum_reduced_levels_Speed not found in file premeshin.txt\n");
				my_amg_manager.maximum_delete_levels_Speed = 0; // нет удаляемых уровней.
				printf("my_amg_manager.maximum_delete_levels_Speed =%lld\n", my_amg_manager.maximum_delete_levels_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("inum_reduced_levels_Pressure", idin)) {
				// Найдено успешно.
				my_amg_manager.maximum_delete_levels_Pressure  = static_cast<integer>(idin);
				//printf("my_amg_manager.maximum_delete_levels_Pressure =%lld\n", my_amg_manager.maximum_delete_levels_Pressure);
			}
			else {
				printf("WARNING!!! inum_reduced_levels_Pressure not found in file premeshin.txt\n");
				my_amg_manager.maximum_delete_levels_Pressure = 0; // нет удаляемых уровней.
				printf("my_amg_manager.maximum_delete_levels_Pressure =%lld\n", my_amg_manager.maximum_delete_levels_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("inum_reduced_levels_Stress", idin)) {
				// Найдено успешно.
				my_amg_manager.maximum_delete_levels_Stress = static_cast<integer>(idin);
				//printf("my_amg_manager.maximum_delete_levels_Stress =%lld\n", my_amg_manager.maximum_delete_levels_Stress);
			}
			else {
				printf("WARNING!!! inum_reduced_levels_Stress not found in file premeshin.txt\n");
				my_amg_manager.maximum_delete_levels_Stress = 0; // нет удаляемых уровней.
				printf("my_amg_manager.maximum_delete_levels_Stress =%lld\n", my_amg_manager.maximum_delete_levels_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// type interpolation procedure:			
			if (imakesource("interpolation", idin)) {
				// Найдено успешно.
				my_amg_manager.number_interpolation_procedure_Temperature = static_cast<integer>(idin);
				//printf("my_amg_manager.number_interpolation_procedure_Temperature =%lld\n", my_amg_manager.number_interpolation_procedure_Temperature);
			}
			else {
				printf("WARNING!!! interpolation not found in file premeshin.txt\n");
				my_amg_manager.number_interpolation_procedure_Temperature = 0; // номер интерполяционной процедуры.
				printf("my_amg_manager.number_interpolation_procedure_Temperature =%lld\n", my_amg_manager.number_interpolation_procedure_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("interpolationSpeed", idin)) {
				// Найдено успешно.
				my_amg_manager.number_interpolation_procedure_Speed = static_cast<integer>(idin);
				//printf("my_amg_manager.number_interpolation_procedure_Speed =%lld\n", my_amg_manager.number_interpolation_procedure_Speed);
			}
			else {
				printf("WARNING!!! interpolationSpeed not found in file premeshin.txt\n");
				my_amg_manager.number_interpolation_procedure_Speed = 0; // номер интерполяционной процедуры.
				printf("my_amg_manager.number_interpolation_procedure_Speed =%lld\n", my_amg_manager.number_interpolation_procedure_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("interpolationPressure", idin)) {
				// Найдено успешно.
				my_amg_manager.number_interpolation_procedure_Pressure = static_cast<integer>(idin);
				//printf("my_amg_manager.number_interpolation_procedure_Pressure =%lld\n", my_amg_manager.number_interpolation_procedure_Pressure);
			}
			else {
				printf("WARNING!!! interpolationPressure not found in file premeshin.txt\n");
				my_amg_manager.number_interpolation_procedure_Pressure = 0; // номер интерполяционной процедуры.
				printf("my_amg_manager.number_interpolation_procedure_Pressure =%lld\n", my_amg_manager.number_interpolation_procedure_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("interpollationStress", idin)) {
				// Найдено успешно.
				my_amg_manager.number_interpolation_procedure_Stress  = static_cast<integer>(idin);
				//printf("my_amg_manager.number_interpolation_procedure_Stress =%lld\n", my_amg_manager.number_interpolation_procedure_Stress);
			}
			else {
				printf("WARNING!!! interpollationStress not found in file premeshin.txt\n");
				my_amg_manager.number_interpolation_procedure_Stress = 0; // номер интерполяционной процедуры.
				printf("my_amg_manager.number_interpolation_procedure_Stress =%lld\n", my_amg_manager.number_interpolation_procedure_Stress);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("CFalgorithmandDataStruct_Temperature", idin)) {
				// Найдено успешно.
				
				switch (idin) {
				case 0: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
					break;
				case 1: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				case 2: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
					break;
				case 3: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
					break;
				case 4: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
					break;
				case 5: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
					break;
				case 6: my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
					break;
				default:
					my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				}
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Temperature =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Temperature);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Temperature not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Temperature = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP; //  3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Temperature =RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("CFalgorithmandDataStruct_Speed", idin)) {
				// Найдено успешно.
				
				switch (idin) {
				case 0: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
					break;
				case 1: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				case 2: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
					break;
				case 3: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
					break;
				case 4: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
					break;
				case 5: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
					break;
				case 6: my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
					break;
				default:
					my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				}
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Speed =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Speed);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Speed not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Speed = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP; //  3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Speed =RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("CFalgorithmandDataStruct_Pressure", idin)) {
				// Найдено успешно.
				
				switch (idin) {
				case 0: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
					break;
				case 1: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				case 2: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
					break;
				case 3: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
					break;
				case 4: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
					break;
				case 5: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
					break;
				case 6: my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
					break;
				default:
					my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				}

				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Pressure =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Pressure);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Pressure not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Pressure = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP; //  3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Pressure =RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP\n");
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("CFalgorithmandDataStruct_Stress", idin)) {
				// Найдено успешно.
				
				switch (idin) {
				case 0: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE;
					break;
				case 1: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				case 2: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP;
					break;
				case 3: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP;
					break;
				case 4: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE;
					break;
				case 5: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP;
					break;
				case 6: my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE;
					break;
				default:
					my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE;
					break;
				}
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Stress =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Stress);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Stress not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Stress = RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP; // 3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Stress =RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP\n");
				if (bSTOP_Reading) system("pause");
			}

			// Устарело:
			//my_amg_manager.itypemodifyinterpol = din;
			//my_amg_manager.baglomeration_with_consistency_scaling = din;
			
			if (imakesource("modinterpol", idin)) {
				// Найдено успешно.
				my_amg_manager.bdiagonal_dominant = static_cast<integer>(idin);
				//printf("my_amg_manager.bdiagonal_dominant =%lld\n", my_amg_manager.bdiagonal_dominant);
			}
			else {
				printf("WARNING!!! modinterpol not found in file premeshin.txt\n");
				my_amg_manager.bdiagonal_dominant = 1; // diagonal dominant.
				printf("my_amg_manager.bdiagonal_dominant =%lld\n", my_amg_manager.bdiagonal_dominant);
				if (bSTOP_Reading) system("pause");
			}

			// 12,01,2018 заглушка, данные параметры более не используются в коде солвера.
			//fscanf_s(fp, "%lld", &din);
			//my_amg_manager.inumberadaptpass = din;


			// 23.02.2018
			// print matrix portrait
			if (imakesource("TemperatureMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bTemperatureMatrixPortrait = false;
					}
					else {
						my_amg_manager.bTemperatureMatrixPortrait = true;
					}
					//printf("my_amg_manager.bTemperatureMatrixPortrait =%lld\n", my_amg_manager.bTemperatureMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bTemperatureMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bTemperatureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
				}
			}
			else {
				printf("WARNING!!! TemperatureMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bTemperatureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
				printf("my_amg_manager.bTemperatureMatrixPortrait =NO_PRINT\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("SpeedMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bSpeedMatrixPortrait = false;
					}
					else {
						my_amg_manager.bSpeedMatrixPortrait = true;
					}
					//printf("my_amg_manager.bSpeedMatrixPortrait =%lld\n", my_amg_manager.bSpeedMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bSpeedMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bSpeedMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
				}
			}
			else {
				printf("WARNING!!! SpeedMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bSpeedMatrixPortrait = false; //  false - NO_PRINT, true - PRINT.
				printf("my_amg_manager.bSpeedMatrixPortrait =NO_PRINT\n");
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("PressureMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bPressureMatrixPortrait = false;
					}
					else {
						my_amg_manager.bPressureMatrixPortrait = true;
					}
					//printf("my_amg_manager.bPressureMatrixPortrait =%lld\n", my_amg_manager.bPressureMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bPressureMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bPressureMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
				}
			}
			else {
				printf("WARNING!!! PressureMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bPressureMatrixPortrait = false; //  false - NO_PRINT, true - PRINT.
				printf("my_amg_manager.bPressureMatrixPortrait =NO_PRINT\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("StressMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bStressMatrixPortrait = false;
					}
					else {
						my_amg_manager.bStressMatrixPortrait = true;
					}
					//printf("my_amg_manager.bStressMatrixPortrait =%lld\n", my_amg_manager.bStressMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bStressMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bStressMatrixPortrait = false; // false - NO_PRINT, true - PRINT.
				}
			}
			else {
				printf("WARNING!!! StressMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bStressMatrixPortrait = false; //  false - NO_PRINT, true - PRINT.
				printf("my_amg_manager.bStressMatrixPortrait =NO_PRINT\n");
				if (bSTOP_Reading) system("pause");
			}

			// truncation of interpolation:
			if (imakesource("truncationT", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.itruncation_interpolation_Temperature = static_cast<integer>(idin);
					//printf("my_amg_manager.itruncation_interpolation_Temperature =%lld\n", my_amg_manager.itruncation_interpolation_Temperature);
				}
				else {
					printf("my_amg_manager.itruncation_interpolation_Temperature must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.itruncation_interpolation_Temperature = 0; // 0 - OFF, 1 - ON.
				}
			}
			else {
				printf("WARNING!!! truncationT not found in file premeshin.txt\n");
				my_amg_manager.itruncation_interpolation_Temperature = 0; // 0 - OFF, 1 - ON.
				printf("my_amg_manager.itruncation_interpolation_Temperature =%lld\n", my_amg_manager.itruncation_interpolation_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("truncationSpeed", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.itruncation_interpolation_Speed = static_cast<integer>(idin);
					//printf("my_amg_manager.itruncation_interpolation_Speed =%lld\n", my_amg_manager.itruncation_interpolation_Speed);
				}
				else {
					printf("my_amg_manager.itruncation_interpolation_Speed must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.itruncation_interpolation_Speed = 0; // 0 - OFF, 1 - ON.
				}
			}
			else {
				printf("WARNING!!! truncationSpeed not found in file premeshin.txt\n");
				my_amg_manager.itruncation_interpolation_Speed = 0; // 0 - OFF, 1 - ON.
				printf("my_amg_manager.itruncation_interpolation_Speed =%lld\n", my_amg_manager.itruncation_interpolation_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("truncationPressure", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.itruncation_interpolation_Pressure = static_cast<integer>(idin);
					//printf("my_amg_manager.itruncation_interpolation_Pressure =%lld\n", my_amg_manager.itruncation_interpolation_Pressure);
				}
				else {
					printf("my_amg_manager.itruncation_interpolation_Pressure must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.itruncation_interpolation_Pressure = 0; // 0 - OFF, 1 - ON.
				}
			}
			else {
				printf("WARNING!!! truncationPressure not found in file premeshin.txt\n");
				my_amg_manager.itruncation_interpolation_Pressure = 0; // 0 - OFF, 1 - ON.
				printf("my_amg_manager.itruncation_interpolation_Pressure =%lld\n", my_amg_manager.itruncation_interpolation_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("truncationStress", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.itruncation_interpolation_Stress = static_cast<integer>(idin);
					//printf("my_amg_manager.itruncation_interpolation_Stress =%lld\n", my_amg_manager.itruncation_interpolation_Stress);
				}
				else {
					printf("my_amg_manager.itruncation_interpolation_Stress must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.itruncation_interpolation_Stress = 0; // 0 - OFF, 1 - ON.
				}
			}
			else {
				printf("WARNING!!! truncationStress not found in file premeshin.txt\n");
				my_amg_manager.itruncation_interpolation_Stress = 0; // 0 - OFF, 1 - ON.
				printf("my_amg_manager.itruncation_interpolation_Stress =%lld\n", my_amg_manager.itruncation_interpolation_Stress);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("truncation_Temperature_drop_tolerance", fin)) {
				// Найдено успешно.
				if ((fin >= 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Temperature = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.truncation_interpolation_Temperature must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.truncation_interpolation_Temperature = 0.2; //  default 0.2 - truncation of interpolation.
				}
			}
			else {
				printf("WARNING!!! truncation_Temperature_drop_tolerance not found in file premeshin.txt\n");
				my_amg_manager.truncation_interpolation_Temperature = 0.2; // default 0.2 - truncation of interpolation.
				printf("my_amg_manager.truncation_interpolation_Temperature =%e\n", my_amg_manager.truncation_interpolation_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("truncation_Speed_drop_tolerance", fin)) {
				// Найдено успешно.
				if ((fin >= 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Speed = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.truncation_interpolation_Speed must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.truncation_interpolation_Speed = 0.2; //  default 0.2 - truncation of interpolation.
				}
			}
			else {
				printf("WARNING!!! truncation_Speed_drop_tolerance not found in file premeshin.txt\n");
				my_amg_manager.truncation_interpolation_Speed = 0.2; // default 0.2 - truncation of interpolation.
				printf("my_amg_manager.truncation_interpolation_Speed =%e\n", my_amg_manager.truncation_interpolation_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("truncation_Pressure_drop_tolerance", fin)) {
				// Найдено успешно.
				if ((fin >= 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Pressure = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.truncation_interpolation_Pressure must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.truncation_interpolation_Pressure = 0.2; //  default 0.2 - truncation of interpolation.
				}
			}
			else {
				printf("WARNING!!! truncation_Pressure_drop_tolerance not found in file premeshin.txt\n");
				my_amg_manager.truncation_interpolation_Pressure = 0.2; // default 0.2 - truncation of interpolation.
				printf("my_amg_manager.truncation_interpolation_Pressure  =%e\n", my_amg_manager.truncation_interpolation_Pressure);
				if (bSTOP_Reading) system("pause");
			}
			
			if (fmakesource("truncation_Stress_drop_tolerance", fin)) {
				// Найдено успешно.
				if ((fin >= 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Stress = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.truncation_interpolation_Stress must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.truncation_interpolation_Stress = 0.2; //  default 0.2 - truncation of interpolation.
				}
			}
			else {
				printf("WARNING!!! truncation_Stress_drop_tolerance not found in file premeshin.txt\n");
				my_amg_manager.truncation_interpolation_Stress = 0.2; // default 0.2 - truncation of interpolation.
				printf("my_amg_manager.truncation_interpolation_Stress =%e\n", my_amg_manager.truncation_interpolation_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// number nFinnest sweeps:
			if (imakesource("nFinnest", idin)) {
				// Найдено успешно.
				my_amg_manager.nFinnest_Temperature = static_cast<integer>(idin);
				//printf("my_amg_manager.nFinnest_Temperature =%lld\n", my_amg_manager.nFinnest_Temperature);
			}
			else {
				printf("WARNING!!! nFinnest not found in file premeshin.txt\n");
				my_amg_manager.nFinnest_Temperature = 2; // Количество сглаживаний на самой подробной сетке.
				printf("my_amg_manager.nFinnest_Temperature =%lld\n", my_amg_manager.nFinnest_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("nFinnestSpeed", idin)) {
				// Найдено успешно.
				my_amg_manager.nFinnest_Speed = static_cast<integer>(idin);
				//printf("my_amg_manager.nFinnest_Speed =%lld\n", my_amg_manager.nFinnest_Speed);
			}
			else {
				printf("WARNING!!! nFinnestSpeed not found in file premeshin.txt\n");
				my_amg_manager.nFinnest_Speed = 2; // Количество сглаживаний на самой подробной сетке.
				printf("my_amg_manager.nFinnest_Speed =%lld\n", my_amg_manager.nFinnest_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("nFinnestPressure", idin)) {
				// Найдено успешно.
				my_amg_manager.nFinnest_Pressure = static_cast<integer>(idin);
				//printf("my_amg_manager.nFinnest_Pressure =%lld\n", my_amg_manager.nFinnest_Pressure);
			}
			else {
				printf("WARNING!!! nFinnestPressure not found in file premeshin.txt\n");
				my_amg_manager.nFinnest_Pressure = 2; // Количество сглаживаний на самой подробной сетке.
				printf("my_amg_manager.nFinnest_Pressure =%lld\n", my_amg_manager.nFinnest_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("nFinnestStress", idin)) {
				// Найдено успешно.
				my_amg_manager.nFinnest_Stress = static_cast<integer>(idin);
				//printf("my_amg_manager.nFinnest_Stress =%lld\n", my_amg_manager.nFinnest_Stress);
			}
			else {
				printf("WARNING!!! nFinnestStress not found in file premeshin.txt\n");
				my_amg_manager.nFinnest_Stress = 2; // Количество сглаживаний на самой подробной сетке.
				printf("my_amg_manager.nFinnest_Stress =%lld\n", my_amg_manager.nFinnest_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// number presweeps:
			if (imakesource("numberpresmothers", idin)) {
				// Найдено успешно.
				my_amg_manager.nu1_Temperature = static_cast<integer>(idin);
				//printf("my_amg_manager.nu1_Temperature =%lld\n", my_amg_manager.nu1_Temperature);
			}
			else {
				printf("WARNING!!! numberpresmothers not found in file premeshin.txt\n");
				my_amg_manager.nu1_Temperature = 1; // Количество предварительных сглаживаний.
				printf("my_amg_manager.nu1_Temperature =%lld\n", my_amg_manager.nu1_Temperature);
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("numberpresmothersSpeed", idin)) {
				// Найдено успешно.
				my_amg_manager.nu1_Speed = static_cast<integer>(idin);
				//printf("my_amg_manager.nu1_Speed =%lld\n", my_amg_manager.nu1_Speed);
			}
			else {
				printf("WARNING!!! numberpresmothersSpeed not found in file premeshin.txt\n");
				my_amg_manager.nu1_Speed = 1; // Количество предварительных сглаживаний.
				printf("my_amg_manager.nu1_Speed =%lld\n", my_amg_manager.nu1_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("numberpresmothersPressure", idin)) {
				// Найдено успешно.
				my_amg_manager.nu1_Pressure = static_cast<integer>(idin);
				//printf("my_amg_manager.nu1_Pressure =%lld\n", my_amg_manager.nu1_Pressure);
			}
			else {
				printf("WARNING!!! numberpresmothersPressure not found in file premeshin.txt\n");
				my_amg_manager.nu1_Pressure = 1; // Количество предварительных сглаживаний.
				printf("my_amg_manager.nu1_Pressure =%lld\n", my_amg_manager.nu1_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("numberpresmoothersStress", idin)) {
				// Найдено успешно.
				my_amg_manager.nu1_Stress = static_cast<integer>(idin);
				//printf("my_amg_manager.nu1_Stress =%lld\n", my_amg_manager.nu1_Stress);
			}
			else {
				printf("WARNING!!! numberpresmoothersStress not found in file premeshin.txt\n");
				my_amg_manager.nu1_Stress = 1; // Количество предварительных сглаживаний.
				printf("my_amg_manager.nu1_Stress =%lld\n", my_amg_manager.nu1_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// number postsweeps:
			if (imakesource("numberpostsweeps", idin)) {
				// Найдено успешно.
				my_amg_manager.nu2_Temperature = static_cast<integer>(idin);
				//printf("my_amg_manager.nu2_Temperature =%lld\n", my_amg_manager.nu2_Temperature);
			}
			else {
				printf("WARNING!!! numberpostsweeps not found in file premeshin.txt\n");
				my_amg_manager.nu2_Temperature = 2; // Количество пост сглаживаний.
				printf("my_amg_manager.nu2_Temperature =%lld\n", my_amg_manager.nu2_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("numberpostsweepsSpeed", idin)) {
				// Найдено успешно.
				my_amg_manager.nu2_Speed = static_cast<integer>(idin);
				//printf("my_amg_manager.nu2_Speed =%lld\n", my_amg_manager.nu2_Speed);
			}
			else {
				printf("WARNING!!! numberpostsweepsSpeed not found in file premeshin.txt\n");
				my_amg_manager.nu2_Speed = 2; // Количество пост сглаживаний.
				printf("my_amg_manager.nu2_Speed =%lld\n", my_amg_manager.nu2_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("numberpostsweepsPressure", idin)) {
				// Найдено успешно.
				my_amg_manager.nu2_Pressure = static_cast<integer>(idin);
				//printf("my_amg_manager.nu2_Pressure =%lld\n", my_amg_manager.nu2_Pressure);
			}
			else {
				printf("WARNING!!! numberpostsweepsPressure not found in file premeshin.txt\n");
				my_amg_manager.nu2_Pressure = 2; // Количество пост сглаживаний.
				printf("my_amg_manager.nu2_Pressure =%lld\n", my_amg_manager.nu2_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("numberpostsweepsStress", idin)) {
				// Найдено успешно.
				my_amg_manager.nu2_Stress = static_cast<integer>(idin);
				//printf("my_amg_manager.nu2_Stress =%lld\n", my_amg_manager.nu2_Stress);
			}
			else {
				printf("WARNING!!! numberpostsweepsStress not found in file premeshin.txt\n");
				my_amg_manager.nu2_Stress = 2; // Количество пост сглаживаний.
				printf("my_amg_manager.nu2_Stress =%lld\n", my_amg_manager.nu2_Stress);
				if (bSTOP_Reading) system("pause");
			}


			// memory size:
		    if (imakesource("memorysize", idin)) {
				// Найдено успешно.
				my_amg_manager.memory_size_Temperature = static_cast<integer>(idin);
				//printf("my_amg_manager.memory_size_Temperature =%lld\n", my_amg_manager.memory_size_Temperature);
			}
			else {
				printf("WARNING!!! memorysize not found in file premeshin.txt\n");
				my_amg_manager.memory_size_Temperature = 13; //Количество потребляемой памяти в размерах матрицы А.
				printf("my_amg_manager.memory_size_Temperature =%lld\n", my_amg_manager.memory_size_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("memorysizeSpeed", idin)) {
				// Найдено успешно.
				my_amg_manager.memory_size_Speed = static_cast<integer>(idin);
				//printf("my_amg_manager.memory_size_Speed =%lld\n", my_amg_manager.memory_size_Speed);
			}
			else {
				printf("WARNING!!! memorysizeSpeed not found in file premeshin.txt\n");
				my_amg_manager.memory_size_Speed = 13; // Количество потребляемой памяти в размерах матрицы А.
				printf("my_amg_manager.memory_size_Speed =%lld\n", my_amg_manager.memory_size_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("memorysizePressure", idin)) {
				// Найдено успешно.
				my_amg_manager.memory_size_Pressure = static_cast<integer>(idin);
				//printf("my_amg_manager.memory_size_Pressure =%lld\n", my_amg_manager.memory_size_Pressure);
			}
			else {
				printf("WARNING!!! memorysizePressure not found in file premeshin.txt\n");
				my_amg_manager.memory_size_Pressure = 13; // Количество потребляемой памяти в размерах матрицы А.
				printf("my_amg_manager.memory_size_Pressure =%lld\n", my_amg_manager.memory_size_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("memorysizeStress", idin)) {
				// Найдено успешно.
				my_amg_manager.memory_size_Stress = static_cast<integer>(idin);
				//printf("my_amg_manager.memory_size_Stress =%lld\n", my_amg_manager.memory_size_Stress);
			}
			else {
				printf("WARNING!!! memorysizeStress not found in file premeshin.txt\n");
				my_amg_manager.memory_size_Stress = 22; // Количество потребляемой памяти в размерах матрицы А.
				printf("my_amg_manager.memory_size_Stress =%lld\n", my_amg_manager.memory_size_Stress);
				if (bSTOP_Reading) system("pause");
			}


			// Параметр верхней релаксации в сглаживателе
			// Код считывания полностью удален 24.09.2019.
			

			// использовать ли ilu2 smoother.
			if (imakesource("smoothertypeTemperature", idin)) {
				// Найдено успешно.
				din = static_cast<integer>(idin);
				my_amg_manager.ilu2_smoother_Temperature = din;
				my_amg_manager.b_gmresTemp = false;
				my_amg_manager.b_ChebyshevSmootherTemp = false;
				my_amg_manager.b_spai0Temp = false;

				/*
				if (din == 3) {
					din = 0;
					my_amg_manager.ilu2_smoother_Temperature = din;
					my_amg_manager.b_gmresTemp = true;
				}
				*/

				//printf("my_amg_manager.ilu2_smoother_Temperature =%lld\n", my_amg_manager.ilu2_smoother_Temperature);
			}
			else {
				printf("WARNING!!! smoothertypeTemperature not found in file premeshin.txt\n");
				my_amg_manager.ilu2_smoother_Temperature = 4; // Сглаживатель Якоби.
				my_amg_manager.b_gmresTemp = false;
				my_amg_manager.b_ChebyshevSmootherTemp = false;
				my_amg_manager.b_spai0Temp = false;
				printf("my_amg_manager.ilu2_smoother_Temperature =%lld\n", my_amg_manager.ilu2_smoother_Temperature);
				if (bSTOP_Reading) system("pause");
			}
					   		
			if (imakesource("smoothertypeSpeed", idin)) {
				// Найдено успешно.
				din = static_cast<integer>(idin);
				my_amg_manager.ilu2_smoother_Speed = din;
				my_amg_manager.b_gmresSpeed = false;	
				my_amg_manager.b_ChebyshevSmootherSpeed = false;
				my_amg_manager.b_spai0Speed = false;
				/*
				if (din == 3) {
					din = 0;
					my_amg_manager.ilu2_smoother_Speed = din;
					my_amg_manager.b_gmresSpeed = true;
				}
				*/

				//printf("my_amg_manager.ilu2_smoother_Speed =%lld\n", my_amg_manager.ilu2_smoother_Speed);
			}
			else {
				printf("WARNING!!! smoothertypeSpeed not found in file premeshin.txt\n");
				my_amg_manager.ilu2_smoother_Speed = 0; // Сглаживатель Якоби.
				my_amg_manager.b_gmresSpeed = false;
				my_amg_manager.b_ChebyshevSmootherSpeed = false;
				my_amg_manager.b_spai0Speed = false;
				printf("my_amg_manager.ilu2_smoother_Speed =%lld\n", my_amg_manager.ilu2_smoother_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("smoothertypePressure", idin)) {
				// Найдено успешно.
				din = static_cast<integer>(idin);
				my_amg_manager.ilu2_smoother_Pressure = din;
				my_amg_manager.b_gmresPressure = false;		
				my_amg_manager.b_ChebyshevSmootherPressure = false;
				my_amg_manager.b_spai0Pressure = false;
				/*
				if (din == 3) {
					din = 0;
					my_amg_manager.ilu2_smoother_Pressure = din;
					my_amg_manager.b_gmresPressure = true;
				}
				*/

				//printf("my_amg_manager.ilu2_smoother_Pressure =%lld\n", my_amg_manager.ilu2_smoother_Pressure);
			}
			else {
				printf("WARNING!!! smoothertypePressure not found in file premeshin.txt\n");
				my_amg_manager.ilu2_smoother_Pressure = 0; // Сглаживатель Якоби.
				my_amg_manager.b_gmresPressure = false;
				my_amg_manager.b_ChebyshevSmootherPressure = false;
				my_amg_manager.b_spai0Pressure = false;
				printf("my_amg_manager.ilu2_smoother_Pressure =%lld\n", my_amg_manager.ilu2_smoother_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("smoothertypeStress", idin)) {
				// Найдено успешно.
				din = static_cast<integer>(idin);
				my_amg_manager.ilu2_smoother_Stress = din;
				my_amg_manager.b_gmresStress = false;	
				my_amg_manager.b_ChebyshevSmootherStress = false;
				my_amg_manager.b_spai0Stress = false;
				/*
				if (din == 3) {
					din = 0;
					my_amg_manager.ilu2_smoother_Stress = din;
					my_amg_manager.b_gmresStress = true;
				}
				*/
				//printf("my_amg_manager.ilu2_smoother_Stress =%lld\n", my_amg_manager.ilu2_smoother_Stress);
			}
			else {
				printf("WARNING!!! smoothertypeStress not found in file premeshin.txt\n");
				my_amg_manager.ilu2_smoother_Stress = 0; // Сглаживатель Якоби.
				my_amg_manager.b_gmresStress = false;
				my_amg_manager.b_ChebyshevSmootherStress = false;
				my_amg_manager.b_spai0Stress = false;
				printf("my_amg_manager.ilu2_smoother_Stress =%lld\n", my_amg_manager.ilu2_smoother_Stress);
				if (bSTOP_Reading) system("pause");
			}

			my_amg_manager.gold_const_Temperature = return_gold_const(my_amg_manager.ilu2_smoother_Temperature);
			my_amg_manager.gold_const_Speed = return_gold_const(my_amg_manager.ilu2_smoother_Speed);
			my_amg_manager.gold_const_Pressure = return_gold_const(my_amg_manager.ilu2_smoother_Pressure);
			my_amg_manager.gold_const_Stress = return_gold_const(my_amg_manager.ilu2_smoother_Stress);


			if (imakesource("Chebyshev_degree", idin)) {
				// Найдено успешно.
				din = static_cast<integer>(idin);
				my_amg_manager.Chebyshev_degree = din;
				
				if (!((din == 5) || (din == 8) || (din == 16) || (din == 25) || (din == 32) || (din == 64) || (din == 128) || (din == 256) || (din == 512))) {
					std::cout << "Chebyshev degree undefined. Chebyshev degree mast be equal: 5, 8, 16, 25, 32, 64, 128, 256, 512 in this programm. Chebyshev degree == " << din << std::endl;
					system("pause");
					my_amg_manager.Chebyshev_degree = 5;
				}

				//printf("my_amg_manager.Chebyshev_degree =%lld\n", my_amg_manager.Chebyshev_degree);
			}
			else {
				printf("WARNING!!!Chebyshev_degree not found in file premeshin.txt\n");
				my_amg_manager.Chebyshev_degree = 5;
				
				printf("my_amg_manager.Chebyshev_degree =%lld\n", my_amg_manager.Chebyshev_degree);
				if (bSTOP_Reading) system("pause");
			}

			// strength threshold:
			
			if (fmakesource("threshold", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.theta_Temperature = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.theta_Temperature must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.theta_Temperature = 0.23; //  threshold default 0.23.
				}
			}
			else {
				printf("WARNING!!! threshold not found in file premeshin.txt\n");
				my_amg_manager.theta_Temperature = 0.23; // threshold default 0.23.
				printf("my_amg_manager.theta_Temperature =%e\n", my_amg_manager.theta_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("thresholdSpeed", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.theta_Speed = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.theta_Speed must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.theta_Speed = 0.23; //  threshold default 0.23.
				}
			}
			else {
				printf("WARNING!!! thresholdSpeed not found in file premeshin.txt\n");
				my_amg_manager.theta_Speed = 0.23; // threshold default 0.23.
				printf("my_amg_manager.theta_Speed =%e\n", my_amg_manager.theta_Speed);
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("thresholdPressure", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.theta_Pressure = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.theta_Pressure must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.theta_Pressure = 0.23; //  threshold default 0.23.
				}
			}
			else {
				printf("WARNING!!! thresholdPressure not found in file premeshin.txt\n");
				my_amg_manager.theta_Pressure = 0.23; // threshold default 0.23.
				printf("my_amg_manager.theta_Pressure =%e\n", my_amg_manager.theta_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("thresholdStress", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.theta_Stress = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.theta_Stress must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.theta_Stress = 0.23; //  threshold default 0.23.
				}
			}
			else {
				printf("WARNING!!! thresholdStress not found in file premeshin.txt\n");
				my_amg_manager.theta_Stress = 0.23; // threshold default 0.23.
				printf("my_amg_manager.theta_Stress =%e\n", my_amg_manager.theta_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// magic threshold:
			// magic <=> F_to_F
			if (fmakesource("magicT", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.F_to_F_Temperature = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.F_to_F_Temperature must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.F_to_F_Temperature = 0.4; //  F_to_F default 0.4.
				}
			}
			else {
				printf("WARNING!!! magicT not found in file premeshin.txt\n");
				my_amg_manager.F_to_F_Temperature = 0.4; // F_to_F default 0.4.
				printf("my_amg_manager.F_to_F_Temperature =%e\n", my_amg_manager.F_to_F_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("magicSpeed", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.F_to_F_Speed = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.F_to_F_Speed must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.F_to_F_Speed = 0.4; //  F_to_F default 0.4.
				}
			}
			else {
				printf("WARNING!!! magicSpeed not found in file premeshin.txt\n");
				my_amg_manager.F_to_F_Speed = 0.4; // F_to_F default 0.4.
				printf("my_amg_manager.F_to_F_Speed =%e\n", my_amg_manager.F_to_F_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("magicPressure", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.F_to_F_Pressure = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.F_to_F_Pressure must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.F_to_F_Pressure = 0.4; //  F_to_F default 0.4.
				}
			}
			else {
				printf("WARNING!!! magicPressure not found in file premeshin.txt\n");
				my_amg_manager.F_to_F_Pressure = 0.4; // F_to_F default 0.4.
				printf("my_amg_manager.F_to_F_Pressure =%e\n", my_amg_manager.F_to_F_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("magicStress", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.F_to_F_Stress = static_cast<doublereal>(fin);
				}
				else {
					printf("my_amg_manager.F_to_F_Stress must be  <1 and >0. current_value=%e\n", fin);
					system("pause");
					my_amg_manager.F_to_F_Stress = 0.4; //  F_to_F default 0.4.
				}
			}
			else {
				printf("WARNING!!! magicStress not found in file premeshin.txt\n");
				my_amg_manager.F_to_F_Stress = 0.4; // F_to_F default 0.4.
				printf("my_amg_manager.F_to_F_Stress  =%e\n", my_amg_manager.F_to_F_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// AMG Splitting (coarsening)
			// Способ построения C/F разбиения: 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
			// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
			if (imakesource("coarseningTemp", idin)) {
				// Найдено успешно.			

				switch (idin) {
				case 0: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
					break;
				case 1: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
					break;
				case 2: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
					break;
				case 3: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
					break;
				case 4: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
					break;
				case 5: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
					break;
				case 6: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
					break;
				case 7: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
					break;
				case 8: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				case 9: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
					break;
				case 10: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
					break;
				case 11: my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
					break;
				default:
					my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				}

				//printf("my_amg_manager.icoarseningTemp =%lld\n", my_amg_manager.icoarseningTemp);
			}
			else {
				printf("WARNING!!! coarseningTemp not found in file premeshin.txt\n");
				my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION; // 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
				printf("my_amg_manager.icoarseningTemp = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("coarseningSpeed", idin)) {
				// Найдено успешно.
				switch (idin) {
				case 0: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
					break;
				case 1: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
					break;
				case 2: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
					break;
				case 3: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
					break;
				case 4: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
					break;
				case 5: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
					break;
				case 6: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
					break;
				case 7: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
					break;
				case 8: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				case 9: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
					break;
				case 10: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
					break;
				case 11: my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
					break;
				default:
					my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				}

				//printf("my_amg_manager.icoarseningSpeed =%lld\n", my_amg_manager.icoarseningSpeed);
			}
			else {
				printf("WARNING!!! coarseningSpeed not found in file premeshin.txt\n");
				my_amg_manager.icoarseningSpeed = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION; // 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
				printf("my_amg_manager.icoarseningSpeed =MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("coarseningPressure", idin)) {
				// Найдено успешно.
				switch (idin) {
				case 0: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
					break;
				case 1: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
					break;
				case 2: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
					break;
				case 3: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
					break;
				case 4: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
					break;
				case 5: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
					break;
				case 6: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
					break;
				case 7: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
					break;
				case 8: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				case 9: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
					break;
				case 10: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
					break;
				case 11: my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
					break;
				default:
					my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				}

				//printf("my_amg_manager.icoarseningPressure =%lld\n", my_amg_manager.icoarseningPressure);
			}
			else {
				printf("WARNING!!! coarseningPressure not found in file premeshin.txt\n");
				my_amg_manager.icoarseningPressure = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION; // 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
				printf("my_amg_manager.icoarseningPressure =MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("coarseningStress", idin)) {
				// Найдено успешно.
				switch (idin) {
				case 0: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
					break;
				case 1: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
					break;
				case 2: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
					break;
				case 3: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
					break;
				case 4: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION;
					break;
				case 5: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION;
					break;
				case 6: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION;
					break;
				case 7: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION;
					break;
				case 8: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				case 9: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS;
					break;
				case 10: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2;
					break;
				case 11: my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full;
					break;
				default:
					my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS;
					break;
				}
				//printf("my_amg_manager.icoarseningStress =%lld\n", my_amg_manager.icoarseningStress);
			}
			else {
				printf("WARNING!!! coarseningStress not found in file premeshin.txt\n");
				my_amg_manager.icoarseningStress = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION; // 0 - standart, 1 - RS 2, 2 - standart Strong Transpose, 3 - RS2 Strong Transpose.
				printf("my_amg_manager.icoarseningStress =MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION\n");
				if (bSTOP_Reading) system("pause");
			}

			// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
			// Если din==1, Stabilization BiCGStab.
			// 8.01.2017 Метод Хенка ван дер Ворста BiCGStab 1992
			// предобусловленный алгебраическим многосеточным методом.
			// 9.01.2018 Если din==2, FGMRes предобусловленный алгебраическим многосеточным методом.
			if (imakesource("StabilizationTemp", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) || (idin == 2) || (idin == 3)) {
					my_amg_manager.istabilizationTemp = static_cast<integer>(idin);
					//printf("my_amg_manager.istabilizationTemp =%lld\n", my_amg_manager.istabilizationTemp);
				}
				else {
					printf("my_amg_manager.istabilizationTemp must be equal 3, 2, 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.istabilizationTemp = 3; //  0 - none, 1 - BiCGStab, 2 - FGMRes, 3 - for Non Linear problem.
				}
			}
			else {
				printf("WARNING!!! StabilizationTemp not found in file premeshin.txt\n");
				my_amg_manager.istabilizationTemp = 0; // 0 - none
				printf("my_amg_manager.istabilizationTemp =%lld\n", my_amg_manager.istabilizationTemp);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("StabilizationSpeed", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) || (idin == 2)) {
					my_amg_manager.istabilizationSpeed = static_cast<integer>(idin);
					//printf("my_amg_manager.istabilizationSpeed =%lld\n", my_amg_manager.istabilizationSpeed);
				}
				else {
					printf("my_amg_manager.istabilizationSpeed must be equal 2, 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.istabilizationSpeed = 0; //  0 - none, 1 - BiCGStab, 2 - FGMRes.
				}
			}
			else {
				printf("WARNING!!! StabilizationSpeed not found in file premeshin.txt\n");
				my_amg_manager.istabilizationSpeed = 0; // 0 - none
				printf("my_amg_manager.istabilizationSpeed =%lld\n", my_amg_manager.istabilizationSpeed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("StabilizationPressure", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) || (idin == 2) || (idin == 3)) {
					my_amg_manager.istabilizationPressure = static_cast<integer>(idin);
					//printf("my_amg_manager.istabilizationPressure =%lld\n", my_amg_manager.istabilizationPressure);
				}
				else {
					printf("my_amg_manager.istabilizationPressure must be equal 3, 2, 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.istabilizationPressure = 0; //  0 - none, 1 - BiCGStab, 2 - FGMRes.
				}
			}
			else {
				printf("WARNING!!! StabilizationPressure not found in file premeshin.txt\n");
				my_amg_manager.istabilizationPressure = 0; // 0 - none
				printf("my_amg_manager.istabilizationPressure =%lld\n", my_amg_manager.istabilizationPressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("StabilizationStress", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) || (idin == 2) || (idin == 3)) {
					my_amg_manager.istabilizationStress = static_cast<integer>(idin);
					//printf("my_amg_manager.istabilizationStress =%lld\n", my_amg_manager.istabilizationStress);
				}
				else {
					printf("my_amg_manager.istabilizationStress must be equal 3, 2, 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.istabilizationStress = 0; //  0 - none, 1 - BiCGStab, 2 - FGMRes.
				}
			}
			else {
				printf("WARNING!!! StabilizationStress not found in file premeshin.txt\n");
				my_amg_manager.istabilizationStress = 0; // 0 - none, 1 - BiCGStab, 2 - FGMRes.
				printf("my_amg_manager.istabilizationStress =%lld\n", my_amg_manager.istabilizationStress);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("Patch_number", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 7)) {
					my_amg_manager.ipatch_number = static_cast<integer>(idin);
					//printf("my_amg_manager.ipatch_number =%lld\n", my_amg_manager.ipatch_number);
				}
				else {
					printf("my_amg_manager.ipatch_number must be equal 7 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.ipatch_number = 7; // 7 - Strong Transpose, 0 - патч не применяется.
				}
			}
			else {
				printf("WARNING!!! Patch_number not found in file premeshin.txt\n");
				my_amg_manager.ipatch_number = 7; // 7 -Strong Transpose, 0 - патч не применяется.
				printf("my_amg_manager.ipatch_number =%lld\n", my_amg_manager.ipatch_number);
				if (bSTOP_Reading) system("pause");
			}

			// Печать лога на консоль.
		    if (imakesource("printlogTemperature", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.iprint_log_Temperature = static_cast<integer>(idin);
					//printf("my_amg_manager.iprint_log_Temperature =%lld\n", my_amg_manager.iprint_log_Temperature);
				}
				else {
					printf("my_amg_manager.iprint_log_Temperature must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.iprint_log_Temperature = 0; // 1 - печатать лог на консоль.
				}
			}
			else {
				printf("WARNING!!! printlogTemperature not found in file premeshin.txt\n");
				my_amg_manager.iprint_log_Temperature = 0; // 1 - печатать лог на консоль.
				printf("my_amg_manager.iprint_log_Temperature =%lld\n", my_amg_manager.iprint_log_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("printlogSpeed", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.iprint_log_Speed = static_cast<integer>(idin);
					//printf("my_amg_manager.iprint_log_Speed =%lld\n", my_amg_manager.iprint_log_Speed);
				}
				else {
					printf("my_amg_manager.iprint_log_Speed must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.iprint_log_Speed = 0; // 1 - печатать лог на консоль.
				}
			}
			else {
				printf("WARNING!!! printlogSpeed not found in file premeshin.txt\n");
				my_amg_manager.iprint_log_Speed = 0; // 1 - печатать лог на консоль.
				printf("my_amg_manager.iprint_log_Speed =%lld\n", my_amg_manager.iprint_log_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("printlogPressure", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.iprint_log_Pressure = static_cast<integer>(idin);
					//printf("my_amg_manager.iprint_log_Pressure =%lld\n", my_amg_manager.iprint_log_Pressure);
				}
				else {
					printf("my_amg_manager.iprint_log_Pressure must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.iprint_log_Pressure = 0; // 1 - печатать лог на консоль.
				}
			}
			else {
				printf("WARNING!!! printlogPressure not found in file premeshin.txt\n");
				my_amg_manager.iprint_log_Pressure = 0; // 1 - печатать лог на консоль.
				printf("my_amg_manager.iprint_log_Pressure =%lld\n", my_amg_manager.iprint_log_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("printlogStress", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.iprint_log_Stress = static_cast<integer>(idin);
					//printf("my_amg_manager.iprint_log_Stress =%lld\n", my_amg_manager.iprint_log_Stress);
				}
				else {
					printf("my_amg_manager.iprint_log_Stress must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.iprint_log_Stress = 0; // 1 - печатать лог на консоль.
				}
			}
			else {
				printf("WARNING!!! printlogStress not found in file premeshin.txt\n");
				my_amg_manager.iprint_log_Stress = 0; // 1 - печатать лог на консоль.
				printf("my_amg_manager.iprint_log_Stress =%lld\n", my_amg_manager.iprint_log_Stress);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bthreshold_Temperature_auto", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bthreshold_Temperature_auto = false;
					}
					else {
						my_amg_manager.bthreshold_Temperature_auto = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bthreshold_Temperature_auto = false\n");
					// } else {
					// printf("my_amg_manager.bthreshold_Temperature_auto = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bthreshold_Temperature_auto must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
					my_amg_manager.bthreshold_Temperature_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				}
			}
			else {
				printf("WARNING!!! bthreshold_Temperature_auto not found in file premeshin.txt\n");
				// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
				my_amg_manager.bthreshold_Temperature_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				printf("my_amg_manager.bthreshold_Temperature_auto =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bthreshold_Speed_auto", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bthreshold_Speed_auto = false;
					}
					else {
						my_amg_manager.bthreshold_Speed_auto = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bthreshold_Speed_auto = false\n");
					// } else {
					// printf("my_amg_manager.bthreshold_Speed_auto = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bthreshold_Speed_auto must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
					my_amg_manager.bthreshold_Speed_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				}
			}
			else {
				printf("WARNING!!! bthreshold_Speed_auto not found in file premeshin.txt\n");
				// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
				my_amg_manager.bthreshold_Speed_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				printf("my_amg_manager.bthreshold_Speed_auto =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bthreshold_Pressure_auto", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bthreshold_Pressure_auto = false;
					}
					else {
						my_amg_manager.bthreshold_Pressure_auto = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bthreshold_Pressure_auto = false\n");
					// } else {
					// printf("my_amg_manager.bthreshold_Pressure_auto = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bthreshold_Pressure_auto must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
					my_amg_manager.bthreshold_Pressure_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				}
			}
			else {
				printf("WARNING!!! bthreshold_Pressure_auto not found in file premeshin.txt\n");
				// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
				my_amg_manager.bthreshold_Pressure_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				printf("my_amg_manager.bthreshold_Pressure_auto =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bthreshold_Stress_auto", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bthreshold_Stress_auto = false;
					}
					else {
						my_amg_manager.bthreshold_Stress_auto = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bthreshold_Stress_auto = false\n");
					// } else {
					// printf("my_amg_manager.bthreshold_Stress_auto = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bthreshold_Stress_auto must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
					my_amg_manager.bthreshold_Stress_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				}
			}
			else {
				printf("WARNING!!! bthreshold_Stress_auto not found in file premeshin.txt\n");
				// 0 - Значение пользователя для threshold одно на всех уровнях выставляемое вручную.
				my_amg_manager.bthreshold_Stress_auto = false; // 1 - автоматически выставляемые значения threshold на каждом уровне сетки.
				printf("my_amg_manager.bthreshold_Stress_auto =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bcf_reorder_Temperature", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bcf_reorder_Temperature = false;
					}
					else {
						my_amg_manager.bcf_reorder_Temperature = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bcf_reorder_Temperature = false\n");
					// } else {
					// printf("my_amg_manager.bcf_reorder_Temperature = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bcf_reorder_Temperature must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - CF упорядочивание при итерировании не используется.
					my_amg_manager.bcf_reorder_Temperature = false; // 1 - при итерировании используется CF упорядочивание.
				}
			}
			else {
				printf("WARNING!!! bcf_reorder_Temperature not found in file premeshin.txt\n");
				// 0 - CF упорядочивание при итерировании не используется.
				my_amg_manager.bcf_reorder_Temperature = false; // 1 - при итерировании используется CF упорядочивание.
				printf("my_amg_manager.bcf_reorder_Temperature =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bcf_reorder_Speed", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bcf_reorder_Speed = false;
					}
					else {
						my_amg_manager.bcf_reorder_Speed = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bcf_reorder_Speed = false\n");
					// } else {
					// printf("my_amg_manager.bcf_reorder_Speed = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bcf_reorder_Speed must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - CF упорядочивание при итерировании не используется.
					my_amg_manager.bcf_reorder_Speed = false; // 1 - при итерировании используется CF упорядочивание.
				}
			}
			else {
				printf("WARNING!!! bcf_reorder_Speed not found in file premeshin.txt\n");
				// 0 - CF упорядочивание при итерировании не используется.
				my_amg_manager.bcf_reorder_Speed = false; // 1 - при итерировании используется CF упорядочивание.
				printf("my_amg_manager.bcf_reorder_Speed =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bcf_reorder_Pressure", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bcf_reorder_Pressure = false;
					}
					else {
						my_amg_manager.bcf_reorder_Pressure = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bcf_reorder_Pressure = false\n");
					// } else {
					// printf("my_amg_manager.bcf_reorder_Pressure = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bcf_reorder_Pressure must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - CF упорядочивание при итерировании не используется.
					my_amg_manager.bcf_reorder_Pressure = false; // 1 - при итерировании используется CF упорядочивание.
				}
			}
			else {
				printf("WARNING!!! bcf_reorder_Pressure not found in file premeshin.txt\n");
				// 0 - CF упорядочивание при итерировании не используется.
				my_amg_manager.bcf_reorder_Pressure = false; // 1 - при итерировании используется CF упорядочивание.
				printf("my_amg_manager.bcf_reorder_Pressure =false\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("bcf_reorder_Stress", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 0) {
						my_amg_manager.bcf_reorder_Stress = false;
					}
					else {
						my_amg_manager.bcf_reorder_Stress = true;
					}
					// if (idin == 0) {
					//printf("my_amg_manager.bcf_reorder_Stress = false\n");
					// } else {
					// printf("my_amg_manager.bcf_reorder_Stress = true\n");
					//} 
				}
				else {
					printf("my_amg_manager.bcf_reorder_Stress must be equal 1 or 0. current_value=(bool)(%d)\n", idin);
					system("pause");
					// 0 - CF упорядочивание при итерировании не используется.
					my_amg_manager.bcf_reorder_Stress = false; // 1 - при итерировании используется CF упорядочивание.
				}
			}
			else {
				printf("WARNING!!! bcf_reorder_Stress not found in file premeshin.txt\n");
				// 0 - CF упорядочивание при итерировании не используется.
				my_amg_manager.bcf_reorder_Stress = false; // 1 - при итерировании используется CF упорядочивание.
				printf("my_amg_manager.bcf_reorder_Stress =false\n");
				if (bSTOP_Reading) system("pause");
			}



			if (imakesource("amgcl_smoother", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 11)) {
					my_amg_manager.amgcl_smoother = static_cast<integer>(idin);
					//printf("my_amg_manager.amgcl_smoother =%lld\n", my_amg_manager.amgcl_smoother);
				}
				else {
					printf("my_amg_manager.amgcl_smoother must be 0<= value <=11. current_value=%d\n", idin);
					system("pause");
					// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped-jacobi.
					// 4 - spai1; 5 - chebyshev; 6 - ilu1 (iluk,k=1); 7 - ilu2 (iluk,k=2);
					// 8 - iluk, k=4; 9 - iluk, k=6; 10 - iluk, k=8; 11 - iluk, k=10.
					my_amg_manager.amgcl_smoother = 0; 
				}
			}
			else {
				printf("WARNING!!! amgcl_smoother not found in file premeshin.txt\n");
				// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped-jacobi.
                // 4 - spai1; 5 - chebyshev; 6 - ilu1 (iluk,k=1); 7 - ilu2 (iluk,k=2);
				// 8 - iluk, k=4; 9 - iluk, k=6; 10 - iluk, k=8; 11 - iluk, k=10.
				my_amg_manager.amgcl_smoother = 0; 
				printf("my_amg_manager.amgcl_smoother =%lld\n", my_amg_manager.amgcl_smoother);
				printf("0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped-jacobi;\n");
				printf("4 - spai1; 5 - chebyshev; 6 - ilu1 (iluk,k=1); 7 - ilu2 (iluk,k=2).\n");
				printf("8 - iluk, k=4; 9 - iluk, k=6; 10 - iluk, k=8; 11 - iluk, k=10.\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("amgcl_selector", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 1)) {
					my_amg_manager.amgcl_selector = static_cast<integer>(idin);
					//printf("my_amg_manager.amgcl_selector =%lld\n", my_amg_manager.amgcl_selector);
				}
				else {
					printf("my_amg_manager.amgcl_selector must be 0<= value <=1. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.amgcl_selector = 1; // 0 - John W. Ruge и K. Stueben (analog amg1r5); 1 - smoother aggregation.
				}
			}
			else {
				printf("WARNING!!! amgcl_selector not found in file premeshin.txt\n");
				my_amg_manager.amgcl_selector = 1; // 0 - John W. Ruge и K. Stueben  (analog amg1r5); 1 - smoother aggregation.
				printf("my_amg_manager.amgcl_selector =%lld\n", my_amg_manager.amgcl_selector);
				printf(" 0 - Ruge-Stueben (analog amg1r5); 1 - smoother aggregation. \n");
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("amgcl_iterator", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 1)) {
					switch (idin) {
					case 0: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
						break;
					case 1: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::FGMRes;
						break;
					default: my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab;
						break;
					}
					//printf("my_amg_manager.amgcl_iterator =%lld\n", my_amg_manager.amgcl_iterator);
				}
				else {
					printf("my_amg_manager.amgcl_iterator must be 0<= value <=1. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab; // 0 - BiCGStab; 1 - FGMRes.
				}
			}
			else {
				printf("WARNING!!! amgcl_iterator not found in file premeshin.txt\n");
				my_amg_manager.amgcl_iterator = AMGCL_ITERATOR_ALG::BiCGStab; // 0 - BiCGStab; 1 - FGMRes.
				printf("my_amg_manager.amgcl_iterator =AMGCL_ITERATOR_ALG::BiCGStab\n");
				printf("0 - BiCGStab; 1 - FGMRes.\n");
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("lfil", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 7)) {
					my_amg_manager.lfil = static_cast<integer>(idin);
					//printf("my_amg_manager.lfil =%lld\n", my_amg_manager.lfil);
				}
				else {
					printf("my_amg_manager.lfil must be 0<= value <=7. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.lfil = 2; // lfil=1. заполнение в неполном LU разложение.
				}
			}
			else {
				printf("WARNING!!! lfil not found in file premeshin.txt\n");
				my_amg_manager.lfil = 1; // lfil=1. заполнение в неполном LU разложение.
				printf("my_amg_manager.lfil =%lld\n", my_amg_manager.lfil);
				if (bSTOP_Reading) system("pause");
			}

			
			if (imakesource("reconstruct", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					ireconstruction_free_construct_alloc = static_cast<integer>(idin);
					//printf("ireconstruction_free_construct_alloc =%lld\n", ireconstruction_free_construct_alloc);
				}
				else {
					printf("ireconstruction_free_construct_alloc must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					ireconstruction_free_construct_alloc = 1; // 1 - перестраивать все структуры при решении экономия оперативной памяти.
				}
			}
			else {
				printf("WARNING!!! reconstruct not found in file premeshin.txt\n");
				ireconstruction_free_construct_alloc = 1; // 1 - перестраивать все структуры при решении экономия оперативной памяти.
				printf("ireconstruction_free_construct_alloc =%lld\n", ireconstruction_free_construct_alloc);
				if (bSTOP_Reading) system("pause");
			}

		

			
			if (imakesource("AnimationFields", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					ianimation_write_on = static_cast<integer>(idin);
					//printf("ianimation_write_on =%lld\n", ianimation_write_on);
				}
				else {
					printf("ianimation_write_on must be equal 1 or 0. current_value=%d\n",idin);
					system("pause");
					ianimation_write_on = 0; // 0 - не писать анимацию.
				}
			}
			else {
				printf("WARNING!!! AnimationFields not found in file premeshin.txt\n");
				ianimation_write_on = 0; // 0 - не писать анимацию.
				printf("ianimation_write_on =%lld\n", ianimation_write_on);
				if (bSTOP_Reading) system("pause");
			}


			// выделение оперативной памяти.
			gtdps = new TEMP_DEP_POWER[ltdp];
			matlist = new TPROP[lmatmax];
			if (b == nullptr) {
				b = new BLOCK[lb];
				for (integer i_7 = 0; i_7 < lb; i_7++) {
					b[i_7].temp_Sc = nullptr;
					b[i_7].arr_Sc = nullptr;
					b[i_7].g.hi = nullptr;
					b[i_7].g.xi = nullptr;
					b[i_7].g.yi = nullptr;
					b[i_7].g.zi = nullptr;
					b[i_7].g.root_CAD_STL = nullptr;
				}
			}
			else {
				for (integer i_7 = 0; i_7 < lb; i_7++) {

					b[i_7].g.clear_CAD_STL();

					if (b[i_7].temp_Sc != nullptr) {
						delete[] b[i_7].temp_Sc;
						b[i_7].temp_Sc = nullptr;
					}
					if (b[i_7].arr_Sc != nullptr) {
						delete[] b[i_7].arr_Sc;
						b[i_7].arr_Sc = nullptr;
					}
					if (b[i_7].g.hi != nullptr) {
						delete[] b[i_7].g.hi;
						b[i_7].g.hi = nullptr;
					}
					if (b[i_7].g.xi != nullptr) {
						delete[] b[i_7].g.xi;
						b[i_7].g.xi = nullptr;
					}
					if (b[i_7].g.yi != nullptr) {
						delete[] b[i_7].g.yi;
						b[i_7].g.yi = nullptr;
					}
					if (b[i_7].g.zi != nullptr) {
						delete[] b[i_7].g.zi;
						b[i_7].g.zi = nullptr;
					}
				}
				delete[] b;
				b = new BLOCK[lb];
				for (integer i_7 = 0; i_7 < lb; i_7++) {
					b[i_7].temp_Sc = nullptr;
					b[i_7].arr_Sc = nullptr;
					b[i_7].g.hi = nullptr;
					b[i_7].g.xi = nullptr;
					b[i_7].g.yi = nullptr;
					b[i_7].g.zi = nullptr;
					b[i_7].g.root_CAD_STL = nullptr;
				}
			}
			//if (s == nullptr) {
				//s = new SOURCE[ls];
			//}
			//else {
				delete[] s;
				s = new SOURCE[ls];
			//}
			//if (w == nullptr) {
				//w = new WALL[lw+1]; // +1 это гран условие по умолчанию, заглушка, чтобы избежать неопределенности.
			//}
			//else {
				delete[] w;
				w = new WALL[lw+1];// +1 это гран условие по умолчанию, заглушка, чтобы избежать неопределенности.
			//}
			int i = 0; // счётчик цикла for

			for (i = 0; i < ltdp; ++i) {
				// считывание имён файлов.
				gtdps[i].sname = new char[100]; // выделение памяти
				fscanf_s(fp, "%s", gtdps[i].sname, 100);
				//printf("%s",gtdps[i].sname);
				//system("PAUSE");
				// построение таблицы в памяти.
				my_read_power_table(gtdps[i].sname, gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
				printf("printeger table\n"); // debug
				my_print_table(gtdps[i].intemp, gtdps[i].inoffset_drain, gtdps[i].rtemp, gtdps[i].roffset_drain, gtdps[i].rpower_table);
				printf("Please, press any key to continue...\n");
				//system("PAUSE");
				system("pause");

			}

			char buffer[1000];

			// считывание базы материалов
			for (i = 0; i < lmatmax; ++i) {
				// свойства материалов:
				// плотность
			
				char name0[1000] = "matherial";
				
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				
				strcat_s(name0,"_rho");
				float ffin = 0.0;

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].rho = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].rho = 1.0; // плотность.
					printf("matlist[i].rho =%e\n", matlist[i].rho);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "n_cp");
				// теплоёмкость при постоянном давлении
				//fscanf(fp, "%f", &fin);
				//matlist[i].cp = fin;

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].n_cp = static_cast<integer>(idin);
					//printf("matlist[i].n_cp =%lld\n",matlist[i].n_cp );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].n_cp = 1; // .
					printf("matlist[i].n_cp =%lld\n", matlist[i].n_cp);
					if (bSTOP_Reading) system("pause");
				}

				matlist[i].arr_cp = nullptr;
				matlist[i].temp_cp = nullptr;
				matlist[i].arr_cp = new float[matlist[i].n_cp];
				matlist[i].temp_cp = new float[matlist[i].n_cp];
				if (matlist[i].temp_cp == nullptr) {
					printf("problem memory allocation for temp_cp\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_cp == nullptr) {
					printf("problem memory allocation for arr_cp\n");
					system("pause");
					exit(1);
				}
				for (int i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
					// Температура в C.
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "temp_cp");

					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);
					
					float ffin = 0.0;

					if (fmakesource_float_version(name0, ffin)) {
						// Найдено успешно.
						matlist[i].temp_cp[i_4] = (ffin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						matlist[i].temp_cp[i_4] = 30.0; // .
						printf(" matlist[%d].temp_cp[%d]=%e\n", i, i_4, matlist[i].temp_cp[i_4]);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "arr_cp");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);
					
					if (fmakesource_float_version(name0, ffin)) {
						// Найдено успешно.
						matlist[i].arr_cp[i_4] = (ffin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
						matlist[i].arr_cp[i_4] = 700.0; // .
						printf("matlist[%d].arr_cp[%d]  =%e\n",i, i_4, matlist[i].arr_cp[i_4]);
						if (bSTOP_Reading) system("pause");
					}

				}
				// теплопроводность

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "n_lam");

				//fscanf(fp, "%f", &fin);
				//matlist[i].lam = fin;
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].n_lam = static_cast<integer>(idin);
					//printf("matlist[%lld].n_lam =%lld\n", i, matlist[i].n_lam);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].n_lam = 1; // Количество различных значений теплопропроводности при её зависимости от температуры.
					printf("matlist[%d].n_lam =%lld\n",i, matlist[i].n_lam);
					if (bSTOP_Reading) system("pause");
				}

				matlist[i].arr_lam = nullptr;
				matlist[i].temp_lam = nullptr;
				matlist[i].arr_lam = new float[matlist[i].n_lam];
				matlist[i].temp_lam = new float[matlist[i].n_lam];
				if (matlist[i].temp_lam == nullptr) {
					printf("problem memory allocation for temp_lam\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_lam == nullptr) {
					printf("problem memory allocation for arr_lam\n");
					system("pause");
					exit(1);
				}
				for (int i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "temp_lam");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);

					float ffin = 0.0;

					// Температура в C.
					if (fmakesource_float_version(name0, ffin)) {
						// Найдено успешно.
						matlist[i].temp_lam[i_4] = (ffin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
						matlist[i].temp_lam[i_4] = 30.0; // Температура при которой задана теплопроводность.
						printf(" matlist[%d].temp_lam[%d]=%e\n", i, i_4, matlist[i].temp_lam[i_4]);
						if (bSTOP_Reading) system("pause");
					}

					
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "arr_lam");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);


					if (fmakesource_float_version(name0, ffin)) {
						// Найдено успешно.
						matlist[i].arr_lam[i_4] = (ffin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
						matlist[i].arr_lam[i_4] = 4.0; // ЭЧЭС default.
						printf("matlist[%d].arr_lam[%d] =%e\n", i, i_4, matlist[i].arr_lam[i_4]);
						if (bSTOP_Reading) system("pause");
					}
				}
				// ортотропность теплопроводности:
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_lam_x");
				

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_x = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].orthotropy_multiplyer_x = 1.0; // no orthotropy .
					printf("matlist[%d].orthotropy_multiplyer_x =%e\n",i, matlist[i].orthotropy_multiplyer_x);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_lam_y");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_y = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].orthotropy_multiplyer_y = 1.0; // no orthotropy .
					printf("matlist[%d].orthotropy_multiplyer_y =%e\n", i, matlist[i].orthotropy_multiplyer_y);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_lam_z");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_z = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].orthotropy_multiplyer_z = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_z=%e\n", i, matlist[i].orthotropy_multiplyer_z);
					if (bSTOP_Reading) system("pause");
				}

				// 28.08.2020.
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Linear_expansion_coefficient_x");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_x_beta_t_solid = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_x_beta_t_solid = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_x_beta_t_solid=%e\n", i, matlist[i].orthotropy_multiplyer_x_beta_t_solid);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Linear_expansion_coefficient_y");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_y_beta_t_solid = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_y_beta_t_solid = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_y_beta_t_solid=%e\n", i, matlist[i].orthotropy_multiplyer_y_beta_t_solid);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Linear_expansion_coefficient_z");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_z_beta_t_solid = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_z_beta_t_solid = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_z_beta_t_solid=%e\n", i, matlist[i].orthotropy_multiplyer_z_beta_t_solid);
					if (bSTOP_Reading) system("pause");
				}

				// Ортотропность модуля Юнга.
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Young_Module_x");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_x_Young_Module = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_x_Young_Module = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_x_Young_Module=%e\n", i, matlist[i].orthotropy_multiplyer_x_Young_Module);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Young_Module_y");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_y_Young_Module = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_y_Young_Module = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_y_Young_Module=%e\n", i, matlist[i].orthotropy_multiplyer_y_Young_Module);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Young_Module_z");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_z_Young_Module = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_z_Young_Module = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_z_Young_Module=%e\n", i, matlist[i].orthotropy_multiplyer_z_Young_Module);
					if (bSTOP_Reading) system("pause");
				}

				// mult_Poisson_ratio_xy
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Poisson_ratio_xy");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_Poisson_ratio_xy = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_Poisson_ratio_xy = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_Poisson_ratio_xy=%e\n", i, matlist[i].orthotropy_multiplyer_Poisson_ratio_xy);
					if (bSTOP_Reading) system("pause");
				}

				// mult_Poisson_ratio_xz
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Poisson_ratio_xz");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_Poisson_ratio_xz = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_Poisson_ratio_xz = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_Poisson_ratio_xz=%e\n", i, matlist[i].orthotropy_multiplyer_Poisson_ratio_xz);
					if (bSTOP_Reading) system("pause");
				}

				// mult_Poisson_ratio_yz
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Poisson_ratio_yz");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_Poisson_ratio_yz = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_Poisson_ratio_yz = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_Poisson_ratio_yz=%e\n", i, matlist[i].orthotropy_multiplyer_Poisson_ratio_yz);
					if (bSTOP_Reading) system("pause");
				}

				// mult_Poisson_ratio_yx
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Poisson_ratio_yx");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_Poisson_ratio_yx = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_Poisson_ratio_yx = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_Poisson_ratio_yx=%e\n", i, matlist[i].orthotropy_multiplyer_Poisson_ratio_yx);
					if (bSTOP_Reading) system("pause");
				}

				// mult_Poisson_ratio_zx
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Poisson_ratio_zx");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_Poisson_ratio_zx = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_Poisson_ratio_zx = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_Poisson_ratio_zx=%e\n", i, matlist[i].orthotropy_multiplyer_Poisson_ratio_zx);
					if (bSTOP_Reading) system("pause");
				}

				// mult_Poisson_ratio_zy
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_Poisson_ratio_zy");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_Poisson_ratio_zy = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].orthotropy_multiplyer_Poisson_ratio_zy = 1.0; // no orthotropy .
					printf(" matlist[%d].orthotropy_multiplyer_Poisson_ratio_zy=%e\n", i, matlist[i].orthotropy_multiplyer_Poisson_ratio_zy);
					if (bSTOP_Reading) system("pause");
				}

				//bShearModuleActive
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "bShearModuleActive");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					if (idin == 1) {
						matlist[i].bActive_ShearModule = true;
					}
					else {
						matlist[i].bActive_ShearModule = false;
					}
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].bActive_ShearModule = false; // no orthotropy .
					printf(" matlist[%d].bActive_ShearModule=false\n", i);
					if (bSTOP_Reading) system("pause");
				}

				//ShearModuleGxy
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "ShearModuleGxy");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].ShearModule_xy = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].ShearModule_xy = 1.0; // no orthotropy .
					printf(" matlist[%d].ShearModule_xy=%e\n", i, matlist[i].ShearModule_xy);
					if (bSTOP_Reading) system("pause");
				}

				//ShearModuleGyz
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "ShearModuleGyz");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].ShearModule_yz = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].ShearModule_yz = 1.0; // no orthotropy .
					printf(" matlist[%d].ShearModule_yz=%e\n", i, matlist[i].ShearModule_yz);
					if (bSTOP_Reading) system("pause");
				}

				//ShearModuleGxz
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "ShearModuleGxz");

				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].ShearModule_xz = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].ShearModule_xz = 1.0; // no orthotropy .
					printf(" matlist[%d].ShearModule_xz=%e\n", i, matlist[i].ShearModule_xz);
					if (bSTOP_Reading) system("pause");
				}


				// 5.08.2017.
				// Коэффициенты для задачи упругости.
				// Модуль Юнга и коэффициент Пуассона.
				// В.Н.Сидоров, В.В. Вершинин Метод конечных элементов в расчёте сооружений.
				

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "n_Poisson_ratio");
				// Коэффициент Пуассона
				//fscanf(fp, "%f", &fin);
				//matlist[i].arr_Poisson_ratio[0]= fin;

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].n_Poisson_ratio = static_cast<integer>(idin);
					//printf("matlist[i].n_Poisson_ratio =%lld\n",matlist[i].n_Poisson_ratio);

					if (matlist[i].arr_Poisson_ratio != nullptr) {
						delete[] matlist[i].arr_Poisson_ratio;
						matlist[i].arr_Poisson_ratio = nullptr;
					}
					if (matlist[i].temp_Poisson_ratio != nullptr) {
						delete[] matlist[i].temp_Poisson_ratio;
						matlist[i].temp_Poisson_ratio = nullptr;
					}
					
					
					if (matlist[i].n_Poisson_ratio > 0) {
						matlist[i].arr_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
						matlist[i].temp_Poisson_ratio = new float[matlist[i].n_Poisson_ratio];
					}
					if (matlist[i].temp_Poisson_ratio == nullptr) {
						printf("problem memory allocation for temp_Poisson_ratio\n");
						system("pause");
						exit(1);
					}
					if (matlist[i].arr_Poisson_ratio == nullptr) {
						printf("problem memory allocation for arr_Poisson_ratio\n");
						system("pause");
						exit(1);
					}
					for (int i_4 = 0; i_4 < matlist[i].n_Poisson_ratio; i_4++) {
						// Температура в C.
						name0[0] = '\0'; strcat_s(name0, "matherial");
						buffer[0] = '\0';
						_itoa_s(i, buffer, 10);
						strcat_s(name0, buffer);
						strcat_s(name0, "temp_Poisson_ratio");

						buffer[0] = '\0'; _itoa_s(i_4, buffer, 10); strcat_s(name0, buffer);

						if (fmakesource_float_version(name0, ffin)) {
							// Найдено успешно.
							matlist[i].temp_Poisson_ratio[i_4] = (ffin);
						}
						else {

							printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
							matlist[i].temp_Poisson_ratio[i_4] = 30.0; // .
							printf(" matlist[%d].temp_Poisson_ratio[%d]=%e\n", i, i_4, matlist[i].temp_Poisson_ratio[i_4]);
							if (bSTOP_Reading) system("pause");
						}

						name0[0] = '\0'; strcat_s(name0, "matherial");
						buffer[0] = '\0';
						_itoa_s(i, buffer, 10);
						strcat_s(name0, buffer);
						strcat_s(name0, "arr_Poisson_ratio");
						buffer[0] = '\0'; _itoa_s(i_4, buffer, 10); strcat_s(name0, buffer);

						if (fmakesource_float_version(name0, ffin)) {
							// Найдено успешно.
							matlist[i].arr_Poisson_ratio[i_4] = (ffin);
						}
						else {

							printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
							matlist[i].arr_Poisson_ratio[i_4] = 0.3f; // .
							printf("matlist[%d].arr_Poisson_ratio[%d]  =%e\n", i, i_4, matlist[i].arr_Poisson_ratio[i_4]);
							if (bSTOP_Reading) system("pause");

						}

					}

				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].n_Poisson_ratio = 1; // .
					printf("matlist[i].n_Poisson_ratio =%lld\n", matlist[i].n_Poisson_ratio);
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "Poisson_ratio");

					// Poisson ratio

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].n_Poisson_ratio = 1;
						matlist[i].temp_Poisson_ratio = new float[1];
						matlist[i].temp_Poisson_ratio[0] = 25.0;
						matlist[i].arr_Poisson_ratio = new float[1];
						matlist[i].arr_Poisson_ratio[0] = static_cast <float>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						matlist[i].n_Poisson_ratio = 1;
						matlist[i].temp_Poisson_ratio = new float[1];
						matlist[i].temp_Poisson_ratio[0] = 25.0;
						matlist[i].arr_Poisson_ratio = new float[1];
						matlist[i].arr_Poisson_ratio[0] = 0.3f; // default 0.3.
						printf(" matlist[%d].arr_Poisson_ratio =%e\n", i, matlist[i].arr_Poisson_ratio[0]);
						if (bSTOP_Reading) system("pause");
					}
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "n_Young_Module");
				// Модуль Юнга
				//fscanf(fp, "%f", &fin);
				//matlist[i].arr_Young_Module[0]= fin*1e9;

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].n_YoungModule = static_cast<integer>(idin);
					//printf("matlist[i].n_YoungModule =%lld\n",matlist[i].n_YoungModule);

					if (matlist[i].arr_Young_Module != nullptr) {
						delete[] matlist[i].arr_Young_Module;
						matlist[i].arr_Young_Module = nullptr;
					}
					if (matlist[i].temp_Young_Module != nullptr) {
						delete[] matlist[i].temp_Young_Module;
						matlist[i].temp_Young_Module = nullptr;
					}					
					
					matlist[i].arr_Young_Module = new float[matlist[i].n_YoungModule];
					matlist[i].temp_Young_Module = new float[matlist[i].n_YoungModule];
					if (matlist[i].temp_Young_Module == nullptr) {
						printf("problem memory allocation for temp_Young_Module\n");
						system("pause");
						exit(1);
					}
					if (matlist[i].arr_Young_Module == nullptr) {
						printf("problem memory allocation for arr_Young_Module\n");
						system("pause");
						exit(1);
					}
					for (int i_4 = 0; i_4 < matlist[i].n_YoungModule; i_4++) {
						// Температура в C.
						name0[0] = '\0'; strcat_s(name0, "matherial");
						buffer[0] = '\0';
						_itoa_s(i, buffer, 10);
						strcat_s(name0, buffer);
						strcat_s(name0, "temp_Young_Module");

						buffer[0] = '\0'; _itoa_s(i_4, buffer, 10); strcat_s(name0, buffer);

						if (fmakesource_float_version(name0, ffin)) {
							// Найдено успешно.
							matlist[i].temp_Young_Module[i_4] = (ffin);
						}
						else {
							
								printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
								matlist[i].temp_Young_Module[i_4] = 30.0; // .
								printf(" matlist[%d].temp_Young_Module[%d]=%e\n", i, i_4, matlist[i].temp_Young_Module[i_4]);
								if (bSTOP_Reading) system("pause");
						}

						name0[0] = '\0'; strcat_s(name0, "matherial");
						buffer[0] = '\0';
						_itoa_s(i, buffer, 10);
						strcat_s(name0, buffer);
						strcat_s(name0, "arr_Young_Module");
						buffer[0] = '\0'; _itoa_s(i_4, buffer, 10); strcat_s(name0, buffer);

						if (fmakesource(name0, fin)) {
							// Найдено успешно.
							matlist[i].arr_Young_Module[i_4] = static_cast <float>((fin)*1.0e+9);
						}
						else {
							
								printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
								matlist[i].arr_Young_Module[i_4] = static_cast<float>(200.0e+9); // .
								printf("matlist[%d].arr_Young_Module[%d]  =%e\n", i, i_4, matlist[i].arr_Young_Module[i_4]);
								if (bSTOP_Reading) system("pause");
							
						}

					}

				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].n_YoungModule = 1; // .
					printf("matlist[i].n_YoungModule =%lld\n", matlist[i].n_YoungModule);
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "Young_Module");

					// Young_Module*1E+9

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].n_YoungModule = 1;
						matlist[i].temp_Young_Module = new float[1];
						matlist[i].temp_Young_Module[0] = 25.0;
						matlist[i].arr_Young_Module = new float[1];
						matlist[i].arr_Young_Module[0] = static_cast <float>((fin)* 1E+9);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						matlist[i].n_YoungModule = 1;
						matlist[i].temp_Young_Module  = new float[1];
						matlist[i].temp_Young_Module[0] = 25.0;
						matlist[i].arr_Young_Module  = new float[1];
						matlist[i].arr_Young_Module[0] = static_cast < float>(200.0E+9); // default 200GPa.
						printf(" matlist[%d].Young_Module =%e\n", i, matlist[i].arr_Young_Module[0]);
						if (bSTOP_Reading) system("pause");
					}
				}
								
				

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "n_Linear_expansion_coefficient");
				// теплоёмкость при постоянном давлении
				//fscanf(fp, "%f", &fin);
				//matlist[i].cp = fin;

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].n_beta_t_solid = static_cast<integer>(idin);
					//printf("matlist[i].n_beta_t_solid =%lld\n",matlist[i].n_beta_t_solid );

					if (matlist[i].arr_beta_t_solid != nullptr) {
						delete[] matlist[i].arr_beta_t_solid;
						matlist[i].arr_beta_t_solid = nullptr;
					}
					if (matlist[i].temp_beta_t_solid != nullptr) {
						delete[] matlist[i].temp_beta_t_solid;
						matlist[i].temp_beta_t_solid = nullptr;
					}
					
					
					matlist[i].arr_beta_t_solid = new float[matlist[i].n_beta_t_solid];
					matlist[i].temp_beta_t_solid = new float[matlist[i].n_beta_t_solid];
					if (matlist[i].temp_beta_t_solid == nullptr) {
						printf("problem memory allocation for temp_beta_t_solid\n");
						system("pause");
						exit(1);
					}
					if (matlist[i].arr_beta_t_solid == nullptr) {
						printf("problem memory allocation for arr_beta_t_solid\n");
						system("pause");
						exit(1);
					}
					for (int i_4 = 0; i_4 < matlist[i].n_beta_t_solid; i_4++) {
						// Температура в C.
						name0[0] = '\0'; strcat_s(name0, "matherial");
						buffer[0] = '\0';
						_itoa_s(i, buffer, 10);
						strcat_s(name0, buffer);
						strcat_s(name0, "temp_Linear_expansion_coefficient");

						buffer[0] = '\0'; _itoa_s(i_4, buffer, 10); strcat_s(name0, buffer);

						if (fmakesource_float_version(name0, ffin)) {
							// Найдено успешно.
							matlist[i].temp_beta_t_solid[i_4] = (ffin);
						}
						else {
							
								printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
								matlist[i].temp_beta_t_solid[i_4] = 30.0; // .
								printf(" matlist[%d].temp_beta_t_solid[%d]=%e\n", i, i_4, matlist[i].temp_beta_t_solid[i_4]);
								if (bSTOP_Reading) system("pause");
							
						}

						name0[0] = '\0'; strcat_s(name0, "matherial");
						buffer[0] = '\0';
						_itoa_s(i, buffer, 10);
						strcat_s(name0, buffer);
						strcat_s(name0, "arr_Linear_expansion_coefficient");
						buffer[0] = '\0'; _itoa_s(i_4, buffer, 10); strcat_s(name0, buffer);

						if (fmakesource(name0, fin)) {
							// Найдено успешно.
							matlist[i].arr_beta_t_solid[i_4] = static_cast<float>((fin)*1.0e-6);
						}
						else {
							
								printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
								matlist[i].arr_beta_t_solid[i_4] = static_cast < float>(1.0e-6); // .
								printf("matlist[%d].arr_beta_t_solid[%d]  =%e\n", i, i_4, matlist[i].arr_beta_t_solid[i_4]);
								if (bSTOP_Reading) system("pause");
							
						}

					}


				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].n_beta_t_solid = 1; // .
					printf("matlist[i].n_beta_t_solid =%lld\n", matlist[i].n_beta_t_solid);
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "Linear_expansion_coefficient");

					// beta_t_solid*1E-6

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].n_beta_t_solid = 1;
						matlist[i].temp_beta_t_solid = new float[1];
						matlist[i].temp_beta_t_solid[0] = 25.0;
						matlist[i].arr_beta_t_solid = new float[1];
						matlist[i].arr_beta_t_solid[0] = static_cast<float>((fin)* 1E-6);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						matlist[i].n_beta_t_solid = 1;
						matlist[i].temp_beta_t_solid = new float[1];
						matlist[i].temp_beta_t_solid[0] = 25.0;
						matlist[i].arr_beta_t_solid = new float[1];
						matlist[i].arr_beta_t_solid[0] = static_cast < float>(1.0E-6); // default 1.0e-6.
						printf(" matlist[%d].beta_t_solid =%e\n", i, matlist[i].arr_beta_t_solid[0]);
						if (bSTOP_Reading) system("pause");
					}					
				}

				


				// Коэффициенты Ламе не используются в коде.
				
				// коэффициент динамической вязкости
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mu");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].mu = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].mu = static_cast < float>(1.7894e-5); // air.
					printf("matlist[%d].mu =%e\n", i, matlist[i].mu);
					if (bSTOP_Reading) system("pause");
				}

				// коэффициент линейного температурного расширения
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "beta_t");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].beta_t = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].beta_t = 0.003331f; // air.
					printf(" matlist[%d].beta_t =%e\n", i, matlist[i].beta_t);
					if (bSTOP_Reading) system("pause");
				}

				// признак библиотечности материала
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "blibmat");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].blibmat = static_cast<integer>(idin);
					//printf(" matlist[%lld].blibmat=%lld\n", i, matlist[i].blibmat);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].blibmat = 0; // не библиотечный материал.
					printf(" matlist[%d].blibmat=%lld\n",i, matlist[i].blibmat);
					if (bSTOP_Reading) system("pause");
				}
				// номер материала в библиотеке
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ilibident");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].ilibident = static_cast<integer>(idin);
					//printf("matlist[%lld].ilibident =%lld\n", i, matlist[i].ilibident);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].ilibident = 0; // .
					printf(" matlist[%d].ilibident=%lld\n", i, matlist[i].ilibident);
					if (bSTOP_Reading) system("pause");
				}
				//printf("blibmat=%d ilibident=%d\n", matlist[i].blibmat, matlist[i].ilibident);
				//system("pause");

				// для каждого КО данного блока,
				// если только он не перекрывается другими блоками
				// может быть использовано приближение 
				// Обербека-Буссинеска с соответствующей опорной температурой Tref.
				name0[0] = '\0'; strcat_s(name0, "matherial");
				strcat_s(name0, buffer);
				strcat_s(name0, "bBoussinesq");
				
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din= static_cast<integer>(idin);
					switch (din) {
					case 0: matlist[i].bBussineskApproach = false; break;
					case 1: matlist[i].bBussineskApproach = true; break;
					default: matlist[i].bBussineskApproach = false; break;
					}
					//printf("matlist[%lld].bBussineskApproach =%lld\n",i, din);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].bBussineskApproach = true; // Приближение Обербека-Буссинеска.
					printf("matlist[%d].bBussineskApproach=true\n",i);
					if (bSTOP_Reading) system("pause");
				}
				// номер закона для зависимости динамической вязкости от напряжения сдвига
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ilawmu");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].ilawmu = static_cast<integer>(idin);
					//printf("matlist[%lld].ilawmu =%lld\n",i,  matlist[i].ilawmu);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].ilawmu = 0; // .
					printf(" matlist[%d].ilawmu =%lld\n",i, matlist[i].ilawmu);
					if (bSTOP_Reading) system("pause");
				}
				// минимальное значение динамической вязкости
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "mumin");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].mumin = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].mumin = 0.0; // .
					printf("matlist[%d].mumin =%e\n", i, matlist[i].mumin);
					if (bSTOP_Reading) system("pause");
				}

				// максимальное значение динамической вязкости
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "mumax");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].mumax = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].mumax = 1.0; // .
					printf("matlist[%d].mumax =%e\n", i, matlist[i].mumax);
					if (bSTOP_Reading) system("pause");
				}
				// параметры модельных законов для зависимости вязкости от напряжения сдвига
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Amu");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].Amu = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].Amu = 0.0; // .
					printf("matlist[%d].Amu =%e\n", i, matlist[i].Amu);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Bmu");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].Bmu = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].Bmu = 0.0; // .
					printf(" matlist[%d].Bmu =%e\n", i, matlist[i].Bmu);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Cmu");
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].Cmu = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].Cmu = 0.0; // .
					printf("matlist[%d].Cmu  =%e\n",i, matlist[i].Cmu );
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "degreennmu");
				// показатель степени
				
				if (fmakesource_float_version(name0, ffin)) {
					// Найдено успешно.
					matlist[i].degreennmu = (ffin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].degreennmu = 1.0; // .
					printf("matlist[%d].degreennmu =%e\n", i, matlist[i].degreennmu);
					if (bSTOP_Reading) system("pause");
				}

				// печать считанных значений на консоль
				//printf("%e %e %e %e %e\n", matlist[i].rho, matlist[i].cp, matlist[i].lam, matlist[i].mu, matlist[i].beta_t);
				if (0) {
					printf("HEAT_CAPACITY\n");
					printf("t_C HEAT_CAPACITY\n");
					for (int i_4 = 0; i_4 < matlist[i].n_cp; i_4++) {
						printf("%e %e\n", matlist[i].temp_cp[i_4], matlist[i].arr_cp[i_4]);
					}
					printf("lam\n");
					printf("t_C lam\n");
					for (int i_4 = 0; i_4 < matlist[i].n_lam; i_4++) {
						printf("%e %e\n", matlist[i].temp_lam[i_4], matlist[i].arr_lam[i_4]);
					}
					printf("%e %e %e\n", matlist[i].rho, matlist[i].mu, matlist[i].beta_t);
					printf("%lld %lld %lld\n", matlist[i].blibmat, matlist[i].ilibident, matlist[i].ilawmu); // bBoussinesq не печатается
					printf("%e %e %e %e %e %e\n", matlist[i].mumin, matlist[i].mumax, matlist[i].Amu, matlist[i].Bmu, matlist[i].Cmu, matlist[i].degreennmu);
				}
			}

			integer iblock_target = 0;

			// считывание блоков
			for (i = 0; i < lb; ++i) {

				char name0[1000] = "body";
				buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
				strcat_s(name0, "name");

				if (smakesource(name0, b[i].name)) {
					// Найдено успешно.

					//printf(" b[%lld].name=%s\n", i, b[i].name);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].name[0] = '\0'; // уникальное текстовое объёмного геометрического тела.
					printf(" b[%d].name=unknown\n", i);
					//if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iunion");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].iunion_id  = static_cast<integer>(idin); // 0 == Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					//printf(" b[%lld].iunion_id=%lld\n",i, b[i].iunion_id);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].iunion_id = 0; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					printf(" b[%d].iunion_id=%lld\n",i, b[i].iunion_id);
					if (bSTOP_Reading) system("pause");
				}


				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "igeometry_type");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din= static_cast<integer>(idin);
					CHECK_TYPE_GEOM(din);
					b[i].g.itypegeom = din;
					//printf(" b[%lld].g.itypegeom=%lld\n",i,b[i].g.itypegeom );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.itypegeom = 0; // 0 - Prism, 1 - Cylinder, 2 - Polygon, 3 - CAD_STL.
					printf(" b[%d].g.itypegeom=%lld\n", i, b[i].g.itypegeom);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bvisible");
								
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din == 1) {
						b[i].bvisible = true;
					}
					else {
						b[i].bvisible = false;
					}
					//printf("b[%lld].bvisible =%lld\n",i,din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].bvisible = true; // видимый блок.
					printf(" b[%d].bvisible=true\n",i );
					if (bSTOP_Reading) system("pause");
				}

				// геометрия
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.xS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.xS = scale*0.0; // .
					printf(" b[%d].g.xS=%e\n",i, b[i].g.xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.yS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.yS = scale * 0.0; // .
					printf(" b[%d].g.yS=%e\n",i, b[i].g.yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");
								
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.zS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.zS = scale * 0.0; // .
					printf(" b[%d].g.zS =%e\n", i, b[i].g.zS);
					if (bSTOP_Reading) system("pause");
				}
	
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.xE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.xE = scale * 0.0; // .
					printf("b[%d].g.xE =%e\n",i, b[i].g.xE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.yE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.yE = scale * 0.0; // .
					printf("b[%d].g.yE =%e\n",i, b[i].g.yE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");
			
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.zE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.zE = scale * 0.0; // .
					printf("b[%d].g.zE =%e\n",i, b[i].g.zE);
					if (bSTOP_Reading) system("pause");
				}

				// swap
				if (b[i].g.xS > b[i].g.xE) {
					dbuf = b[i].g.xS;
					b[i].g.xS = b[i].g.xE;
					b[i].g.xE = dbuf;
				}
				if (b[i].g.yS > b[i].g.yE) {
					dbuf = b[i].g.yS;
					b[i].g.yS = b[i].g.yE;
					b[i].g.yE = dbuf;
				}
				if (b[i].g.zS > b[i].g.zE) {
					dbuf = b[i].g.zS;
					b[i].g.zS = b[i].g.zE;
					b[i].g.zE = dbuf;
				}

				// Cylinder
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iPlane");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].g.iPlane = static_cast<integer>(idin);
					//printf("b[%lld].g.iPlane =%lld\n",i,b[i].g.iPlane );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.iPlane = 0; // .
					printf("b[%d].g.iPlane =%lld\n", i, b[i].g.iPlane);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xC");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.xC = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.xC = scale * 0.0; // .
					printf("b[%d].g.xC =%e\n",i, b[i].g.xC);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yC");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.yC = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.yC = scale * 0.0; // .
					printf("b[%d].g.yC =%e\n", i, b[i].g.yC);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zC");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.zC = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.zC = scale * 0.0; // .
					printf("b[%d].g.zC =%e\n", i, b[i].g.zC);
					if (bSTOP_Reading) system("pause");
				}
			
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Hcyl");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.Hcyl = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.Hcyl = scale * 0.0; // .
					printf("b[%d].g.Hcyl =%e\n",i, b[i].g.Hcyl);
					if (bSTOP_Reading) system("pause");
				}
				
				if (b[i].g.Hcyl <= 0.0) {
					// ликвидируем отрицательную высоту цилиндра.
					switch (b[i].g.iPlane) {
					case XY_PLANE:
						b[i].g.zC += b[i].g.Hcyl;
						break;
					case XZ_PLANE:
						b[i].g.yC += b[i].g.Hcyl;
						break;
					case YZ_PLANE:
						b[i].g.xC += b[i].g.Hcyl;
						break;
					}
					b[i].g.Hcyl = fabs(b[i].g.Hcyl);
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "R_out_cyl");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.R_out_cyl = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.R_out_cyl = scale * 0.0; // .
					printf("b[%d].g.R_out_cyl =%e\n",i, b[i].g.R_out_cyl);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "R_in_cyl");
								 
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.R_in_cyl = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.R_in_cyl = scale * 0.0; // .
					printf("b[%d].g.R_in_cyl =%e\n",i, b[i].g.R_in_cyl);
					if (bSTOP_Reading) system("pause");
				}

				if (b[i].g.itypegeom == CYLINDER) {
					// patch 10.01.2020
					// Для правильной генерации сеточных линий xpos, ypos, zpos.
					switch (b[i].g.iPlane) {
					case XY_PLANE: 
						b[i].g.xS = b[i].g.xC - b[i].g.R_out_cyl;
						b[i].g.xE = b[i].g.xC + b[i].g.R_out_cyl;
						b[i].g.yS = b[i].g.yC - b[i].g.R_out_cyl;
						b[i].g.yE = b[i].g.yC + b[i].g.R_out_cyl;
						b[i].g.zS = b[i].g.zC;
						b[i].g.zE = b[i].g.zC + b[i].g.Hcyl;
						break;
					case XZ_PLANE:
						b[i].g.xS = b[i].g.xC - b[i].g.R_out_cyl;
						b[i].g.xE = b[i].g.xC + b[i].g.R_out_cyl;
						b[i].g.zS = b[i].g.zC - b[i].g.R_out_cyl;
						b[i].g.zE = b[i].g.zC + b[i].g.R_out_cyl;
						b[i].g.yS = b[i].g.yC;
						b[i].g.yE = b[i].g.yC + b[i].g.Hcyl;
						break;
					case YZ_PLANE:
						b[i].g.zS = b[i].g.zC - b[i].g.R_out_cyl;
						b[i].g.zE = b[i].g.zC + b[i].g.R_out_cyl;
						b[i].g.yS = b[i].g.yC - b[i].g.R_out_cyl;
						b[i].g.yE = b[i].g.yC + b[i].g.R_out_cyl;
						b[i].g.xS = b[i].g.xC;
						b[i].g.xE = b[i].g.xC + b[i].g.Hcyl;
						break;
					}
				}

				// Polygon
				//printf("Polygon\n");
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iPlane_obj2");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].g.iPlane_obj2 = static_cast<integer>(idin);
					//printf("b[%lld].g.iPlane_obj2 =%lld\n", i, b[i].g.iPlane_obj2);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.iPlane_obj2 = 0; // .
					printf("b[%d].g.iPlane_obj2 =%lld\n", i, b[i].g.iPlane_obj2);
					if (bSTOP_Reading) system("pause");
				}
#if doubleintprecision == 1
				//printf("iPlane_obj2=%lld\n", b[i].g.iPlane_obj2);
#else
				//printf("iPlane_obj2=%d\n", b[i].g.iPlane_obj2);
#endif			
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "nsizei");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].g.nsizei = static_cast<integer>(idin);
					//printf(" b[%lld].g.nsizei=%lld\n",i, b[i].g.nsizei );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.nsizei = 3; // .
					printf("b[%d].g.nsizei =%lld\n", i, b[i].g.nsizei);
					if (bSTOP_Reading) system("pause");
				}
#if doubleintprecision == 1
				//printf("nsizei=%lld\n", b[i].g.nsizei);
#else
				//printf("nsizei=%d\n", b[i].g.nsizei);
#endif			
				delete[] b[i].g.hi;
				delete[] b[i].g.xi;
				delete[] b[i].g.yi;
				delete[] b[i].g.zi;
				b[i].g.hi = new doublereal[b[i].g.nsizei];
				b[i].g.xi = new doublereal[b[i].g.nsizei];
				b[i].g.yi = new doublereal[b[i].g.nsizei];
				b[i].g.zi = new doublereal[b[i].g.nsizei];
				for (int i73 = 0; i73 < b[i].g.nsizei; i73++) {
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "hi");

					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.hi[i73] = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.hi[i73] = scale * 0.0; // .
						printf(" b[%d].g.hi[%d]=%e\n", i,i73,  b[i].g.hi[i73]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "xi");
					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.xi[i73] = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.xi[i73] = scale * 0.0; // .
						printf("b[%d].g.xi[%d] =%e\n",i, i73, b[i].g.xi[i73]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "yi");
					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.yi[i73] = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.yi[i73] = scale * 0.0; // .
						printf("b[%d].g.yi[%d] =%e\n", i, i73, b[i].g.yi[i73]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "zi");
					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.zi[i73] = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.zi[i73] = scale * 0.0; // .
						printf(" b[%d].g.zi[%d]=%e\n",i,i73, b[i].g.zi[i73]);
						if (bSTOP_Reading) system("pause");
					}

					if (b[i].g.hi[i73] < 0.0) {
						// ликвидируем отрицательную высоту цилиндра.
						switch (b[i].g.iPlane_obj2) {
						case XY_PLANE:
							b[i].g.zi[i73] += b[i].g.hi[i73];
							break;
						case XZ_PLANE:
							b[i].g.yi[i73] += b[i].g.hi[i73];
							break;
						case YZ_PLANE:
							b[i].g.xi[i73] += b[i].g.hi[i73];
							break;
						}
						b[i].g.hi[i73] = fabs(b[i].g.hi[i73]);
					}

#if doubleintprecision == 1
					//printf("%lld h=%e x=%e y=%e z=%e",i73, b[i].g.hi[i73], b[i].g.xi[i73], b[i].g.yi[i73], b[i].g.zi[i73]);
#else
					//printf("%d h=%e x=%e y=%e z=%e", i73, b[i].g.hi[i73], b[i].g.xi[i73], b[i].g.yi[i73], b[i].g.zi[i73]);
#endif
				//system("PAUSE");
				}

				if ((b[i].g.itypegeom == POLYGON) && (b[i].g.nsizei > 0)) {
					// Мы заносим в окаймляющую призму границы полигона чтобы ускорить 
					// поиски в inmodel_temp быстро отсекая ячейки вне границ полигона.
					doublereal xmin53 = 1.0e30;
					doublereal ymin53 = 1.0e30;
					doublereal zmin53 = 1.0e30;
					doublereal xmax53 = -1.0e30;
					doublereal ymax53 = -1.0e30;
					doublereal zmax53 = -1.0e30;
					for (integer i73 = 0; i73 < b[i].g.nsizei; i73++) {
						if (b[i].g.xi[i73] > xmax53) xmax53 = b[i].g.xi[i73];
						if (b[i].g.yi[i73] > ymax53) ymax53 = b[i].g.yi[i73];
						if (b[i].g.zi[i73] > zmax53) zmax53 = b[i].g.zi[i73];
						if (b[i].g.xi[i73] < xmin53) xmin53 = b[i].g.xi[i73];
						if (b[i].g.yi[i73] < ymin53) ymin53 = b[i].g.yi[i73];
						if (b[i].g.zi[i73] < zmin53) zmin53 = b[i].g.zi[i73];
					}
					switch (b[i].g.iPlane_obj2) {
					case XY_PLANE:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmin53 + b[i].g.hi[0];
						break;
					case XZ_PLANE:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymin53 + b[i].g.hi[0];
						break;
					case YZ_PLANE:
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.xS = xmin53;
						b[i].g.xE = xmin53 + b[i].g.hi[0];
						break;
					}
				}

				if (b[i].g.itypegeom == CYLINDER) {
					// Cylinder
					
					//printf("%e %e %e %e %e %e %e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE, b[i].g.xC, b[i].g.yC, b[i].g.zC, b[i].g.Hcyl, b[i].g.R_out_cyl, b[i].g.R_in_cyl);
					//system("PAUSE");
					
				}

				if (b[i].g.itypegeom == CAD_STL) {
					// CAD_STL

					//printf("i=%d CAD\n", i);
					//system("pause");

					// Считывание бинарного STL файла.
					b[i].g.ReadSTL_binary(lb, b[0].g);

				}

				doublereal vol_poly_memo=0.0;
				if (b[i].g.itypegeom == POLYGON) {
					vol_poly_memo = Volume_polygon(b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2);
				}

				if ((i != 0)&&((b[i].g.itypegeom == PRISM)||
					((b[i].g.itypegeom == CYLINDER)&&(fabs(b[i].g.R_in_cyl)<1.0e-30)) ||
					((b[i].g.itypegeom == POLYGON)&&(b[i].g.nsizei == 4)))) 
				{
					// 15.11.2020
					// Преобразование в CAD STL для тестирования.

					//b[i].g.PRISM2CAD_STL_forTEST();
				}


				// Ввод предполагается корректным.
				// emissivity
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissW");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissW = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissW = 0.0; // изоляция.
					printf("b[%d].radiation.emissW =%e\n", i, b[i].radiation.emissW);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissE");				
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissE = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissE = 0.0; // изоляция.
					printf("b[%d].radiation.emissE =%e\n", i, b[i].radiation.emissE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissS = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissS = 0.0; // изоляция.
					printf("b[%d].radiation.emissS =%e\n", i, b[i].radiation.emissS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissN");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissN = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissN = 0.0; // изоляция.
					printf("b[%d].radiation.emissN  =%e\n", i, b[i].radiation.emissN);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissB");			
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissB = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissB = 0.0; // изоляция.
					printf("b[%d].radiation.emissB =%e\n", i, b[i].radiation.emissB);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissT");			
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissT = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissT = 0.0; // изоляция.
					printf("b[%d].radiation.emissT =%e\n",i, b[i].radiation.emissT);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "binternalRadiation");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din == 0) {
						// Блок не является вакуумным промежутком.
						b[i].radiation.binternalRadiation = false;
					}
					else {
						// блок является вакуумным промежутком.
						b[i].radiation.binternalRadiation = true;
						if (bvacuumPrism) {
							bdouble_vacuum_PRISM = true;
						}
						bvacuumPrism = true;
						// Вычисляем View Factors.
						calculate_view_factors(b[i]);
					}
					//printf(" b[%lld].radiation.binternalRadiation=%lld\n",i,din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.binternalRadiation = false; // Блок не является вакуумным промежутком.
					printf("b[%d].radiation.binternalRadiation =false\n", i);
					if (bSTOP_Reading) system("pause");
				}
				b[i].radiation.nodelistW = nullptr;
				b[i].radiation.nodelistE = nullptr;
				b[i].radiation.nodelistS = nullptr;
				b[i].radiation.nodelistN = nullptr;
				b[i].radiation.nodelistB = nullptr;
				b[i].radiation.nodelistT = nullptr;
				b[i].radiation.nodelistWsize = 0;
				b[i].radiation.nodelistEsize = 0;
				b[i].radiation.nodelistSsize = 0;
				b[i].radiation.nodelistNsize = 0;
				b[i].radiation.nodelistBsize = 0;
				b[i].radiation.nodelistTsize = 0;

				// идентификатор материала в базе материалов
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "imatid");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].imatid = static_cast<integer>(idin);
					//printf(" b[%lld].imatid=%lld\n", i, b[i].imatid);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].imatid = 0; // идентификатор материала блока.
					printf(" b[%d].imatid=%lld\n", i, b[i].imatid);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bCylinderFixed");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					// bCylinderFixed
					if (din == 1) {
						b[i].CylinderFixed = true;
					}
					else {
						b[i].CylinderFixed = false;
					}
					//printf("b[%lld].CylinderFixed =%lld\n", i, din);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].CylinderFixed = false; // .
					printf("b[%d].CylinderFixed = false\n", i);
					if (bSTOP_Reading) system("pause");
				}

				// мощность тепловыделения
				//fscanf_s(fp, "%f", &fin);
				//b[i].Sc = fin;
				// 19 november 2016 температурно зависимая мощность тепловыделения.
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "n_power");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].n_Sc = static_cast<integer>(idin);
					//printf("b[%lld].n_Sc =%lld\n",i,b[i].n_Sc );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].n_Sc = 1; // .
					printf(" b[%d].n_Sc=%lld\n", i, b[i].n_Sc);
					if (bSTOP_Reading) system("pause");
				}

				//b[i].arr_Sc = nullptr;
				//b[i].temp_Sc = nullptr;
				delete[] b[i].arr_Sc;
				delete[] b[i].temp_Sc;
				b[i].arr_Sc = new doublereal[b[i].n_Sc];
				b[i].temp_Sc = new doublereal[b[i].n_Sc];
				if (b[i].temp_Sc == nullptr) {
					printf("problem memory allocation for temp_Sc\n");
					system("pause");
					exit(1);
				}
				if (b[i].arr_Sc == nullptr) {
					printf("problem memory allocation for arr_Sc\n");
					system("pause");
					exit(1);
				}
				// Объём полигона.
				doublereal vol_poly;
				if (b[i].g.itypegeom == POLYGON) {
					vol_poly = Volume_polygon(b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2);
				}
				if (b[i].g.itypegeom == CAD_STL)
				{
					// Объём CAD STL геометрии
					vol_poly = b[i].g.volume_CAD_STL();
					doublereal vol_prism = (b[i].g.xE - b[i].g.xS) * (b[i].g.yE - b[i].g.yS) * (b[i].g.zE - b[i].g.zS);
					if (vol_poly_memo > 1.0e-30) {
						vol_prism = vol_poly_memo;
						if (fabs((vol_poly / vol_prism) - 1.0) < 0.1) {
							iblock_target++;
						}
						char ch = 'n';
						if (b[i].g.CAD_is_PRISM()) {
							ch = 'y';
						}
						std::wcout << ch << "id= " << i << "CAD *.stl obj " << b[i].name << "cad vol = " << vol_poly << "polygon vol =" << vol_prism << "  ratio = " << (vol_poly / vol_prism) << std::endl;
					}
					else {
						if (fabs((vol_poly / vol_prism) - 1.0) < 0.1)
						{
							iblock_target++;
						}
						char ch = 'n';
						if (b[i].g.CAD_is_PRISM()) {
							ch = 'y';
						}
						std::wcout << ch << "id= " << i << "CAD *.stl obj " << b[i].name << "cad vol = " << vol_poly << "prism vol =" << vol_prism << "  ratio = " << (vol_poly / vol_prism) << std::endl;
					}
				}

				for (int i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
					// Температура в C.
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "temp_power");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].temp_Sc[i_4] = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].temp_Sc[i_4] = 30.0; // .
						printf("b[%d].temp_Sc[%d] =%e\n", i, i_4, b[i].temp_Sc[i_4]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "arr_power");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);
										
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						if (fin != fin) {
							b[i].arr_Sc[i_4] = 0.0;
						}
						else {

							if (b[i].g.itypegeom == POLYGON) {
								// Для полигона передается из интерфейса просто мощность, а не удельная мощность.
								// Т.к. интерфейс не содержит функцию расчёта объёма полигона.
								// Для единообразия здесь мощность преобразуется в удельную мощность.

								if (vol_poly > 1.0e-30) {

									b[i].arr_Sc[i_4] = static_cast<doublereal>(fin) / vol_poly;

								}
								else {
									printf("error zero volume in polygon %s...\n", b[i].name);
									system("PAUSE");
									exit(1);
								}
							}
							else {
								b[i].arr_Sc[i_4] = static_cast<doublereal>(fin);
							}
						}
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].arr_Sc[i_4] = 0.0; // .
						printf("b[%d].arr_Sc[%d] =%e\n",i,i_4, b[i].arr_Sc[i_4]);
						if (bSTOP_Reading) system("pause");
					}
				}

				// debug
				//if (fabs(b[i].Sc)>1.0e-30) {
				//printf("%e\n", b[i].Sc);
				//system("PAUSE");
				//}
				// стиль зависимости мощности тепловыделения в блоке от времени.
				// 0 - не зависит от времени, 1 - square wave зависимость, 
				// 2 - square wave #2 зависимость, 3 - hot cold режим,
				// 4 - piecewise const.
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ipower_time_depend");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					switch (idin) {
					case 0: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
						break;
					case 1: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE;
						break;
					case 2: b[i].ipower_time_depend = POWER_TIME_DEPEND::SQUARE_WAVE2;
						break;
					case 3: b[i].ipower_time_depend = POWER_TIME_DEPEND::HOT_COLD;
						break;
					case 4: b[i].ipower_time_depend = POWER_TIME_DEPEND::PIECEWISE_CONST;
						break;
					default: b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER;
						break;
					}
					//printf("b[%lld].ipower_time_depend =%lld\n", i, b[i].ipower_time_depend);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].ipower_time_depend = POWER_TIME_DEPEND::CONST_POWER; // не зависит от времени.
					printf("b[%d].ipower_time_depend =POWER_TIME_DEPEND::CONST\n", i);
					if (bSTOP_Reading) system("pause");
				}
				// тип блока
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "itype");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					
					switch (idin) {
					case 1: b[i].itype = PHYSICS_TYPE_IN_BODY::SOLID;
						break;
					case 2: b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
						break;
					case 3: b[i].itype = PHYSICS_TYPE_IN_BODY::FLUID;
						break;
					default:
						b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW;
						break;
					}
					//printf("b[%lld].itype  =%lld\n",i,b[i].itype  );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].itype = PHYSICS_TYPE_IN_BODY::HOLLOW; // 1 - SOLID, 2 - HOLLOW, 3 - FLUID.
					printf(" b[%d].itype =PHYSICS_TYPE_IN_BODY::HOLLOW\n", i);
					if (bSTOP_Reading) system("pause");
				}

				// печать считанных значений на консоль
				//printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				//printf("%lld %lld\n", b[i].imatid, b[i].itype);
				//printf("temperature depend power\n");
				//printf("t_C power_W\n");
				//for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					//printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
				//}
			}

			std::cout << "iblock_target=" << iblock_target << " lb=" << lb << std::endl;
			//system("pause");
			// считывание источников тепла
			for (i = 0; i < ls; ++i) {

				char name0[1000] = "source";
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "name");

				if (smakesource(name0, s[i].name)) {
					// Найдено успешно.

					//printf(" s[%lld].name=%s\n", i, s[i].name);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].name[0] = '\0'; // уникальное текстовое имя источника тепла.
					printf(" s[%d].name=unknown\n", i);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
				strcat_s(name0, "iunion");			
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					s[i].iunion_id = static_cast<integer>(idin);
					//printf(" s[%lld].iunion_id=%lld\n",i,s[i].iunion_id );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].iunion_id = 0; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					printf(" s[%d].iunion_id=%lld\n",i, s[i].iunion_id);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Power");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].power = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].power = 0.0; // .
					printf(" s[%d].power=%e\n",i, s[i].power);
					if (bSTOP_Reading) system("pause");
				}


				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "itempdep");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din == 0) {
						// const
						// задаётся постоянное значение рассеиваемой в тепло мощности.
						s[i].power_multiplyer = 1.0;
						s[i].bgarber_depend = false;
					}
					else if (din == 1) {
						// рассеиваемая в тепло мощность задаётся таблично
						// в зависимости от максимальной температуры и рабочего 
						// значения смещения стока.
						s[i].bgarber_depend = true;
						s[i].power_multiplyer = s[i].power;
						// мощность будет вычислена при ambient Temperature несколькими строчками позже.
					}
					//printf(" s[%lld].bgarber_depend=%lld\n", i, din);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					// const
					// задаётся постоянное значение рассеиваемой в тепло мощности.
					s[i].power_multiplyer = 1.0;
					s[i].bgarber_depend = false; 
					printf("s[%d].bgarber_depend = false\n",i );
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "id_table");

				// уникальный номер таблицы.
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					s[i].igarber_depend = static_cast<integer>(idin);
					//printf(" s[%lld].igarber_depend=%lld\n", i, s[i].igarber_depend);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].igarber_depend = 0; // .
					printf("s[%d].igarber_depend =%lld\n", i, s[i].igarber_depend);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "operatingoffsetdrain");

				
				 // рабочее значение смещения стока.
				//printf("offset drain is %e\n",s[i].roperation_offset_drain);
				//system("PAUSE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].roperation_offset_drain = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].roperation_offset_drain = 28.0; // .
					printf("s[%d].roperation_offset_drain =%e\n",i, s[i].roperation_offset_drain);
					if (bSTOP_Reading) system("pause");
				}

				bool bsplinereadOk = true;
				if (s[i].bgarber_depend) {
					s[i].power = my_splain_interpol_power_table(gtdps[s[i].igarber_depend].intemp,
						gtdps[s[i].igarber_depend].inoffset_drain,
						gtdps[s[i].igarber_depend].rtemp,
						gtdps[s[i].igarber_depend].roffset_drain,
						gtdps[s[i].igarber_depend].rpower_table,
						operatingtemperature,
						s[i].roperation_offset_drain);
					if (bsplinereadOk) {
						// одиночная тестовая проверка сплайновой аппроксимации.
						printf("single test validation spline approximation...\n");
						printf("calculate initial power=%e\n", s[i].power);
						printf("please, press any key to continue...");
						// system("PAUSE");
						system("pause");

						bsplinereadOk = false;
					}
					s[i].power *= s[i].power_multiplyer; // домножение на корректирующий множитель.
				}


				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iPlane");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					s[i].iPlane = static_cast<integer>(idin);
					//printf(" s[%lld].iPlane=%lld\n",i,s[i].iPlane );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].iPlane = 0; // .
					printf("s[%d].iPlane =%lld\n",i, s[i].iPlane);
					if (bSTOP_Reading) system("pause");
				}
				// геометрия
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.xS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.xS = scale * 0.0; // .
					printf(" s[%d].g.xS=%e\n",i, s[i].g.xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.yS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.yS = scale * 0.0; // .
					printf("s[%d].g.yS  =%e\n", i, s[i].g.yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.zS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.zS = scale * 0.0; // .
					printf("s[%d].g.zS =%e\n", i, s[i].g.zS);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.xE  = scale *  static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.xE = scale *  0.0; // .
					printf(" s[%d].g.xE=%e\n",i, s[i].g.xE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.yE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.yE = scale *  0.0; // .
					printf(" s[%d].g.yE =%e\n",i, s[i].g.yE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.zE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.zE = scale * 0.0; // .
					printf(" s[%d].g.zE=%e\n",i, s[i].g.zE);
					if (bSTOP_Reading) system("pause");
				}


				// swap
				if (s[i].g.xS > s[i].g.xE) {
					dbuf = s[i].g.xS;
					s[i].g.xS = s[i].g.xE;
					s[i].g.xE = dbuf;
				}
				if (s[i].g.yS > s[i].g.yE) {
					dbuf = s[i].g.yS;
					s[i].g.yS = s[i].g.yE;
					s[i].g.yE = dbuf;
				}
				if (s[i].g.zS > s[i].g.zE) {
					dbuf = s[i].g.zS;
					s[i].g.zS = s[i].g.zE;
					s[i].g.zE = dbuf;
				}
				switch (s[i].iPlane) {
				case XY_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ_PLANE: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ_PLANE: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				//printf("source %e %lld %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
			}

			// 19.01.2021
		// Стенка заглушка с условием прилипания по скорости(нулевой вектор скорости), 
		// нулевым тепловым потоком по температуре, 
		// свободной free границе в механике без приложенной силы.
			w[lw].iunion_id = 0; // cabinet
			w[lw].ifamily = WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY;
			w[lw].Tamb = 0.0;
			w[lw].emissivity = 0.0;
			w[lw].film_coefficient = 0.0;
			w[lw].ViewFactor = 1.0;
			w[lw].hf = 0.0;
			w[lw].bsymmetry = false;
			w[lw].bpressure = false;
			w[lw].bopening = false;
			w[lw].Vx = 0.0;
			w[lw].Vy = 0.0;
			w[lw].Vz = 0.0;
			w[lw].P = 0.0;
			w[lw].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
			w[lw].xForce = 0.0;
			w[lw].yForce = 0.0;
			w[lw].zForce = 0.0;
			w[lw].iPlane = XY_PLANE;
			w[lw].g.xS = 0.0;
			w[lw].g.yS = 0.0;
			w[lw].g.zS = 0.0;
			w[lw].g.xE = 0.0;
			w[lw].g.yE = 0.0;
			w[lw].g.zE = 0.0;

			// считывание твёрдых стенок
			for (i = 0; i < lw; ++i) {

				char name0[1000] = "wall";
				buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
				strcat_s(name0, "name");

				if (smakesource(name0, w[i].name)) {
					// Найдено успешно.
					
					//printf(" w[%lld].name=%s\n", i, w[i].name);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].name[0] = '\0'; // уникальное текстовое имя твёрдой стенки.
					printf(" w[%d].name=unknown\n", i);
					//if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iunion");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					w[i].iunion_id = static_cast<integer>(idin);
					//printf(" w[%lld].iunion_id=%lld\n", i, w[i].iunion_id);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].iunion_id = 0; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					printf(" w[%d].iunion_id=%lld\n",i, w[i].iunion_id);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "family");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					switch (idin) {
						case 1 : w[i].ifamily= WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY; break;
						case 2: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
						case 3: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY; break;
						case 4: w[i].ifamily= WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY; break;
						default: w[i].ifamily= WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; break;
					}
					//w[i].ifamily = static_cast<integer>(idin);
					din= static_cast<integer>(idin); // for switch() см. далее.
					//printf("w[%lld].ifamily =%lld\n", i, w[i].ifamily);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].ifamily = WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY; // Однородное условие Неймана.
					din=2;//NEIMAN_FAMILY
					printf("w[%d].ifamily =%d\n", i, w[i].ifamily);
					if (bSTOP_Reading) system("pause");
				}

				switch (din) {
				case 1:  
					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "Tamb");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].Tamb = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%d].Tamb =%e\n", i, w[i].Tamb);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "heat_transfer_coefficient_vs_emissivity");

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "ViewFactor");

					// Stefan Bolcman
					// termostability wall
					w[i].emissivity = 0.0;
					w[i].ViewFactor = 1.0; // фактор видимости.
					w[i].film_coefficient = 0.0;					

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "HF");
					
					w[i].hf = 0.0;
					
					break; // первого рода
				case 2:  
					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "Tamb");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].Tamb = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%d].Tamb =%e\n", i, w[i].Tamb);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "heat_transfer_coefficient_vs_emissivity");

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "ViewFactor");

					// Stefan Bolcman
					// adiabatic wall
					w[i].emissivity = 0.0;
					w[i].ViewFactor = 1.0; // фактор видимости.
					w[i].film_coefficient = 0.0;
					

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "HF");

					
					w[i].hf = 0.0;
					
					break; // однородное условие Неймана
				case 3:  
					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "Tamb");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].Tamb = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%d].Tamb =%e\n", i, w[i].Tamb);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "heat_transfer_coefficient_vs_emissivity");

					// Stefan Bolcman
					 // Newton-Richman condition, film coefficient.
					w[i].emissivity = 0.0;
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						printf("film_coefficient=%e\n", fin);
						w[i].film_coefficient = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].film_coefficient = 0.0; // .
						printf(" w[%d].film_coefficient=%e\n",i, w[i].film_coefficient);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "ViewFactor");

					w[i].ViewFactor = 1.0; // фактор видимости.

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "HF");

					
					w[i].hf = 0.0;
					
					break; // Ньютон-Рихман.
				case 4:  
					
					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "Tamb");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].Tamb = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%d].Tamb =%e\n", i, w[i].Tamb);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "heat_transfer_coefficient_vs_emissivity");

					// Stefan Bolcman
					// Stefan - Bolcman condition
					w[i].film_coefficient = 0.0;
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].emissivity = static_cast<doublereal>(fin); // Коэффициент излучения==поглощения.
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].emissivity = 0.0; // Коэффициент излучения==поглощения.
						printf("w[%d].emissivity =%e\n",i, w[i].emissivity);
						if (bSTOP_Reading) system("pause");
					}


					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "ViewFactor");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].ViewFactor = static_cast<doublereal>(fin);
						w[i].ViewFactor = fmax(0.0, fmin(w[i].ViewFactor, 1.0));
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].ViewFactor  = 1.0; // фактор видимости.
						printf("w[%d].ViewFactor =%e\n", i, w[i].ViewFactor);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "HF");
										
					w[i].hf = 0.0;
					
					break; // Стефан-Больцман.
				default:
					printf("error: wall unlnown boundary condition type.\n");
					system("PAUSE");
					break;
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bsymmetry");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din == 1) w[i].bsymmetry = true; else w[i].bsymmetry = false;
					//printf("w[%lld].bsymmetry =%lld\n",i, din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].bsymmetry = false; // .
					printf(" w[%d].bsymmetry = false\n", i);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bpressure");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din == 1) w[i].bpressure = true; else w[i].bpressure = false;
					//printf(" w[%lld].bpressure=%lld\n", i,din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].bpressure = false; // .
					printf("w[%d].bpressure = false\n",i );
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bopening");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din == 1) w[i].bopening = true; else w[i].bopening = false;
					//printf("  w[i].bopening =%lld\n",i, din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].bopening = false; // .
					printf(" w[%d].bopening = false\n",i );
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Vx");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].Vx = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].Vx = 0.0; // .
					printf(" w[%d].Vx =%e\n", i, w[i].Vx);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Vy");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].Vy = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].Vy = 0.0; // .
					printf("w[%d].Vy =%e\n", i, w[i].Vy);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Vz");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].Vz = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].Vz = 0.0; // .
					printf("w[%d].Vz =%e\n", i, w[i].Vz);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "P");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].P = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].P  = 0.0; // .
					printf("w[%d].P =%e\n",i, w[i].P);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ithermal_stress_boundary_condition");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);

					if (din >= 0 && din < 11) {
						// 0- FREE
						
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
						switch (din) {
						case 0: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE;
							break;
						case 1: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FIXIT;
							break;
						case 2: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FIXIT;
							break;
						case 3: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FIXIT;
							break;
						case 4: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XY_FIXIT;
							break;
						case 5: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::XZ_FIXIT;
							break;
						case 6: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::YZ_FIXIT;
							break;
						case 7: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::ALL_FIXIT;
							break;
						case 8: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::X_FORCE;
							break;
						case 9: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Y_FORCE;
							break;
						case 10: w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::Z_FORCE;
							break;
						default:
							w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // Free all
							break;
						}
					}
					else {
						printf("error: unknown ithermal_Stress_boundary_condition\n");
						printf("ithermal_Stress_boundary_condition=%lld\n", din);
						system("PAUSE");
						w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // Free all
					}
					//printf("w[%lld].ithermal_Stress_boundary_condition =%lld\n",i,din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].ithermal_Stress_boundary_condition = THERMAL_STRESS_BOUNDARY_CONDITION::FREE; // 0 - FREE.
					printf(" w[%d].ithermal_Stress_boundary_condition=FREE BOUNDARY\n",i );
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xForce");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].xForce = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].xForce = 0.0; // .
					printf("w[%d].xForce =%e\n",i, w[i].xForce);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yForce");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].yForce = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].yForce = 0.0; // .
					printf("w[%d].yForce=%e\n", i, w[i].yForce);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zForce");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].zForce = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].zForce = 0.0; // .
					printf(" w[%d].zForce=%e\n",i, w[i].zForce);
					if (bSTOP_Reading) system("pause");
				}				
				//printf("Force Fx=%e Fy=%e Fz=%e\n", w[i].xForce, w[i].yForce, w[i].zForce);
				//system("PAUSE");

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iPlane");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					w[i].iPlane = static_cast<integer>(idin);
					//printf("w[%lld].iPlane =%lld\n", i, w[i].iPlane);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].iPlane = 0; // .
					printf(" w[%d].iPlane=%lld\n", i, w[i].iPlane);
					if (bSTOP_Reading) system("pause");
				}
				// геометрия
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.xS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.xS = scale * 0.0; // .
					printf("w[%d].g.xS =%e\n", i, w[i].g.xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.yS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.yS = scale * 0.0; // .
					printf(" w[%d].g.yS=%e\n", i, w[i].g.yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.zS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.zS = scale * 0.0; // .
					printf(" w[%d].g.zS  =%e\n", i, w[i].g.zS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.xE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.xE = scale * 0.0; // .
					printf("w[%d].g.xE  =%e\n", i, w[i].g.xE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.yE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.yE = scale * 0.0; // .
					printf("w[%d].g.yE =%e\n", i, w[i].g.yE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.zE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.zE = scale * 0.0; // .
					printf("w[%d].g.zE =%e\n", i, w[i].g.zE);
					if (bSTOP_Reading) system("pause");
				}
				// swap
				if (w[i].g.xS > w[i].g.xE) {
					dbuf = w[i].g.xS;
					w[i].g.xS = w[i].g.xE;
					w[i].g.xE = dbuf;
				}
				if (w[i].g.yS > w[i].g.yE) {
					dbuf = w[i].g.yS;
					w[i].g.yS = w[i].g.yE;
					w[i].g.yE = dbuf;
				}
				if (w[i].g.zS > w[i].g.zE) {
					dbuf = w[i].g.zS;
					w[i].g.zS = w[i].g.zE;
					w[i].g.zE = dbuf;
				}
				//printf("wall %lld %e %e %lld %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
			}


			// АСЕМБЛЕСЫ.		
			int lu_loc;
			if (imakesource("iu_number", idin)) {
				// Найдено успешно.
				lu = static_cast<integer>(idin);

				lu_loc = lu;// Количество активных асемблесов.
				// Так как могут быть многие неактивные асемблесы то номер асемблеса может быть значительно выше.
				lu = 0;
				// Поиск максимального номера юниона который активен. 23.10.2021
				for (int iuid = 99; iuid >= 0; iuid--) {


					char name0[1000] = "assembles";
					buffer[0] = '\0'; _itoa_s(iuid, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "xS");

				
					if (fmakesource(name0, fin)) {
						lu = iuid+1;
						break;
					}			
					
				}
				//printf(" =%lld\n", );
			}
			else {
				printf("WARNING!!! iu_number not found in file premeshin.txt\n");
				lu = 0; // .
				printf(" lu=%d\n",lu );
				if (bSTOP_Reading) system("pause");
			}
			if (lu == 0) {
				my_union = nullptr;
			}
			else {
				my_union = new UNION[lu];
				// инициализация.
				for (i = 0; i < lu; ++i) {
					my_union[i].f = nullptr;
					my_union[i].xpos = nullptr;
					my_union[i].ypos = nullptr;
					my_union[i].zpos = nullptr;
					my_union[i].xposadd = nullptr;
					my_union[i].yposadd = nullptr;
					my_union[i].zposadd = nullptr;
					my_union[i].iunion_parent = -1; // Cabinet
					my_union[i].iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER; // 2 - CoarseMeshGen
					my_union[i].inxadd = -1;
					my_union[i].inyadd = -1;
					my_union[i].inzadd = -1;
					my_union[i].flow_interior = 0;
					my_union[i].active = false;
				}
			}

			i = 0;
			int lu_found = -1;
			for (int i_72 = 0; i_72 < lu_loc; i_72++) {

				// Так как могут быть многие неактивные асемблесы то номер асемблеса может быть значительно выше.
				i = lu_found + 1;
				// Поиск максимального номера юниона который активен. 23.10.2021
				for (int iuid = i; iuid <= 99; iuid++) {


					char name0[1000] = "assembles";
					buffer[0] = '\0'; _itoa_s(iuid, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "xS");


					if (fmakesource(name0, fin)) {
						lu_found = iuid;
						break;
					}

				}

				// он обязательно будет найден так как количество активных юнионов заранее известно и прописано в файле premeshin.
				i = lu_found;

				my_union[i].active = true;


				char name0[1000] = "assembles";
				buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
				strcat_s(name0, "name");

				if (smakesource(name0, my_union[i].name)) {
					// Найдено успешно.

					//printf(" my_union[%lld].name=%s\n", i, my_union[i].name);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].name[0] = '\0'; // уникальное текстовое имя твёрдой стенки.
					printf(" my_union[%d].name=unknown\n", i);
					//if (bSTOP_Reading) system("pause");
				}




				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
				strcat_s(name0, "iunionparent");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = static_cast<integer>(idin);
					if (din >=-1) my_union[i].iunion_parent = din; else my_union[i].iunion_parent = -1; // Cabinet
					//printf("my_union[%lld].iunion_parent =%lld\n",i, din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].iunion_parent = -1; // .
					printf(" my_union[%d].iunion_parent = -1 (Cabinet)\n", i);
					if (bSTOP_Reading) system("pause");
				}


				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].xS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].xS = scale * 0.0; // .
					printf(" my_union[%d].xS=%e\n", i, my_union[i].xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].xE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].xE = scale * 0.0; // .
					printf(" my_union[%d].xE=%e\n",i, my_union[i].xE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].yS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].yS = scale * 0.0; // .
					printf(" my_union[%d].yS=%e\n", i, my_union[i].yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].yE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].yE = scale * 0.0; // .
					printf("my_union[%d].yE =%e\n", i, my_union[i].yE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].zS = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].zS = scale * 0.0; // .
					printf(" my_union[%d].zS=%e\n",i, my_union[i].zS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].zE = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].zE = scale * 0.0; // .
					printf(" my_union[%d].zE=%e\n", i, my_union[i].zE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "identifire");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].id  = static_cast<integer>(idin);
					//printf(" my_union[%lld].id=%lld\n", i, my_union[i].id);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].id = static_cast<integer>(i+1); // Уникальный идентификатор АССЕМБЛЕСА.
					printf(" my_union[%d].id=%lld\n",i, my_union[i].id);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "inx");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].inx = static_cast<integer>(idin);
					//printf(" my_union[%lld].inx=%lld\n", i, my_union[i].inx);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].inx = 23; // .
					printf(" my_union[%d].inx=%d\n", i, my_union[i].inx);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iny");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].iny = static_cast<integer>(idin);
					//printf(" my_union[%lld].iny=%lld\n", i, my_union[i].iny);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].iny = 23; // .
					printf("my_union[%d].iny =%d\n", i, my_union[i].iny);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "inz");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].inz = static_cast<integer>(idin);
					//printf(" my_union[%lld].inz =%lld\n", i, my_union[i].inz );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].inz = 23; // .
					printf(" my_union[%d].inz =%d\n", i, my_union[i].inz);
					if (bSTOP_Reading) system("pause");
				}
			}

			// считывание информации о наборе решаемых уравнений			
			
			if (imakesource("egddata_itemper", idin)) {
				// Найдено успешно.
				// 0 - none;
				// 1 - метод контрольного объема;
				// 2 - метод конечных элементов.
				if ((idin >= 0) && (idin <= 2)) {
					eqin.itemper = static_cast<integer>(idin);
				}
				else {
					eqin.itemper = 1; // Решаем уравнение теплопередачи методом контрольного объема.
				}
				//printf(" eqin.itemper=%lld\n",eqin.itemper );
			}
			else {
				printf("WARNING!!! egddata_itemper not found in file premeshin.txt\n");
				eqin.itemper = 1; // Решаем уравнение теплопередачи методом контрольного объема.
				printf(" eqin.itemper=%lld\n", eqin.itemper);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("egddata_imaxflD", idin)) {
				// Найдено успешно.
				eqin.imaxflD = static_cast<integer>(idin);
				//printf(" eqin.imaxflD=%lld\n", eqin.imaxflD);
			}
			else {
				printf("WARNING!!! egddata_imaxflD not found in file premeshin.txt\n");
				eqin.imaxflD = 1; // По умолчанию у нас всегда есть одна общая гидродинамическая область.
				printf(" eqin.imaxflD=%lld\n", eqin.imaxflD);
				if (bSTOP_Reading) system("pause");
			}
			//printf("itemper=%lld eqin.imaxflD=%lld\n", eqin.itemper, eqin.imaxflD);
			//system("PAUSE");
			if (eqin.imaxflD == 0) {
				eqin.fluidinfo = nullptr;
			}
			else
			{
				// выделение оперативной памяти
				if (eqin.fluidinfo != nullptr) {
					delete eqin.fluidinfo;
					eqin.fluidinfo = nullptr;
				}
				eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
				for (i = 0; i < eqin.imaxflD; ++i) {
					// Считывание координат опорной точки
					char name0[1000] = "egddata";
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "xc");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].xc = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].xc = scale * 0.0; // .
						printf(" eqin.fluidinfo[%d].xc=%e\n", i, eqin.fluidinfo[i].xc);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "yc");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].yc = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].yc = scale * 0.0; // .
						printf(" eqin.fluidinfo[%d].yc=%e\n", i, eqin.fluidinfo[i].yc);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "zc");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].zc = scale * static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].zc = scale * 0.0; // .
						printf(" eqin.fluidinfo[%d].zc=%e\n", i, eqin.fluidinfo[i].zc);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "iflow");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].iflow = static_cast<integer>(idin);
						//printf("eqin.fluidinfo[%lld].iflow =%lld\n", i, eqin.fluidinfo[i].iflow);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].iflow = 0; // .
						printf(" eqin.fluidinfo[%d].iflow=%lld\n", i, eqin.fluidinfo[i].iflow);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "iflowregime");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						switch (idin) {
						case 0: eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
							break;
						case 1:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::TURBULENT;
							break;
						default:  eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR;
							break;
						}
						
						//printf(" eqin.fluidinfo[%lld].iflowregime=%lld\n",i,eqin.fluidinfo[i].iflowregime );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].iflowregime = FLOW_REGIME::LAMINAR; // Laminar.
						printf(" eqin.fluidinfo[%d].iflowregime=FLOW_REGIME::LAMINAR\n", i);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "iturbmodel");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						// Режим течения: ламинарный или конкретная модель турбулентности.
						switch (idin) {
						case 0: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
							break;
						case 1: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::SMAGORINSKY;
							break;
						case 2: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RNG_LES;
							break;
						case 3: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_SPALART_ALLMARES;
							break;
						case 4: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_MENTER_SST;
							break;
						case 5: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_STANDART_K_EPS;
							break;
						case 6: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::RANS_LANGTRY_MENTOR_SST;
							break;
						default: eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD;
							break;
						}
						//printf(" eqin.fluidinfo[%lld].iturbmodel=%lld\n", i, eqin.fluidinfo[i].iturbmodel);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].iturbmodel = TURBULENT_MODEL::ZEROEQMOD; // Zero Equation Turbulence Model.
						printf(" eqin.fluidinfo[%d].iturbmodel=VISCOSITY_MODEL::LAMINAR\n", i);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "SmagConst");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].Cs = static_cast<doublereal>(fin);// постоянная Смагоринского.
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].Cs = 0.151; // default 0.151.
						printf(" eqin.fluidinfo[%d].Cs=%e\n", i, eqin.fluidinfo[i].Cs);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "iDynamicStressGermano");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din = static_cast<integer>(idin);
						// учёт динамического определения квадрата постоянной Смагоринского.
						// включает Dynamic Subgrid Scale Model Германо 1991 года.
						if (din == 1) {
							eqin.fluidinfo[i].bDynamic_Stress = true;
						}
						else {
							eqin.fluidinfo[i].bDynamic_Stress = false;
						}
						//printf("eqin.fluidinfo[%lld].bDynamic_Stress  =%lld\n",i, din );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bDynamic_Stress = false; // .
						printf(" eqin.fluidinfo[%d].bDynamic_Stress = false\n", i);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "iLimitersCs");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din = static_cast<integer>(idin);
						// включает ограничение сверху и снизу на возможные значения постоянной Смагоринского.
						if (din == 1) {
							eqin.fluidinfo[i].bLimiters_Cs = true;
						}
						else {
							eqin.fluidinfo[i].bLimiters_Cs = false;
						}
						//printf("eqin.fluidinfo[%lld].bLimiters_Cs =%lld\n", i, din);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bLimiters_Cs = false; // .
						printf(" eqin.fluidinfo[%d].bLimiters_Cs = false\n", i);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "minCs");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].rminCs = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rminCs = -1.0e20; // минимальное возможное значение постоянной Смагоринского.
						printf("eqin.fluidinfo[%d].rminCs =%e\n", i, eqin.fluidinfo[i].rminCs);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "maxCs");


					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].rmaxCs = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rmaxCs = 1.0e23; // максимальное возможное значение постоянной Смагоринского.
						printf(" eqin.fluidinfo[%d].rmaxCs=%e\n", i, eqin.fluidinfo[i].rmaxCs);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "itypeFiltrGermano");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].itypeFiltrGermano = static_cast<integer>(idin);
						//printf(" eqin.fluidinfo[%lld].itypeFiltrGermano=%lld\n",i,eqin.fluidinfo[i].itypeFiltrGermano );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].itypeFiltrGermano = 2; // тип фильтра в модели Германо 1991 года. 
						printf(" eqin.fluidinfo[%d].itypeFiltrGermano=%lld\n", i, eqin.fluidinfo[i].itypeFiltrGermano);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "roughness");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].roughness = 1.0e-6 * static_cast<doublereal>(fin); // шероховатость стенки в м.
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].roughness = 1.0e-6 * 10.0; // шероховатость стенки в м.
						printf(" eqin.fluidinfo[%d].roughness=%e\n", i, eqin.fluidinfo[i].roughness);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "rRimult");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						// множитель корректирующий турбулентное число Ричардсона.
						eqin.fluidinfo[i].rRi_mult = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rRi_mult = 1.0; // множитель корректирующий турбулентное число Ричардсона.
						printf(" eqin.fluidinfo[%d].rRi_mult=%e\n", i, eqin.fluidinfo[i].rRi_mult);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "rSelectiveAngle");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].rSelectiveAngle = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rSelectiveAngle = 15.0; // пороговое значение угла в модели Selective Smagorinsky.
						printf(" eqin.fluidinfo[%d].rSelectiveAngle=%e\n", i, eqin.fluidinfo[i].rSelectiveAngle);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "ipowerroughness");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].ipowerroughness = static_cast<integer>(idin);
						//printf(" eqin.fluidinfo[%lld].ipowerroughness=%lld\n", i, eqin.fluidinfo[i].ipowerroughness);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].ipowerroughness = 2; // показатель степени в модели учёта шероховатости.
						printf(" eqin.fluidinfo[%d].ipowerroughness=%lld\n", i, eqin.fluidinfo[i].ipowerroughness);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "itypefiltr");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = static_cast<integer>(idin);
						//printf(" eqin.fluidinfo[%lld].itypeSelectiveSmagorinsky_filtr=%lld\n", i, eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = 2; // тип фильтра в модели Selective Smagorinsky.
						printf(" eqin.fluidinfo[%d].itypeSelectiveSmagorinsky_filtr=%lld\n", i, eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "bfdelta");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].bfdelta = static_cast<integer>(idin);
						//printf(" eqin.fluidinfo[%lld].bfdelta=%lld\n", i, eqin.fluidinfo[i].bfdelta);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bfdelta = 1; // учёт неравномерности сетки.
						printf(" eqin.fluidinfo[%d].bfdelta=%lld\n", i, eqin.fluidinfo[i].bfdelta);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "bSmagorinsky_Lilly");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].bSmagorinskyLilly = static_cast<integer>(idin);
						//printf(" eqin.fluidinfo[%lld].bSmagorinskyLilly=%lld\n", i, eqin.fluidinfo[i].bSmagorinskyLilly);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bSmagorinskyLilly = 1; // модель Смагоринского-Лиллу.
						printf(" eqin.fluidinfo[%d].bSmagorinskyLilly=%lld\n", i, eqin.fluidinfo[i].bSmagorinskyLilly);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "bsurface_roughness");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].bsurface_roughness = static_cast<integer>(idin);
						//printf(" eqin.fluidinfo[%lld].bsurface_roughness=%lld\n", i, eqin.fluidinfo[i].bsurface_roughness);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bsurface_roughness = 0; // учёт шероховатости стенки.
						printf(" eqin.fluidinfo[%d].bsurface_roughness=%lld\n", i, eqin.fluidinfo[i].bsurface_roughness);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "bSwirlamendment");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din = static_cast<integer>(idin);
						// учёт течений с кривизной линий тока.
						if (din == 1) {
							eqin.fluidinfo[i].bSwirlAmendment = true;
						}
						else {
							eqin.fluidinfo[i].bSwirlAmendment = false;
						}
						//printf(" eqin.fluidinfo[%lld].bSwirlAmendment=%lld\n",i,din );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bSwirlAmendment = true; // .
						printf(" eqin.fluidinfo[%d].bSwirlAmendment = true\n", i);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i, buffer, 10); strcat_s(name0, buffer);
					strcat_s(name0, "bSelectiveSmagorinsky");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din = static_cast<integer>(idin);
						// учёт избирательности в модели Смагоринского
						if (din == 1) {
							eqin.fluidinfo[i].bSelectiveSmagorinsky = true;
						}
						else {
							eqin.fluidinfo[i].bSelectiveSmagorinsky = false;
						}
						//printf(" eqin.fluidinfo[%lld].bSelectiveSmagorinsky=%lld\n",i,din );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bSelectiveSmagorinsky = false; // .
						printf(" eqin.fluidinfo[%d].bSelectiveSmagorinsky = false\n", i);
						if (bSTOP_Reading) system("pause");
					}

				}
				// Параметры преобразователя картинок для отчетов.
				// 5.01.2018

				if (fmakesource("report_min", fin)) {
					// Найдено успешно.
					pfpir.fminimum = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! report_min not found in file premeshin.txt\n");
					pfpir.fminimum = scale * 0.0; // .
					printf(" pfpir.fminimum =%e\n", pfpir.fminimum);
					if (bSTOP_Reading) system("pause");
				}

				if (fmakesource("report_max", fin)) {
					// Найдено успешно.
					pfpir.fmaximum = scale * static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! report_max not found in file premeshin.txt\n");
					pfpir.fmaximum = scale * 0.0; // .
					printf("pfpir.fmaximum =%e\n", pfpir.fmaximum);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("report_directional", idin)) {
					// Найдено успешно.
					switch (idin) {
					case 0: pfpir.idir = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; break;
					case 1: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					case 2: pfpir.idir = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL; break;
					default: pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; break;
					}
					//pfpir.idir = static_cast<integer>(idin);
					//printf(" pfpir.idir=%lld\n", pfpir.idir);
				}
				else {
					printf("WARNING!!! report_directional not found in file premeshin.txt\n");
					pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL;
					printf(" pfpir.idir=%d\n", pfpir.idir);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("amg1r6_checker", idin)) {
					// Найдено успешно.
					AMG1R6_LABEL = static_cast<integer>(idin);
					//printf(" AMG1R6_LABEL=%lld\n",AMG1R6_LABEL );
				}
				else {
					printf("WARNING!!! amg1r6_checker not found in file premeshin.txt\n");
					AMG1R6_LABEL = 0; // default amg1r5.
					printf(" AMG1R6_LABEL=%lld\n", AMG1R6_LABEL);
					if (bSTOP_Reading) system("pause");
				}


				if (imakesource("nrd", idin)) {
					// Найдено успешно.
					if (idin > 9999) {
						b_iluk_amg1r5_LABEL_D = true;
						nrd_LABEL = static_cast<integer>(idin-10000);
					}
					else {
						nrd_LABEL = static_cast<integer>(idin);
					}
					//printf(" nrd_LABEL=%lld\n",nrd_LABEL );
				}
				else {
					printf("WARNING!!! nrd not found in file premeshin.txt\n");
					nrd_LABEL = 1131; // default C/F relaxation.
					printf(" nrd_LABEL=%lld\n", nrd_LABEL);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("nru", idin)) {
					// Найдено успешно.
					if (idin > 9999) {
						b_iluk_amg1r5_LABEL_U = true;
						nru_LABEL = static_cast<integer>(idin - 10000);
					}
					else {
						nru_LABEL = static_cast<integer>(idin);
					}
					//printf(" nru_LABEL=%lld\n",nru_LABEL );
				}
				else {
					printf("WARNING!!! nru not found in file premeshin.txt\n");
					nru_LABEL = 1131; // default C/F relaxation.
					printf(" nru_LABEL=%lld\n", nru_LABEL);
					if (bSTOP_Reading) system("pause");
				}

				if (fmakesource("ecg2", fin)) {
					// Найдено успешно.
					ecg2_LABEL = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! ecg2 not found in file premeshin.txt\n");
					ecg2_LABEL = 0.25; // .
					printf("ecg2 or strong threshold amg1r5 =%e\n", ecg2_LABEL);
					if (bSTOP_Reading) system("pause");
				}

				if (fmakesource("ewt2", fin)) {
					// Найдено успешно.
					ewt2_LABEL = static_cast<doublereal>(fin);
				}
				else {
					printf("WARNING!!! ewt2 not found in file premeshin.txt\n");
					ewt2_LABEL = 0.35; // .
					printf("ewt2 or F to F amg1r5 =%e\n", ewt2_LABEL);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("number_processors", idin)) {
					// Найдено успешно.
					number_processors_global_var = (int)(idin);
					//printf(" number_processors_global_var=%lld\n",number_processors_global_var );
				}
				else {
					printf("WARNING!!! number_processors not found in file premeshin.txt\n");
					number_processors_global_var = 1; // default 1 thread.
					printf(" number_processors=%d\n", number_processors_global_var);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("gpu_id", idin)) {
					// Найдено успешно.
					idevice_Tesla = (int)(idin);
					//printf(" idevice_Tesla=%lld\n",idevice_Tesla );
				}
				else {
					printf("WARNING!!! idevice_Tesla not found in file premeshin.txt\n");
					idevice_Tesla = 0; // default gpu id = 0.
					printf(" idevice_Tesla=%d\n", idevice_Tesla);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("number_iterations_SIMPLE_algorithm", idin)) {
					// Найдено успешно.
					if (idin>=0) {
						number_iteration_SIMPLE_algorithm = (unsigned int)(idin);
					}
					else {
						number_iteration_SIMPLE_algorithm =0;
						std::cout<<"WARNING : number_iteration_SIMPLE_algorithm =0;"<<std::endl;
					}
					//printf(" number_iterations_SIMPLE_algorithm=%lld\n",number_iteration_SIMPLE_algorithm );
				}
				else {
					printf("WARNING!!! number_iterations_SIMPLE_algorithm not found in file premeshin.txt\n");
					number_iteration_SIMPLE_algorithm = 0; // default 1000 iterations.
					printf(" number_iterations_SIMPLE_algorithm=%u\n", number_iteration_SIMPLE_algorithm);
					if (bSTOP_Reading) system("pause");
				}
				
				if (imakesource("stabilization_amg1r5_algorithm", idin)) {
					// Найдено успешно.
					switch (idin) {
					case 0: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::NONE_only_amg1r5;
						break;
					case 1: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
						break;
					case 2: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::FGMRes_plus_amg1r5;
						break;
					case 3: stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5;
						break;
					default:
						stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5;
						break;
					}
					//printf(" stabilization_amg1r5_algorithm=%lld\n",din );
				}
				else {
					printf("WARNING!!! stabilization_amg1r5_algorithm not found in file premeshin.txt\n");
					stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5; // default BiCGStab+amg1r5.
					printf(" stabilization_amg1r5_algorithm= BiCGStab_plus_amg1r5\n");
					if (bSTOP_Reading) system("pause");
				}


				if (imakesource("ilevelMeshCabinet", idin)) {
					// Найдено успешно.

					BonLevelDrobim = idin;

					//printf(" ilevelMeshCabinet=%lld\n",idin );
					//system("pause");
				}
				else {
					//printf("WARNING!!! ilevelMeshCabinet not found in file premeshin.txt\n");
					BonLevelDrobim = -1;
					//printf(" ilevelMeshCabinet= BonLevelDrobim\n");
					//if (bSTOP_Reading) system("pause");
				}

				if (imakesource("isizearrAdaptRegion", idin)) {
					// Найдено успешно.

					isizearrAdaptRegion = idin;

					//printf(" ilevelMeshCabinet=%lld\n",idin );
				}
				else {
					
					isizearrAdaptRegion = 0;
					
				}

				if (isizearrAdaptRegion > 0) {

					AdaptRegion = new ADAPTREGION[isizearrAdaptRegion];

					for (int i_63 = 0; i_63 < isizearrAdaptRegion; ++i_63) {

						// Адаптация сетки региона номер 1.
						if (i_63 == 0) {

							if (fmakesource("AdaptRegionxS1", fin)) {
								// Найдено успешно.
								AdaptRegion[0].xS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionxS1 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionxE1", fin)) {
								// Найдено успешно.
								AdaptRegion[0].xE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionxE1 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionyS1", fin)) {
								// Найдено успешно.
								AdaptRegion[0].yS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionyS1 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionyE1", fin)) {
								// Найдено успешно.
								AdaptRegion[0].yE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionyE1 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionzS1", fin)) {
								// Найдено успешно.
								AdaptRegion[0].zS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionzS1 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionzE1", fin)) {
								// Найдено успешно.
								AdaptRegion[0].zE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionzE1 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (imakesource("ilevelMeshAdaptRegion1", idin)) {
								// Найдено успешно.

								AdaptRegion[0].ilevelMeshAdaptRegion = idin;

								//printf(" ilevelMeshAdaptRegion1=%lld\n",idin );
							}
							else {
								//printf("WARNING!!! ilevelMeshAdaptRegion1 not found in file premeshin.txt\n");
								AdaptRegion[0].ilevelMeshAdaptRegion = -1;
								//printf(" ilevelMeshAdaptRegion1 = AdaptRegion[0].ilevelMeshAdaptRegion\n");
								//if (bSTOP_Reading) system("pause");
							}

						}
					
					
						if (i_63 == 1) {

							if (fmakesource("AdaptRegionxS2", fin)) {
								// Найдено успешно.
								AdaptRegion[1].xS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionxS2 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionxE2", fin)) {
								// Найдено успешно.
								AdaptRegion[1].xE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionxE2 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionyS2", fin)) {
								// Найдено успешно.
								AdaptRegion[1].yS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionyS2 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionyE2", fin)) {
								// Найдено успешно.
								AdaptRegion[1].yE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionyE2 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionzS2", fin)) {
								// Найдено успешно.
								AdaptRegion[1].zS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionzS2 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionzE2", fin)) {
								// Найдено успешно.
								AdaptRegion[1].zE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionzE2 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (imakesource("ilevelMeshAdaptRegion2", idin)) {
								// Найдено успешно.

								AdaptRegion[1].ilevelMeshAdaptRegion = idin;

								//printf(" ilevelMeshAdaptRegion2=%lld\n",idin );
							}
							else {
								//printf("WARNING!!! ilevelMeshAdaptRegion2 not found in file premeshin.txt\n");
								AdaptRegion[1].ilevelMeshAdaptRegion = -1;
								//printf(" ilevelMeshAdaptRegion2 = AdaptRegion[1].ilevelMeshAdaptRegion\n");
								//if (bSTOP_Reading) system("pause");
							}

						}
					
					
						if (i_63 == 2) {

							if (fmakesource("AdaptRegionxS3", fin)) {
								// Найдено успешно.
								AdaptRegion[2].xS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionxS3 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionxE3", fin)) {
								// Найдено успешно.
								AdaptRegion[2].xE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionxE3 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionyS3", fin)) {
								// Найдено успешно.
								AdaptRegion[2].yS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionyS3 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionyE3", fin)) {
								// Найдено успешно.
								AdaptRegion[2].yE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionyE3 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionzS3", fin)) {
								// Найдено успешно.
								AdaptRegion[2].zS = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionzS3 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (fmakesource("AdaptRegionzE3", fin)) {
								// Найдено успешно.
								AdaptRegion[2].zE = scale * fin;
							}
							else {
								std::cout << "Fatal error not found AdaptRegionzE3 in premeshin.txt\n";
								system("pause");
								exit(1);
							}

							if (imakesource("ilevelMeshAdaptRegion3", idin)) {
								// Найдено успешно.

								AdaptRegion[2].ilevelMeshAdaptRegion = idin;

								//printf(" ilevelMeshAdaptRegion3=%lld\n",idin );
							}
							else {
								//printf("WARNING!!! ilevelMeshAdaptRegion3 not found in file premeshin.txt\n");
								AdaptRegion[2].ilevelMeshAdaptRegion = -1;
								//printf(" ilevelMeshAdaptRegion3 = AdaptRegion[2].ilevelMeshAdaptRegion\n");
								//if (bSTOP_Reading) system("pause");
							}

						}

                    }

				}

				{
					char name0[1000] = "egddata";
					name0[0] = '\0'; strcat_s(name0, "free_debug_parametr1");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						// Свободный параметр для выполнения отладки передаваемый из интерфейса внутрь солвера.
						free_debug_parametr1 = static_cast<doublereal>(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						free_debug_parametr1 = 0.0; // Свободный параметр для выполнения отладки передаваемый из интерфейса внутрь солвера.
						printf(" free_debug_parametr1=%e\n", free_debug_parametr1);
						if (bSTOP_Reading) system("pause");
					}
				}

			}


			integer ilb_p = 0;// Количество блоков внутри которых задана тепловая мощность.
			// Для правильного учёта нужно использовать функцию whot_is_block, а для этого нужна сетка.
			// Потому что блоки перекрываются и мощность надо учитывать с учётом экранирования.
			// Здесь это никоим образом не учитывается.


			int lb_ic = 0; // Количество пропущенных объектов которые полностью перекрыты
			// другими телами.
			// Список объектов перекрытых другими телами.
			int* b_list_delete_in_scan = new int[lb];
			for (int i_1 = 1; i_1 < lb; ++i_1) {
				for (int i_2 = 0; i_2 < i_1; ++i_2) {
					if (b[i_1].g.itypegeom==PRISM && b[i_2].g.itypegeom == PRISM &&
						(b[i_1].g.xS <= b[i_2].g.xS) && (b[i_1].g.xE >= b[i_2].g.xE) &&
						(b[i_1].g.yS <= b[i_2].g.yS) && (b[i_1].g.yE >= b[i_2].g.yE) &&
						(b[i_1].g.zS <= b[i_2].g.zS) && (b[i_1].g.zE >= b[i_2].g.zE))
					{
						if (lb_ic < lb) {
							b_list_delete_in_scan[lb_ic] = i_2;
							lb_ic++;
						}
					}
					// здесь еще надо предусмотреть различные сочетания призм и цилиндров.
					// Но все равно 100% корректность даст только сетка и функция whot_is_block.
					// Два цилиндра.
					
					if ((b[i_1].g.itypegeom == CYLINDER) && (b[i_2].g.itypegeom == CYLINDER) &&
						(b[i_1].g.iPlane == b[i_2].g.iPlane)&&(fabs(b[i_1].g.R_in_cyl)<1.0e-30))
					{
						doublereal x1, y1, z1;
						doublereal x2, y2, z2;
						doublereal x3, y3, z3;
						doublereal x4, y4, z4;

						switch (b[i_1].g.iPlane)
						{
						case XY_PLANE:
							// находим окаймляющий четырёхугольник.
							x1 = b[i_2].g.xC + b[i_2].g.R_out_cyl;
							y1 = b[i_2].g.yC + b[i_2].g.R_out_cyl;
							x2 = b[i_2].g.xC - b[i_2].g.R_out_cyl;
							y2 = b[i_2].g.yC + b[i_2].g.R_out_cyl;
							x3 = b[i_2].g.xC + b[i_2].g.R_out_cyl;
							y3 = b[i_2].g.yC - b[i_2].g.R_out_cyl;
							x4 = b[i_2].g.xC - b[i_2].g.R_out_cyl;
							y4 = b[i_2].g.yC - b[i_2].g.R_out_cyl;
							if ((b[i_1].g.zC <= b[i_2].g.zC) && (b[i_1].g.zC + b[i_1].g.Hcyl >= b[i_2].g.zC + b[i_2].g.Hcyl) &&
								(sqrt((x1 - b[i_1].g.xC)*(x1 - b[i_1].g.xC) + (y1 - b[i_1].g.yC)*(y1 - b[i_1].g.yC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((x2 - b[i_1].g.xC)*(x2 - b[i_1].g.xC) + (y2 - b[i_1].g.yC)*(y2 - b[i_1].g.yC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((x3 - b[i_1].g.xC)*(x3 - b[i_1].g.xC) + (y3 - b[i_1].g.yC)*(y3 - b[i_1].g.yC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((x4 - b[i_1].g.xC)*(x4 - b[i_1].g.xC) + (y4 - b[i_1].g.yC)*(y4 - b[i_1].g.yC)) <= b[i_1].g.R_out_cyl))
							{
								if (lb_ic < lb) {
									b_list_delete_in_scan[lb_ic] = i_2;
									lb_ic++;
								}
							}
							break;
						case XZ_PLANE:
							// находим окаймляющий четырёхугольник.
							x1 = b[i_2].g.xC + b[i_2].g.R_out_cyl;
							z1 = b[i_2].g.zC + b[i_2].g.R_out_cyl;
							x2 = b[i_2].g.xC - b[i_2].g.R_out_cyl;
							z2 = b[i_2].g.zC + b[i_2].g.R_out_cyl;
							x3 = b[i_2].g.xC + b[i_2].g.R_out_cyl;
							z3 = b[i_2].g.zC - b[i_2].g.R_out_cyl;
							x4 = b[i_2].g.xC - b[i_2].g.R_out_cyl;
							z4 = b[i_2].g.zC - b[i_2].g.R_out_cyl;
							if ((b[i_1].g.yC <= b[i_2].g.yC) && (b[i_1].g.yC + b[i_1].g.Hcyl >= b[i_2].g.yC + b[i_2].g.Hcyl) &&
								(sqrt((x1 - b[i_1].g.xC)*(x1 - b[i_1].g.xC) + (z1 - b[i_1].g.zC)*(z1 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((x2 - b[i_1].g.xC)*(x2 - b[i_1].g.xC) + (z2 - b[i_1].g.zC)*(z2 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((x3 - b[i_1].g.xC)*(x3 - b[i_1].g.xC) + (z3 - b[i_1].g.zC)*(z3 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((x4 - b[i_1].g.xC)*(x4 - b[i_1].g.xC) + (z4 - b[i_1].g.zC)*(z4 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl))
							{
								if (lb_ic < lb) {
									b_list_delete_in_scan[lb_ic] = i_2;
									lb_ic++;
								}
							}
							break;
						case YZ_PLANE:
							// находим окаймляющий четырёхугольник.
							y1 = b[i_2].g.yC + b[i_2].g.R_out_cyl;
							z1 = b[i_2].g.zC + b[i_2].g.R_out_cyl;
							y2 = b[i_2].g.yC - b[i_2].g.R_out_cyl;
							z2 = b[i_2].g.zC + b[i_2].g.R_out_cyl;
							y3 = b[i_2].g.yC + b[i_2].g.R_out_cyl;
							z3 = b[i_2].g.zC - b[i_2].g.R_out_cyl;
							y4 = b[i_2].g.yC - b[i_2].g.R_out_cyl;
							z4 = b[i_2].g.zC - b[i_2].g.R_out_cyl;
							if ((b[i_1].g.xC <= b[i_2].g.xC) && (b[i_1].g.xC + b[i_1].g.Hcyl >= b[i_2].g.xC + b[i_2].g.Hcyl) &&
								(sqrt((y1 - b[i_1].g.yC)*(y1 - b[i_1].g.yC) + (z1 - b[i_1].g.zC)*(z1 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((y2 - b[i_1].g.yC)*(y2 - b[i_1].g.yC) + (z2 - b[i_1].g.zC)*(z2 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((y3 - b[i_1].g.yC)*(y3 - b[i_1].g.yC) + (z3 - b[i_1].g.zC)*(z3 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl) &&
								(sqrt((y4 - b[i_1].g.yC)*(y4 - b[i_1].g.yC) + (z4 - b[i_1].g.zC)*(z4 - b[i_1].g.zC)) <= b[i_1].g.R_out_cyl))
							{
								if (lb_ic < lb) {
									b_list_delete_in_scan[lb_ic] = i_2;
									lb_ic++;
								}
							}
							break;
						}
					}
					
				}
			}


			doublereal dpower = 0.0; // Суммарная тепловая мощность в блоках.
			integer ipoly = 0, icyl = 0, iprism = 0, icad_stl = 0;
			integer ihol = 0, isol = 0, iflui = 0;
			for (integer i_1 = 0; i_1 < lb; ++i_1) {
				bool bscalok = true;
				for (integer i_4 = 0; i_4 < lb_ic; i_4++) {
					if (i_1 == b_list_delete_in_scan[i_4]) {
						bscalok = false; // Тело перекрыто другим телом.
					}
				}
				if (bscalok) {
					if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
						ihol++;
					}
					if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						iflui++;
					}
					if (b[i_1].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
						isol++;
					}
					if (b[i_1].g.itypegeom == CAD_STL)
					{
						if (b[i_1].n_Sc > 0) {
							doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
							doublereal vol = b[i_1].g.volume_CAD_STL();
							if (vol < 1.0e-36) {
								printf("ERROR: zero volume in CAD_STL block number %lld\n", i_1);
								system("PAUSE");
								exit(1);
							}
							//if (fabs(b[i_1].arr_Sc[0]) > 0.0) {
							if (pdiss > 0.0) {
								ilb_p++;
								//dpower += b[i_1].arr_Sc[0];
								dpower += pdiss * vol;
								//printf("PRISM id=%lld pdiss=%e vol=%e dpower=%e\n", i_1, pdiss, vol, dpower);
								//system("pause");
							}
						}
						icad_stl++;
					}
					if (b[i_1].g.itypegeom == PRISM) {
						// 0 - PRISM object

						// 13.08.2019 Автомат настройки допусков сетки. 
						// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
						// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
						// Она не может быть больше стороны кабинета деленной на 15.
						if (0.1*fabs(b[i_1].g.xE - b[i_1].g.xS) < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * fabs(b[i_1].g.xE - b[i_1].g.xS);
						if (0.1*fabs(b[i_1].g.yE - b[i_1].g.yS) < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * fabs(b[i_1].g.yE - b[i_1].g.yS);
						if (0.1*fabs(b[i_1].g.zE - b[i_1].g.zS) < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * fabs(b[i_1].g.zE - b[i_1].g.zS);

						if (lb == 1) {// Течение в каверне или тест Дэвиса.
							shorter_length_for_simplificationX_BASIC = 1.0e-10;
							shorter_length_for_simplificationY_BASIC = 1.0e-10;
							shorter_length_for_simplificationZ_BASIC = 1.0e-10;
						}
						else {
							if (shorter_length_for_simplificationX_BASIC > 0.067*(fabs(b[0].g.xE - b[0].g.xS))) shorter_length_for_simplificationX_BASIC = 0.067*(fabs(b[0].g.xE - b[0].g.xS));
							if (shorter_length_for_simplificationY_BASIC > 0.067*(fabs(b[0].g.yE - b[0].g.yS))) shorter_length_for_simplificationY_BASIC = 0.067*(fabs(b[0].g.yE - b[0].g.yS));
							if (shorter_length_for_simplificationZ_BASIC > 0.067*(fabs(b[0].g.zE - b[0].g.zS))) shorter_length_for_simplificationZ_BASIC = 0.067*(fabs(b[0].g.zE - b[0].g.zS));
						}

						iprism++;
						if (b[i_1].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) {
							if (b[i_1].n_Sc > 0) {
								doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
								doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS)*fabs(b[i_1].g.yE - b[i_1].g.yS)*fabs(b[i_1].g.zE - b[i_1].g.zS);
								if (vol < 1.0e-36) {
									printf("xE=%e xS=%e yE=%e yS=%e zE=%e zS=%e\n", b[i_1].g.xE, b[i_1].g.xS, b[i_1].g.yE, b[i_1].g.yS, b[i_1].g.zE, b[i_1].g.zS );
									std::cout << b[i_1].name << std::endl;
									printf("ERROR: zero volume in PRISM block number %lld\n", i_1);
									system("PAUSE");
									exit(1);
								}
								//if (fabs(b[i_1].arr_Sc[0]) > 0.0) {
								if (pdiss > 0.0) {
									ilb_p++;
									//dpower += b[i_1].arr_Sc[0];
									dpower += pdiss * vol;
									//printf("PRISM id=%lld pdiss=%e vol=%e dpower=%e\n", i_1, pdiss, vol, dpower);
									//system("pause");
								}
							}
						}
					}
					if (b[i_1].g.itypegeom == CYLINDER) {
						// Cylinder

						// На тот случай если геометрия состоит только из цилиндров.
						// 13.08.2019 Автомат настройки допусков сетки. 
						// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
						// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
						switch (b[i_1].g.iPlane) {
						case XY_PLANE:
							if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.Hcyl;
							if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
							if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
							if (b[i_1].g.R_in_cyl > 1.0e-36) {
								// Если внутренний радиус существует (задавался пользователем).
								if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
								if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
							}
							break;
						case XZ_PLANE:
							if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.Hcyl;
							if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
							if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
							if (b[i_1].g.R_in_cyl > 1.0e-36) {
								// Если внутренний радиус существует (задавался пользователем).
								if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
								if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
							}
							break;
						case YZ_PLANE:
							if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.Hcyl;
							if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
							if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
							if (b[i_1].g.R_in_cyl > 1.0e-36) {
								// Если внутренний радиус существует (задавался пользователем).
								if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
								if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
							}
							break;
						}


						icyl++;
						if (b[i_1].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) {
							if (b[i_1].n_Sc > 0) {
								doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
								doublereal vol = 0.0;
								vol = b[i_1].g.Hcyl*M_PI*(b[i_1].g.R_out_cyl*b[i_1].g.R_out_cyl - b[i_1].g.R_in_cyl*b[i_1].g.R_in_cyl);
								if (vol < 1.0e-36) {
									printf("ERROR: zero volume in CYLINDER block number %lld\n", i_1);
									system("PAUSE");
									exit(1);
								}
								if (pdiss > 0.0) {
									ilb_p++;

									dpower += pdiss * vol;
									//printf("Cylinder id==%lld pdiss=%e vol=%e dpower=%e\n", i_1, pdiss, vol, dpower);
									//system("pause");
									//printf("ERROR: non zero power in cylinder object.\n");
									//system("PAUSE");
									//exit(1);
								}
							}
						}
					}
					if (b[i_1].g.itypegeom == POLYGON) {
						// Polygon
						ipoly++;
						if (b[i_1].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) {
							if (b[i_1].n_Sc > 0) {
								doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);

								// Объём полигона.
								doublereal vol_poly = Volume_polygon(b[i_1].g.nsizei, b[i_1].g.xi, b[i_1].g.yi, b[i_1].g.zi, b[i_1].g.hi, b[i_1].g.iPlane_obj2);
								if (vol_poly < 1.0e-36) {
									printf("ERROR: zero volume in POLYGON block number %lld\n", i_1);
									system("PAUSE");
									exit(1);
								}
								
								if (pdiss > 0.0) {
									ilb_p++;

									dpower += pdiss*vol_poly;

									//printf("ERROR: non zero power in polygon object.\n");
									//system("PAUSE");
									//exit(1);
								}
							}
						}
					}
				}// bscanlok
            }

			delete[] b_list_delete_in_scan;
			b_list_delete_in_scan = nullptr;

			doublereal dsoupow = 0.0; // интегральная тепловая мощность плоских источников тепла.
			for (integer i_1 = 0; i_1 < ls; ++i_1) {
				dsoupow += s[i_1].power;
			}


			printf("Apriory quick model statistics:\n");
			std::cout << "number of thermal power blocks lb_p=" << ilb_p << std::endl;
			printf("Blocks integral power =%e W\n", dpower);
			std::cout << "number of sources ls=" << ls << std::endl;
			printf("Sources integral power = %e W\n", dsoupow);
			printf("Full total power = %e W\n", dpower + dsoupow);
			// Запоминаем полное тепловыделение в модели.
			d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = dpower + dsoupow;
			std::cout << "number of blocks lb=" << lb << std::endl;
			std::cout << "PRISMS = " << iprism << ", CYLINDERS = " << icyl << ", POLYGONS = " << ipoly << ", CAD_STL = " << icad_stl << "\n";
			std::cout <<  "SOLID: " <<  isol << std::endl;
			std::cout <<  "HOLLOW: " <<  ihol << std::endl;
			std::cout <<  "FLUID: " <<  iflui << std::endl;
			std::cout << "number of walls lw=" <<  lw << std::endl;
			std::cout << "number of units lu=" <<  lu << std::endl;

			if (fp != NULL) {
				fclose(fp); // закрытие файла
				fp = NULL;
			}
			
		}
	}

#endif
	printf("OK. \n");
} // premeshin_new


// считывание параметров из 
// входного файла premeshin.txt
void premeshin(const char *fname, integer &lmatmax, int &lb, int &ls, int &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w,
	doublereal &dgx, doublereal &dgy, doublereal &dgz, int &inx, int &iny, int &inz, doublereal &operatingtemperature,
	integer &ltdp, TEMP_DEP_POWER* &gtdps, int &lu, UNION* &my_union)
{

	FILE* fp = nullptr;
	

#ifdef MINGW_COMPILLER
	int err1 = 0;
	fp = fopen64(fname, "r");
	if (fp != nullptr) {
		err1 = 0;
	}
	else {
		err1 = 1;
	}
#else
	errno_t err1;
	err1 = fopen_s(&fp, fname, "r");
#endif

	if (err1 != 0) {
		printf("Error!!! No input File premeshin.txt.\n");
		printf("You must use the graphical user interface\n");
		printf("AliceMesh_v0_45.exe in Delphi, which will \n");
		printf("prepare the model and write the premeshin.txt\n");
		printf("file in the required format.\n");
		system("pause");
		// Если файла premeshin.txt нет то мы не можем ничего обрабатывать
		// т.к. элементарно отсутствуют входные данные. Поэтому мы выходим из 
		// приложения.
		exit(1);
	}
	else
	{
		if (fp != NULL) {

			
			fclose(fp);
			fp = NULL;
			

			int idin = 0;

			if (imakesource("lb", idin)) {
				// Найдено успешно.
				lb = static_cast<integer>(idin);
				//printf("lb =%lld\n", lb);

				if ((lb > 0)&&(lb < 100000)) {

					


					// начало 21.08.2019 - окончание 23.08.2019
			        // Сдержит простейший парсер.
					bool bSTOP_Reading = true;
					premeshin_new(fname, lmatmax, lb, ls, lw, matlist, b, s, w, dgx, dgy, dgz, inx, iny, inz,
						operatingtemperature, ltdp, gtdps, lu, my_union, bSTOP_Reading);
				}
				else {
					printf("ERROR!!! your model incorrect...\n");
					printf("Number of blocks is equal = %d \n",lb);
					system("pause");
					exit(1);
				}
			}
			else {
				// Если файл записан в старом формате то поддерживается возможность его быстрого считывания.
				// Начиная с 23.08.2019 старый формат начинает считаться устаревшим.				


				// Стабильная и быстрая версия.
			    premeshin_old(fname, lmatmax, lb, ls, lw, matlist, b, s, w, dgx, dgy, dgz, inx, iny, inz,
				              operatingtemperature, ltdp, gtdps, lu, my_union);
			}			
		}
		else {
			printf("Problem open File premeshin.txt \n");
			system("pause");
		}
	}
}

#endif
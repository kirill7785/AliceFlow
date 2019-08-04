// файл inputlaplas.cpp
// содержит ввод структур данных  
// построенных графической чертилкой.

#pragma once
#ifndef _INPUTLAPLAS_CPP_
#define _INPUTLAPLAS_CPP_ 1


#include <math.h> // математические функции
//#include <string.h> // функции обработки строк
#include "power_temperature_depend.cpp" // сплайновая интерполляция таблично заданной функции (зависимость мощности от температуры).

// передаваемый из файла параметр отвечающий за подробность расчётной сетки.
doublereal etalon_max_size_ratio = 2.0;
// параметр отвечающий за качество расчётной сетки.
doublereal etalon_max_size_ratio2 = 30.0; // 30 оптимум по документации FlowVision.
bool bsnap_TO_global = true; // snap to grid

// плоскости
const unsigned char XY = 1;
const unsigned char XZ = 2;
const unsigned char YZ = 3;

// тип блока
const unsigned char SOLID = 1;
const unsigned char HOLLOW = 2;
const unsigned char FLUID = 3;

// Геометрическая форма блока
const unsigned char PRISM = 0;
const unsigned char CYLINDER = 1;
const unsigned char POLYGON = 2;

// точка в трёхмерном пространстве
typedef struct TPOinteger {
	doublereal x=0.0, y=0.0, z=0.0;
} TOCHKA;

// только для enumerate_volume_improved
// и для uniformsimplemeshgen.cpp
// См. также заголовочный файл constr_struct.cpp
typedef struct TBlock_indexes {
	integer iL=-1, iR=-2, jL=-1, jR=-2, kL=-1, kR=-2;
} Block_indexes;

// геометрическое описание
typedef struct TGEOM {
	integer itypegeom= PRISM; // 0 - Prism, 1 - Cylinder, 2 - Polygon
	// Prism
	doublereal xS=0.0, yS=0.0, zS=0.0; // координаты начала объекта
	doublereal xE=0.0, yE=0.0, zE=0.0; // координаты конца объекта
	// Cylinder
	integer iPlane=-2; // плоскость в которой лежит нижнее основание цилиндра.
	doublereal xC=0.0, yC=0.0, zC=0.0, Hcyl=0.0, R_out_cyl=0.0, R_in_cyl=0.0;
	// Polygon
	integer iPlane_obj2=-2; // плоскость в которой лежит нижнее основание полигона.
	integer nsizei=-2; // Количество опорных точек образующих выпуклый полигон.
	doublereal *hi=NULL, *xi=NULL, *yi=NULL, *zi=NULL;
} GEOM;

// Проверка типа геометрии при вводе.
// CHECK_TYPE_GEOM(din);
void CHECK_TYPE_GEOM(integer itypegeom) {
	switch (itypegeom) {
	case PRISM: break;
	case CYLINDER: break;
	case POLYGON: break;
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
	// cp - теплоёмкость,
	// lam - теплопроводность,
	// mu - динамическая вязкость,
	// beta_t - коэффициент линейного температурного расширения.
	//doublereal rho, cp, lam;
	doublereal rho=-1.0e30;
	// 16_11_2016
	// температура передаётся в градусах Цельсия.
	// температурно - зависимые теплопроводность и теплоёмкость.
	integer n_lam = -1, n_cp = -1; // число точек для линейной интерполляции.
	// arr_lam=f(temp_lam);
	// arr_cp=f(temp_cp);
	doublereal *temp_lam=NULL, *temp_cp=NULL, *arr_lam=NULL, *arr_cp=NULL;
	// Ортотропность теплопроводности :
	// позволяет моделировать тепловые трубы и материалы с ортотропной теплопроводностью :
	// CuMoCu, SiC и текстолитовые платы.
	// теплопроводность в заданном координатном направлении домножается на этот множитель.
	doublereal orthotropy_multiplyer_x=1.0, orthotropy_multiplyer_y = 1.0, orthotropy_multiplyer_z = 1.0;
	doublereal mu = -1.0e30, beta_t = -1.0e30, beta_t_solid = -1.0e30;// beta_t_solid для Механики.
	// Следующие два параметра относятся к внутренней библиотеке 
	// материалов:
	integer blibmat=-3; // 0 - материал пользователя и используются 
	// постоянные параметры предложенные выше. 1 - библиотечный материал.
	integer ilibident=-3; // идентификатор библиотечного материала.
	// идентификатор: 
	// 0 - текучая среда которой нет в библиотеке, 
	// 100 - твёрдое тело которого нет в библиотеке.
	// В случае если идентификатор принимает значения 0 и 100 то
	// используются постоянные параметры заданные пользователем, т.к.
	// в этом случае blibmat - обязательно равен 0.
	// 1-99 - библиотечная текучая среда: жидкости, газы и жидкие металлы.
	// 101-<infinity библиотечное твёрдое тело.
	bool bBussineskApproach=true; // нужно ли использовать приближение Обербека-Буссинеска

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
	integer ilawmu=-3; // номер закона

	doublereal mumin = -1.0e30, mumax = -1.0e30; // ограничители динамической вязкости
	doublereal Amu = -1.0e30, Bmu = -1.0e30, Cmu = -1.0e30, degreennmu = -1.0e30; // константы моделей для зависимости вязкости от напряжения сдвига.

	// 8.4.2017 
	// Параметры прочности материала.
	doublereal mu_Lame = -1.0e30, lambda_Lame = -1.0e30; // Коэффициенты Лямэ.

} PROP;

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
	doublereal MPI = 3.1415926;
	doublereal a = ((1.0 + h*h)*(1.0 + w*w)) / (1.0 + h*h + w*w);
	doublereal b = (w*w*(1.0 + h*h + w*w)) / ((1.0 + w*w)*(h*h + w*w));
	doublereal c = (h*h*(1.0 + h*h + w*w)) / ((1.0 + h*h)*(h*h + w*w));
	doublereal Fi_j = (1.0 / (MPI*w))*(h*atan(1.0 / h) + w*atan(1.0 / w) -
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
	doublereal MPI = 3.1415926;
	doublereal x_tilda1 = sqrt(1.0 + x_tilda*x_tilda);
	doublereal y_tilda1 = sqrt(1.0 + y_tilda*y_tilda);
	doublereal Fi_j = (1.0 / (MPI*x_tilda*y_tilda))*(log((x_tilda1*x_tilda1*y_tilda1*y_tilda1) / (x_tilda1*x_tilda1 + y_tilda1*y_tilda1 - 1.0)) +
		2.0*x_tilda*(y_tilda1*atan(x_tilda / y_tilda1) - atan(x_tilda)) + 2.0*y_tilda*(x_tilda1*atan(y_tilda / x_tilda1) - atan(y_tilda)));
	return Fi_j;
} //view_factor_aligned_parallel_rectangles

// Для вычисления усреднённой температуры на грани нужен
// список узлов грани (два внутрених соседа для внутренней грани) и 
// (граничный узел и ближайший к нему соседний внутренний для узлов на границе расчётной области).
typedef struct TMY_PAIR {
	integer node1, node21, node22, node23, node24; // если отсутствует то стоит -1.
	doublereal  dS1=0.0, dS2=0.0, dS3=0.0, dS4=0.0; // площадь ячейки.
	// 20 сентября 2016.
	// dS должно быть позже полностью удалено.
	//doublereal dS;
} MY_PAIR;

// Параметр нижней релаксации для вычисления плотностей радиационных потоков.
const doublereal alpha_Radiation_block = 0.1;
// Постоянная Стефана-Больцмана для вычисления плотностей радиационных потоков.
const doublereal sigma_Radiation_block = 5.670367e-8; // Wxm!-2xK!-4.
// фиксированное число итераций для нахождения плотностей радиационых потоков.
const integer maxiter_Radiation_block = 1000; // вообще хватает и 200, здесь с запасом.


typedef struct TBLOCKRADIATION {
	// emissivity:
	doublereal emissW = -1.0e30, emissE = -1.0e30, emissS = -1.0e30, emissN = -1.0e30, emissB = -1.0e30, emissT = -1.0e30;
	bool binternalRadiation=false;
	// View Factors for Prism Object.
	doublereal FWE = -1.0e30, FWS = -1.0e30, FWN = -1.0e30, FWB = -1.0e30, FWT = -1.0e30;
	doublereal FEW = -1.0e30, FES = -1.0e30, FEN = -1.0e30, FEB = -1.0e30, FET = -1.0e30;
	doublereal FSW = -1.0e30, FSE = -1.0e30, FSN = -1.0e30, FSB = -1.0e30, FST = -1.0e30;
	doublereal FNW = -1.0e30, FNE = -1.0e30, FNS = -1.0e30, FNB = -1.0e30, FNT = -1.0e30;
	doublereal FBW = -1.0e30, FBE = -1.0e30, FBS = -1.0e30, FBN = -1.0e30, FBT = -1.0e30;
	doublereal FTW = -1.0e30, FTE = -1.0e30, FTS = -1.0e30, FTN = -1.0e30, FTB = -1.0e30;
	// Температуры в Кельвинах на гранях Prism Object:
	// среднее арифметическое температуры на грани Prism Object.
	doublereal TempW = -1.0e30, TempE = -1.0e30, TempS = -1.0e30, TempN = -1.0e30, TempB = -1.0e30, TempT = -1.0e30;
	// Плотности радиационных тепловых потоков на гранях Prism Object.
	doublereal JW=0.0, JE=0.0, JS = 0.0, JN = 0.0, JB = 0.0, JT = 0.0;
	

	// список узлов на каждой из граней и их количество:
	MY_PAIR* nodelistW=NULL;
	integer nodelistWsize=-4;
	MY_PAIR* nodelistE=NULL;
	integer nodelistEsize = -4;
	MY_PAIR* nodelistS=NULL;
	integer nodelistSsize = -4;
	MY_PAIR* nodelistN=NULL;
	integer nodelistNsize = -4;
	MY_PAIR* nodelistB=NULL;
	integer nodelistBsize = -4;
	MY_PAIR* nodelistT=NULL;
	integer nodelistTsize = -4;
} BLOCKRADIATION;

// блок
typedef struct TBLOCK {
	integer itype = -4; // тип SOLID, HOLLOW или FLUID.
	GEOM g;
	integer imatid = -4; // идентификатор материала в библиотеке
	//doublereal Sc; // мощность тепловыделения на единицу объёма
	// 19 11 2016 Температурно зависимая мощность тепловыделения.
	integer n_Sc = -4;
	doublereal *arr_Sc=NULL, *temp_Sc=NULL;
    // стиль зависимости мощности тепловыделения от времени.
	integer ipower_time_depend = -4; // 0 - не зависит, 1 - square wave зависимость.
	// Всё относящееся к теплообмену излучением:
	// излучательные способности поверхностей, модель вакуумного промежутка.
	BLOCKRADIATION radiation;
	bool bvisible=true; // Виден ли блок при экспорте в техплот.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id = -4;
	// Фиксировать боковую стенку цилиндра ?
	bool CylinderFixed=false;
} BLOCK;



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
doublereal get_lam(integer nsize_lam, doublereal* &temp_lam, doublereal* &arr_lam, doublereal t_C) {
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
				return (af*t_C+bf);
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

  // Возвращает значение теплоёмкости при температуре t_C в градусах.
doublereal get_cp(integer nsize_cp, doublereal* &temp_cp, doublereal* &arr_cp, doublereal t_C) {
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
					return (af*t_C + bf);
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
// по заданным значениям :
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
		for (integer k = 0; k <= maxiter_Radiation_block; k++)
		{
			// Возможно недоэтерированность это очень хорошо а то иначе мы получим глобальную расходимость.

			// ВНИМАНИЕ !!! b.radiation.Temp в Кельвинах.

			JW_new = alpha_Radiation_block*((sigma_Radiation_block*b.radiation.TempW*b.radiation.TempW*b.radiation.TempW*b.radiation.TempW + ((1.0 - b.radiation.emissW) / b.radiation.emissW)*
				(b.radiation.FWE*b.radiation.JE + b.radiation.FWS*b.radiation.JS + b.radiation.FWN*b.radiation.JN + b.radiation.FWB*b.radiation.JB +
				b.radiation.FWT*b.radiation.JT)) / apW) + (1.0 - alpha_Radiation_block)*b.radiation.JW;

			JE_new = alpha_Radiation_block*((sigma_Radiation_block*b.radiation.TempE*b.radiation.TempE*b.radiation.TempE*b.radiation.TempE + ((1.0 - b.radiation.emissE) / b.radiation.emissE)*
				(b.radiation.FEW*b.radiation.JW + b.radiation.FES*b.radiation.JS + b.radiation.FEN*b.radiation.JN + b.radiation.FEB*b.radiation.JB +
				b.radiation.FET*b.radiation.JT)) / apE) + (1.0 - alpha_Radiation_block)*b.radiation.JE;

			JS_new = alpha_Radiation_block*((sigma_Radiation_block*b.radiation.TempS*b.radiation.TempS*b.radiation.TempS*b.radiation.TempS + ((1.0 - b.radiation.emissS) / b.radiation.emissS)*
				(b.radiation.FSW*b.radiation.JW + b.radiation.FSE*b.radiation.JE + b.radiation.FSN*b.radiation.JN + b.radiation.FSB*b.radiation.JB +
				b.radiation.FST*b.radiation.JT)) / apS) + (1.0 - alpha_Radiation_block)*b.radiation.JS;

			JN_new = alpha_Radiation_block*((sigma_Radiation_block*b.radiation.TempN*b.radiation.TempN*b.radiation.TempN*b.radiation.TempN + ((1.0 - b.radiation.emissN) / b.radiation.emissN)*
				(b.radiation.FNW*b.radiation.JW + b.radiation.FNE*b.radiation.JE + b.radiation.FNS*b.radiation.JS + b.radiation.FNB*b.radiation.JB +
				b.radiation.FNT*b.radiation.JT)) / apN) + (1.0 - alpha_Radiation_block)*b.radiation.JN;

			JB_new = alpha_Radiation_block*((sigma_Radiation_block*b.radiation.TempB*b.radiation.TempB*b.radiation.TempB*b.radiation.TempB + ((1.0 - b.radiation.emissB) / b.radiation.emissB)*
				(b.radiation.FBW*b.radiation.JW + b.radiation.FBE*b.radiation.JE + b.radiation.FBS*b.radiation.JS + b.radiation.FBN*b.radiation.JN +
				b.radiation.FBT*b.radiation.JT)) / apB) + (1.0 - alpha_Radiation_block)*b.radiation.JB;

			JT_new = alpha_Radiation_block*((sigma_Radiation_block*b.radiation.TempT*b.radiation.TempT*b.radiation.TempT*b.radiation.TempT + ((1.0 - b.radiation.emissT) / b.radiation.emissT)*
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

		printf("TW=%e TE=%e TS=%e TN=%e TB=%e TT=%e\n", b.radiation.TempW, b.radiation.TempE, b.radiation.TempS, b.radiation.TempN, b.radiation.TempB, b.radiation.TempT);
		printf("JW=%e JE=%e JS=%e JN=%e JB=%e JT=%e\n", b.radiation.JW, b.radiation.JE, b.radiation.JS, b.radiation.JN, b.radiation.JB, b.radiation.JT);
	//	system("PAUSE");

#ifdef MINGW_COMPILLER
		fprintf(fp_radiation_log, "%e %e ", b.radiation.TempB, b.radiation.TempT);
#else
		fprintf_s(fp_radiation_log, "%e %e ", b.radiation.TempB, b.radiation.TempT);
#endif			
		
	}
} // calculation_density_radiation_heat_flux

// Мощность рассеиваемая в тепло зависит от
// температуры и смещения стока.
// Информация о зависимости хранится в таблично заданной 
// функция. Таблица взята из результатов расчёта в програме
// Г.З. Гарбера.
typedef struct TTEMP_DEP_POWER {
	char* sname=NULL; // имя текстового файла содержащего табличную функцию.
	
	integer intemp=0; // количество различных дискретных значений температуры.
	integer inoffset_drain=0; // количество различных дискретных хначений смещения стока.
	
	// таблица: левый столбец температуры,
	// верхняя строка без первого значения - смещения стока.
	// Остальные значения в таблице соответствуют рассеиваемой мощности.
	doublereal* rtemp=NULL; // значения температуры
	doublereal* roffset_drain=NULL; // значения смещения стока
	doublereal** rpower_table=NULL; // таблица мощностей.

} TEMP_DEP_POWER;

// источник тепла
typedef struct TSOURCE {
	// мощность источника тепла
	// в случае если мощность зависит от максимальной температуры и смещения стока то
	// это значение автоматически пересчитывается перед сборкой матрицы, 
	// так что в коде сборки граничного условия ничего менять не нужно.
	doublereal power=0.0; 
	// Мощность может зависеть от температуры и смещения стока.
	// Зависимость задаётся таблично на основе расчётов в программе Г.З. Гарбера.
	bool bgarber_depend=false; // задана ли зависимость рассеиваемой мощности от температуры и смещения стока.
	// В случае bgarber_depend == true, значение power является коэффициентом на который домножается 
	// истинное табличное значение мощности.
	integer igarber_depend=-1; // уникальный номер табличной зависимости. (нумерация начинается с нуля).
	doublereal roperation_offset_drain=28.0; // рабочее значение смещения стока.
	doublereal power_multiplyer=1.0; // на эту константу домножается реальная мощность, это удобно для корректирования
	// табличной мощности например когда нужно считать четверть транзистора.

	doublereal square=0.0; // площадь источника тепла
	integer iPlane=XZ; // плоскость в которой лежит источник тепла
	// 1 - XY, 2 - XZ, 3 - YZ.
	GEOM g; // границы объекта.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id=0;
} SOURCE;

// стенка (идеальный теплооотвод)
typedef struct TWALL {
	// 1 - Дирихле, 2 - однородное условие Неймана.
	integer ifamily=2; // род краевого условия по температуре
	doublereal Tamb=20.0, hf=0.0; // значение температуры на идеальном теплоотводе и тепловой поток.
	doublereal emissivity=0.8, film_coefficient=3.0;
	bool bsymmetry = false, bpressure = false, bopening=false;
	doublereal Vx=0.0, Vy=0.0, Vz=0.0, P=0.0;
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
	integer ithermal_Stress_boundary_condition=0;
	doublereal xForce=0.0, yForce=0.0, zForce=0.0; // Три компоненты силы, приложенные к узлу в Ньютонах.
    integer iPlane=XZ; // плоскость в которой лежит стенка
	// 1 - XY, 2 - XZ, 3 - YZ.
	GEOM g; // границы объекта.
	// Фиксированная граница с нулевым смещением для
	// расчёта прочности конструкции.
	//bool bfixboundary; // true фиксированная, false свободная.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id=0;
} WALL;


// Пересекаются ли два отрезка 02.08.2019
bool b_is_intersect(doublereal x1, doublereal y1,
	doublereal x2, doublereal y2, doublereal x3, doublereal y3, 
	doublereal x4, doublereal y4)
{
	bool bintersect = true;
	doublereal Ua, Ub, numerator_a, numerator_b, denominator;

	denominator = (y4 - y3)*(x1 - x2) - (x4 - x3)*(y1 - y2);
    
	if (fabs(denominator) < 1.0e-40) {
		if (fabs((x1*y2 - x2 * y1)*(x4 - x3) - (x3*y4 - x4 * y3)*(x2 - x1)) <= 1.0e-40
			&& fabs((x1*y2 - x2 * y1)*(y4 - y3) - (x3*y4 - x4 * y3)*(y2 - y1)) <= 1.0e-40)
		{
			
			if ((fabs(x1 - x2) < 1.0e-40) && (fabs(x1 - x3) < 1.0e-40) && (fabs(x1 - x4) < 1.0e-40)) {
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
					bintersect = true; // пересекаются
				}
			}
			else if ((fabs(y1 - y2) < 1.0e-40) && (fabs(y1 - y3) < 1.0e-40) && (fabs(y1 - y4) < 1.0e-40)) {
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
					bintersect = true; // пересекаются
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
				bintersect = true; // пересекаются
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
		
		if (Ua > 1.0e-40 && Ua < 1.0-1.0e-40 && Ub > 1.0e-40 && Ub < 1.0-1.0e-40) {
			printf("denominator=%e\n", denominator);
			printf("Ua=%e Ub=%e\n", Ua, Ub);
			bintersect = true;
		}
		else {
			bintersect = false;
		}
	}

	return bintersect;
}

// Проверяет если ли выход за пределы кабинета
// среди блоков, стенок и источников тепла. 02.08.2019.
void BODY_CHECK(BLOCK* &b, integer lb, WALL* &w, integer lw, SOURCE* &s, integer ls) {
	bool bOk = true;
	for (integer i = lb - 1; i >= 0; i--) {
		// Проверка типа геометрии.
		CHECK_TYPE_GEOM(b[i].g.itypegeom);
		if (b[i].g.itypegeom == PRISM) {
			// Ошибка не число.
			if (b[i].g.xS != b[i].g.xS) {
				printf("Error body[%lld] xS NOT VALUE=%e.\n", i, b[i].g.xS);
				system("pause");
				bOk = false;
			}
			if (b[i].g.xE != b[i].g.xE) {
				printf("Error body[%lld] xE NOT VALUE=%e.\n", i, b[i].g.xE);
				system("pause");
				bOk = false;
			}
			if (b[i].g.yS != b[i].g.yS) {
				printf("Error body[%lld] yS NOT VALUE=%e.\n", i, b[i].g.yS);
				system("pause");
				bOk = false;
			}
			if (b[i].g.yE != b[i].g.yE) {
				printf("Error body[%lld] yE NOT VALUE=%e.\n", i, b[i].g.yE);
				system("pause");
				bOk = false;
			}
			if (b[i].g.zS != b[i].g.zS) {
				printf("Error body[%lld] zS NOT VALUE=%e.\n", i, b[i].g.zS);
				system("pause");
				bOk = false;
			}
			if (b[i].g.zE != b[i].g.zE) {
				printf("Error body[%lld] zE NOT VALUE=%e.\n", i, b[i].g.zE);
				system("pause");
				bOk = false;
			}
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (b[i].g.xS >= b[i].g.xE) {
				printf("body[%lld].xS=%e >= body[%lld].xE=%e\n", i, b[i].g.xS, i, b[i].g.xE);
				system("pause");
				doublereal temp = b[i].g.xE;
				b[i].g.xE = b[i].g.xS;
				b[i].g.xS = temp;
				bOk = false;
			}
			if (b[i].g.yS >= b[i].g.yE) {
				printf("body[%lld].yS=%e >= body[%lld].yE=%e\n", i, b[i].g.yS, i, b[i].g.yE);
				system("pause");
				doublereal temp = b[i].g.yE;
				b[i].g.yE = b[i].g.yS;
				b[i].g.yS = temp;
				bOk = false;
			}
			if (b[i].g.zS >= b[i].g.zE) {
				printf("body[%lld].zS=%e >= body[%lld].zE=%e\n", i, b[i].g.zS, i, b[i].g.zE);
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
					printf("b[%lld].g.xE=%e > cabinet.xE=%e.\n", i, b[i].g.xE, b[0].g.xE);
					system("pause");
					b[i].g.xE = b[0].g.xE; // исправление.
					bOk = false;
				}
				if (b[i].g.yE > b[0].g.yE) {
					printf("ERROR. Your model is incorrect.\n");
					printf("b[%lld].g.yE=%e > cabinet.yE=%e.\n", i, b[i].g.yE, b[0].g.yE);
					system("pause");
					b[i].g.yE = b[0].g.yE; // исправление.
					bOk = false;
				}
				if (b[i].g.zE > b[0].g.zE) {
					printf("ERROR. Your model is incorrect.\n");
					printf("b[%lld].g.zE=%e > cabinet.zE=%e.\n", i, b[i].g.zE, b[0].g.zE);
					system("pause");
					b[i].g.zE = b[0].g.zE; // исправление.
					bOk = false;
				}
				if (b[i].g.xS < b[0].g.xS) {
					printf("ERROR. Your model is incorrect.\n");
					printf("b[%lld].g.xS=%e < cabinet.xS=%e.\n", i, b[i].g.xS, b[0].g.xS);
					system("pause");
					b[i].g.xS = b[0].g.xS; // исправление.
					bOk = false;
				}
				if (b[i].g.yS < b[0].g.yS) {
					printf("ERROR. Your model is incorrect.\n");
					printf("b[%lld].g.yS=%e < cabinet.yS=%e.\n", i, b[i].g.yS, b[0].g.yS);
					system("pause");
					b[i].g.yS = b[0].g.yS; // исправление.
					bOk = false;
				}
				if (b[i].g.zS < b[0].g.zS) {
					printf("ERROR. Your model is incorrect.\n");
					printf("b[%lld].g.zS=%e < cabinet.zS=%e.\n", i, b[i].g.zS, b[0].g.zS);
					system("pause");
					b[i].g.zS = b[0].g.zS; // исправление.
					bOk = false;
				}
			}
		}
		if (b[i].g.itypegeom == CYLINDER) {
			if (b[i].g.xC != b[i].g.xC) {
				printf("Error body[%lld] xC NOT VALUE=%e.\n",i, b[i].g.xC);
				system("pause");
				bOk = false;
			}
			if (b[i].g.yC != b[i].g.yC) {
				printf("Error body[%lld] yC NOT VALUE=%e.\n", i, b[i].g.yC);
				system("pause");
				bOk = false;
			}
			if (b[i].g.zC != b[i].g.zC) {
				printf("Error body[%lld] zC NOT VALUE=%e.\n", i, b[i].g.zC);
				system("pause");
				bOk = false;
			}
			if (b[i].g.Hcyl != b[i].g.Hcyl) {
				printf("Error body[%lld] Hcyl NOT VALUE=%e.\n", i, b[i].g.Hcyl);
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_out_cyl != b[i].g.R_out_cyl) {
				printf("Error body[%lld] R_out_cyl NOT VALUE=%e.\n", i, b[i].g.R_out_cyl);
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl != b[i].g.R_in_cyl) {
				printf("Error body[%lld] R_in_cyl NOT VALUE=%e.\n", i, b[i].g.R_in_cyl);
				system("pause");
				bOk = false;
			}
			if (b[i].g.Hcyl <= 0.0) {
				printf("Error: negative or zero Hcyl=%e of cylinder body[%lld]\n", b[i].g.Hcyl, i);
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_out_cyl <= 0.0) {
				printf("Error: negative or zero R_out_cyl=%e of cylinder body[%lld]\n", b[i].g.R_out_cyl, i);
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl < 0.0) {
				printf("Error: negative R_in_cyl=%e of cylinder body[%lld]\n", b[i].g.R_in_cyl, i);
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl > b[i].g.R_out_cyl) {
				printf("Error: b[%lld].g.R_in_cyl=%e > b[%lld].g.R_out_cyl=%e\n", i, b[i].g.R_in_cyl, i, b[i].g.R_out_cyl);
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
				case XY:
					xS_box_cyl -= b[i].g.R_out_cyl;
					xE_box_cyl += b[i].g.R_out_cyl;
					yS_box_cyl -= b[i].g.R_out_cyl;
					yE_box_cyl += b[i].g.R_out_cyl;
					zE_box_cyl += b[i].g.Hcyl;
					break;
				case XZ:
					xS_box_cyl -= b[i].g.R_out_cyl;
					xE_box_cyl += b[i].g.R_out_cyl;
					zS_box_cyl -= b[i].g.R_out_cyl;
					zE_box_cyl += b[i].g.R_out_cyl;
					yE_box_cyl += b[i].g.Hcyl;
					break;
				case YZ:
					zS_box_cyl -= b[i].g.R_out_cyl;
					zE_box_cyl += b[i].g.R_out_cyl;
					yS_box_cyl -= b[i].g.R_out_cyl;
					yE_box_cyl += b[i].g.R_out_cyl;
					xE_box_cyl += b[i].g.Hcyl;
					break;
				}

				if (xE_box_cyl > b[0].g.xE) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					printf("body[%lld] xE_box_cyl=%e > cabinet.xE=%e.\n", i, xE_box_cyl, b[0].g.xE);
					system("pause");

					bOk = false;
				}
				if (yE_box_cyl > b[0].g.yE) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					printf("body[%lld] yE_box_cyl=%e > cabinet.yE=%e.\n", i, yE_box_cyl, b[0].g.yE);
					system("pause");

					bOk = false;
				}
				if (zE_box_cyl > b[0].g.zE) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					printf("body[%lld] zE_box_cyl=%e > cabinet.zE=%e.\n", i, zE_box_cyl, b[0].g.zE);
					system("pause");

					bOk = false;
				}
				if (xS_box_cyl < b[0].g.xS) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					printf("body[%lld] xS_box_cyl=%e < cabinet.xS=%e.\n", i, xS_box_cyl, b[0].g.xS);
					system("pause");

					bOk = false;
				}
				if (yS_box_cyl < b[0].g.yS) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					printf("body[%lld] yS_box_cyl=%e < cabinet.yS=%e.\n", i, yS_box_cyl, b[0].g.yS);
					system("pause");

					bOk = false;
				}
				if (zS_box_cyl < b[0].g.zS) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					printf("body[%lld] zS_box_cyl=%e < cabinet.zS=%e.\n", i, zS_box_cyl, b[0].g.zS);
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
					printf("Error number points on polygon very small: b[%lld].g.nsizei=%lld<2\n",i, b[i].g.nsizei);
					system("pause");
					bOk = false;
				}

				if (b[i].g.xi == NULL) {
					printf("Error: no memory allocate for b[%lld].g.xi array.\n",i);
					system("pause");
					bOk = false;
				}

				if (b[i].g.yi == NULL) {
					printf("Error: no memory allocate for b[%lld].g.yi array.\n", i);
					system("pause");
					bOk = false;
				}

				if (b[i].g.zi == NULL) {
					printf("Error: no memory allocate for b[%lld].g.zi array.\n", i);
					system("pause");
					bOk = false;
				}

				if (b[i].g.hi == NULL) {
					printf("Error: no memory allocate for b[%lld].g.hi array.\n", i);
					system("pause");
					bOk = false;
				}

				for (integer ip = 0; ip < b[i].g.nsizei; ip++) {
					if (b[i].g.xi[ip] != b[i].g.xi[ip]) {
						printf("Error body[%lld] xi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.xi[ip]);
						system("pause");
						bOk = false;
					}
					if (b[i].g.yi[ip] != b[i].g.yi[ip]) {
						printf("Error body[%lld] yi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.yi[ip]);
						system("pause");
						bOk = false;
					}
					if (b[i].g.zi[ip] != b[i].g.zi[ip]) {
						printf("Error body[%lld] zi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.zi[ip]);
						system("pause");
						bOk = false;
					}
					if (b[i].g.hi[ip] != b[i].g.hi[ip]) {
						printf("Error body[%lld] hi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.hi[ip]);
						system("pause");
						bOk = false;
					}
				}

				integer i_t1, i_t2, i_t3, it_c=0;
				switch (b[i].g.iPlane_obj2) {
				case XY: 
					for (integer ip = 0; ip < b[i].g.nsizei; ip++) {

						// Детектирование ОШИБКИ самопересечения отрезков в полигоне.
						// Мы проверяем каждый с каждым
						it_c = 0;
						for (integer ip_t = 2; ip_t < b[i].g.nsizei; ip_t++) {
							i_t1 = ip +1;
							if (i_t1 >= b[i].g.nsizei) i_t1 -= b[i].g.nsizei;
							i_t2 = ip + it_c+2;
							if (i_t2 >= b[i].g.nsizei) i_t2 -= b[i].g.nsizei;
							i_t3 = ip + it_c+3;
							if (i_t3 >= b[i].g.nsizei) i_t3 -= b[i].g.nsizei;

							// Мы попарно роверяем есть ли самопересечения каждый с каждым.
							if (b_is_intersect(b[i].g.xi[ip], b[i].g.yi[ip], b[i].g.xi[i_t1], b[i].g.yi[i_t1],
								b[i].g.xi[i_t2], b[i].g.yi[i_t2], b[i].g.xi[i_t3], b[i].g.yi[i_t3])) {
								// ОШИБКА: обнаружено самопересечение отрезков в полигоне.
								printf("Error in POLYGON : Self-intersection of segments in body[%lld]\n",i);
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
							if (fabs(b[i].g.zi[ip]- b[i].g.zi[0]) > 1.0e-40) {
								printf("Error: NON planar polygon. body[%lld] \n",i);
								printf("Diagnostic: fabs(b[%lld].g.zi[%lld]- b[%lld].g.zi[0]) > 1.0e-40\n",i,ip,i);
								system("pause");
								bOk = false;
							}
							if (fabs(b[i].g.hi[ip] - b[i].g.hi[0]) > 1.0e-40) {
								printf("Error: NON planar polygon. body[%lld] \n", i);
								printf("Diagnostic: fabs(b[%lld].g.hi[%lld]- b[%lld].g.hi[0]) > 1.0e-40\n", i, ip, i);
								system("pause");
								bOk = false;
							}
						}
						if (b[i].g.hi[ip] <= 0.0) {
							printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							system("pause");
							bOk = false;
						}
						zS_box_polygon = b[i].g.zi[ip];
						zE_box_polygon = b[i].g.zi[ip]+ b[i].g.hi[ip];
					}
					break;
				case XZ: 
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

							// Мы попарно роверяем есть ли самопересечения каждый с каждым.
							if (b_is_intersect(b[i].g.xi[ip], b[i].g.zi[ip], b[i].g.xi[i_t1], b[i].g.zi[i_t1],
								b[i].g.xi[i_t2], b[i].g.zi[i_t2], b[i].g.xi[i_t3], b[i].g.zi[i_t3])) {
								// ОШИБКА: обнаружено самопересечение отрезков в полигоне.
								printf("Error in POLYGON : Self-intersection of segments in body[%lld]\n", i);
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
							if (fabs(b[i].g.yi[ip] - b[i].g.yi[0]) > 1.0e-40) {
								printf("Error: NON planar polygon. body[%lld] \n", i);
								printf("Diagnostic: fabs(b[%lld].g.yi[%lld]- b[%lld].g.yi[0]) > 1.0e-40\n", i, ip, i);
								system("pause");
								bOk = false;
							}
							if (fabs(b[i].g.hi[ip] - b[i].g.hi[0]) > 1.0e-40) {
								printf("Error: NON planar polygon. body[%lld] \n", i);
								printf("Diagnostic: fabs(b[%lld].g.hi[%lld]- b[%lld].g.hi[0]) > 1.0e-40\n", i, ip, i);
								system("pause");
								bOk = false;
							}
						}
						if (b[i].g.hi[ip] <= 0.0) {
							printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							system("pause");
							bOk = false;
						}
						yS_box_polygon = b[i].g.yi[ip];
						yE_box_polygon = b[i].g.yi[ip] + b[i].g.hi[ip];
					}
					break;
				case YZ:
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

							// Мы попарно роверяем есть ли самопересечения каждый с каждым.
							if (b_is_intersect(b[i].g.yi[ip], b[i].g.zi[ip], b[i].g.yi[i_t1], b[i].g.zi[i_t1],
								b[i].g.yi[i_t2], b[i].g.zi[i_t2], b[i].g.yi[i_t3], b[i].g.zi[i_t3])) {
								// ОШИБКА: обнаружено самопересечение отрезков в полигоне.
								printf("Error in POLYGON : Self-intersection of segments in body[%lld]\n", i);
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
							if (fabs(b[i].g.xi[ip] - b[i].g.xi[0]) > 1.0e-40) {
								printf("Error: NON planar polygon. body[%lld] \n", i);
								printf("Diagnostic: fabs(b[%lld].g.xi[%lld]- b[%lld].g.xi[0]) > 1.0e-40\n", i, ip, i);
								system("pause");
								bOk = false;
							}
							if (fabs(b[i].g.hi[ip] - b[i].g.hi[0]) > 1.0e-40) {
								printf("Error: NON planar polygon. body[%lld] \n", i);
								printf("Diagnostic: fabs(b[%lld].g.hi[%lld]- b[%lld].g.hi[0]) > 1.0e-40\n", i, ip, i);
								system("pause");
								bOk = false;
							}
						}
						if (b[i].g.hi[ip] <= 0.0) {
							printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							system("pause");
							bOk = false;
						}
						xS_box_polygon = b[i].g.xi[ip];
						xE_box_polygon = b[i].g.xi[ip] + b[i].g.hi[ip];
					}
					break;
				}

				if (xE_box_polygon > b[0].g.xE) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					printf("body[%lld] xE_box_cyl=%e > cabinet.xE=%e.\n", i, xE_box_polygon, b[0].g.xE);
					system("pause");

					bOk = false;
				}
				if (yE_box_polygon > b[0].g.yE) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					printf("body[%lld] yE_box_cyl=%e > cabinet.yE=%e.\n", i, yE_box_polygon, b[0].g.yE);
					system("pause");

					bOk = false;
				}
				if (zE_box_polygon > b[0].g.zE) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					printf("body[%lld] zE_box_cyl=%e > cabinet.zE=%e.\n", i, zE_box_polygon, b[0].g.zE);
					system("pause");

					bOk = false;
				}
				if (xS_box_polygon < b[0].g.xS) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					printf("body[%lld] xS_box_cyl=%e < cabinet.xS=%e.\n", i, xS_box_polygon, b[0].g.xS);
					system("pause");

					bOk = false;
				}
				if (yS_box_polygon < b[0].g.yS) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					printf("body[%lld] yS_box_cyl=%e < cabinet.yS=%e.\n", i, yS_box_polygon, b[0].g.yS);
					system("pause");

					bOk = false;
				}
				if (zS_box_polygon < b[0].g.zS) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					printf("body[%lld] zS_box_cyl=%e < cabinet.zS=%e.\n", i, zS_box_polygon, b[0].g.zS);
					system("pause");

					bOk = false;
				}

			}

		}
    }

	for (integer i = 0; i < lw; i++) {

		// Ошибка не число.
		if (w[i].g.xS != w[i].g.xS) {
			printf("Error wall[%lld] xS NOT VALUE=%e.\n", i, w[i].g.xS);
			system("pause");
			bOk = false;
		}
		if (w[i].g.xE != w[i].g.xE) {
			printf("Error wall[%lld] xE NOT VALUE=%e.\n", i, w[i].g.xE);
			system("pause");
			bOk = false;
		}
		if (w[i].g.yS != w[i].g.yS) {
			printf("Error wall[%lld] yS NOT VALUE=%e.\n", i, w[i].g.yS);
			system("pause");
			bOk = false;
		}
		if (w[i].g.yE != w[i].g.yE) {
			printf("Error wall[%lld] yE NOT VALUE=%e.\n", i, w[i].g.yE);
			system("pause");
			bOk = false;
		}
		if (w[i].g.zS != w[i].g.zS) {
			printf("Error wall[%lld] zS NOT VALUE=%e.\n", i, w[i].g.zS);
			system("pause");
			bOk = false;
		}
		if (w[i].g.zE != w[i].g.zE) {
			printf("Error wall[%lld] zE NOT VALUE=%e.\n", i, w[i].g.zE);
			system("pause");
			bOk = false;
		}

		switch (w[i].iPlane) {
		case XY :
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.xS >= w[i].g.xE) {
				printf("wall[%lld].xS=%e >= wall[%lld].xE=%e\n", i, w[i].g.xS, i, w[i].g.xE);
				system("pause");
				doublereal temp = w[i].g.xE;
				w[i].g.xE = w[i].g.xS;
				w[i].g.xS = temp;
				bOk = false;
			}
			if (w[i].g.yS >= w[i].g.yE) {
				printf("wall[%lld].yS=%e >= wall[%lld].yE=%e\n", i, w[i].g.yS, i, w[i].g.yE);
				system("pause");
				doublereal temp = w[i].g.yE;
				w[i].g.yE = w[i].g.yS;
				w[i].g.yS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.zS - w[i].g.zE) > 1.0e-40) {
				printf("non union iso position on plane wall[%lld]: zS=%e zE=%e\n",i, w[i].g.zS, w[i].g.zE);
				system("pause");
				w[i].g.zS = w[i].g.zE;
				bOk = false;
			}
			break;
		case XZ :
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.xS >= w[i].g.xE) {
				printf("wall[%lld].xS=%e >= wall[%lld].xE=%e\n", i, w[i].g.xS, i, w[i].g.xE);
				system("pause");
				doublereal temp = w[i].g.xE;
				w[i].g.xE = w[i].g.xS;
				w[i].g.xS = temp;
				bOk = false;
			}
			if (w[i].g.zS >= w[i].g.zE) {
				printf("wall[%lld].zS=%e >= wall[%lld].zE=%e\n", i, w[i].g.zS, i, w[i].g.zE);
				system("pause");
				doublereal temp = w[i].g.zE;
				w[i].g.zE = w[i].g.zS;
				w[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.yS - w[i].g.yE) > 1.0e-40) {
				printf("non union iso position on plane wall[%lld]: yS=%e yE=%e\n", i, w[i].g.yS, w[i].g.yE);
				system("pause");
				w[i].g.yS=w[i].g.yE;
				bOk = false;
			}
			break;
		case YZ :
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.yS >= w[i].g.yE) {
				printf("wall[%lld].yS=%e >= wall[%lld].yE=%e\n", i, w[i].g.yS, i, w[i].g.yE);
				system("pause");
				doublereal temp = w[i].g.yE;
				w[i].g.yE = w[i].g.yS;
				w[i].g.yS = temp;
				bOk = false;
			}
			if (w[i].g.zS >= w[i].g.zE) {
				printf("wall[%lld].zS=%e >= wall[%lld].zE=%e\n", i, w[i].g.zS, i, w[i].g.zE);
				system("pause");
				doublereal temp = w[i].g.zE;
				w[i].g.zE = w[i].g.zS;
				w[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.xS - w[i].g.xE) > 1.0e-40) {
				printf("non union iso position on plane wall[%lld]: xS=%e xE=%e\n", i, w[i].g.xS, w[i].g.xE);
				system("pause");
				w[i].g.xS = w[i].g.xE;
				bOk = false;
			}
			break;
		}
		if (w[i].g.xE > b[0].g.xE) {
			printf("ERROR. Your model is incorrect.\n");
			printf("wall[%lld].g.xE=%e > cabinet.xE=%e.\n", i, w[i].g.xE, b[0].g.xE);
			system("pause");
			w[i].g.xE = b[0].g.xE; // исправление.
			bOk = false;
		}
		if (w[i].g.yE > b[0].g.yE) {
			printf("ERROR. Your model is incorrect.\n");
			printf("wall[%lld].g.yE=%e > cabinet.yE=%e.\n", i, w[i].g.yE, b[0].g.yE);
			system("pause");
			w[i].g.yE = b[0].g.yE; // исправление.
			bOk = false;
		}
		if (w[i].g.zE > b[0].g.zE) {
			printf("ERROR. Your model is incorrect.\n");
			printf("wall[%lld].g.zE=%e > cabinet.zE=%e.\n", i, w[i].g.zE, b[0].g.zE);
			system("pause");
			w[i].g.zE = b[0].g.zE; // исправление.
			bOk = false;
		}
		if (w[i].g.xS < b[0].g.xS) {
			printf("ERROR. Your model is incorrect.\n");
			printf("wall[%lld].g.xS=%e < cabinet.xS=%e.\n", i, w[i].g.xS, b[0].g.xS);
			system("pause");
			w[i].g.xS = b[0].g.xS; // исправление.
			bOk = false;
		}
		if (w[i].g.yS < b[0].g.yS) {
			printf("ERROR. Your model is incorrect.\n");
			printf("wall[%lld].g.yS=%e < cabinet.yS=%e.\n", i, w[i].g.yS, b[0].g.yS);
			system("pause");
			w[i].g.yS = b[0].g.yS; // исправление.
			bOk = false;
		}
		if (w[i].g.zS < b[0].g.zS) {
			printf("ERROR. Your model is incorrect.\n");
			printf("wall[%lld].g.zS=%e < cabinet.zS=%e.\n", i, w[i].g.zS, b[0].g.zS);
			system("pause");
			w[i].g.zS = b[0].g.zS; // исправление.
			bOk = false;
		}
	}

	for (integer i = 0; i < ls; i++) {

		// Ошибка не число.
		if (s[i].g.xS != s[i].g.xS) {
			printf("Error source[%lld] xS NOT VALUE=%e.\n", i, s[i].g.xS);
			system("pause");
			bOk = false;
		}
		if (s[i].g.xE != s[i].g.xE) {
			printf("Error source[%lld] xE NOT VALUE=%e.\n", i, s[i].g.xE);
			system("pause");
			bOk = false;
		}
		if (s[i].g.yS != s[i].g.yS) {
			printf("Error source[%lld] yS NOT VALUE=%e.\n", i, s[i].g.yS);
			system("pause");
			bOk = false;
		}
		if (s[i].g.yE != s[i].g.yE) {
			printf("Error source[%lld] yE NOT VALUE=%e.\n", i, s[i].g.yE);
			system("pause");
			bOk = false;
		}
		if (s[i].g.zS != s[i].g.zS) {
			printf("Error source[%lld] zS NOT VALUE=%e.\n", i, s[i].g.zS);
			system("pause");
			bOk = false;
		}
		if (s[i].g.zE != s[i].g.zE) {
			printf("Error source[%lld] zE NOT VALUE=%e.\n", i, s[i].g.zE);
			system("pause");
			bOk = false;
		}

		switch (s[i].iPlane) {
		case XY:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.xS >= s[i].g.xE) {
				printf("source[%lld].xS=%e >= source[%lld].xE=%e\n", i, s[i].g.xS, i, s[i].g.xE);
				system("pause");
				doublereal temp = s[i].g.xE;
				s[i].g.xE = s[i].g.xS;
				s[i].g.xS = temp;
				bOk = false;
			}
			if (s[i].g.yS >= s[i].g.yE) {
				printf("source[%lld].yS=%e >= source[%lld].yE=%e\n", i, s[i].g.yS, i, s[i].g.yE);
				system("pause");
				doublereal temp = s[i].g.yE;
				s[i].g.yE = s[i].g.yS;
				s[i].g.yS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.zS - s[i].g.zE) > 1.0e-40) {
				printf("non union iso position on plane source[%lld]: zS=%e zE=%e\n", i, s[i].g.zS, s[i].g.zE);
				system("pause");
				s[i].g.zS = s[i].g.zE;
				bOk = false;
			}
			break;
		case XZ:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.xS >= s[i].g.xE) {
				printf("source[%lld].xS=%e >= source[%lld].xE=%e\n", i, s[i].g.xS, i, s[i].g.xE);
				system("pause");
				doublereal temp = s[i].g.xE;
				s[i].g.xE = s[i].g.xS;
				s[i].g.xS = temp;
				bOk = false;
			}
			if (s[i].g.zS >= s[i].g.zE) {
				printf("source[%lld].zS=%e >= source[%lld].zE=%e\n", i, s[i].g.zS, i, s[i].g.zE);
				system("pause");
				doublereal temp = s[i].g.zE;
				s[i].g.zE = s[i].g.zS;
				s[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.yS - s[i].g.yE) > 1.0e-40) {
				printf("non union iso position on plane source[%lld]: yS=%e yE=%e\n", i, s[i].g.yS, s[i].g.yE);
				system("pause");
				s[i].g.yS = s[i].g.yE;
				bOk = false;
			}
			break;
		case YZ:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.yS >= s[i].g.yE) {
				printf("source[%lld].yS=%e >= source[%lld].yE=%e\n", i, s[i].g.yS, i, s[i].g.yE);
				system("pause");
				doublereal temp = s[i].g.yE;
				s[i].g.yE = s[i].g.yS;
				s[i].g.yS = temp;
				bOk = false;
			}
			if (s[i].g.zS >= s[i].g.zE) {
				printf("source[%lld].zS=%e >= source[%lld].zE=%e\n", i, s[i].g.zS, i, s[i].g.zE);
				system("pause");
				doublereal temp = s[i].g.zE;
				s[i].g.zE = s[i].g.zS;
				s[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.xS - s[i].g.xE) > 1.0e-40) {
				printf("non union iso position on plane source[%lld]: xS=%e xE=%e\n", i, s[i].g.xS, s[i].g.xE);
				system("pause");
				s[i].g.xS = s[i].g.xE;
				bOk = false;
			}
			break;
		}
		if (s[i].g.xE > b[0].g.xE) {
			printf("ERROR. Your model is incorrect.\n");
			printf("source[%lld].g.xE=%e > cabinet.xE=%e.\n", i, s[i].g.xE, b[0].g.xE);
			system("pause");
			s[i].g.xE = b[0].g.xE; // исправление.
			bOk = false;
		}
		if (s[i].g.yE > b[0].g.yE) {
			printf("ERROR. Your model is incorrect.\n");
			printf("source[%lld].g.yE=%e > cabinet.yE=%e.\n", i, s[i].g.yE, b[0].g.yE);
			system("pause");
			s[i].g.yE = b[0].g.yE; // исправление.
			bOk = false;
		}
		if (s[i].g.zE > b[0].g.zE) {
			printf("ERROR. Your model is incorrect.\n");
			printf("source[%lld].g.zE=%e > cabinet.zE=%e.\n", i, s[i].g.zE, b[0].g.zE);
			system("pause");
			s[i].g.zE = b[0].g.zE; // исправление.
			bOk = false;
		}
		if (s[i].g.xS < b[0].g.xS) {
			printf("ERROR. Your model is incorrect.\n");
			printf("source[%lld].g.xS=%e < cabinet.xS=%e.\n", i, s[i].g.xS, b[0].g.xS);
			system("pause");
			s[i].g.xS = b[0].g.xS; // исправление.
			bOk = false;
		}
		if (s[i].g.yS < b[0].g.yS) {
			printf("ERROR. Your model is incorrect.\n");
			printf("source[%lld].g.yS=%e < cabinet.yS=%e.\n", i, s[i].g.yS, b[0].g.yS);
			system("pause");
			s[i].g.yS = b[0].g.yS; // исправление.
			bOk = false;
		}
		if (s[i].g.zS < b[0].g.zS) {
			printf("ERROR. Your model is incorrect.\n");
			printf("source[%lld].g.zS=%e < cabinet.zS=%e.\n", i, s[i].g.zS, b[0].g.zS);
			system("pause");
			s[i].g.zS = b[0].g.zS; // исправление.
			bOk = false;
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


// Параметры которые используются для 
// настройки модели Смагоринского.
typedef struct TSMAGORINSKYINFO {
	doublereal Cs; // константа Смагоринского.
	bool bDynamic_Stress; // Для определения константы Смагоринского Cs используется динамическая модель Германо.
	bool bLimiters_Cs; // ограничивать ли постоянную Смагоринского ?
	doublereal minCs, maxCs; // минимальное и максимальное значения константы Смагоринского.
	integer itypeFiltrGermano; // тип фильтра который используется для осреднения в модели Германо.
	doublereal roughness; // значение шероховатости на твёрдой неподвижной стенке.
						  // показатель степени для учёта шероховатости на стенке.
	integer ipowerroughness; // может принимать значения только 1 или 2.
	bool bfdelta; // использовать ли поправку дающую улучшение на неравномерной сетке ?
	bool bSmagorinsky_Lilly; // использовать ли модель Смагоринского-Лиллу ?
	bool bsurface_roughness; // использовать ли поправку учитывающую шероховатость стенки ?
	bool bSelectiveSmagorinsky; // использовать ли Selective Smagorinsky Model ?
	integer itypeFILTRSelectiveSmagorinsky; // тип фильтра который используется для осреднения в модели Selective Smagorinsky.
	doublereal SSangle; // угол между вихрем и осреднённым вихрем в модели Selective Smagorinsky;
	bool bRichardsonCorrect; // использовать ли поправку связанную с числом Ричардсона для течений с кривизной линий тока.
	doublereal rRichardsonMultiplyer; // коэффициент в поправочной формуле связанной с кривизной линий тока.
} SMAGORINSKYINFO;

// Для полилинейного метода:
typedef struct TNODELR {
	integer id=-1; // идентификатор внутреннего узла
	struct TNODELR *next = NULL; // ссылка на следующий узел или NULL
} NODELR;

typedef struct TNODELR_BASE {
	integer ilineid=-1; // идентификатор сеточной линии
	integer iN=-1; // количество узловых точек включая граничные
	struct TNODELR *root = NULL; // корень новой сеточной линии
	struct TNODELR_BASE *next = NULL; // указатель на следующую сеточную линию или NULL

									  // Особый случай :
									  // обработка плоского бесконечно 
									  // тонкого источника тепла.
									  // Переменные равны истине если связь с источником 
									  // прервана (в случае если воздух граничит с плоским источником).
	bool bNeimanStart=false;
	bool bNeimanEnd=false;
} NODELR_BASE;


typedef struct TBOUND {
	// граничный узел и узлы вглубь расчётной области
	integer iB=-1, iI=-1, iII=-1;
	integer iI1=-1, iI2=-1;
	integer Norm=-10000; // внутренняя нормаль
				  // marker boundary:
				  /*
				  * О маркере границы:
				  * if (MCB < ls) значит источник с номером MCB;
				  * else if (MCB < ls+lw) значит стенка c номером MCB-ls;
				  * else значение по умолчанию. MCB==ls+lw.
				  */
	integer MCB=-1; // marker boundary

				 // соседи на границе области,
				 // их позиции нужны в матрице:
				 // На этих позициях будут стоять нули.
				 // Если позиция отсутствует, то стоит -1.
				 // Всего возможно не более 4 позиций, но 
				 // здесь сделано больше из-за 
				 // алгоритмического удобства.
				 // Позиции соответствуют сторонам света.
	integer iW[6] = {-1, -1, -1, -1, -1, -1 };

	// Для внутреннего источника указывает
	// на сосендний внутренний узел с обратной стороны
	// источника противоположной текущей стороне, или на -1
	// если источник на границе hollow блока.
	// TODO 6 мая 2016.

	// для теплообмена излучением хранит emissivity.
	doublereal emissivity=0.8; // излучательная способность границы.

						   // Площадь грани.
						   // Это новое поле введённое лишь 20 сентября 2016.
	doublereal dS=0.0;

	// координаты центра грани.
	TOCHKA p_c;

} BOUND;


// информация о жидкой зоне
typedef struct TFLOWINFO {
	// опорная точка
	doublereal xc=1.0e30, yc=1.0e30, zc=1.0e30;
	// нужно ли расчитывать поле течения
	integer iflow=0;
	// режим течения : 0 - ламинарный, 1 - турбулентный
	integer iflowregime=0;
	// модель турбулентности
	// 0 - Zero Equation Turbulence Model (RANS).
	// 1 - модель Смагоринского (LES). (как опция модели, модель Германо 1991 (LES)).
	// 2 - Spalart-Allmares (RANS). 19.04.2019
	integer iturbmodel=0;
	// параметры модели Смагоринского:
	doublereal Cs; // постоянная Смагоринского.
	bool bDynamic_Stress = false; // если == true то включает динамическую модель Германо для определения квадрата постоянной Смагоринского.
	bool bLimiters_Cs = false; // включает ограничение на постоянную Смагоринского : ограничение на минимальное и максимальное значения.
	doublereal rminCs, rmaxCs; // минимальное и максимальное ограничение для константы Смагоринского.
	integer itypeFiltrGermano; // тестовый фильтр который используется для осреднения в модели Германо 1991 года.
	// поправка на неравномерной сетке,
	// модель Смагоринского-Лиллу,
	// учёт шероховатости стенки.
	integer bfdelta, bSmagorinskyLilly, bsurface_roughness;
	doublereal roughness; // шероховатость стенки.
	integer ipowerroughness; // показатель степени в модели учёта шероховатости стенки.
	doublereal rRi_mult; // корректирующий множитель подправляющий турбулентное число Ричардсона.
	doublereal rSelectiveAngle; // пороговое значение угла в модели Selective Smagorinsky.
	integer itypeSelectiveSmagorinsky_filtr; // тип фильтра в модели Смагоринского.
	// поправка на течения с кривизной линий тока,
	// избирательная модель Смагоринского.
	bool bSwirlAmendment = false, bSelectiveSmagorinsky=false;

	// Spalart-Allmares (RANS). 19.04.2019
	// Эмпирические константы, определяющие SA модель турбулентности.
	doublereal sigma_nu = 2.0 / 3.0;
	doublereal c_b1 = 0.1355;
	doublereal c_b2 = 0.622;
	doublereal karman = 0.41;
	doublereal c_w1 = (c_b1/(karman*karman))+((1.0+c_b2)/sigma_nu);
	doublereal c_w2 = 0.3;
	doublereal c_w3 = 2.0; 
	doublereal c_nu1 = 7.1;
	doublereal C_t3 = 1.2;
	doublereal C_t4 = 0.5;

} FLOWINFO;

// вся информация о полном наборе решаемых уравнений
typedef struct TEQUATIONINFO {
	// нужно ли решать уравнение теплопроводности
	// 0 - ненужно, 1 - нужно.
	integer itemper=1;
	// максимальное количество изолированных по жидкости
	// жидких зон (FLUID interior)
	integer imaxflD=0;
	// база данных с информацией о жидких зонах.
	FLOWINFO* fluidinfo=NULL;
} EQUATIONINFO;

EQUATIONINFO eqin; // информация о наборе решаемых уравнений

struct Tdatabase {
	doublereal *x = NULL, *y = NULL, *z = NULL; // координаты узлов.
	integer maxelm=0;
	integer** nvtxcell = NULL;
	integer ncell=0;
	// связь теплопередачи с гидродинамикой.
	integer **ptr = NULL;// для тестирования алгебраического многосеточного метода
};

// одна строка матрицы СЛАУ в 3D варианте
// для внутреннего КО.
typedef struct Tequation3D {
	doublereal ap=0.0, ae=0.0, an=0.0, aw=0.0, as=0.0, at=0.0, ab=0.0, b=0.0;
	integer iP = -1, iE = -1, iN = -1, iT = -1, iW = -1, iS = -1, iB = -1;
	// Расширение структуры данных под АЛИС сетку.
	// На АЛИС сетке шаблон имеет переменное число связей и 
	// их число колеблется от семиточечного шаблона до
	// 25 точечного максимально.
	// На АЛИС сетке даже чисто диффузионная матрица не имеет диагонального преобладания и 
	// для решения СЛАУ с такой матрицей нужны специальные "робастые" методы типа BiCGStab + ILU2.
	// Признак существования дополнительных связей:
	// true если коэффициент существует.
	bool bE2=false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;
	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;
	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;
	// Значение матричного коэффициента:
	doublereal ae2 = 0.0, an2 = 0.0, aw2 = 0.0, as2 = 0.0, at2 = 0.0, ab2 = 0.0;
	doublereal ae3 = 0.0, an3 = 0.0, aw3 = 0.0, as3 = 0.0, at3 = 0.0, ab3 = 0.0;
	doublereal ae4 = 0.0, an4 = 0.0, aw4 = 0.0, as4 = 0.0, at4 = 0.0, ab4 = 0.0;
	// индексация для дополнительных связей.
	integer iE2=-1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1;
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1;
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1;
	
} equation3D;

// одна строка матрицы СЛАУ в 3D варианте
// для граничного КО.
typedef struct Tequation3D_bon {
	doublereal aw = 0.0, ai = 0.0, b = 0.0; // wall, internal, правая часть
	integer iW = -1, iI = -1;

	// соседи на границе области,
	// их позиции нужны в матрице:
	// На этих позициях будут стоять нули.
	integer iW1 = -1, iW2 = -1, iW3 = -1, iW4 = -1;
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
	bool free_temper_level1=true;
	// флаг отвечающий за освобождение памяти
	// второго уровня. Когда матрица СЛАУ перезаписывается
	// в формат SIMPLESPARSE. Исходную матрицу в формате 
	// equation3D можно убрать из оперативной памяти компьютера.
	bool free_temper_level2=true;

	integer maxnod=0; // максимальный номер узла (размерность массива)
					// pa[0..maxnod-1];
	TOCHKA* pa = NULL; // координаты узлов сетки принадлежащие расчётной области

	integer maxelm=0; // число ненулевых контрольных объёмов
					// nvtx[0..7][0..maxelm-1]
	integer **nvtx = NULL; // список узлов для каждого элемента (ненулевого контрольного объёма)
						   // sosedi[0..11][0..maxelm-1]
						   //integer **sosedi; // соседние контрольные объёмы для каждого внутреннего контрольного объёма
						   // AliceMesh
	ALICE_PARTITION **sosedi = NULL;// соседние контрольные объёмы для каждого внутреннего контрольного объёма
	integer maxbound=0; // число граничных узлов
	integer maxp=0; // maxp == maxelm + maxbound;
				  // sosedb[0..maxbound-1];
	BOUND* sosedb = NULL; // граничные узлы расчётной области 
						  // для всех граничных КО. Равно истине если имеем дело с
						  // границей строго внутри расчётной области причём на ней 
						  // расположен именно плоский бесконечно тонкий источник и по 
						  // одну его сторону расположена жидкость а по другую твёрдое тело.
	bool* binternalsource = NULL;

	// какому блоку принадлежит внутренний КО
	integer* whot_is_block = NULL;

	integer **ptr = NULL; // Связь с гидродинамикой.

						  // для графической визуализации
	integer ncell=0; // количество связей для контрольных объёмов.
	integer **nvtxcell = NULL; // связи для контрольных объёмов.

							   // Для АЛИС сетки хранит номер уровня ячейки, это 
							   // потребуется при сборке матрицы.
							   // ilevel_alice[0..maxelm-1].
	integer *ilevel_alice = NULL;


	doublereal *potent = NULL; // массив узловых потенциалов (искомых функций)
	doublereal **total_deformation = NULL; // Полная деформация.

										   // Свойства материалов разделяются для  
										   // внутренних и для граничных КО.
										   // сначала prop[0..2][0..maxelm-1]
	doublereal **prop = NULL; // свойства материалов для внутренних КО.
							  // сначала prop_b[0..2][0..maxbound-1]
	doublereal **prop_b = NULL; // свойства для граничных КО.


	doublereal *Sc = NULL; // объёмное тепловыделение приходящееся на один выбранный контрольный объём.
	integer *ipower_time_depend = NULL; // закон зависимости мощности тепловыделения от времени.

	doublereal alpha=1.0; // параметр релаксации для температуры
	equation3D *slau = NULL; // коэффициенты матрицы СЛАУ для внутренних КО
	equation3D_bon *slau_bon = NULL; // коэффициенты матрицы СЛАУ для граничных КО

									 // для полилинейного метода :
									 // полилинейный метод рекомендован проф. Минесотского университета С. Патанкаром.
									 // полилинейный метод обладает фирменной особенностью - за первые несколько итераций 
									 // невязка падает на  несколько порядков. Это говорит о том что полилинейный метод 
									 // может быть использован как предобуславливатель в алгоритме Ван-Дер-Ворста - BiCGStab.
	NODELR_BASE *rootWE = NULL;
	NODELR_BASE *rootSN = NULL;
	NODELR_BASE *rootBT = NULL;

	integer iWE=-1000, iSN = -1000, iBT = -1000; // число сеточных линий вдоль каждого из направлений.

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
	doublereal resLR1sk=0.0; // O(h!3)

						 // Копия сеточных размеров для восстановления данных.
	integer inx_copy = 0, iny_copy = 0, inz_copy=0;
	doublereal operatingtemperature_copy=0.0;

	doublereal *xpos_copy = NULL, *ypos_copy = NULL, *zpos_copy = NULL;

	// 9 августа 2015.
	// Для распараллеливания 
	integer *ifrontregulationgl = NULL;
	integer *ibackregulationgl = NULL; // обратное преобразование.

	doublereal operatingtemperature = 20.0;

	Tdatabase database;

}  TEMPER;

typedef struct TFLOW {
	integer maxnod=0; // максимальный номер узла (размерность массива)
	integer maxelm = 0; // число внутренних контрольных объёмов
					// nvtx[0..7][0..maxelm-1]
	integer **nvtx = NULL; // список узлов для каждого внутреннего элемента (контрольного объёма)

						   // pa[0..maxnod-1]
	TOCHKA* pa = NULL; // координаты узлов сетки принадлежащие расчётной области

					   // sosedi[0..11][0..maxelm-1]
					   //integer **sosedi; // соседние контрольные объёмы для каждого внутреннего КО
					   // Для ALICEMESH сетки.
	ALICE_PARTITION **sosedi = NULL;// соседние контрольные объёмы для каждого внутреннего КО
	integer maxbound = 0; // число граничных КО
	integer maxp = 0; // число уравнений
				  // sosedb[0..maxbound-1];
	BOUND* sosedb = NULL; // граничные узлы расчётной области

	integer *ptr = NULL; // Связь с теплопроводностью


						 // какому блоку принадлежит внутренний КО
	integer* whot_is_block = NULL;

	// potent[iVar][0..maxp-1]
	doublereal **potent = NULL; // массив узловых потенциалов (искомых функций)
								// prop[0..2][0..maxelm-1]
	doublereal **prop = NULL; // свойства материалов
							  // prop_b[0..2][0..maxbound-1]
	doublereal **prop_b = NULL; // свойства материалов для граничных КО

	doublereal *alpha = NULL; // параметры нижней релаксации
	equation3D **slau = NULL; // коэффициенты матрицы СЛАУ для внутренних КО.
	equation3D_bon **slau_bon = NULL; // коэффициенты матрицы СЛАУ для граничных КО
									  // для реализации монотонизатора Рхи-Чоу требуется хранить диагональные коэффициенты.
	doublereal **diag_coef = NULL;
	doublereal OpTemp=0.0; // Operating Temperature


	bool bactive=true; // нужно-ли рассчитывать поле течения
	bool bPressureFix = false; // нужно ли фиксировать давление в одной точке
	bool bLR1free=false; // нужно ли применять плавающий полилинейный солвер (он показывает более быструю сходимость).

				   // для полилинейного метода :
				   // полилинейный метод рекомендован проф. Минесотского университета С. Патанкаром.
				   // полилинейный метод обладает фирменной особенностью - за первые несколько итераций 
				   // невязка падает на  несколько порядков. Это говорит о том что полилинейный метод 
				   // может быть использован как предобуславливатель в алгоритме Ван-Дер-Ворста - BiCGStab.

	integer iWE=-1000, iSN = -1000, iBT = -1000; // число сеточных линий вдоль каждого из направлений.
	integer** iN = NULL; //iN[3][max3(iWE,iSN,iBT)];
	integer*** id = NULL; //id[3][max3(iWE,iSN,iBT)][max(iN)]; 

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
	integer iflowregime= 0; // default LAMINAR
						 // Кратчайшее расстояние до ближайшей стенки [0..maxelm-1]
	doublereal* rdistWall = NULL; // расстояние до ближайшей твёрдой стенки.
								  // Толщина пограничного слоя в формуле Эскудиера
	doublereal rdistWallmax=0.0;
	// S инвариант тензора скоростей-деформаций
	doublereal* SInvariantStrainRateTensor = NULL; // [0..maxelm+maxbound-1]; // инициализируется нулём.

												   // массовый поток через грани КО :
	doublereal** mf = NULL;

	SMAGORINSKYINFO smaginfo; // параметры модели Смагоринского.

							  // выходная невязка для поправки давления
							  // согласованная с точностью аппроксимации уравнения.
	doublereal resICCG=0.0; // O(h!2)
	doublereal resLR1sk=0.0; // O(h!3)

						 // Для правильной работы mass balance для естественно конвективных задач
						 // нужен идентификатор разных по связности гидродинамических областей для
						 // внутренних КО.
	integer *icolor_different_fluid_domain = NULL; // Разные гидродинамические подобласти имеют разные цвета.

												   // 9 августа 2015.
												   // Для распараллеливания 
	integer *ifrontregulationgl = NULL;
	integer *ibackregulationgl = NULL; // обратное преобразование.

} FLOW;




// Объединение
// Для блочно структурированной расчётной сетки.
// Начало разработки 25.04.2018.
typedef struct TUNION {
	// id передаётся из интерфейса.
	integer id; // Уникальный номер объединения.

				// Из кабинета юнион видится как Hollow блок
				// в виде прмой прямоугольной призмы.
				// размеры передаются из интерфейса.
	doublereal xS, xE, yS, yE, zS, zE;

	// Внутренняя сетка union.
	// Union является кабинетом для своих внутренних блоков.
	// для внутреннего пользования.
	doublereal *xpos = NULL, *ypos = NULL, *zpos = NULL;
	doublereal *xposadd = NULL, *yposadd = NULL, *zposadd = NULL;

	// Для внутренней сетки (размерности).
	// Передаётся из интерфейса.
	integer inx, iny, inz;
	// для внутреннего пользования.
	integer inxadd = -1, inyadd = -1, inzadd = -1;

	//  Тип сеточного генератора
	// передаётся из интерфейса.
	integer iswitchMeshGenerator = 2; // 2 - CoarseMeshGen

									  // Локальные объявления.
	TEMPER t;
	integer flow_interior; // Суммарное число FLUID зон
	FLOW* f = NULL;
} UNION;



// считывание параметров из 
// входного файла premeshin.txt
void premeshin(const char *fname, integer &lmatmax, integer &lb, integer &ls, integer &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w, 
	           doublereal &dgx, doublereal &dgy, doublereal &dgz, integer &inx, integer &iny, integer &inz, doublereal &operatingtemperature, 
			   integer &ltdp, TEMP_DEP_POWER* &gtdps, integer &lu, UNION* &my_union) {


#ifdef MINGW_COMPILLER
	// eqin - информация о наборе решаемых уравнений.

	// dgx, dgy, dgz - вектор силы тяжести.
	// inx, iny, inz - количество точек по каждой из осей.

	FILE *fp;
	errno_t err1=0;
	fp=fopen64(fname, "r");
	if (err1 != 0) {
		printf("No input File premeshin.txt \n");
		//system("PAUSE");
		system("pause");

	}
	else
	{
		if (fp != NULL) {
			float fin = 0.0;
			integer din = 0;
			doublereal scale = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания

			fscanf(fp, "%d", &din);
			ionly_solid_visible=din;
			printf("ionly_solid_visible =%d\n", ionly_solid_visible);
			fscanf(fp, "%f", &fin);
			scale = fin;
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
			idirectional_for_XY_Plot = din;

			fscanf(fp, "%f", &fin);
			etalon_max_size_ratio = fin; // подробность расчётной сетки.

			fscanf(fp, "%f", &fin);
			etalon_max_size_ratio2 = fin; // Критерий качества расчётной сетки на основе FlowVision.

			fscanf(fp, "%d", &din);
			switch (din) {
			case 0: bsnap_TO_global = false;  break;
			default : bsnap_TO_global = true;  break;
			}
			

			fscanf(fp, "%d", &din);
			iswitchsolveramg_vs_BiCGstab_plus_ILU2 = din; // Выбор решающего устройства : либо amg1r5 либо BiCGStab+ILU2.

			fscanf(fp, "%d", &din);
			iswitchsolveramg_vs_BiCGstab_plus_ILU6 = din; // Выбор решающего устройства : либо РУМБА0.14 либо BiCGStab+ILU6.

			fscanf(fp, "%d", &din);
			if (din == 1) {
				// SIMPLEC algorithm.
				iSIMPLE_alg = SIMPLEC_Van_Doormal_and_Raithby;
			}
			else {
				// SIMPLE algorithm 1972.
				iSIMPLE_alg = SIMPLE_Carretto;
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
			default: iFLOWScheme = 2; break; // UDS самая стабильная схема.
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
			default: iTEMPScheme = 2; break; // UDS самая стабильная схема.
			}


			// Выбор сеточного генератора.
			fscanf(fp, "%d", &din);
			iswitchMeshGenerator = din;


			fscanf(fp, "%d", &din);
			steady_or_unsteady_global_determinant = 2;
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 5) || (din == 6)|| (din == 7) || (din == 8) || (din == 9)) {
				// 0 - thermal only steady state calculation,
				// 1 - thermal only unsteady calculation,
				// 2 - mesh generator only.
				// 3 - fluid dynamic steady state.
				// 5 - Static Structural (Thermal solver #2)
				// 6 - Thermal Stress
				// 7 - Unsteady thermal solver #2
				// 8 - Visualisation only
				// 9 - cfd unsteady fluid dynamic.
				steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady -1 calculation.
			}
			else {
				printf("error input parametr steady or unsteady calculation\n");
				system("PAUSE");
				exit(1);
			}

			fscanf(fp, "%d", &din);
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3)) {
				glTSL.id_law = din;
			}
			else {
				printf("error input parametr timestep law\n");
				system("PAUSE");
				exit(1);
			}
			fscanf(fp, "%f", &fin);
			glTSL.Factor_a_for_Linear = fin; // Factor_a
			if ((fin<=0.0)||(fin>=1.0)) {
				printf("error input parametr timestep law Factor a\n");
				system("PAUSE");
				exit(1);
            }
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law tau must be strongly positive\n");
				system("PAUSE");
				exit(1);
            }
			glTSL.tau = fin; // длительность импульса.
			fscanf(fp, "%f", &fin);
			glTSL.Q = fin;  // Скважность.
			// Параметры импульсного режима для темы АППАРАТ.
			fscanf(fp, "%f", &fin);
			if ((fin<=0.0)||(fin>=1.0)) {
				printf("error input parametr timestep law APPARAT multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law APPARAT tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;			
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law APPARAT tau2 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau2 = fin;
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law APPARAT tau_pause must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau_pause = fin;
			fscanf(fp, "%d", &din);
			glTSL.n_cycle = din;
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law APPARAT Period must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.T_all = fin;
			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl <= 0.0) {
				printf("error in parameters Square Wave APPARAT time step law.\n");
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
			adiabatic_vs_heat_transfer_coeff = din;  // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
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
			itype_ALICE_Mesh=din;
			fscanf(fp, "%d", &din);
			my_amg_manager.m_restart = din;
			// classical algebraic multigrid parameters:
			// only for my_agregat_amg.cu.
			fscanf(fp, "%d", &din);
			my_amg_manager.imySortAlgorithm = din;
			fscanf(fp, "%d", &din);
			//my_amg_manager.maximum_levels = din;
			my_amg_manager.maximum_delete_levels_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Stress = din;

			// type interpolation procedure :
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
			my_amg_manager.iCFalgorithm_and_data_structure_Temperature = din;// 3-Treap.
			fscanf(fp, "%d", &din);
			my_amg_manager.iCFalgorithm_and_data_structure_Speed = din;// 3-Treap.
			fscanf(fp, "%d", &din);
			my_amg_manager.iCFalgorithm_and_data_structure_Pressure = din;// 3-Treap.
			fscanf(fp, "%d", &din);
			my_amg_manager.iCFalgorithm_and_data_structure_Stress = din;// 3-Treap.

			fscanf(fp, "%d", &din);
			//my_amg_manager.itypemodifyinterpol = din;
			//my_amg_manager.baglomeration_with_consistency_scaling = din;
			my_amg_manager.bdiagonal_dominant = din;
			fscanf(fp, "%d", &din);
			//my_amg_manager.inumberadaptpass = din;

			// 23.02.2018
			// print matrix portrait
			fscanf(fp, "%d", &din);
			my_amg_manager.bTemperatureMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
			fscanf(fp, "%d", &din);
			my_amg_manager.bSpeedMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
			fscanf(fp, "%d", &din);
			my_amg_manager.bPressureMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
			fscanf(fp, "%d", &din);
			my_amg_manager.bStressMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.

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

			// number nFinnest sweeps :
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

			// number postsweeps :
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.nu2_Stress = din;

			// memory size :
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Temperature = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Speed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Pressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.memory_size_Stress = din;


			// Параметр верхней релаксации в сглаживателе.
			fscanf(fp, "%f", &fin);
			my_amg_manager.gold_const_Temperature = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.gold_const_Speed = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.gold_const_Pressure = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.gold_const_Stress = fin;

			// использовать ли ilu2 smoother.
			fscanf(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Temperature = din;
				my_amg_manager.b_gmresTemp = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Temperature = din;
			}
			fscanf(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Speed = din;
				my_amg_manager.b_gmresSpeed = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Speed = din;
			}			
			fscanf(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Pressure = din;
				my_amg_manager.b_gmresPressure = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Pressure = din;
			}			
			fscanf(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Stress = din;
				my_amg_manager.b_gmresStress = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Stress = din;
			}
			

			// strength threshold :
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Temperature = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Speed = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Pressure = fin;
			fscanf(fp, "%f", &fin);
			my_amg_manager.theta_Stress = fin;

			// magic threshold :
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
			// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
			// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
			fscanf(fp, "%d", &din);
			my_amg_manager.icoarseningTemp = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.icoarseningSpeed = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.icoarseningPressure = din;
			fscanf(fp, "%d", &din);
			my_amg_manager.icoarseningStress = din;

			// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
			// Если din==1, Stabilization BiCGStab.
			// 8.01.2017 Метод ван дер Ворста BiCGStab 
			// предобусловленный алгебраичесеким многосеточным методом.
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
			w = new WALL[lw];
			integer i = 0; // счётчик цикла for

			for (i = 0; i < ltdp; i++) {
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
			for (i = 0; i < lmatmax; i++) {
				// свойства материалов:
				// плотность
				fscanf(fp, "%f", &fin);
				matlist[i].rho = fin;
				// теплоёмкость при постоянном давлении
				//fscanf(fp, "%f", &fin);
				//matlist[i].cp = fin;
				fscanf(fp, "%d", &din);
				matlist[i].n_cp = din;
				matlist[i].arr_cp = NULL;
				matlist[i].temp_cp = NULL;
				matlist[i].arr_cp = new doublereal[matlist[i].n_cp];
				matlist[i].temp_cp= new doublereal[matlist[i].n_cp];
				if (matlist[i].temp_cp == NULL) {
					printf("problem memory allocation for temp_cp\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_cp == NULL) {
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
				matlist[i].arr_lam = NULL;
				matlist[i].temp_lam = NULL;
				matlist[i].arr_lam = new doublereal[matlist[i].n_lam];
				matlist[i].temp_lam = new doublereal[matlist[i].n_lam];
				if (matlist[i].temp_lam == NULL) {
					printf("problem memory allocation for temp_lam\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_lam == NULL) {
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
				// ортотропность теплопроводности :
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x=fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y=fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z=fin;
				// 5.08.2017.
				// Коэффициенты для задачи упругости.
				// Модуль Юнга и коэффициент Пуассона.
				doublereal Poissonratio = 0.154;
				doublereal Youngmodule = 217.5e9;
				fscanf(fp, "%f", &fin);
				Poissonratio = fin;
				fscanf(fp, "%f", &fin);
				Youngmodule = fin*1e9;
				fscanf(fp, "%f", &fin);
				matlist[i].beta_t_solid = fin*1E-6;
				// Коэффициенты Лямэ.
				doublereal E1_koef = Youngmodule / (1.0- Poissonratio*Poissonratio);
				doublereal nu1_koef = Poissonratio / (1.0- Poissonratio);
				matlist[i].mu_Lame = E1_koef/(2.0*(1.0+ nu1_koef));
				matlist[i].lambda_Lame = (E1_koef*nu1_koef)/(1.0- nu1_koef*nu1_koef);
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
				printf("%e %e %e %e %e\n", matlist[i].rho, matlist[i].arr_cp[0], matlist[i].arr_lam[0], matlist[i].mu, matlist[i].beta_t);
				printf("%d %d %d\n", matlist[i].blibmat, matlist[i].ilibident, matlist[i].ilawmu); // bBoussinesq не печатается
				printf("%e %e %e %e %e %e\n", matlist[i].mumin, matlist[i].mumax, matlist[i].Amu, matlist[i].Bmu, matlist[i].Cmu, matlist[i].degreennmu);
			}

			// считывание блоков
			for (i = 0; i<lb; i++) {

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
					case XY:
						b[i].g.zC += b[i].g.Hcyl;
						break;
					case XZ:
						b[i].g.yC += b[i].g.Hcyl;
						break;
					case YZ:
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
						case XY:
							b[i].g.zi[i73] += b[i].g.hi[i73];
							break;
						case XZ:
							b[i].g.yi[i73] += b[i].g.hi[i73];
							break;
						case YZ:
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
					case XY:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmin53 + b[i].g.hi[0];
						break;
					case XZ:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymin53 + b[i].g.hi[0];
						break;
					case YZ:
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.xS = xmin53;
						b[i].g.xE = xmin53 + b[i].g.hi[0];
						break;
					}
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
				b[i].radiation.nodelistW = NULL;
				b[i].radiation.nodelistE = NULL;
				b[i].radiation.nodelistS = NULL;
				b[i].radiation.nodelistN = NULL;
				b[i].radiation.nodelistB = NULL;
				b[i].radiation.nodelistT = NULL;
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
				b[i].arr_Sc = NULL;
				b[i].temp_Sc = NULL;
				b[i].arr_Sc = new doublereal[b[i].n_Sc];
				b[i].temp_Sc = new doublereal[b[i].n_Sc];
				if (b[i].temp_Sc == NULL) {
					printf("problem memory allocation for temp_Sc\n");
					system("pause");
					exit(1);
				}
				if (b[i].arr_Sc == NULL) {
					printf("problem memory allocation for arr_Sc\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
					// Температура в C.
					fscanf(fp, "%f", &fin);
					b[i].temp_Sc[i_4] = fin;
					fscanf(fp, "%f", &fin);
					if (fin != fin) {
						b[i].arr_Sc[i_4] = 0.0;
					}
					else {
						b[i].arr_Sc[i_4] = fin;
					}
				}

                // стиль зависимости мощности тепловыделения в блоке от времени.
				fscanf(fp, "%d", &din);
				b[i].ipower_time_depend=din;
				// тип блока
				fscanf(fp, "%d", &din);
				b[i].itype = din;

				// печать считанных значений на консоль
				printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				printf("%d %d %d\n", b[i].imatid,  b[i].itype, b[i].ipower_time_depend);
				printf("temperature depend power\n");
				printf("t_C power_W\n");
				for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
				}
			}

			// считывание источников тепла
			for (i = 0; i < ls; i++) {

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
						printf("calculate initial power=%e\n", s[i].power);
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
				case XY: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				printf("%e %d %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
			}

			// считывание твёрдых стенок



			for (i = 0; i < lw; i++) {

				fscanf(fp, "%d", &din);
				w[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				fscanf(fp, "%d", &din);
				w[i].ifamily = din;
				switch (din) {
				case 1:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// termostability wall
					w[i].emissivity=0.0;
					w[i].film_coefficient=0.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // первого рода
				case 2:  fscanf(fp, "%f", &fin);
					w[i].Tamb = 0.0;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// adiabatic wall
					w[i].emissivity=0.0;
					w[i].film_coefficient=0.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // однородное условие Неймана
				case 3:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// Newton-Richman condition, film coefficient.
					w[i].emissivity=0.0;
					w[i].film_coefficient=fin;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Ньютон-Рихман.
				case 4:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// Stefan - Bolcman condition
					w[i].emissivity=fin;
					w[i].film_coefficient=0.0;
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
				w[i].ithermal_Stress_boundary_condition = din;
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
				printf("%d %e %e %d %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
			}


			// АСЕМБЛЕСЫ.
			fscanf(fp, "%d", &din);
			lu = din;
			if (lu == 0) {
				my_union = NULL;
			}
			else {
				my_union = new UNION[lu];
				// инициализация.
				for (i = 0; i < lu; i++) {
					my_union[i].f = NULL;
					my_union[i].xpos = NULL;
					my_union[i].ypos = NULL;
					my_union[i].zpos = NULL;
					my_union[i].xposadd = NULL;
					my_union[i].yposadd = NULL;
					my_union[i].zposadd = NULL;
					my_union[i].iswitchMeshGenerator = 2; // 2 - CoarseMeshGen
					my_union[i].inxadd = -1;
					my_union[i].inyadd = -1;
					my_union[i].inzadd = -1;
					my_union[i].flow_interior = 0;
				}
			}
			for (i = 0; i < lu; i++) {
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
				eqin.fluidinfo = NULL;
			}
			else
			{
				// выделение оперативной памяти
				if (eqin.fluidinfo != NULL) {
					delete eqin.fluidinfo;
					eqin.fluidinfo = NULL;
				}
				eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
				for (i = 0; i < eqin.imaxflD; i++) {
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
					eqin.fluidinfo[i].iflowregime = din;
					fscanf(fp, "%d", &din);
					eqin.fluidinfo[i].iturbmodel = din;
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
					pfpir.idir = din;

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					AMG1R6_LABEL = din;

				}
			}

			fclose(fp); // закрытие файла
		}
	
	}
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
		printf("No input File premeshin.txt \n");
		//system("PAUSE");
		system("pause");

	}
	else
	{
		if (fp != NULL) {
			float fin = 0.0;
			integer din = 0;
			doublereal scale = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания



			fscanf_s(fp, "%d", &din);
			ionly_solid_visible = din;
			printf("ionly_solid_visible =%d\n", ionly_solid_visible);
			fscanf_s(fp, "%f", &fin);
			scale = fin;
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
			idirectional_for_XY_Plot = din;

			fscanf_s(fp, "%f", &fin);
			etalon_max_size_ratio = fin; // подробность расчётной сетки.

			fscanf_s(fp, "%f", &fin);
			etalon_max_size_ratio2 = fin; // Критерий качества расчётной сетки на основе FlowVision.

			//printf("etalon_max_size_ratio=%e etalon_max_size_ratio2=%e\n", etalon_max_size_ratio, etalon_max_size_ratio2);
			//system("PAUSE");

			fscanf_s(fp, "%d", &din);
			switch (din) {
			case 0: bsnap_TO_global = false;  break;
			default: bsnap_TO_global = true;  break;
			}

			fscanf_s(fp, "%d", &din);
			iswitchsolveramg_vs_BiCGstab_plus_ILU2 = din; // Выбор решающего устройства : либо amg1r5 либо BiCGStab+ILU2.

			fscanf_s(fp, "%d", &din);
			iswitchsolveramg_vs_BiCGstab_plus_ILU6 = din; // Выбор решающего устройства : либо РУМБА0.14 либо BiCGStab+ILU6.

			fscanf_s(fp, "%d", &din);
			if (din == 1) {
				// SIMPLEC algorithm.
				iSIMPLE_alg = SIMPLEC_Van_Doormal_and_Raithby;
			}
			else {
				// SIMPLE algorithm 1972.
				iSIMPLE_alg = SIMPLE_Carretto;
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
			default: iFLOWScheme = 2; break; // UDS самая стабильная схема.
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
			default: iTEMPScheme = 2; break; // UDS самая стабильная схема.
			}



			// Выбор сеточного генератора.
			fscanf_s(fp, "%d", &din);
			iswitchMeshGenerator = din;

			fscanf_s(fp, "%d", &din);
			steady_or_unsteady_global_determinant = 2;
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 5) || (din == 6) || (din == 7) || (din == 8) || (din == 9)) {
				// 0 - thermal only steady state calculation,
				// 1 - thermal only unsteady calculation,
				// 2 - mesh generator only.
				// 3 - fluid dynamic steady state.
				// 5 - Static Structural (Thermal solver #2)
				// 6 - Thermal Stress
				// 7 - Unsteady thermal solver #2
				// 8 - Visualisation only
				// 9 - cfd unsteady fluid dynamic.
				steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady -1 calculation.
			}
			else {
				printf("error input parametr steady or unsteady calculation\n");
				system("PAUSE");
				exit(1);
			}

			fscanf_s(fp, "%d", &din);
			if ((din == 0) || (din == 1) || (din == 2) || (din == 3)) {
				glTSL.id_law = din;
			}
			else {
				printf("error input parametr timestep law\n");
				system("PAUSE");
				exit(1);
			}
			fscanf_s(fp, "%f", &fin);
			if ((fin <= 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law Factor a\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.Factor_a_for_Linear = fin; // Factor_a
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law tau must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau = fin; // длительность импульса.
			fscanf_s(fp, "%f", &fin);
			glTSL.Q = fin;  // Скважность.
			// Параметры импульсного режима для темы АППАРАТ.
			fscanf_s(fp, "%f", &fin);
			if ((fin <= 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law APPARAT multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law APPARAT tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law APPARAT tau2 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau2 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law APPARAT tau_pause must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau_pause = fin;
			fscanf_s(fp, "%d", &din);
			glTSL.n_cycle = din;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law APPARAT Period must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.T_all = fin;
			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle * (2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl <= 0.0) {
				printf("error in parameters Square Wave APPARAT time step law.\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}

			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
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
			adiabatic_vs_heat_transfer_coeff = din; // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
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
			itype_ALICE_Mesh = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.m_restart = din;
			// classical algebraic multigrid parameters:
			// only for my_agregat_amg.cu.
			// only for my_agregat_amg.cu.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.imySortAlgorithm = din;
			fscanf_s(fp, "%d", &din);
			//my_amg_manager.maximum_levels = din;
			my_amg_manager.maximum_delete_levels_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.maximum_delete_levels_Stress = din;

			// type interpolation procedure :
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
			my_amg_manager.iCFalgorithm_and_data_structure_Temperature = din;// 3-Treap.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iCFalgorithm_and_data_structure_Speed = din;// 3-Treap.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iCFalgorithm_and_data_structure_Pressure = din;// 3-Treap.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.iCFalgorithm_and_data_structure_Stress = din;// 3-Treap.

			fscanf_s(fp, "%d", &din);
			//my_amg_manager.itypemodifyinterpol = din;
			//my_amg_manager.baglomeration_with_consistency_scaling = din;
			my_amg_manager.bdiagonal_dominant = din;
			fscanf_s(fp, "%d", &din);
			//my_amg_manager.inumberadaptpass = din;


			// 23.02.2018
			// print matrix portrait
			fscanf(fp, "%d", &din);
			my_amg_manager.bTemperatureMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
			fscanf(fp, "%d", &din);
			my_amg_manager.bSpeedMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
			fscanf(fp, "%d", &din);
			my_amg_manager.bPressureMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
			fscanf(fp, "%d", &din);
			my_amg_manager.bStressMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.

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

			// number nFinnest sweeps :
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

			// number postsweeps :
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.nu2_Stress = din;

			// memory size :
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Temperature = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Speed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Pressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.memory_size_Stress = din;


			// Параметр верхней релаксации в сглаживателе.
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.gold_const_Temperature = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.gold_const_Speed = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.gold_const_Pressure = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.gold_const_Stress = fin;

			// использовать ли ilu2 smoother.
			// использовать ли ilu2 smoother.
			fscanf_s(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Temperature = din;
				my_amg_manager.b_gmresTemp = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Temperature = din;
			}
			fscanf_s(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Speed = din;
				my_amg_manager.b_gmresSpeed = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Speed = din;
			}
			fscanf_s(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Pressure = din;
				my_amg_manager.b_gmresPressure = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Pressure = din;
			}
			fscanf_s(fp, "%d", &din);
			if (din == 3) {
				din = 0;
				my_amg_manager.ilu2_smoother_Stress = din;
				my_amg_manager.b_gmresStress = true;
			}
			else {
				my_amg_manager.ilu2_smoother_Stress = din;
			}


			// strength threshold :
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Temperature = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Speed = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Pressure = fin;
			fscanf_s(fp, "%f", &fin);
			my_amg_manager.theta_Stress = fin;

			// magic threshold :
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
			// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
			// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
			fscanf_s(fp, "%d", &din);
			my_amg_manager.icoarseningTemp = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.icoarseningSpeed = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.icoarseningPressure = din;
			fscanf_s(fp, "%d", &din);
			my_amg_manager.icoarseningStress = din;

			// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
			// Если din==1, Stabilization BiCGStab.
			// 8.01.2017 Метод ван дер Ворста BiCGStab 
			// предобусловленный алгебраичесеким многосеточным методом.
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

			fscanf(fp, "%d", &din);
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
			w = new WALL[lw];
			integer i = 0; // счётчик цикла for

			for (i = 0; i < ltdp; i++) {
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
			for (i = 0; i < lmatmax; i++) {
				// свойства материалов:
				// плотность
				fscanf_s(fp, "%f", &fin);
				matlist[i].rho = fin;
				// теплоёмкость при постоянном давлении
				//fscanf(fp, "%f", &fin);
				//matlist[i].cp = fin;
				fscanf_s(fp, "%d", &din);
				matlist[i].n_cp = din;
				matlist[i].arr_cp = NULL;
				matlist[i].temp_cp = NULL;
				matlist[i].arr_cp = new doublereal[matlist[i].n_cp];
				matlist[i].temp_cp = new doublereal[matlist[i].n_cp];
				if (matlist[i].temp_cp == NULL) {
					printf("problem memory allocation for temp_cp\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_cp == NULL) {
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
				matlist[i].arr_lam = NULL;
				matlist[i].temp_lam = NULL;
				matlist[i].arr_lam = new doublereal[matlist[i].n_lam];
				matlist[i].temp_lam = new doublereal[matlist[i].n_lam];
				if (matlist[i].temp_lam == NULL) {
					printf("problem memory allocation for temp_lam\n");
					system("pause");
					exit(1);
				}
				if (matlist[i].arr_lam == NULL) {
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
				// ортотропность теплопроводности :
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_x = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y = fin;
				fscanf_s(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z = fin;
				// 5.08.2017.
				// Коэффициенты для задачи упругости.
				// Модуль Юнга и коэффициент Пуассона.
				doublereal Poissonratio = 0.154;
				doublereal Youngmodule = 217.5e9;
				fscanf_s(fp, "%f", &fin);
				Poissonratio = fin;
				fscanf_s(fp, "%f", &fin);
				Youngmodule = fin * 1e9;
				fscanf_s(fp, "%f", &fin);
				// beta_t_solid*1E-6
				matlist[i].beta_t_solid = fin * 1E-6;
				// Коэффициенты Лямэ.
				doublereal E1_koef = Youngmodule / (1.0 - Poissonratio * Poissonratio);
				doublereal nu1_koef = Poissonratio / (1.0 - Poissonratio);
				matlist[i].mu_Lame = E1_koef / (2.0 * (1.0 + nu1_koef));
				matlist[i].lambda_Lame = (E1_koef * nu1_koef) / (1.0 - nu1_koef * nu1_koef);
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
			for (i = 0; i < lb; i++) {

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
					case XY:
						b[i].g.zC += b[i].g.Hcyl;
						break;
					case XZ:
						b[i].g.yC += b[i].g.Hcyl;
						break;
					case YZ:
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
						case XY:
							b[i].g.zi[i73] += b[i].g.hi[i73];
							break;
						case XZ:
							b[i].g.yi[i73] += b[i].g.hi[i73];
							break;
						case YZ:
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
					case XY:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmin53 + b[i].g.hi[0];
						break;
					case XZ:
						b[i].g.xS = xmin53;
						b[i].g.xE = xmax53;
						b[i].g.zS = zmin53;
						b[i].g.zE = zmax53;
						b[i].g.yS = ymin53;
						b[i].g.yE = ymin53 + b[i].g.hi[0];
						break;
					case YZ:
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
				b[i].radiation.nodelistW = NULL;
				b[i].radiation.nodelistE = NULL;
				b[i].radiation.nodelistS = NULL;
				b[i].radiation.nodelistN = NULL;
				b[i].radiation.nodelistB = NULL;
				b[i].radiation.nodelistT = NULL;
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
				b[i].arr_Sc = NULL;
				b[i].temp_Sc = NULL;
				b[i].arr_Sc = new doublereal[b[i].n_Sc];
				b[i].temp_Sc = new doublereal[b[i].n_Sc];
				if (b[i].temp_Sc == NULL) {
					printf("problem memory allocation for temp_Sc\n");
					system("pause");
					exit(1);
				}
				if (b[i].arr_Sc == NULL) {
					printf("problem memory allocation for arr_Sc\n");
					system("pause");
					exit(1);
				}
				for (integer i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
					// Температура в C.
					fscanf_s(fp, "%f", &fin);
					b[i].temp_Sc[i_4] = fin;
					fscanf_s(fp, "%f", &fin);
					if (fin != fin) {
						b[i].arr_Sc[i_4] = 0.0;
					}
					else {
						b[i].arr_Sc[i_4] = fin;
					}
				}

				// debug
				//if (fabs(b[i].Sc)>1.0e-30) {
					//printf("%e\n", b[i].Sc);
					//system("PAUSE");
				//}
				// стиль зависимости мощности тепловыделения в блоке от времени.
				// 0 - не зависит от времени, 1 - square wave зависимость, 
				// 2 - square wave apparat зависимость.
				fscanf_s(fp, "%d", &din);
				b[i].ipower_time_depend = din;
				// тип блока
				fscanf_s(fp, "%d", &din);
				b[i].itype = din;

				// печать считанных значений на консоль
				printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				printf("%d %d %d\n", b[i].imatid, b[i].itype, b[i].ipower_time_depend);
				printf("temperature depend power\n");
				printf("t_C power_W\n");
				for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
				}
			}

			// считывание источников тепла
			for (i = 0; i < ls; i++) {

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
				case XY: s[i].square = fabs(s[i].g.xE - s[i].g.xS) * fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ: s[i].square = fabs(s[i].g.xE - s[i].g.xS) * fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ: s[i].square = fabs(s[i].g.yE - s[i].g.yS) * fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				printf("%e %d %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
			}

			// считывание твёрдых стенок



			for (i = 0; i < lw; i++) {

				fscanf_s(fp, "%d", &din);
				w[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

				fscanf_s(fp, "%d", &din);
				w[i].ifamily = din;
				switch (din) {
				case 1:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf_s(fp, "%f", &fin); // Stefan Bolcman
					// termostability wall
					w[i].emissivity = 0.0;
					w[i].film_coefficient = 0.0;
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // первого рода
				case 2:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = 0.0;
					fscanf_s(fp, "%f", &fin);  // Stefan Bolcman
					// adiabatic wall
					w[i].emissivity = 0.0;
					w[i].film_coefficient = 0.0;
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // однородное условие Неймана
				case 3:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf_s(fp, "%f", &fin); // Stefan Bolcman
					// Newton-Richman condition, film coefficient.
					w[i].emissivity = 0.0;
					w[i].film_coefficient = fin;
					fscanf_s(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Ньютон-Рихман.
				case 4:  fscanf_s(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf_s(fp, "%f", &fin); // Stefan Bolcman
					// Stefan - Bolcman condition
					w[i].emissivity = fin;
					w[i].film_coefficient = 0.0;
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
				fscanf(fp, "%d", &din);
				w[i].ithermal_Stress_boundary_condition = din;
				fscanf(fp, "%f", &fin);
				w[i].xForce = fin;
				fscanf(fp, "%f", &fin);
				w[i].yForce = fin;
				fscanf(fp, "%f", &fin);
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
			fscanf_s(fp, "%d", &din);
			lu = din;
			if (lu == 0) {
				my_union = NULL;
			}
			else {
				my_union = new UNION[lu];
				// инициализация.
				for (i = 0; i < lu; i++) {
					my_union[i].f = NULL;
					my_union[i].xpos = NULL;
					my_union[i].ypos = NULL;
					my_union[i].zpos = NULL;
					my_union[i].xposadd = NULL;
					my_union[i].yposadd = NULL;
					my_union[i].zposadd = NULL;
					my_union[i].iswitchMeshGenerator = 2; // 2 - CoarseMeshGen
					my_union[i].inxadd = -1;
					my_union[i].inyadd = -1;
					my_union[i].inzadd = -1;
					my_union[i].flow_interior = 0;
				}
			}
			for (i = 0; i < lu; i++) {
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
			//printf("eqin.imaxflD= %d\n", eqin.imaxflD); getchar();
			if (eqin.imaxflD == 0) {
				eqin.fluidinfo = NULL;
			}
			else
			{
				// выделение оперативной памяти
				if (eqin.fluidinfo != NULL) {
					delete eqin.fluidinfo;
					eqin.fluidinfo = NULL;
				}
				eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
				for (i = 0; i < eqin.imaxflD; i++) {
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
					eqin.fluidinfo[i].iflowregime = din;
					fscanf_s(fp, "%d", &din);
					eqin.fluidinfo[i].iturbmodel = din;
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
					pfpir.idir = din;

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					AMG1R6_LABEL = din;

				}
			}

			fclose(fp); // закрытие файла
		}
	}


#else
// eqin - информация о наборе решаемых уравнений.

// dgx, dgy, dgz - вектор силы тяжести.
// inx, iny, inz - количество точек по каждой из осей.

FILE* fp;
errno_t err1;
if ((err1 = fopen_s(&fp, fname, "r")) != 0) {
	printf("No input File premeshin.txt \n");
	//system("PAUSE");
	system("pause");

}
else
{
	if (fp != NULL) {
		float fin = 0.0;
		integer din = 0;
		doublereal scale = 1.0;
		doublereal dbuf; // для упорядочивания в порядке возрастания



		fscanf_s(fp, "%lld", &din);
		ionly_solid_visible = din;
		printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
		fscanf_s(fp, "%f", &fin);
		scale = fin;
		fscanf_s(fp, "%lld", &din);
		lmatmax = din;
		fscanf_s(fp, "%lld", &din);
		lb = din;
		fscanf_s(fp, "%lld", &din);
		ls = din;
		fscanf_s(fp, "%lld", &din);
		lw = din;
		fscanf_s(fp, "%lld", &din);
		ltdp = din; // количество уникальных данных с табличными данными по зависимости расеиваемой мощности от температуры.
		

					// Считываем значение вектора силы тяжести:
		fscanf_s(fp, "%f", &fin);
		dgx = fin;
		fscanf_s(fp, "%f", &fin);
		dgy = fin;
		fscanf_s(fp, "%f", &fin);
		dgz = fin;

		// считываем количество точек на каждой координатной оси
		fscanf_s(fp, "%lld", &din);
		inx = din;
		fscanf_s(fp, "%lld", &din);
		iny = din;
		fscanf_s(fp, "%lld", &din);
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
		Tochka_position_X0_for_XY_Plot = scale*fin;
		fscanf_s(fp, "%f", &fin);
		Tochka_position_Y0_for_XY_Plot = scale*fin;
		fscanf_s(fp, "%f", &fin);
		Tochka_position_Z0_for_XY_Plot = scale*fin;
		// Направление линии совпадает с направлением 
		// одной из осей декартовой прямоугольной системы координат:
		// 0 - Ox, 1 - Oy, 2 - Oz.
		fscanf_s(fp, "%lld", &din);
		idirectional_for_XY_Plot = din;

		fscanf_s(fp, "%f", &fin);
		etalon_max_size_ratio = fin; // подробность расчётной сетки.

		fscanf_s(fp, "%f", &fin);
		etalon_max_size_ratio2 = fin; // Критерий качества расчётной сетки на основе FlowVision.

		fscanf_s(fp, "%lld", &din);
		switch (din) {
		case 0: bsnap_TO_global = false;  break;
		default: bsnap_TO_global = true;  break;
		}

		fscanf_s(fp, "%lld", &din);
		iswitchsolveramg_vs_BiCGstab_plus_ILU2 = din; // Выбор решающего устройства : либо amg1r5 либо BiCGStab+ILU2.

		fscanf_s(fp, "%lld", &din);
		iswitchsolveramg_vs_BiCGstab_plus_ILU6 = din; // Выбор решающего устройства : либо РУМБА0.14 либо BiCGStab+ILU6.

		fscanf_s(fp, "%lld", &din);
		if (din == 1) {
			// SIMPLEC algorithm.
			iSIMPLE_alg = SIMPLEC_Van_Doormal_and_Raithby;
		}
		else {
			// SIMPLE algorithm 1972.
			iSIMPLE_alg = SIMPLE_Carretto;
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
		default: iFLOWScheme = 2; break; // UDS самая стабильная схема.
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
		default: iTEMPScheme = 2; break; // UDS самая стабильная схема.
		}



		// Выбор сеточного генератора.
		fscanf_s(fp, "%lld", &din);
		iswitchMeshGenerator = din;

		fscanf_s(fp, "%lld", &din);
		steady_or_unsteady_global_determinant = 2;
		if ((din == 0) || (din == 1) || (din == 2) || (din == 3) || (din == 5) || (din == 6) || (din == 7) || (din == 8) || (din == 9)) {
			// 0 - thermal only steady state calculation,
			// 1 - thermal only unsteady calculation,
			// 2 - mesh generator only.
			// 3 - fluid dynamic steady state.
			// 5 - Static Structural (Thermal solver #2)
			// 6 - Thermal Stress
			// 7 - Unsteady thermal solver #2
			// 8 - Visualisation only
			// 9 - cfd unsteady fluid dynamic.
			steady_or_unsteady_global_determinant = din; // thermal only: steady  - 0, or unsteady -1 calculation.
			//printf("steady_or_unsteady_global_determinant =%lld\n",din);
			//system("PAUSE");
		}
		else {
			printf("error input parametr steady or unsteady calculation\n");
			system("PAUSE");
			exit(1);
		}

		fscanf_s(fp, "%lld", &din);
		if ((din == 0) || (din == 1) || (din == 2) || (din == 3)) {
			glTSL.id_law = din;
		}
		else {
			printf("error input parametr timestep law\n");
			system("PAUSE");
			exit(1);
		}
		fscanf_s(fp, "%f", &fin);
		if ((fin <= 0.0) || (fin >= 1.0)) {
			printf("error input parametr timestep law Factor a\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.Factor_a_for_Linear = fin; // Factor_a
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law tau must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau = fin; // длительность импульса.
		fscanf_s(fp, "%f", &fin);
		glTSL.Q = fin;  // Скважность.
						// Параметры импульсного режима для темы АППАРАТ.
		fscanf_s(fp, "%f", &fin);
		if ((fin <= 0.0) || (fin >= 1.0)) {
			printf("error input parametr timestep law APPARAT multiplyer\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.m1 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law APPARAT tau1 must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau1 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law APPARAT tau2 must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau2 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law APPARAT tau_pause must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau_pause = fin;
		fscanf_s(fp, "%lld", &din);
		glTSL.n_cycle = din;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law APPARAT Period must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.T_all = fin;
		doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
		if (t_pause_gl <= 0.0) {
			printf("error in parameters Square Wave APPARAT time step law.\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}

		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
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
		adiabatic_vs_heat_transfer_coeff = din; // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
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
		itype_ALICE_Mesh = din;

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.m_restart = din;
		// classical algebraic multigrid parameters:
		// only for my_agregat_amg.cu.
		// only for my_agregat_amg.cu.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.imySortAlgorithm = din;
		fscanf_s(fp, "%lld", &din);
		//my_amg_manager.maximum_levels = din;
		my_amg_manager.maximum_delete_levels_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.maximum_delete_levels_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.maximum_delete_levels_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.maximum_delete_levels_Stress = din;

		// type interpolation procedure :
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
		my_amg_manager.iCFalgorithm_and_data_structure_Temperature = din;// 3-Treap.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iCFalgorithm_and_data_structure_Speed = din;// 3-Treap.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iCFalgorithm_and_data_structure_Pressure = din;// 3-Treap.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.iCFalgorithm_and_data_structure_Stress = din;// 3-Treap.

		fscanf_s(fp, "%lld", &din);
		//my_amg_manager.itypemodifyinterpol = din;
		//my_amg_manager.baglomeration_with_consistency_scaling = din;
		my_amg_manager.bdiagonal_dominant = din;
		fscanf_s(fp, "%lld", &din);
		//my_amg_manager.inumberadaptpass = din;


		// 23.02.2018
		// print matrix portrait
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.bTemperatureMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.bSpeedMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.bPressureMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.bStressMatrixPortrait = din; // 0 - NO_PRINT, 1 - PRINT.

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

		// number nFinnest sweeps :
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

		// number postsweeps :
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.nu2_Stress = din;

		// memory size :
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Temperature = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Speed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Pressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.memory_size_Stress = din;


		// Параметр верхней релаксации в сглаживателе.
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.gold_const_Temperature = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.gold_const_Speed = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.gold_const_Pressure = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.gold_const_Stress = fin;

		// использовать ли ilu2 smoother.
		// использовать ли ilu2 smoother.
		fscanf_s(fp, "%lld", &din);
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Temperature = din;
			my_amg_manager.b_gmresTemp = true;
		}
		else {
			my_amg_manager.ilu2_smoother_Temperature = din;
		}
		fscanf_s(fp, "%lld", &din);
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Speed = din;
			my_amg_manager.b_gmresSpeed = true;
		}
		else {
			my_amg_manager.ilu2_smoother_Speed = din;
		}
		fscanf_s(fp, "%lld", &din);
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Pressure = din;
			my_amg_manager.b_gmresPressure = true;
		}
		else {
			my_amg_manager.ilu2_smoother_Pressure = din;
		}
		fscanf_s(fp, "%lld", &din);
		if (din == 3) {
			din = 0;
			my_amg_manager.ilu2_smoother_Stress = din;
			my_amg_manager.b_gmresStress = true;
		}
		else {
			my_amg_manager.ilu2_smoother_Stress = din;
		}

		// strength threshold :
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Temperature = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Speed = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Pressure = fin;
		fscanf_s(fp, "%f", &fin);
		my_amg_manager.theta_Stress = fin;

		// magic threshold :
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
		// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
		// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.icoarseningTemp = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.icoarseningSpeed = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.icoarseningPressure = din;
		fscanf_s(fp, "%lld", &din);
		my_amg_manager.icoarseningStress = din;

		// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
		// Если din==1, Stabilization BiCGStab.
		// 8.01.2017 Метод ван дер Ворста BiCGStab 
		// предобусловленный алгебраичесеким многосеточным методом.
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
		w = new WALL[lw];
		integer i = 0; // счётчик цикла for

		for (i = 0; i < ltdp; i++) {
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
		for (i = 0; i < lmatmax; i++) {
			// свойства материалов:
			// плотность
			fscanf_s(fp, "%f", &fin);
			matlist[i].rho = fin;
			// теплоёмкость при постоянном давлении
			//fscanf(fp, "%f", &fin);
			//matlist[i].cp = fin;
			fscanf_s(fp, "%lld", &din);
			matlist[i].n_cp = din;
			matlist[i].arr_cp = NULL;
			matlist[i].temp_cp = NULL;
			matlist[i].arr_cp = new doublereal[matlist[i].n_cp];
			matlist[i].temp_cp = new doublereal[matlist[i].n_cp];
			if (matlist[i].temp_cp == NULL) {
				printf("problem memory allocation for temp_cp\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_cp == NULL) {
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
			matlist[i].arr_lam = NULL;
			matlist[i].temp_lam = NULL;
			matlist[i].arr_lam = new doublereal[matlist[i].n_lam];
			matlist[i].temp_lam = new doublereal[matlist[i].n_lam];
			if (matlist[i].temp_lam == NULL) {
				printf("problem memory allocation for temp_lam\n");
				system("pause");
				exit(1);
			}
			if (matlist[i].arr_lam == NULL) {
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
			// ортотропность теплопроводности :
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_x = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_y = fin;
			fscanf_s(fp, "%f", &fin);
			matlist[i].orthotropy_multiplyer_z = fin;
			// 5.08.2017.
			// Коэффициенты для задачи упругости.
			// Модуль Юнга и коэффициент Пуассона.
			// Алюминий стр. 232 В.Н.Сидоров, В.В. Вершинин Метод конечных элементов в расчёте сооружений.
			//doublereal Poissonratio = 0.33;
			//doublereal Youngmodule = 71.7e9;
			// steel
			// Tensile Yield Strength 4.12E+8
			// Compressive Yield Strength 4.12E+8
			// Tensile Ultimate Strength 5.7E+8
			//doublereal Poissonratio = 0.154;
			//doublereal Youngmodule = 217.5E9;
			doublereal Poissonratio = 0.3;
			doublereal Youngmodule = 200.0E9;
			fscanf_s(fp, "%f", &fin);
			Poissonratio = fin;
			fscanf_s(fp, "%f", &fin);
			Youngmodule = fin*1e9;
			fscanf_s(fp, "%f", &fin);
			// beta_t_solid*1E-6
			matlist[i].beta_t_solid = fin*1E-6;
			// Коэффициенты Лямэ.
			//doublereal E1_koef = Youngmodule / (1.0 - Poissonratio*Poissonratio);
			//doublereal nu1_koef = Poissonratio / (1.0 - Poissonratio);
			//matlist[i].mu_Lame = E1_koef / (2.0*(1.0 + nu1_koef));
			//matlist[i].lambda_Lame = (E1_koef*nu1_koef) / (1.0 - nu1_koef*nu1_koef);
			// стр. 25 В.Н.Сидоров, В.В. Вершинин Метод конечных элементов в расчёте сооружений.
			//+ 19.10.2018 проверено.
			matlist[i].mu_Lame = Youngmodule / (2.0*(1.0+ Poissonratio));
		    matlist[i].lambda_Lame = (Poissonratio*Youngmodule) / ((1.0+ Poissonratio)*(1.0-2.0*Poissonratio));
			//printf("E=%e N/m^2 mu=%e lambda=%e\n", Youngmodule, matlist[i].mu_Lame, matlist[i].lambda_Lame);
			//system("PAUSE");
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
		for (i = 0; i<lb; i++) {

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
				case XY:
					b[i].g.zC += b[i].g.Hcyl;
					break;
				case XZ:
					b[i].g.yC += b[i].g.Hcyl;
					break;
				case YZ:
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
					case XY:
						b[i].g.zi[i73] += b[i].g.hi[i73];
						break;
					case XZ:
						b[i].g.yi[i73] += b[i].g.hi[i73];
						break;
					case YZ:
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
				case XY:
					b[i].g.xS = xmin53;
					b[i].g.xE = xmax53;
					b[i].g.yS = ymin53;
					b[i].g.yE = ymax53;
					b[i].g.zS = zmin53;
					b[i].g.zE = zmin53 + b[i].g.hi[0];
					break;
				case XZ:
					b[i].g.xS = xmin53;
					b[i].g.xE = xmax53;
					b[i].g.zS = zmin53;
					b[i].g.zE = zmax53;
					b[i].g.yS = ymin53;
					b[i].g.yE = ymin53 + b[i].g.hi[0];
					break;
				case YZ:
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
			b[i].radiation.nodelistW = NULL;
			b[i].radiation.nodelistE = NULL;
			b[i].radiation.nodelistS = NULL;
			b[i].radiation.nodelistN = NULL;
			b[i].radiation.nodelistB = NULL;
			b[i].radiation.nodelistT = NULL;
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
			b[i].arr_Sc = NULL;
			b[i].temp_Sc = NULL;
			b[i].arr_Sc = new doublereal[b[i].n_Sc];
			b[i].temp_Sc = new doublereal[b[i].n_Sc];
			if (b[i].temp_Sc == NULL) {
				printf("problem memory allocation for temp_Sc\n");
				system("pause");
				exit(1);
			}
			if (b[i].arr_Sc == NULL) {
				printf("problem memory allocation for arr_Sc\n");
				system("pause");
				exit(1);
			}
			for (integer i_4 = 0; i_4 < b[i].n_Sc; i_4++) {
				// Температура в C.
				fscanf_s(fp, "%f", &fin);
				b[i].temp_Sc[i_4] = fin;
				fscanf_s(fp, "%f", &fin);
				if (fin != fin) {
					b[i].arr_Sc[i_4] = 0.0;
				}
				else {
					b[i].arr_Sc[i_4] = fin;
				}
			}

			// debug
			//if (fabs(b[i].Sc)>1.0e-30) {
			//printf("%e\n", b[i].Sc);
			//system("PAUSE");
			//}
			// стиль зависимости мощности тепловыделения в блоке от времени.
			// 0 - не зависит от времени, 1 - square wave зависимость, 
			// 2 - square wave apparat зависимость.
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			b[i].ipower_time_depend = din;
			// тип блока
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			b[i].itype = din;

			// печать считанных значений на консоль
			//printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
			//printf("%lld %lld %lld\n", b[i].imatid, b[i].itype, b[i].ipower_time_depend);
			//printf("temperature depend power\n");
			//printf("t_C power_W\n");
			//for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
				//printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
			//}
		}

		// считывание источников тепла
		for (i = 0; i < ls; i++) {

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
			case XY: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
			case XZ: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
			case YZ: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
			default: break;
			}
			//printf("source %e %lld %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
		}

		// считывание твёрдых стенок



		for (i = 0; i < lw; i++) {

			fscanf_s(fp, "%lld", &din);
			w[i].iunion_id = din; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.

#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			w[i].ifamily = din;
			switch (din) {
			case 1:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = fin;
				fscanf_s(fp, "%f", &fin); // Stefan Bolcman
										  // termostability wall
				w[i].emissivity = 0.0;
				w[i].film_coefficient = 0.0;
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // первого рода
			case 2:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = 0.0;
				fscanf_s(fp, "%f", &fin);  // Stefan Bolcman
										   // adiabatic wall
				w[i].emissivity = 0.0;
				w[i].film_coefficient = 0.0;
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // однородное условие Неймана
			case 3:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = fin;
				fscanf_s(fp, "%f", &fin); // Stefan Bolcman
										  // Newton-Richman condition, film coefficient.
				w[i].emissivity = 0.0;
				w[i].film_coefficient = fin;
				fscanf_s(fp, "%f", &fin);
				w[i].hf = 0.0;
				break; // Ньютон-Рихман.
			case 4:  fscanf_s(fp, "%f", &fin);
				w[i].Tamb = fin;
				fscanf_s(fp, "%f", &fin); // Stefan Bolcman
										  // Stefan - Bolcman condition
				w[i].emissivity = fin;
				w[i].film_coefficient = 0.0;
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
				w[i].ithermal_Stress_boundary_condition = din;
			}
			else {
				printf("error: unknown ithermal_Stress_boundary_condition\n");
				printf("ithermal_Stress_boundary_condition=%lld\n", din);
				system("PAUSE");
				w[i].ithermal_Stress_boundary_condition = 0; // Free all
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
#if doubleintprecision == 1
		fscanf_s(fp, "%lld", &din);
#else
		fscanf_s(fp, "%d", &din);
#endif
		lu = din;
		if (lu == 0) {
			my_union = NULL;
		}
		else {
			my_union = new UNION[lu];
			// инициализация.
			for (i = 0; i < lu; i++) {
				my_union[i].f = NULL;
				my_union[i].xpos = NULL;
				my_union[i].ypos = NULL;
				my_union[i].zpos = NULL;
				my_union[i].xposadd = NULL;
				my_union[i].yposadd = NULL;
				my_union[i].zposadd = NULL;
				my_union[i].iswitchMeshGenerator = 2; // 2 - CoarseMeshGen
				my_union[i].inxadd = -1;
				my_union[i].inyadd = -1;
				my_union[i].inzadd = -1;
				my_union[i].flow_interior = 0;
			}
		}
		for (i = 0; i < lu; i++) {
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
			my_union[i].inx = din;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].iny = din;
#if doubleintprecision == 1
			fscanf_s(fp, "%lld", &din);
#else
			fscanf_s(fp, "%d", &din);
#endif
			my_union[i].inz = din;

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
			eqin.fluidinfo = NULL;
		}
		else
		{
			// выделение оперативной памяти
			if (eqin.fluidinfo != NULL) {
				delete eqin.fluidinfo;
				eqin.fluidinfo = NULL;
			}
			eqin.fluidinfo = new FLOWINFO[eqin.imaxflD];
			for (i = 0; i < eqin.imaxflD; i++) {
				// Считывание координат опорной точки
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].xc = scale*fin;
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].yc = scale*fin;
				fscanf_s(fp, "%f", &fin);
				eqin.fluidinfo[i].zc = scale*fin;
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
				eqin.fluidinfo[i].iflowregime = din;
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				eqin.fluidinfo[i].iturbmodel = din;
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
				eqin.fluidinfo[i].roughness = 1.0e-6*fin; // шероховатость стенки в м.
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
				pfpir.fminimum = (doublereal)(scale * fin);
				fscanf_s(fp, "%f", &fin);
				pfpir.fmaximum = (doublereal)(scale * fin);
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				pfpir.idir = din;

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				AMG1R6_LABEL = din;

			}
		}


		integer ilb_p = 0;// Количество блоков внутри которых задана тепловая мощность.
		doublereal dpower = 0.0; // Суммарная тепловая мощность в блоках.
		integer ipoly = 0, icyl = 0, iprism = 0;
		integer ihol = 0, isol = 0, iflui = 0;
		for (integer i_1 = 0; i_1 < lb; i_1++) {
			if (b[i_1].itype == HOLLOW) {
				ihol++;
			}
			if (b[i_1].itype == FLUID) {
				iflui++;
			}
			if (b[i_1].itype == SOLID) {
				isol++;
			}
			if (b[i_1].g.itypegeom == PRISM) {
				// 0 - PRISM object
				iprism++;
				if (b[i_1].n_Sc > 0) {
					doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
					doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS)*fabs(b[i_1].g.yE - b[i_1].g.yS)*fabs(b[i_1].g.zE - b[i_1].g.zS);
					if (vol < 1.0e-40) {
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
				icyl++;
				if (b[i_1].n_Sc > 0) {
					doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
					doublereal vol = 0.0;
					const doublereal MPI0 = 3.1415926;
					vol = b[i_1].g.Hcyl*MPI0*(b[i_1].g.R_out_cyl*b[i_1].g.R_out_cyl - b[i_1].g.R_in_cyl*b[i_1].g.R_in_cyl);				
					if (vol < 1.0e-40) {
						printf("ERROR: zero volume in CYLINDER block number %lld\n", i_1);
						system("PAUSE");
						exit(1);
					}
					if (pdiss > 0.0) {
						ilb_p++;

						dpower += pdiss*vol;
						//printf("ERROR : non zero power in cylinder object.\n");
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
					if (pdiss > 0.0) {
						ilb_p++;
						printf("ERROR : non zero power in polygon object.\n");
						system("PAUSE");
						exit(1);
					}
				}
			}
		}

		doublereal dsoupow = 0.0; // интегральная тепловая мощность плоских источников тепла.
		for (integer i_1 = 0; i_1 < ls; i_1++) {
			dsoupow += s[i_1].power;
		}


		printf("Apriory quick model statistics:\n");
		printf("number of thermal power blocks lb_p=%lld\n", ilb_p);
		printf("Blocks integral power =%e W\n", dpower);
		printf("number of sources ls=%lld\n",ls);
		printf("Sources integral power = %e W\n", dsoupow);
		printf("Full total power = %e W\n", dpower + dsoupow);
		printf("number of blocks lb=%lld\n",lb);
		printf("PRISMS = %lld, CYLINDERS = %lld, POLYGONS = %lld\n", iprism, icyl, ipoly);
		printf("SOLID: %lld\n", isol);
		printf("HOLLOW: %lld\n", ihol);
		printf("FLUID: %lld\n", iflui);
		printf("number of walls lw=%lld\n",lw);
		printf("number of units lu=%lld\n", lu);

		fclose(fp); // закрытие файла
	}
}

#endif
#endif
	printf("OK. \n");
} // premeshin

#endif
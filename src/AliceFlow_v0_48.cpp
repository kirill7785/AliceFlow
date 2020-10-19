// AliceFlow_v0_48.cpp
// 24.06.2020 - 16 июля 2020 графовый температурный солвер - network_T solver.
// 2 июня 2020 нелинейный температурный солвер на базе amg1r5.
// 27 февраля 2020 PMIS aggregator, FF1 amg interpolation. Acomp=3.0;
// Как в магистерской диссертации Джефри Баттлера 2006. 
// Авторы агрегатора PMIS Янг и де Стерк 2004год.
// 12 января 2020 iluk сглаживатель для amg1r5 алгоритма.
// 10 октября 2019 запрограммировал модель турбулентности K-Omega SST Ментера (RANS). 
// 2 октября 2019 запрограммировал модель турбулентности Спаларта Аллмареса (RANS).
// 4 августа 2019; 03 ноября 2019 откомпилировал компилятором gcc (g++ 9.1). GNU project. 
// 7-8 мая 2019 присоединил алгебраический многосеточный метод amgcl Дениса Демидова.
// 6 апреля 2019 откомпилирована в visual studio community edition 2019 (open source).
// 25.03.2019 Начал использовать PVS-Studio 6.0
// 19 марта(03) 2019 заработала гидродинамика на АЛИС сетках.
// 6.05.2018 LINK: fatal error LNK1102: недостаточно памяти 2015 VS community.
// Выход: компиляция проекта с опцией /bigobj
// Подсветка синтаксиса для cuda:
// Tools->Options->Text Editor->File Extension Ввести cu и нажать кнопку add.
// 9 июля 2017 переход на 64 битные целые int64_t.
// 15 апреля 2017 откомпилирована в vs community edition 2017 (open source).
// 1 октября 2016 откомпилировал на nvidia cuda 8.0. 
// АЛИС сетки введены в строй для теплопередачи в твёрдом теле. 
// 11 января 2016 года добавил cl_agl_amg_v0_14.
// 15 августа 2015 года. Теперь в Visual Studio 2015.
// AliceFlow_v0_21.cpp
// 15 августа 2015. Теперь в Visual Studio 2013.
// AliceFlow_v0_20.cpp
// 14 августа 2015 Действительно правильное распараллеливание lusol и
// ilu2 decomposition на 2 потока.
// 
// AliceFlow_v0_07.cpp на основе  AliceFlow_v0_06.cpp, но теперь с LES моделью турбулентности.
// Программа не прошла тестирование и содержит ошибки в модели турбулентности Germano.
// 17 апреля 2013 года. Правильное распараллеливание lusol_.
// 1 апреля 2013. Теперь в Visual Studio 2012.
//
// AliceFlow_v0_06.cpp:
// 3D программа AliceFlow_v0_06.cpp наследует свойства AliceFlowv0_05.cpp
// развивая их дальше, наращивая стабильность и функциональность.
// 
// Программа AliceFlow_v0_05.cpp, 
// улучшенный вариант AliceFlow_v0_03.cpp, преобразует 
// файл из формата визуального проектирования
// в формат поступающий на вход солверу.
// В этой же программе генерируется 
// неравномерная структурированная HEX 
// расчётная сетка.
// begin one 17 мая 2011 года.
//
// 3D программа AliceFlowv0_05.cpp 
// осуществляет следующие операции:
// 1. Генератор конечнообъёмной сетки.
// 2. Создание и заполнение всех структур
//    данных необходимых для сборки матрицы СЛАУ.
// 3. Сборка матрицы СЛАУ.
// 4. Решатель СЛАУ.
// 5. Экспорт в визуализатор tecplot360.
// begin two 30 июня 2011 года.
// begin three 14 октября 2011 года. Теперь на Visual Studio 2010.
// begin four 12 марта 2012 года. (по мотивам книги Р. Лафоре - переход на ООП).
//
// Реальные размеры источника тепла в транзисторе типа TGF2023_*
// равны 0.2*120мкм толщиной 100 ангстрем (10-20нм).
 

// закомментировать в случае если сборка приложения 
//осуществляется компилятором gcc (g++ 9.1) от GNU.
//#include "stdafx.h"
//#include "pch.h"

//#define VISUAL_TUDIO_2008_COMPILLER

#ifdef VISUAL_TUDIO_2008_COMPILLER
#define nullptr NULL // для Visual Studio 2008
#endif

// Раскомментировать в случае если сборка приложения 
//осуществляется компилятором g++ (g++ 9.1) от GNU.
// g++ AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt
// gcc_log.txt - файл с диагностическими сообщениями компилятора.
#define MINGW_COMPILLER 1

#ifdef MINGW_COMPILLER

// Типа float 32 бита не хватает - на больших моделях наблюдается расходимость вычислительного процесса.
// Тип float (если сходимость вычислительного процесса имеет место быть) дает в двое большее число итераций
// по сравнению с типом double (12 против 7).
// Тип float128 128 бит считает в 32 раза медленнее чем тип double 64 бита. Он дает сокращение итераций 
// 31 против 38 у типа double. Модели не сходящиеся в double не сходятся и в типе float128. Причина необъяснимая 
// дефектность модели подаваемой на вход солверу.
// По умолчанию рекомендуется использовать только тип double.


/*
// Источник информации
// boost.org/doc/libs/1_65_1/libs/math/example/float128_example.cpp

// Четверная точность. Тип __float128
//#include <quadmath.h>
#include <boost/cstdfloat.hpp> // For float_64_t, float128_t. Must be first include!
//#include <boost/config.hpp>
#include <boost/multiprecision/float128.hpp>
//#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions.hpp> // For gamma function.
#include <boost/math/constants/calculate_constants.hpp> // For constants pi, e ...
#include <typeinfo> //

#include <cmath> // for pow function
*/
/*
C:\AliceFlow_v0_48_GCC>g++ -std=c++11 -fexceptions -std=gnu++17 -g -fext-numeric-literals -II:\modular-boost\libs\math\include -Ii:\modular-boost -c AliceFlow_v0_48.cpp -o obj\Debug\AliceFlow_v0_48.o
C:\AliceFlow_v0_48_GCC>g++ -o bin\Debug\AliceFlow_v0_48.exe obj\Debug\AliceFlow_v0_48.o -lquadmath
*/

//#include <iostream>

//using namespace boost::multiprecision;

//#include <stdlib.h>
//#include <stdio.h>

#endif



// /fp:except option
//#pragma float_control( except, on )

// для std::locale::global(std::locale("en_US.UTF-8"));
// Не работает.
//#include <locale.h>

//using namespace std;

// Проект собран под язык cuda c.

//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"


#include <iostream>

#include <stdio.h>
#include <string.h>


#include <stdlib.h> 
#include <omp.h> // OpenMP
//-Xcompiler "/openmp" 

//#define NOMINMAX // раскомментировать
//#include <stdio.h>
//#include <stdlib.h>
//#include <windows.h> // раскомментировать

// Раскомментировать в случае если ошибки NOT NUMBER отслеживаются.
//#define MY_DEBUG_NOT_NUMBER

unsigned int calculation_main_start_time_global_Depend = 0; // начало счёта мс.
bool CAD_GEOMETRY_OCTREE_MESHGEN = false;
bool SPARSE_MESHER = false;// true;// дополнительное разрежение расчётной сетки в случае амплитрона.

// Одна клетка по Oz
bool b_one_cell_z = false;

// Свободный параметр используемый для отладки 
// программы AliceFlov_v.0.48 и передаваемый из
// программы интерфейса
// AliceMesh_v.0.45.
double free_debug_parametr1 = 0.0;

// Число потоков исполнения. 
// Передаётся из интерфейса пользователя.
// Шина памяти насыщается после 4 потоков.
int number_processors_global_var = 1;

int number_cores() {
	return number_processors_global_var;

	// 4 xeon 2630 v4.
	// Шина памяти насыщается после 4 потоков.
} 

//using namespace System;

// 0 - Проект собран как стандартное консольное приложение windows.
#define _CUDA_IS_ACTIVE_ 0

typedef float real_mix_precision;
// Неизвестных врядли будет больше чем уместится в типе unsigned int,
// однако для прохода по иерархии матриц по прежнему используется тип 
// int64_t.
typedef unsigned int integer_mix_precision;

// Вещественная арифметика.
#define doubleprecision 1
#if doubleprecision == 1
#ifdef MINGW_COMPILLER
// Четверная точность.
//#define doublereal __float128
//#define doublereal  _Quad
//#define doublereal float128
//#define doublereal float // 32 бит
//#define doublereal double // 64 бит
typedef double doublereal;
#else
//#define doublereal  _Quad
//#define doublereal float128
//#define doublereal double // 64 бит
typedef double doublereal;
//#define doublereal long double //double // модель вещественного числа двойной точности
//#define doublereal Decimal // decimal
#endif
#else
//#define doublereal float //float // модель вещественного числа одинарной точности
typedef float doublereal;
#endif

 

#ifdef VISUAL_TUDIO_2008_COMPILLER

doublereal fmin(doublereal a, doublereal b) {
	doublereal r=a;
	if (b<r) r=b;
	return r;
}

doublereal fmax(doublereal a, doublereal b) {
	doublereal r=a;
	if (b>r) r=b;
	return r;
}

template <typename VType>
bool isfinite(VType x) {
	//return std::isfinite(x);

	if (x!=x) {
		return false;
	}
	else {
        return true;
	}
}

#include <math.h>
doublereal expm1(doublereal x) {
	return (exp(x) - 1.0);
}

#endif

#define doubleintprecision 1


#if doubleintprecision == 1
#include <cinttypes> // для типа 64 битного целого числа int64_t

// Внимание !!! при типе int64_t не работают все солверы библиотеки ViennaCL.
// Библиотека Дениса Демидова AMGCL работает и с типом int64_t. 
//#define integer int64_t
typedef int64_t integer;
const int64_t big_FIBO_integer_Value = 9223372036854775803; // для int64t
#else
//#define integer int
typedef int integer;
const integer big_FIBO_integer_Value = 4294967295; // для int
#endif


//const integer parabola_MNK = 0; // производная от параболы построенная по способу наименьших квадратов.
	// Внимание: line_MNK лучше не использовать, т.к. экспериментально установлено, на примере задачи
	// статики tgf2023_01 что сходимости при таком задании производной не наступает, наблюдаются некие автоколебания на протяжении
	// более чем 40 глобальных итераций и более. Возможна даже расходимость.
//const integer line_MNK = 1; // тангенс угла наклона прямой, прямая построена по способу МНК по четырём точкам.
// данный метод также не рекомендуется использовать ! ввиду отсутствия сходимости.
//const integer cubic_parabola = 2; // кубическая парабола по 4 точкам и от неё производная.

//const unsigned int ADIABATIC_WALL_BC = 0; // Адиабатическая стенка.
//const unsigned int NEWTON_RICHMAN_BC = 1; // Нелинейное условие Ньютона-Рихмана.
//const unsigned int STEFAN_BOLCMAN_BC = 2; // Нелинейное граничное условие Стефана-Больцмана.
//const unsigned int MIX_CONDITION_BC = 3; // Смешанное граничное условие Ньютон Рихман + Стефан Больцман.

enum class DEFAULT_CABINET_BOUNDARY_CONDITION { ADIABATIC_WALL_BC = 0, 
	NEWTON_RICHMAN_BC = 1, STEFAN_BOLCMAN_BC = 2, MIX_CONDITION_BC = 3};

//const unsigned int DIRICHLET_FAMILY=1;
//const unsigned int NEIMAN_FAMILY=2;
//const unsigned int NEWTON_RICHMAN_FAMILY=3;
//const unsigned int STEFAN_BOLCMAN_FAMILY=4;

enum class WALL_BOUNDARY_CONDITION { DIRICHLET_FAMILY=1,
	NEIMAN_FAMILY=2, NEWTON_RICHMAN_FAMILY=3, STEFAN_BOLCMAN_FAMILY=4};

enum ORDER_INTERPOLATION { parabola_MNK=0, line_MNK=1, cubic_parabola=2};
enum class ORDER_DERIVATIVE { FIRST_ORDER=1, SECOND_ORDER=2};

enum class INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 {ZERO_INIT = 0, RANDOM_INIT = 1};

// Для функции solve_Thermal(...);
const integer bARRAYrealesation = 1; // На основе отсортированного массива, но вставка в отсортированный массив вызывает проблемы быстродействия.
const integer bAVLrealesation = 2; // На основе АВЛ дерева.

const integer iGLOBAL_RESTART_LIMIT = 2;// 28.05.2020 достаточно двух раз(тестирование данной настройки.). // было 6;
bool bglobal_restart_06_10_2018 = false;
// При стагнации мы всё равно продолжаем но мы признаём что могут быть проблемы при визуализации.
bool bglobal_restart_06_10_2018_stagnation[iGLOBAL_RESTART_LIMIT + 1] = { false,false,false };//{false,false,false, false,false,false, false}; // было 6;
integer iPnodes_count_shadow_memo = 0;

// Запоминаем полное тепловыделение в твердотельной модели, для проверок во время исполнения. 
doublereal d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = 0.0; // Вт
doublereal d_my_optimetric1_6_12_2019 = 0.0;// Глобальная переменная для оптимизации.
doublereal d_my_optimetric2_6_12_2019 = 0.0;// Глобальная переменная для оптимизации.
doublereal d_my_optimetric3_6_12_2019 = 0.0;// Глобальная переменная для оптимизации.

enum class LINE_DIRECTIONAL { X_LINE_DIRECTIONAL=0, 
	Y_LINE_DIRECTIONAL=1, Z_LINE_DIRECTIONAL=2};

// Параметры преобразователя xyplot графиков для отчетов.
// 5.01.2018
typedef struct Tpatcher_for_print_in_report {
	doublereal fminimum, fmaximum;
	LINE_DIRECTIONAL idir; // 0 - X, 1 - Y, 2 - Z.
	Tpatcher_for_print_in_report() {
		fminimum = -1.0e+30;
		fmaximum = 1.0e+30;
		idir= LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; // 0 - X, 1 - Y, 2 - Z.
	}
} Patcher_for_print_in_report;

Patcher_for_print_in_report pfpir; // xyplot графики.

// 9 september 2017.
// делать ли освобождения оперативной памяти и новые построения структур данных.
// Полигоны вызывают проблемы при перестроении сетки, т.к. сетка каждый раз строится по новому, а поле температур 
// остаётся старым. Другое решение данной проблемы состоит из перевыделения поля температур под новую сетку с 
// последующей переинтерполяцией значений температур на новую сетку.
integer ireconstruction_free_construct_alloc = 1; // 0 - off, 1 - on.
// Записывать ли анимацию в текстовый файл
// по окончанию каждого нового шага по времени.
integer ianimation_write_on = 0; // 0 - off, 1 - on.

// Для достижения сходимости мы первые 300 итераций 
// считаем со скоростью vel*rGradual_changes а потом 
// делаем перемасштабирование
// и досчитываем уже со скоростью vel.
// 1.0 - не используется. Не устанавливать 0.1!!!
doublereal rGradual_changes = 1.0; 

// При считывании скоростей из файла load.txt они домножаются
// на множитель my_multiplyer_velocity_load;
doublereal my_multiplyer_velocity_load = 1.0;// 33.3333; // 14.09.2020

// 1 - включить правки amg1r6.f версии; 
// 0 - оставить amg1r5.f в силе.
integer AMG1R6_LABEL = 0; 
integer nrd_LABEL = 1131;
integer nru_LABEL = 1131;
doublereal ecg2_LABEL = 0.25; // strong threshold amg1r5
doublereal ewt2_LABEL = 0.35; // F to F threshold amg1r5
bool b_iluk_amg1r5_LABEL_D = false; // вниз
bool b_iluk_amg1r5_LABEL_U = false; // вверх
// stabilization_amg1r5_algorithm:
// 0 - none(amg1r5); 1 - BiCGStab + amg1r5; 2 - FGMRes+amg1r5;
enum class AMG1R5_OUT_ITERATOR {NONE_only_amg1r5,  BiCGStab_plus_amg1r5,  FGMRes_plus_amg1r5, Non_Linear_amg1r5};
AMG1R5_OUT_ITERATOR stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5; // BiCGStab + amg1r5.

// инициализация компонент скорости константой.
// Единые значения для всей расчётной области.
// initialization value.
doublereal starting_speed_Vx = 0.0;
doublereal starting_speed_Vy = 0.0;
doublereal starting_speed_Vz = 0.0;

// Опорная линия для XY-Plot (variation Plot).
// Мы сохраняем опорную точку, через которую проходит линия, и направление 
// линии вдоль одной из осей декартовой прямоугольной системы координат.
doublereal Tochka_position_X0_for_XY_Plot = 0.0;
doublereal Tochka_position_Y0_for_XY_Plot = 0.0;
doublereal Tochka_position_Z0_for_XY_Plot = 0.0;
LINE_DIRECTIONAL idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; // 0 - Ox axis. 

// При iVar==TEMP && lw==1 выход из солвера можно осуществлять когда максимальная температура между V циклами отличается менее 0.5K.
bool bPhysics_stop = false;
// Особое выделение количества оперативной памяти для ПТБШ.
bool bPhysics_PTBSH_memory = false;
// Решаем только теплопередачу в твёрдом теле:
bool bonly_solid_calculation = false;

// Схемы для аппроксимации конвективного потока на неравномерной сетке.
// 3 августа 2015 схемы стали доступны через GUI пользователя
// в связи с чем объявление идентификаторов вынесено в самое начало кода.
// ограниченные схемы
const unsigned int UNEVEN_MUSCL = 1017;  // van Leer (1977)
const unsigned int UNEVEN_SOUCUP = 1018; // MINMOD
const unsigned int UNEVEN_HLPA = 1019;
const unsigned int UNEVEN_SMART = 1020; // Gaskell and Lau (1988)
const unsigned int UNEVEN_WACEB = 1021;
const unsigned int UNEVEN_SMARTER = 1022;
const unsigned int UNEVEN_STOIC = 1023; // Darwish (1993)
const unsigned int UNEVEN_CLAM = 1024;
const unsigned int UNEVEN_OSHER = 1025; // Chakravarthy and Osher (1983)
const unsigned int UNEVEN_VONOS = 1026;
const unsigned int UNEVEN_LPPA = 1027;
const unsigned int UNEVEN_EXPONENTIAL = 1028;
const unsigned int UNEVEN_SUPER_C = 1029;
const unsigned int UNEVEN_ISNAS = 1030;
const unsigned int UNEVEN_CUBISTA = 1031;
const unsigned int UNEVEN_GAMMA = 1032; // схема с параметром beta_m
const unsigned int UNEVEN_COPLA = 1033; // 1 08 2015
const unsigned int UNEVEN_SECBC = 1034; // 2 08 2015 Yu et al., (2001b) Сингапур, Малазия.
const unsigned int UNEVEN_SGSD = 1035; // 3 08 2015 Li and Tao (2002)

bool bglobal_first_start_radiation = true;

// Если мы решаем нестационарную задачу теплопередачи в твердом теле.
bool bglobal_unsteady_temperature_determinant = false;

// Выбор сеточного генератора:
// simplemeshgen == 0 или unevensimplemeshgen ==1.
// По умолчанию используется simplemeshgen == 0.
enum class CONFORMAL_MESH_GENERATOR_SELECTOR {
	SIMPLEMESHGEN_MESHER=0, UNEVENSIMPLEMESHGEN_MESHER=1,
	COARSEMESHGEN_MESHER=2};
CONFORMAL_MESH_GENERATOR_SELECTOR iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; // обычный сеточный генератор.
// нестационарное или стационарное моделирование.
enum class PHYSICAL_MODEL_SWITCH  {STEADY_TEMPERATURE = 0, UNSTEADY_TEMPERATURE = 1,
	MESHER_ONLY=2, CFD_STEADY=3, 
	STEADY_STATIC_STRUCTURAL=5, STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE=6,
	SECOND_TEMPERATURE_SOLVER=7,
    PREOBRAZOVATEL_FOR_REPORT=8, CFD_UNSTEADY=9, 
	NETWORK_T=10, NETWORK_T_UNSTEADY=11, UNSTEADY_STATIC_STRUCTURAL=12, 
	UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE = 13
};
//  0 - thermal only steady state calculation, // 0 - STEADY_TEMPERATURE.
//  1 - thermal only unsteady calculation, // 1 - UNSTEADY_TEMPERATURE.
//  2 - mesh generator only.
//  3 - fluid dynamic steady state.
//  5 - Static Structural (Thermal solver #2)
//  6 - Thermal Stress
//  7 - Unsteady thermal solver #2
//  8 - Visualisation only
//  9 - cfd unsteady fluid dynamic.
// 10 - NETWORK_T steady state. Графовый метод решения уравнения теплопроводности. Стационарное состояние.
// 11 - NETWORK_T unsteady calculation. Графовый метод решения уравнения теплопроводности. Нестационарное состояние.
// 12 - UNSTEADY STRUCTURAL MECHANICS. Нестационарная механика,
// 13 - UNSTEADY STRUCTURAL MECHANICS AND UNSTEADY TEMPERATURE CALCULATION. Нестационарная механика совместно с нестационарной теплопередачей.
PHYSICAL_MODEL_SWITCH steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE;

// Использовать ли адаптивные локально измельчённые расчётные сетки.
bool b_on_adaptive_local_refinement_mesh = false;
enum class TYPE_ALICE_MESH { ONE_PASS_COARSE_ALICE_MESH=0, MULTI_PASS_MEDIUM_ALICE_MESH=1};
TYPE_ALICE_MESH itype_ALICE_Mesh = TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH;// Тип АЛИС сетки.

// Октава 2 (один из возможных циклов).
typedef struct TPiecewiseConstantTimeStepLaw {
	doublereal time, timestep, m;
	TPiecewiseConstantTimeStepLaw() {
		time=0.0;
		timestep=0.0;
		m=0.0;
	}
} PiecewiseConstantTimeStepLawTimeStepLaw;

// 0 - Linear, 
	// 1 - Square Wave,
	// 2 - Square Wave 2, 
	// 3 - Hot Cold (Нагрев до заданного момента времени, а потом остывание заданное время. Шаг по времени двойной логарифмический.)
	// 4 - Piecewise Constant 20.12.2019 - Закон изменения заданный пользователем.
enum class TIME_STEP_lAW_SELECTOR {
	LINEAR = 0, SQUARE_WAVE = 1, SQUARE_WAVE2 = 2, HOT_COLD = 3, PIECEWISE_CONSTANT = 4
};

// Для задания закона изменения шага по времени и тепловой мощности от времени
// при нестационарном моделировании.
typedef struct TTimeStepLaw
{
	// 0 - Linear, 
	// 1 - Square Wave,
	// 2 - Square Wave 2, 
	// 3 - Hot Cold (Нагрев до заданного момента времени, а потом остывание заданное время. Шаг по времени двойной логарифмический.)
	// 4 - Piecewise Constant 20.12.2019 - Закон изменения заданный пользователем.
	TIME_STEP_lAW_SELECTOR id_law;
	doublereal Factor_a_for_Linear;
	doublereal tau; // длительность импульса для Square Wave
	// 06_03_2017 скважность может быть и дробной.
	doublereal Q; // Скважность для Square Wave.
	// Импульсный режим для альтернативного Square Wave 2.
	doublereal m1, tau1, tau2, tau_pause, T_all;
	integer n_cycle; // 20 Циклов.
	// hot cold reshime (double linear)
	doublereal on_time_double_linear; //3c включено.
	// 4 закон изменения шага по времени.
	integer n_string_PiecewiseConst;
	PiecewiseConstantTimeStepLawTimeStepLaw* table_law_piecewise_constant;
	TTimeStepLaw() {
		// 0 - Linear, 
		// 1 - Square Wave,
		// 2 - Square Wave 2, 
		// 3 - Hot Cold (Евдокимова Н.Л.)
		// 4 - Piecewise Constant 20.12.2019 - Закон изменения заданный пользователем.
		id_law= TIME_STEP_lAW_SELECTOR::LINEAR;
		Factor_a_for_Linear=0.2;
		tau=60.0E-6; // длительность импульса для Square Wave
		// 06_03_2017 скважность может быть и дробной.
		Q=10.0; // Скважность для Square Wave.
		// Импульсный режим для альтернативного Square Wave 2.
		m1=1.0; tau1=0.0; tau2=0.0; tau_pause=0.0; T_all=0.0;
		n_cycle=20; // 20 Циклов.
		// hot cold reshime (double linear)
		on_time_double_linear=3.0; //3c включено.
		// 4 закон изменения шага по времени.
		n_string_PiecewiseConst = 0;
		table_law_piecewise_constant = nullptr;
	}
} TimeStepLaw;

TimeStepLaw glTSL;

// 24 декабря 2016. 
// Для ускорения счёта нелинейных задач в РУМБА 0.14 решателе.
typedef struct TQuickNonlinearBoundaryCondition {
	doublereal emissivity;
	doublereal ViewFactor; // Фактор видимости.
	doublereal Tamb, dS;
	doublereal film_coefficient;
	bool bactive;
	bool bStefanBolcman_q_on;
	bool bNewtonRichman_q_on;
	TQuickNonlinearBoundaryCondition() {
		emissivity=0.8;
		ViewFactor = 1.0; // Фактор видимости.
		Tamb=20.0, dS=0.0;
		film_coefficient=3.0;
		bactive=false;
		bStefanBolcman_q_on = false;
		bNewtonRichman_q_on=false;
	}
} QuickNonlinearBoundaryCondition;

QuickNonlinearBoundaryCondition* qnbc = nullptr;
integer iadd_qnbc_maxelm = 0; // для аддитивного сдвига
bool b_sign_on_nonlinear_bc = false;


// Считаем ли мы SIMPLE алгоритмом.
// Это нужно для более точной настройки невязки для уравнения теплопередачи.
// внутри BiCGStab_internal3 решателя.
bool bSIMPLErun_now_for_temperature = false;
// Количество итераций SIMPLE алгоритма 
// которые задаёт пользователь в интерфейсе.
unsigned int number_iteration_SIMPLE_algorithm = 0; // default - 0
// Это нужно для более точной настройки невязки для уравнения теплопередачи
// при расчёте amg1r5 алгоритмом задач с естественной конвекцией.
bool bSIMPLErun_now_for_natural_convection = false;
// Дополнительная нижняя релаксация для температуры.
doublereal* told_temperature_global_for_HOrelax = nullptr;

/*
для внутренних плоских источников тепла организуется виртуальная грань.
эта грань общая для двух КО располагающихся по-бокам от источника.
принцип единственности должен приводить к тому что теплопроводность грани должна быть единственна
при обработке обоих контрольных объемов примыкающих к данной грани. В качестве теплопроводности
берётся среднее геометрическое.
Устаревший код: плоские бесконечно тонкие источники тепла больше неиспользуются.
Используется объёмное тепловыделение которое учитывается в правой части СЛАУ. 06.01.2020
*/
bool *sourse2Dproblem=nullptr;
doublereal *conductivity2Dinsource=nullptr;

// дополнительная нижняя релаксация.
bool bHORF = false;
bool bdontstartsolver = false;
doublereal* bPamendment_source_old = nullptr;
doublereal* bsource_term_radiation_for_relax = nullptr;
doublereal* b_buffer_correct_source = nullptr;
// Во избежании расходимости РУМБА_0.14 на первой же итерации.
bool bfirst_start_nonlinear_process = true;

// Условие Ньютона-Рихмана по дефолту для температуры.
// 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
DEFAULT_CABINET_BOUNDARY_CONDITION adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC;
// При нелинейном граничном условии на стенке надо делать лишь небольшое число V циклов. 
bool breakRUMBAcalc_for_nonlinear_boundary_condition = false;
bool bvacuumPrism = false; // Наличие вакуумных промежутков.
bool bdouble_vacuum_PRISM = false; // Двойной вакуумный промежуток. Это нужно для временной блокировки Стефана - Больцмана на твёрдой стенке.
bool bBlockStefanBolcman = false; // Если true то блокируем Стефана Больцмана.
doublereal film_coefficient = 0.0; // Коэффициент теплоотдачи.
doublereal operating_temperature_for_film_coeff = 20.0; // Tamb for Newton-Richman condition.
// Если начальное поле температур будет инициализировано значением operating_temperature_for_film_coeff
// То в условии Ньютона-Рихмана мы получим нулевой тепловой поток (однородные условия Неймана) что при
// наличии мощности тепловыделения и отсутствию условий Дирихле приведёт к бесконечным температурам и расходимости
// Вычислительного процесса. Чтобы этого избежать используется переменная blocker_Newton_Richman.
bool blocker_Newton_Richman = true;


// 0 - визуализация и жидкости и твердого тела.
// 1 - визуализация только твёрдого тела.
enum class WHAT_VISIBLE_OPTION {FLUID_AND_SOLID_BODY_VISIBLE=0, ONLY_SOLID_BODY_VISIBLE=1};
WHAT_VISIBLE_OPTION ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;

// переключение между алгебраическим многосеточным методом и алгоритмом Ван дер Ворста BiCGStab+ILU2.
// 0 - алгоритм BiCGStab + ILU2.
// 1 - алгоритм Джона Руге и Клауса Штубена алгебраического многосеточного метода amg1r5 (r6).
// 2 - BiCGStab + ADI (Lr1sk).
// 3 - Gibrid: velocity bicgstab + ilu(lfil), Pressure - РУМБА v0.14.
// 4 - BiCGStab + AINV N.S.Bridson nvidia cusp 0.5.1 library.
// 5 - AMGCL bicgstab+samg Денис Демидов.
// 6 - Nvidia cusp 0.5.1 library BiCGStab +samg.
// 7 - Algebraic Multigrid Румба v0.14.
integer iswitchsolveramg_vs_BiCGstab_plus_ILU2 = 0; // BiCGStab + ILU2.
// 0 - BiCGStab+ILU6, 1 - Direct, 2 - Румба 0.14, 3 - amg1r5, 4 - AMGCL_SECONT_T_SOLVER.
// for Stress Solver
enum class SECOND_T_SOLVER_ID_SWITCH { BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER=0,
	DIRECT_SECOND_T_SOLVER=1, CAMG_RUMBA_v0_14_SECOND_T_SOLVER=2,
	AMG1R5_SECOND_T_SOLVER=3, AMGCL_SECONT_T_SOLVER=4};
SECOND_T_SOLVER_ID_SWITCH iswitchsolveramg_vs_BiCGstab_plus_ILU6 = SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; // BiCGStab + ILU6.

bool bwait = false; // если false, то мы игнорируем getchar().
// Если задать нечто отличное от 1e-10 то программа уходит очень долгий цикл
const doublereal admission = 1.0e-30; //1.0e-10 // для определения совпадения двух вещественных чисел.

unsigned int calculation_vorst_seach_time = 0;

// Если температура канала превысит 
// температуру в 200 градусов Цельсия 
// то прибор Выйдет из строя (сгорит).
// В случае превышения температуры равной TEMPERATURE_FAILURE_DC
// на консоль печатается предупреждающее сообщение и после того
// как пользователь прочтёт сообщение, осуществляется выход из программы.
const doublereal TEMPERATURE_FAILURE_DC = 5000.2;


// Красавин Денис Андреевич Повышение порядка точности аппроксимации
// граничного условия. См. его диссертацию а также книжку С. Патанкара.
// BETA_PRECISION 1.0 4/3=1.333333333 6/5=1.2
const doublereal BETA_PRECISION = 1.0;

// схемы для аппроксимации конвекции-диффузии
const unsigned int CR = 1; // Центрально-разностная
const unsigned int UDS = 2; // Противопоточная первого порядка
const unsigned int COMB = 3; // Комбинированная 
const unsigned int POLY = 4; // Полиномиальная C. Патанкара
const unsigned int EXP = 5; // экспоненциальная схема
const unsigned int BULG = 6; // схема В.К. Булгакова формула (23) из статьи
const unsigned int POW = 7; // показательная

// UDS  см. my_approx_convective.c
unsigned int iFLOWScheme = UDS; // Противопоточная первого порядка
unsigned int iTEMPScheme = UDS; // Противопоточная первого порядка

// включает более быстро сходящийся алгоритм SIMPLEC
// SIMPLEC Van Doormal and Raithby, 1984 год.
// SIMPLEC отличается от SIMPLE только в двух вещах:
// 1. В SIMPLEC не используется нижняя релаксация для давления при коррекции давления, т.е. alphaP=1.0.
// 2. В SIMPLEC псевдо время пропорционально tau ~ alpha/(1-alpha), а в SIMPLE tau ~ alpha. 
// В остальном алгоритмы полностью совпадают.
enum class SIMPLE_CFD_ALGORITHM {SIMPLE_Carretto = 0, SIMPLEC_Van_Doormal_and_Raithby = 1};
//const unsigned int SIMPLE_Carretto = 0; // алгоритм SIMPLE Carretto et al., 1973 используется по умолчанию.
//const unsigned int SIMPLEC_Van_Doormal_and_Raithby = 1; // алгоритм SIMPLEC Van Doormal and Raithby, 1984 год.
SIMPLE_CFD_ALGORITHM iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;// SIMPLE_Carretto SIMPLEC_Van_Doormal_and_Raithby

// Точность решения СЛАУ для всех методов для нормы
// (residual,residual) где () скалярное произведение в R^n.
// значение этой константы должно зависеть от точности аппроксимации 
// дифф. уравнения в часных производных на некоторой сетке по способу
// Метода Контрольного Объёма (FVM). Поэтому эта константа будет переопределена
// согласованным образом в модуле mysolverv0_03.c. Под согласованным понимается то,
// что СЛАУ не имеет смысла решать точнее чем обеспечивает точность дискретного аналога,
// иными словами все последующие итерации значение невязки для которых меньше точности аппроксимации
// на данной сетке работают вхолостую (бесполезны). Это нужно учитывать.
// Для сошедшейся задачи подходит уровень среднеквадратических невязок 1.0e-4
// Источник информации мануал по CFX на русском.
// Если норма Чебышёва то невязка сошедшаяся равна 1.0e-3.
// Источник опять же мануал по CFX.
doublereal dterminatedTResudual = 1e-5; // для МСГ Congruate Gradients а также BiCGStab_internal3.

doublereal globalEndTimeUnsteadyTemperatureCalculation = 1.0; // физическое время при нестационарном моделировании теплопередачи в твёрдом теле. 

// используется для ускорения решения 
// задачи твёрдотельной теплопередачи.
bool bsolid_static_only = false;

const integer inumcore = 2; // число ядер процессора
const bool bparallelizm_old = false;

// Структура границ деления:
typedef struct TPARBOUND {
	integer ileft_start, ileft_finish;
	integer iright_start, iright_finish;
	integer iseparate_start, iseparate_finish;
	bool active; // активность декомпозиции.
	TPARBOUND() {
		ileft_start = -1; ileft_finish = -2;
		iright_start = -1; iright_finish = -2;
		iseparate_start = -1; iseparate_finish = -2;
		active=false; // активность декомпозиции.
	}
} PARBOUND;


// Структура данных используемая для распараллеливания.
typedef struct TPARDATA {
	integer ncore; // 1, 2, 4, 8.
	integer *inumerate;
	// это для ncore==2;
	PARBOUND b0;
	// это для ncore==4;
	PARBOUND b00, b01;
	// это для ncore==8;
	PARBOUND b000, b001, b010, b011;
	TPARDATA() {
		ncore=2; // 1, 2, 4, 8.
		inumerate=nullptr;
	}
} PARDATA;

PARDATA nd;

// используется в Румбе для ускорения счёта.
doublereal* rthdsd_no_radiosity_patch = nullptr;

// При дроблении (если bdroblenie4=true)
// каждая из шести граней иожет граничить с четырьмя соседними ячейками.
typedef struct TALICE_PARTITION {
	bool bdroblenie4;
	integer iNODE1, iNODE2, iNODE3, iNODE4;
	TALICE_PARTITION() {
		bdroblenie4=false;
		iNODE1=-1; iNODE2=-1; iNODE3=-1; iNODE4=-1;
	}
	void define_structural_mesh_neighbour(integer id) {
		bdroblenie4 = false;
		iNODE1 = id; iNODE2 = -1; iNODE3 = -1; iNODE4 = -1;
	}
} ALICE_PARTITION;


// Глобальное сгущение АЛИС сетки к источникам тепла. 
// Это должно увеличить точность нахождения температуры.
const bool b_thermal_source_refinement = true;

#include "adaptive_local_refinement_mesh.cpp" // АЛИС
#include "constr_struct.cpp" // заполнение структур данных TEMPER и FLOW
#include "uniformsimplemeshgen.cpp" // сеточный генератор

#include "my_LR.c" // полилинейный метод

// 8 января 2016.
const bool bvery_big_memory = true; // true нет записи в файл всё храним в оперативной памяти. Это существенно быстрее по скорости.

// Данные определения пределов раньше находились в файле MenterSST.cpp.
const doublereal K_limiter_min = 1.0e-20; //1.0e-14; // 1.0e-14 Fluent limits
// 0.1 тоже подходит и сходимость стала лучше для k по сравнению со случаем 1.0.
// при значении 1.0E-20 расходимость. Значение k на уровне 1e-14 до 1e-10. Октава 1 0.5м/с.
// Справляется только метод сглаженной аггрегации ddemidov AMGCL.
// Вточном решении (сошедшимся) omega не опускалась ниже значения 1.5.
const doublereal Omega_limiter_min = 0.1; // 1.0; иначе будет расходимость. // 1.0e-20 Fluent limits
const doublereal Epsilon_limiter_min = 1.0e-20; // 1.0e-14; значение требует уточнения... 1.0e-20 Fluent limits

UNION* my_union = nullptr; // для объединения.

// Глобальное объявление
TEMPER my_global_temperature_struct;
integer flow_interior = 0; // Суммарное число FLUID зон
FLOW* f = nullptr;




// экспорт картинки в программу tecplot360
#include "my_export_tecplot3.c"

#include "my_material_properties.c" // библиотека реальных свойств материалов


#include "my_linalg.cpp" // самописные функции линейной алгебры
// Для функций: 
// eqsolve_simple_gauss - решает СЛАУ методом исключения Гаусса
// eqsolv_simple_holesskii - решает СЛАУ методом разложения Холесского

// аппроксимация обобщённого уравнения конвекции-диффузии
// на совмещённой сетке
#include "pamendment3.c"

#include "shortest_distance.cpp" // вычисление кратчайшего расстояния до стенки


// Информация о формуле вычисления невязок заимствована из 
// icepak user guide.
typedef struct TFLUENT_RESIDUAL{
	// Данные невязки приводятся на момент конца решения СЛАУ.
	// невязки согласованные с программой FLUENT
	// т.е. вычисляемые по формуле fluent.
	doublereal res_vx; // невязка X скорости
	doublereal res_vy; // невязка Y скорости
	doublereal res_vz; // невязка Z скорости
	doublereal res_no_balance; // несбалансированные источники массы.
	doublereal operating_value_b; // значение несбалансированных источников массы с предыдущей итерации.
	doublereal res_nusha; // Невязка модифицированной кинетической турбулентной вязкости.
	doublereal res_turb_kinetik_energy; // Невязка кинетической энергии турбулентных пульсаций в модели SST K-Omega.
	doublereal res_turb_omega; // Невязка удельной скорости диссипации кинетической энергии турбулентных пульсаций в модели SST K-Omega.
	doublereal res_turb_kinetik_energy_std_ke; // Невязка кинетической энергии турбулентных пульсаций в двухслойной модели k-epsilon
	doublereal res_turb_epsilon; // Невязка скорости диссипации кинетической энергии турбулентных пульсаций в двухслойной модели k-epsilon
	TFLUENT_RESIDUAL() {
		// Данные невязки приводятся на момент конца решения СЛАУ.
		// невязки согласованные с программой FLUENT
		// т.е. вычисляемые по формуле fluent.
		res_vx=1.0; // невязка X скорости
		res_vy=1.0; // невязка Y скорости
		res_vz=1.0; // невязка Z скорости
		res_no_balance=1.0; // несбалансированные источники массы.
		operating_value_b=1.0; // значение несбалансированных источников массы с предыдущей итерации.
		res_nusha=1.0; // Невязка модифицированной кинетической турбулентной вязкости.
		res_turb_kinetik_energy = 1.0; // Невязка кинетической энергии турбулентных пульсаций в модели SST K-Omega.
		res_turb_omega = 1.0; // Невязка удельной скорости диссипации кинетической энергии турбулентных пульсаций в модели SST K-Omega.
		res_turb_kinetik_energy_std_ke = 1.0; // Невязка кинетической энергии турбулентных пульсаций в двухслойной модели k-epsilon
		res_turb_epsilon = 1.0; // Невязка скорости диссипации кинетической энергии турбулентных пульсаций в двухслойной модели k-epsilon
	}
} FLUENT_RESIDUAL;


// соединение решателя
#include "mysolverv0_03.c"


// нестационарный солвер для температуры
// на основе стационарного солвера,
// а также нестационарный солвер для 
// гидродинамики на основе стационарного солвера.
#include "my_unsteady_temperature.cpp"

// Препроцессинг для параллельной обработки.
#include "my_nested_dissection.cpp"

#include <ctime> // для замера времени выполнения.

integer ltdp=0; // количество таблично заданных мощностей от температуры и смещения стока.
TEMP_DEP_POWER* gtdps=nullptr; // Garber temperature depend power sequence. 

// база данных материалов:
integer lmatmax=0; // максимальное число материалов
TPROP* matlist=nullptr; // хранилище базы данных материалов

void check_data(TEMPER t) {
	if (t.potent != nullptr) {
		for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
			if (t.potent[i] != t.potent[i]) {
				std::cout << "t.potent[" << i << "] is " << t.potent[i] << "\n";
				system("pause");
			}
		}
	}
} // check_data

int main_body(char ch_EXPORT_ALICE_ONLY = 'y') {
	//printLOGO();

	//system("PAUSE");	

	// количество блоков, источников и стенок, юнионов.
	integer lb = 0, ls = 0, lw = 0, lu = 0;
	BLOCK* b = nullptr; // список блоков
	SOURCE* s = nullptr; // список источников
	WALL* w = nullptr; // список твёрдых стенок

	// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
	// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
	// при повторном обращении идет обращение к ячейки хеш-таблицы.
	// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
	// Хеш-таблицы для automatic
	// Аллокация оперативной памяти под хеш-таблицы.
	shorter_hash_X = new doublereal[isize_shorter_hash];
	shorter_hash_Y = new doublereal[isize_shorter_hash];
	shorter_hash_Z = new doublereal[isize_shorter_hash];
	bshorter_hash_X = new bool[isize_shorter_hash];
	bshorter_hash_Y = new bool[isize_shorter_hash];
	bshorter_hash_Z = new bool[isize_shorter_hash];


	// Инициализация, показываем всё.
	pfpir.fmaximum = 1.0e+30;
	pfpir.fminimum = -1.0e+30;
	pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL;


	// 22.01.2017 инициализация.
	eqin.fluidinfo = nullptr;
	my_global_temperature_struct.rootBT = nullptr;
	my_global_temperature_struct.rootSN = nullptr;
	my_global_temperature_struct.rootWE = nullptr;

	// 29 10 2016.
	// Инициализация общей памяти в ILU буфере.
	milu_gl_buffer.alu_copy = nullptr;
	milu_gl_buffer.jlu_copy = nullptr;
	milu_gl_buffer.ju_copy = nullptr;

	my_amg_manager_init();

	// Замер времени.
	calculation_main_start_time_global_Depend = 0; // начало счёта мс.
	unsigned int calculation_main_end_time = 0; // окончание счёта мс.
	unsigned int calculation_main_seach_time = 0; // время выполнения участка кода в мс.

	calculation_main_start_time_global_Depend = clock(); // момент начала счёта.

	bool bextendedprint = false; // печать на граничных узлах рассчитанных полей.


	//std::locale::global(std::locale("en_US.UTF-8"));
	//system("mode con cols=166 lines=12000");
	// к примеру для того чтобы поменять цвет шрифта в консоли нужно сделать к примеру так
	//HANDLE hOCol = GetStdHandle(STD_OUTPUT_HANDLE);
	//SetConsoleTextAttribute(hOCol, FOREGROUND_GREEN);
	// Установка атрибутов консоли (белый фон)	
	//SetConsoleTextAttribute(hOCol, BACKGROUND_BLUE |
	//	BACKGROUND_GREEN |
	//	BACKGROUND_RED |
	//	BACKGROUND_INTENSITY);

	//system("cls");
	// Недостаток в том, что белый фон появляется лишь там где напечатан текст.

	// для кабинета.
	// Их нельзя делать глобальными. 
	doublereal* xpos = nullptr, * ypos = nullptr, * zpos = nullptr;
	doublereal* xposadd = nullptr, * yposadd = nullptr, * zposadd = nullptr;


	std::cout << "AliceFlow 3D x64 v0.48\n";
#ifdef _OPENMP 
	omp_set_num_threads(inumcore); // установка числа потоков
#endif

	 
	

		//ilu0_Saadtest();
		//printf("the end Saad ilu0 test\n");
		//system("PAUSE");

		// количество точек по каждой из осей.
		//integer inx=120, iny=64, inz=64;
		integer inx = 30, iny = 30, inz = 30;
		integer inxadd = -1, inyadd = -1, inzadd = -1;
		doublereal dgx = 0.0, dgy = 0.0, dgz = 0.0; // сила тяжести
		doublereal operatingtemperature = 20.0; // Operating Temperature 20.0 Град. С.

		lu = 0;
		// lu, my_union
		loadFromFile();
		premeshin("premeshin.txt", lmatmax, lb, ls, lw, matlist, b, s, w,
			dgx, dgy, dgz, inx, iny, inz, operatingtemperature, ltdp, gtdps, lu, my_union);
		freeStringList();

		// Проверяет если ли выход за пределы кабинета
		// среди блоков, стенок и источников тепла. 02.08.2019.
		// а также асемблесов. 23.07.2020
		BODY_CHECK(b, lb, w, lw, s, ls, my_union, lu);

		init_QSBid(lb, b, w, lw, s, ls); // Для ускоренной работы функции myisblock_id.
		
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY)) {
			// При решении уравнений гидродинамики мы удаляем старый load.txt файл.
			remove("load.txt");
		}
		
		if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT) {
			// Преобразование файла с результатами вычислений.
			// для написания отчетов. 05.01.2018.
			tecplot360patcher_for_print_in_report();
			exit(1);
		}

		integer iCabinetMarker = 0;
		if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER) {
			simplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
		}
		else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER) {
			unevensimplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				dgx, dgy, dgz, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
			// генератор неравномерной сетки с автоматической балансировкой неравномерности.
		}
		else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
			// Я стремился сделать coarse Mesh как в Icepak.
			// Реальные модели чрезвычайно многоэлементно-большеразмерные, а 
			// ресурсы персонального компьютера чрезвычайно слабы, т.к. CPU уперлись в 4ГГц а
			// правильное распараллеливание отдельная большая научная проблема.
			coarsemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
		}
		else {
			switch (iswitchMeshGenerator) {
			case CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined "<< "CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER" << "\n";
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER :
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER" << "\n";
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER" << "\n";
				break;
			default :
				std::cout << "error: your mesh generator is undefined "  << "\n";
				break;
			}
			

			system("pause");
			exit(1);
		}


		// Строим расчётную сетку в объединениях.
		for (integer iu = 0; iu < lu; iu++) {
			my_union[iu].inxadd = -1;
			my_union[iu].inyadd = -1;
			my_union[iu].inzadd = -1;
			my_union[iu].xposadd = nullptr;
			my_union[iu].yposadd = nullptr;
			my_union[iu].zposadd = nullptr;
			my_union[iu].xpos = nullptr;
			my_union[iu].ypos = nullptr;
			my_union[iu].zpos = nullptr;
			integer iup1 = iu + 1;
			switch (my_union[iu].iswitchMeshGenerator) {
			case CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER : simplemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
				my_union[iu].inx, my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd, my_union[iu].inxadd,
				my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER: unevensimplemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos, my_union[iu].inx,
				my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				dgx, dgy, dgz, my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd,
				my_union[iu].inxadd, my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				// генератор неравномерной сетки с автоматической балансировкой неравномерности.
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER: coarsemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
				my_union[iu].inx, my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd, my_union[iu].inxadd,
				my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				break;
			default:
				coarsemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
					my_union[iu].inx, my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
					my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd, my_union[iu].inxadd,
					my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				break;
			}
			// Требуется использовать именно addboundary_rudiment а не addboundary, так как в асемблес нужно 
			// добавить по возможности все подходящие линии кабинета. 14.03.2019.
			// Температура чувствует такое добавление сеточных линий очень ощутимо. 
			// Добавляем сеточный линии глобального кабинета для повышения точности аппроксимации.
			for (integer i76 = 0; i76 <= inx; i76++) {
				// Добавляем глобальные сеточные линии кабинета.
				if ((xpos[i76] >= my_union[iu].xS) && (xpos[i76] <= my_union[iu].xE)) {
					addboundary_rudiment(my_union[iu].xpos, my_union[iu].inx, xpos[i76], YZ_PLANE, b, lb, w, lw, s, ls);
				}
			}
			Sort_method(my_union[iu].xpos, my_union[iu].inx);
			for (integer i76 = 0; i76 <= iny; i76++) {
				// Добавляем глобальные сеточные линии кабинета.
				if ((ypos[i76] >= my_union[iu].yS) && (ypos[i76] <= my_union[iu].yE)) {
					addboundary_rudiment(my_union[iu].ypos, my_union[iu].iny, ypos[i76], XZ_PLANE, b, lb, w, lw, s, ls);
				}
			}
			Sort_method(my_union[iu].ypos, my_union[iu].iny);
			for (integer i76 = 0; i76 <= inz; i76++) {
				// Добавляем глобальные сеточные линии кабинета.
				if ((zpos[i76] >= my_union[iu].zS) && (zpos[i76] <= my_union[iu].zE)) {
					addboundary_rudiment(my_union[iu].zpos, my_union[iu].inz, zpos[i76], XY_PLANE, b, lb, w, lw, s, ls);
				}
			}
			Sort_method(my_union[iu].zpos, my_union[iu].inz);
		}



		if (b_on_adaptive_local_refinement_mesh) {
			// задача 1TT113 22.03.2018:
			// AliceCorse 1266911 узел, в гидродинамической подобласти 93592. 
			// 40мин 33с снятие переходной характеристики до 10000с.
			// Время построения сетки 2мин 6с. 14 уровней иерархии.
			// AliceFine 1092571 узел, в гидродинамической подобласти 8050.
			// AliceFine (повторный вызов) 1335076 узел, в гидродинамической подобласти 93448.
			// Время построения сетки 4мин 25с. 14 уровней иерархии.
			// Время построения сетки (повторный вызов) 6мин 32с. 14 уровней иерархии.



			std::cout << "starting ALICE\n";
			if (0) {
				if (TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH == itype_ALICE_Mesh) {
					// Так делать ни в коем случае нельзя по причине нехватки оперативной памяти.
					doublereal* xpos_copy = nullptr;
					// 10 слишком большое значение константы.
					const integer jiterM = my_amg_manager.nu1_Temperature;
					// десятикратное дробление каждого интервала пополам.
					for (integer jiter = 0; jiter < jiterM; jiter++) {
						xpos_copy = new doublereal[2 * (inx + 1) - 1];
#pragma omp parallel for
						for (integer i74 = 0; i74 < inx; i74++) {
							xpos_copy[2 * i74] = xpos[i74];
							xpos_copy[2 * i74 + 1] = 0.5 * (xpos[i74] + xpos[i74 + 1]);
						}
						xpos_copy[2 * (inx + 1) - 2] = xpos[inx];
						delete[] xpos;
						xpos = nullptr;
						xpos = new doublereal[2 * (inx + 1) - 1];
#pragma omp parallel for
						for (integer i74 = 0; i74 < 2 * (inx + 1) - 1; i74++) {
							xpos[i74] = xpos_copy[i74];
						}
						delete[] xpos_copy;
						xpos_copy = nullptr;
						inx = 2 * (inx + 1) - 2;
					}

					for (integer jiter = 0; jiter < jiterM; jiter++) {
						xpos_copy = new doublereal[2 * (iny + 1) - 1];
#pragma omp parallel for
						for (integer i74 = 0; i74 < iny; i74++) {
							xpos_copy[2 * i74] = ypos[i74];
							xpos_copy[2 * i74 + 1] = 0.5 * (ypos[i74] + ypos[i74 + 1]);
						}
						xpos_copy[2 * (iny + 1) - 2] = ypos[iny];
						delete[] ypos;
						ypos = nullptr;
						ypos = new doublereal[2 * (iny + 1) - 1];
#pragma omp parallel for
						for (integer i74 = 0; i74 < 2 * (iny + 1) - 1; i74++) {
							ypos[i74] = xpos_copy[i74];
						}
						delete[] xpos_copy;
						xpos_copy = nullptr;
						iny = 2 * (iny + 1) - 2;
					}

					for (integer jiter = 0; jiter < jiterM; jiter++) {
						xpos_copy = new doublereal[2 * (inz + 1) - 1];
#pragma omp parallel for
						for (integer i74 = 0; i74 < inz; i74++) {
							xpos_copy[2 * i74] = zpos[i74];
							xpos_copy[2 * i74 + 1] = 0.5 * (zpos[i74] + zpos[i74 + 1]);
						}
						xpos_copy[2 * (inz + 1) - 2] = zpos[inz];
						delete[] zpos;
						zpos = nullptr;
						zpos = new doublereal[2 * (inz + 1) - 1];
#pragma omp parallel for
						for (integer i74 = 0; i74 < 2 * (inz + 1) - 1; i74++) {
							zpos[i74] = xpos_copy[i74];
						}
						delete[] xpos_copy;
						xpos_copy = nullptr;
						inz = 2 * (inz + 1) - 2;
					}
				}
			}

			// Это слишком большая величина на многих промышленных задачах,
			// её нельзя задавать во избежании сбоя в операторе new или malloc.
			integer maxelm_loc = (inx + 1) * (iny + 1) * (inz + 1);
			bool bOkal = alice_mesh(xpos, ypos, zpos, inx, iny, inz, b, lb, lw, w, s, ls, maxelm_loc, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
			//system("PAUSE");



			if (0 || itype_ALICE_Mesh == TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH/*1*/) {
				// Вызываем повторные генерации.

				/*
				Когда было введено ограничение на появление сеточных линий
				при котором близкорасположенные сеточные линии игнорируются
				мы не можем удовлетворить всем критериям AliceMedium сетки
				и она дробиться бесконечно а новые сеточные линии не появляются.
				Для того чтобы избежать этого мы прерываем построение сетки
				в тот момент когда сетка перестает меняться.
				17.08.2019
				*/
				doublereal dSTOP_flag1 = 1.0e+4; // будем сравнивать нормы.
				doublereal dSTOP_flag2 = 1.0e+1; // будем сравнивать нормы.

				while (!bOkal) {
					std::cout << "repeat call ALICE...\n";
					//system("PAUSE");

					/* 3.09.2017
					АЛИС сетку лучше строить когда при дроблении ячейки соответствующая геометрическая длина делится
					ровно пополам. Раньше делился пополам целочисленный индекс в результате чего на нашей существенно неравномерной
					сетке АЛИС ячейки были очень сильно вытянутые. В результате точность аппроксимации чрезвычайно страдала. Теперь когда
					дробиться пополам именно геометрическая длина может не хватать ячеек сетки для балансировки сетки (т.е. будет невозможно
					без добавления новых сеточных линий сделать чтобы уровни соседних ячеек отличались не более чем на 1. Поэтому теперь
					такие невозможные для дробления ячейки ищутся в функции if_disbalnce(...) и для каждой такой ячейки принимается решение
					добавить соответствующую новую сеточную линию. Алгоритм освобождает память и возвращается на исходные позиции и построение
					АЛИС сетки начинается заново только теперь базовая сетка уже содержит недостающие сеточные линии.
					*/


					// Нужно освободить память из под octree дерева и перестроить сетку.
					std::cout << "free octree start...\n";
					//system("PAUSE");
					//system("PAUSE");
					free_octree(oc_global, maxelm_loc);
					delete[] my_ALICE_STACK;
					top_ALICE_STACK = 0;
					std::cout << "free octree end...\n";
					doublereal t_1 = NormaV(xpos, inx + 1);
					doublereal t_2 = NormaV(ypos, iny + 1);
					doublereal t_3 = NormaV(zpos, inz + 1);
					dSTOP_flag2 = sqrt(t_1 * t_1 + t_2 * t_2 + t_3 * t_3);
					std::cout << "comparison mesh " << fabs(dSTOP_flag2 - dSTOP_flag1) << "\n";
					//system("pause");

					// Новое построение расчётной сетки.
					delete[] xpos;
					xpos = nullptr;
					inx = 0;
					delete[] ypos;
					ypos = nullptr;
					iny = 0;
					delete[] zpos;
					zpos = nullptr;
					inz = 0;

					std::cout << "free xpos, ypos, zpos\n";
					//system("PAUSE");

					integer iCabinetMarker = 0;
					if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER) {
						simplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union,
							matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
					}
					else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER) {
						unevensimplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union,
							matlist, dgx, dgy, dgz, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker); // генератор неравномерной сетки с автоматической балансировкой неравномерности.
					}
					else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
						// Я стремился сделать coarse Mesh как в ANSYS Icepak.
						// Реальные модели чрезвычайно многоэлементно-большеразмерные, а 
						// ресурсы персонального компьютера чрезвычайно слабы, т.к. CPU уперлись в 4ГГц а
						// правильное распараллеливание отдельная большая научная проблема.
						coarsemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union,
							matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
					}
					else {

						printf("error: your mesh generator is undefined %d\n", iswitchMeshGenerator);

						system("pause");
						exit(1);
					}
					// Новый запуск АЛИС меширования.

					printf("new construct xpos, ypos, zpos\n");
					//system("PAUSE");

					bOkal = alice_mesh(xpos, ypos, zpos, inx, iny, inz, b, lb, lw, w, s, ls, maxelm_loc, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);

					if (fabs(dSTOP_flag2 - dSTOP_flag1) < 1.0e-20) {
						bOkal = true;
						break;
					}
					dSTOP_flag1 = dSTOP_flag2;

					//system("PAUSE");
				}
			}
			printf("end ALICE\n");
		}

		iCabinetMarker = 0;
		load_TEMPER_and_FLOW(my_global_temperature_struct, f, inx, iny, inz, xpos, ypos, zpos, flow_interior,
			b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
			dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iCabinetMarker);

		my_global_temperature_struct.operatingtemperature = operatingtemperature;


		for (integer iu = 0; iu < lu; iu++) {
			integer iup1 = iu + 1;
			load_TEMPER_and_FLOW(my_union[iu].t, my_union[iu].f,
				my_union[iu].inx, my_union[iu].iny, my_union[iu].inz,
				my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
				my_union[iu].flow_interior,
				b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
				dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iup1);

			my_union[iu].t.operatingtemperature = operatingtemperature;
		}




		// Эти копии данных нужны для полного восстановления данных.
		my_global_temperature_struct.inx_copy = inx;
		my_global_temperature_struct.iny_copy = iny;
		my_global_temperature_struct.inz_copy = inz;
		my_global_temperature_struct.operatingtemperature_copy = operatingtemperature;
		my_global_temperature_struct.xpos_copy = new doublereal[inx + 1];
		my_global_temperature_struct.ypos_copy = new doublereal[iny + 1];
		my_global_temperature_struct.zpos_copy = new doublereal[inz + 1];
		// Данная информация нужна для экономии оперативной памяти,
		// некоторые данные будут выгружены из озу а потом восстановлены 
		// путём вычисления.
#pragma omp parallel for
		for (integer i_7 = 0; i_7 < inx + 1; i_7++) {
			my_global_temperature_struct.xpos_copy[i_7] = xpos[i_7];
		}
#pragma omp parallel for
		for (integer i_7 = 0; i_7 < iny + 1; i_7++) {
			my_global_temperature_struct.ypos_copy[i_7] = ypos[i_7];
		}
#pragma omp parallel for
		for (integer i_7 = 0; i_7 < inz + 1; i_7++) {
			my_global_temperature_struct.zpos_copy[i_7] = zpos[i_7];
		}

		for (integer iu = 0; iu < lu; iu++) {
			// Эти копии данных нужны для полного восстановления данных.
			my_union[iu].t.inx_copy = my_union[iu].inx;
			my_union[iu].t.iny_copy = my_union[iu].iny;
			my_union[iu].t.inz_copy = my_union[iu].inz;
			my_union[iu].t.operatingtemperature_copy = operatingtemperature;
			my_union[iu].t.xpos_copy = new doublereal[my_union[iu].inx + 1];
			my_union[iu].t.ypos_copy = new doublereal[my_union[iu].iny + 1];
			my_union[iu].t.zpos_copy = new doublereal[my_union[iu].inz + 1];
			// Данная информация нужна для экономии оперативной памяти,
			// некоторые данные будут выгружены из озу а потом восстановлены 
			// путём вычисления.
#pragma omp parallel for
			for (integer i_7 = 0; i_7 < my_union[iu].inx + 1; i_7++) {
				my_union[iu].t.xpos_copy[i_7] = my_union[iu].xpos[i_7];
			}
#pragma omp parallel for
			for (integer i_7 = 0; i_7 < my_union[iu].iny + 1; i_7++) {
				my_union[iu].t.ypos_copy[i_7] = my_union[iu].ypos[i_7];
			}
#pragma omp parallel for
			for (integer i_7 = 0; i_7 < my_union[iu].inz + 1; i_7++) {
				my_union[iu].t.zpos_copy[i_7] = my_union[iu].zpos[i_7];
			}
		}

		// Освобождение оперативной памяти из под octree дерева.
		if (b_on_adaptive_local_refinement_mesh) {
			printf("free octree start...\n");
			//system("PAUSE");
			//system("PAUSE");
			integer maxelm_loc = (inx + 1) * (iny + 1) * (inz + 1);
			free_octree(oc_global, maxelm_loc);
			delete[] my_ALICE_STACK;
			top_ALICE_STACK = 0;
			printf("free octree end...\n");
			//system("PAUSE");
			//system("PAUSE");
		}

		if (0) {
			xyplot(f, flow_interior, my_global_temperature_struct);
			printf("after load temper and flow. OK.\n");
			//system("PAUSE"); // debug avtosave
			system("pause");
		}

		// На АЛИС сетке полинейный метод не является работоспособным.
		if (!b_on_adaptive_local_refinement_mesh) {
			// создаёт информацию о сеточных линиях для полилинейного метода LR:
			if (2 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
				// Lr1sk algorithm
				constr_line(f, flow_interior);  // для гидродинамики
			}
			my_global_temperature_struct.rootBT = nullptr;
			my_global_temperature_struct.rootSN = nullptr;
			my_global_temperature_struct.rootWE = nullptr;
			if (2 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
				// Lr1sk algorithm
				constr_line_temp(my_global_temperature_struct, b, lb); // для теплопроводности
				printf("LR preprocessing finish...\n");
			}
		}

		// глобальная память для алгебраического многосеточного метода.

		amgGM.a = nullptr;
		amgGM.f = nullptr;
		amgGM.ia = nullptr;
		amgGM.ig = nullptr;
		amgGM.ja = nullptr;
		amgGM.u = nullptr;
		amgGM.nda = -1;
		amgGM.ndf = -1;
		amgGM.ndia = -1;
		amgGM.ndig = -1;
		amgGM.ndja = -1;
		amgGM.ndu = -1;


		//PARDATA nd;
		nd.ncore = 2; // два ядра.
		// По умолчанию все разбиения неактивны.
		nd.b0.active = false;
		nd.b00.active = false;
		nd.b01.active = false;
		nd.b000.active = false;
		nd.b001.active = false;
		nd.b010.active = false;
		nd.b011.active = false;
		if (0 && (1 == flow_interior)) {
			calc_front(f, f[0], my_global_temperature_struct, flow_interior, ls, lw, w, nd, b, lb, s);
			// разделение выполнено !
			printf("separator compleate...\n");
			//system("PAUSE");
		}



		my_global_temperature_struct.free_temper_level1 = false; // чистая теплопроводность освобождение памяти необходимой для сборки матрицы после успешной сборки.
		my_global_temperature_struct.free_temper_level2 = false; // освобождение памяти под хранение матрицы при перезаписи её в SIMPLESPARSE формат.	

		printf("construction of all structures...\n");
		printf("mesh check start...\n");
		const doublereal d_zalipanie = 1.0e-23;// метров
#if doubleintprecision == 1
		doublereal minimum_gap = 1.0e60;
		for (integer i = 0; i < inx; i++) {
			if ((xpos[i + 1] - xpos[i]) < minimum_gap) minimum_gap = (xpos[i + 1] - xpos[i]);
			if ((xpos[i + 1] - xpos[i]) < d_zalipanie) {
				//printf("error: zalipanie po X: xpos[%lld]=%e xpos[%lld]=%e inx=%lld\n", i, xpos[i], i + 1, xpos[i + 1], inx);
				std::cout << "error: zalipanie po X: xpos[" << i << "]=" << xpos[i] << " xpos[" << i + 1 << "]=" << xpos[i + 1] << " inx=" << inx << std::endl;
				std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
				system("pause");
			}
		}
		std::cout << "minimum gap X=" << minimum_gap << std::endl;
		minimum_gap = 1.0e60;
		for (integer i = 0; i < iny; i++) {
			if ((ypos[i + 1] - ypos[i]) < minimum_gap) minimum_gap = (ypos[i + 1] - ypos[i]);
			if ((ypos[i + 1] - ypos[i]) < d_zalipanie) {
				//printf("error: zalipanie po Y: ypos[%lld]=%e ypos[%lld]=%e iny=%lld\n", i, ypos[i], i + 1, ypos[i + 1], iny);
				std::cout << "error: zalipanie po Y: ypos[" << i << "]=" << ypos[i] << " ypos[" << i + 1 << "]=" << ypos[i + 1] << " iny=" << iny << std::endl;
				std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
				system("pause");
			}
		}
		std::cout << "minimum gap Y=" << minimum_gap << std::endl;
		minimum_gap = 1.0e60;
		for (integer i = 0; i < inz; i++) {
			if ((zpos[i + 1] - zpos[i]) < minimum_gap) minimum_gap = (zpos[i + 1] - zpos[i]);
			if ((zpos[i + 1] - zpos[i]) < d_zalipanie) {
				//printf("error: zalipanie po Z: zpos[%lld]=%e zpos[%lld]=%e inz=%lld\n", i, zpos[i], i + 1, zpos[i + 1], inz);
				std::cout << "error: zalipanie po Z: zpos[" << i << "]=" << zpos[i] << " zpos[" << i + 1 << "]=" << zpos[i + 1] << " inz=" << inz << std::endl;
				std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
				system("pause");
			}
		}
		std::cout << "minimum gap Z=" << minimum_gap << std::endl;
		for (integer iP = 0; iP < my_global_temperature_struct.maxelm; iP++) {
			if ((my_global_temperature_struct.nvtx[0][iP] == 0) || (my_global_temperature_struct.nvtx[1][iP] == 0) || (my_global_temperature_struct.nvtx[2][iP] == 0) || (my_global_temperature_struct.nvtx[3][iP] == 0) || (my_global_temperature_struct.nvtx[4][iP] == 0) || (my_global_temperature_struct.nvtx[5][iP] == 0) || (my_global_temperature_struct.nvtx[6][iP] == 0) || (my_global_temperature_struct.nvtx[7][iP] == 0)) {
				printf("nvtx[%lld]: %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, my_global_temperature_struct.nvtx[0][iP] - 1, my_global_temperature_struct.nvtx[1][iP] - 1, my_global_temperature_struct.nvtx[2][iP] - 1, my_global_temperature_struct.nvtx[3][iP] - 1, my_global_temperature_struct.nvtx[4][iP] - 1, my_global_temperature_struct.nvtx[5][iP] - 1, my_global_temperature_struct.nvtx[6][iP] - 1, my_global_temperature_struct.nvtx[7][iP] - 1);
			}
		}
#else
		doublereal minimum_gap = 1.0e60;
		for (integer i = 0; i < inx; i++) {
			if ((xpos[i + 1] - xpos[i]) < minimum_gap) minimum_gap = (xpos[i + 1] - xpos[i]);
			if ((xpos[i + 1] - xpos[i]) < d_zalipanie) {
				//printf("error: zalipanie po X: xpos[%d]=%e xpos[%d]=%e inx=%d\n", i, xpos[i], i + 1, xpos[i + 1], inx);
				std::cout << "error: zalipanie po X: xpos[" << i << "]=" << xpos[i] << " xpos[" << i + 1 << "]=" << xpos[i + 1] << " inx=" << inx << std::endl;
				std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
				system("pause");
			}
		}
		std::cout << "minimum gap X=" << minimum_gap << std::endl;
		minimum_gap = 1.0e60;
		for (integer i = 0; i < iny; i++) {
			if ((ypos[i + 1] - ypos[i]) < minimum_gap) minimum_gap = (ypos[i + 1] - ypos[i]);
			if ((ypos[i + 1] - ypos[i]) < d_zalipanie) {
				//	printf("error: zalipanie po Y: ypos[%d]=%e ypos[%d]=%e iny=%d\n", i, ypos[i], i + 1, ypos[i + 1], iny);
				std::cout << "error: zalipanie po Y: ypos[" << i << "]=" << ypos[i] << " ypos[" << i + 1 << "]=" << ypos[i + 1] << " iny=" << iny << std::endl;
				std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
				system("pause");
			}
		}
		std::cout << "minimum gap Y=" << minimum_gap << std::endl;
		minimum_gap = 1.0e60;
		for (integer i = 0; i < inz; i++) {
			if ((zpos[i + 1] - zpos[i]) < minimum_gap) minimum_gap = (zpos[i + 1] - zpos[i]);
			if ((zpos[i + 1] - zpos[i]) < d_zalipanie) {
				//		printf("error: zalipanie po Z: zpos[%d]=%e zpos[%d]=%e inz=%d\n", i, zpos[i], i + 1, zpos[i + 1], inz);
				std::cout << "error: zalipanie po Z: zpos[" << i << "]=" << zpos[i] << " zpos[" << i + 1 << "]=" << zpos[i + 1] << " inz=" << inz << std::endl;
				std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
				system("pause");
			}
		}
		std::cout << "minimum gap Z=" << minimum_gap << std::endl;
		for (integer iP = 0; iP < my_global_temperature_struct.maxelm; iP++) {
			if ((my_global_temperature_struct.nvtx[0][iP] == 0) || (my_global_temperature_struct.nvtx[1][iP] == 0) || (my_global_temperature_struct.nvtx[2][iP] == 0) || (my_global_temperature_struct.nvtx[3][iP] == 0) || (my_global_temperature_struct.nvtx[4][iP] == 0) || (my_global_temperature_struct.nvtx[5][iP] == 0) || (my_global_temperature_struct.nvtx[6][iP] == 0) || (my_global_temperature_struct.nvtx[7][iP] == 0)) {
				printf("nvtx[%d]: %d %d %d %d %d %d %d %d \n", iP, my_global_temperature_struct.nvtx[0][iP] - 1, my_global_temperature_struct.nvtx[1][iP] - 1, my_global_temperature_struct.nvtx[2][iP] - 1, my_global_temperature_struct.nvtx[3][iP] - 1, my_global_temperature_struct.nvtx[4][iP] - 1, my_global_temperature_struct.nvtx[5][iP] - 1, my_global_temperature_struct.nvtx[6][iP] - 1, my_global_temperature_struct.nvtx[7][iP] - 1);
			}
		}
#endif


		// Не имеет смысла решать СЛАУ с точностью превышающей порядок аппроксимации.
		// Порядок аппроксимации O(h!2) второго порядка определяется размерами контрольного объёма,
		// т.е. зависит от подробности расчётной сетки.
		for (integer i = 0; i < flow_interior; i++) {
#if doubleintprecision == 1
			printf("FLUID %lld\n", i);
#else
			printf("FLUID %d\n", i);
#endif

			// точность с которой аппроксимировано уравнение для поправки давления.
			f[i].resICCG = rterminate_residual_ICCG_Oh2(f[i]); // O(h!2)
			//printf("residual O(h!2) is equal=%e\n", f[i].resICCG);
			std::cout << "residual O(h!2) is equal=" << f[i].resICCG << std::endl;
			f[i].resLR1sk = rterminate_residual_LR1sk_Oh3(f[i]); // O(h!3)
			//printf("residual O(h!3) is equal=%e\n", f[i].resLR1sk);
			std::cout << "residual O(h!3) is equal=" << f[i].resLR1sk << std::endl;
		}
		printf("TEMPERATURE\n");
		my_global_temperature_struct.resLR1sk = rterminate_residual_LR1sk_temp_Oh3(my_global_temperature_struct); // O(h!3)		
		//printf("temp residual O(h!3) is equal=%e\n", t.resLR1sk);
		std::cout << "temp residual O(h!3) is equal=" << my_global_temperature_struct.resLR1sk << std::endl;
		printf("mesh check.\n");
		if (bwait) {
			//system("PAUSE");
			system("pause");
		}

		sourse2Dproblem = new bool[my_global_temperature_struct.maxbound];
		conductivity2Dinsource = new doublereal[my_global_temperature_struct.maxbound];
		// Нижняя релаксация источникового члена при радиационных потоках.
		bsource_term_radiation_for_relax = new doublereal[my_global_temperature_struct.maxelm];
		for (integer i_init = 0; i_init < my_global_temperature_struct.maxelm; i_init++) bsource_term_radiation_for_relax[i_init] = 0.0;
		b_buffer_correct_source = new doublereal[my_global_temperature_struct.maxelm];


		// невязка continity будет измеряться по отношению к уровню 1e0.
		doublereal* continity_start = nullptr;
		continity_start = new doublereal[flow_interior];
		if (continity_start == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for continity start in main...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}
		for (integer i = 0; i < flow_interior; i++) continity_start[i] = 1.0;

		integer* inumber_iteration_SIMPLE = nullptr;
		inumber_iteration_SIMPLE = new integer[flow_interior];
		if (nullptr == inumber_iteration_SIMPLE) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for inumber_iteration_SIMPLE in main...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}
		for (integer i = 0; i < flow_interior; i++) inumber_iteration_SIMPLE[i] = 0; // начальная итерация алгоритма SIMPLE для каждой FLUID зоны.

		// считывание состояния расчёта из файла для возобновления расчёта
		bool breadOk = false;
		avtoreadvalue(f, my_global_temperature_struct, flow_interior, inumber_iteration_SIMPLE, continity_start, breadOk, b, lb, s, ls, w, lw);
		// Если считывание прошло неуспешно, то breadOk==false и это значит что счёт начнётся заново со значений заданных при инициализации.

		if (b_on_adaptive_local_refinement_mesh) {
			// Инвариант корректности АЛИС сетки.
			printf("the invariant correctness...\n");
			ANES_ALICE_CORRECT(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx);
		}

		// экспорт результата вычисления в программу tecplot360:
		// можно использовать как проверку построенной сетки.
		if (0) {
			exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
			printf("read values. OK.\n");
			//system("PAUSE"); // debug avtosave
			system("pause");
		}



		/*for (integer i=0; i<lw; i++) {
		printf("%e  \n",w[i].Tamb);
		}
		//exporttecplotxy360T_3D(t.maxelm, t.ncell, t.nvtx, t.nvtxcell, t.pa, t.potent);
		exporttecplotxy360T_3D_part2(t.maxelm, t.potent);
		system("PAUSE"); // debug
		*/


		// 29.01.2017
		// if (1 && steady_or_unsteady_global_determinant == MESHER_ONLY)  
		// То мы просто вызываем мешер не вызывая солвера.
		if (1 && ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY)||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T))) {
			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

			// Экспортирует внешнюю поверхность геометрии пользователя в .stl формате.
			// 09.09.2019.
			// Закоментировал 27.05.2020
			//if (steady_or_unsteady_global_determinant == MESHER_ONLY) {
				//export_User_Geom_in_STL_format(my_global_temperature_struct);
			//}

			// Вычисление массы модели.
			massa_cabinet(my_global_temperature_struct, f,  flow_interior,
				b, lb, operatingtemperature,
				matlist);

			calculation_end_time = clock(); // момент окончания счёта.
			calculation_seach_time = calculation_end_time - calculation_start_time;
			unsigned int im = 0, is = 0, ims = 0;
			im = (unsigned int)(calculation_seach_time / 60000); // минуты
			is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
			ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

			// 24.06.2020
			// Графовый метод решения уравнения теплопередачи.
			if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T)) {
				calculate_Network_T(my_global_temperature_struct, b, lb, w, lw, ls, matlist);
			}

			printf("time export to tecplot360 is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);
			//system("pause");
			//exit(1);

			// 1 - solver/solid_static/
			bool bMechanical = false;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
			// Печатает расход жидкости через выходную границу потока.
			// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
			// проходящий через выходную границу потока. 28,10,2019
			//report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);

			// Экспорт сетки в tecplot 360.
			if (1) {
				if (!b_on_adaptive_local_refinement_mesh) {
					if (0 == lu) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
					}
					else {
						//exporttecplotxy360T_3D_part2_assembles(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0, lu, my_union);
						// сетка без разрывов.
						exporttecplot_assembles_mesh(my_global_temperature_struct, lu, my_union);
					}
				}
				else {
					// Экспорт в программу tecplot температуры.
					//С АЛИС сетки.
					ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
				}
			}


		}

		//char ch_EXPORT_ALICE_ONLY = 'y';

		// steady Temperature Finite Volume Method
		if (1 && (steady_or_unsteady_global_determinant== PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE) &&
			(1==eqin.itemper)) {

			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

#pragma omp parallel for
			for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) 
				my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // инициализация.

			// Решаем только теплопередачу в твёрдом теле.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			if (1 == lw) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка:
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
				// Мы инициализируем случайной величиной чтобы избавится от неопределённости при первой сборке СЛАУ.
				//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
			}

			// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
			bsolid_static_only = true;
			bool bcleantemp = false;
			// Условие уже проверялось выше.
			//if (1 == eqin.itemper) 
			{
				bcleantemp = true;
				for (integer i = 0; i < flow_interior; i++) {
					if (1 == eqin.fluidinfo[i].iflow) bcleantemp = false;
				}
				// если bcleantemp==true то мы решаем задачу чистой теплопередачи без учёта конвекции.
			}

			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой теплопроводности в трёхмерной области. 
				printf("solution of pure heat...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//system("PAUSE");
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
				}
				else {
					// Загрузка распределения начальной скорости.
					
					FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
					int err_inicialization_data = 0;
					fp_inicialization_data = fopen64("load.txt", "r");
					if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
					errno_t err_inicialization_data = 0;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif
					if (0 == err_inicialization_data) {
						// Открытие удачно и файл присутствует.
						if (f[0].maxelm > 0) {
							// Только в случае наличия жидких ячеек мы решаем задачу с конвекцией.
							bmyconvective = true;
						}
						fclose(fp_inicialization_data);
					}
					//printf("maxelm flow = %lld\n",f[0].maxelm);
					//system("PAUSE");
				}

				// if (flow_interior>0) bmyconvective=true;
				// массив отладочной информации,
				// конкретно для проверки подхода Рхи-Чоу 1983
				doublereal** rhie_chow = nullptr;
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

				{
					integer* color = nullptr;
					integer dist_max = 3;
					calculate_color_for_temperature(color, my_global_temperature_struct,inx,xpos);

					// если flow_interior == 0 то f[0] просто формальный параметр  
					solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
						rhie_chow,
						b, lb, s, ls, w, lw,
						dbeta, flow_interior,
						bmyconvective, nullptr, 0.001, 0.001,
						false,
						matlist, 0,
						bprintmessage,
						gtdps, ltdp, 1.0, // последний параметр равный 1.0 означает что мощность подаётся.
						m, nullptr, // скорость с предыдущего временного слоя. 
						nullptr,
						lu, my_union, color, dist_max); // массовый поток через границу с предыдущего временного слоя.
					delete[] color;
				}

				while ((bglobal_restart_06_10_2018)) {

					// Расчётная сетка была перестроена глобально. Требуется перевыделить память.
					if (my_global_temperature_struct.potent != nullptr) {
						delete[] my_global_temperature_struct.potent;
						my_global_temperature_struct.potent = nullptr;
						my_global_temperature_struct.potent = new doublereal[my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound];
						for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // инициализация.
						if (bsource_term_radiation_for_relax != nullptr) {
							delete[] bsource_term_radiation_for_relax;
							bsource_term_radiation_for_relax = nullptr;
							bsource_term_radiation_for_relax = new doublereal[my_global_temperature_struct.maxelm];
						}
						if (sourse2Dproblem != nullptr) {
							delete[] sourse2Dproblem;
							sourse2Dproblem = nullptr;
							sourse2Dproblem = new bool[my_global_temperature_struct.maxbound];
						}
						if (conductivity2Dinsource != nullptr) {
							delete[] conductivity2Dinsource;
							conductivity2Dinsource = nullptr;
							conductivity2Dinsource = new doublereal[my_global_temperature_struct.maxbound];
						}

						if (rthdsd_no_radiosity_patch != nullptr) {
							delete[]	rthdsd_no_radiosity_patch;
							rthdsd_no_radiosity_patch = nullptr;
						}

						if (b_buffer_correct_source != nullptr) {
							delete[] b_buffer_correct_source;
							b_buffer_correct_source = nullptr;
							b_buffer_correct_source = new doublereal[my_global_temperature_struct.maxelm];
						}

						if (my_global_temperature_struct.slau != nullptr) {
							delete[] my_global_temperature_struct.slau;
							my_global_temperature_struct.slau = nullptr;
							my_global_temperature_struct.slau = new equation3D[my_global_temperature_struct.maxelm]; // коэффициенты матрицы СЛАУ для внутренних КО.
							if (my_global_temperature_struct.slau == nullptr) {
								// недостаточно памяти на данном оборудовании.
								printf("Problem: not enough memory on your equipment for slau temperature constr struct...\n");
								printf("Please any key to exit...\n");
								//system("PAUSE");
								system("pause");
								exit(1);
							}
						}

						if (my_global_temperature_struct.slau_bon != nullptr) {
							delete[] my_global_temperature_struct.slau_bon;
							my_global_temperature_struct.slau_bon = nullptr;
							my_global_temperature_struct.slau_bon = new equation3D_bon[my_global_temperature_struct.maxbound]; // коэффициенты матрицы СЛАУ для граничных КО
							if (my_global_temperature_struct.slau_bon == nullptr) {
								// недостаточно памяти на данном оборудовании.
								printf("Problem: not enough memory on your equipment for slau boundary temperature constr struct...\n");
								printf("Please any key to exit...\n");
								//system("PAUSE");
								system("pause");
								exit(1);
							}
						}
					}
					//bglobal_restart_06_10_2018 = false;

					integer* color = nullptr;
					integer dist_max = 3;
					calculate_color_for_temperature(color, my_global_temperature_struct,inx,xpos);


					// если flow_interior == 0 то f[0] просто формальный параметр  
					solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
						rhie_chow,
						b, lb, s, ls, w, lw,
						dbeta, flow_interior,
						bmyconvective, nullptr, 0.001, 0.001,
						false,
						matlist, 0,
						bprintmessage,
						gtdps, ltdp, 1.0, m,
						nullptr, // скорость с предыдущего временного слоя. 
						nullptr,
						lu, my_union, color, dist_max); // массовый поток через границу с предыдущего временного слоя.
									   // последний параметр равный 1.0 означает что мощность подаётся.
					delete[] color;

				}

				// Вычисляет среднюю температуру в жидких ячейках.
				doublereal Tavg_fluid = 0.0;
				doublereal ic_avg_Temp_fluid = 0.0;

#pragma omp parallel for reduction(+ : Tavg_fluid, ic_avg_Temp_fluid)
				for (integer i_7 = 0; i_7 < f[0].maxelm; i_7++) {
					if ((f[0].ptr[i_7] >= 0) && (f[0].ptr[i_7] < my_global_temperature_struct.maxelm)) {
						Tavg_fluid += my_global_temperature_struct.potent[f[0].ptr[i_7]];
						ic_avg_Temp_fluid += 1.0;
					}
				}
				if (ic_avg_Temp_fluid > 1.0e-30) {
					Tavg_fluid /= ic_avg_Temp_fluid;
					printf("average fluid temperature is %e\n", Tavg_fluid);
				}
				else {
					//printf("no fluid cell\n");
				}

				// Вычисление массы модели.
				massa_cabinet(my_global_temperature_struct, f, flow_interior,
					b, lb, operatingtemperature,
					matlist);

				// 10.10.2017
				// Построение двумерного графика вдоль структуры.
				xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);
				//printf("graphics writing sucseful\n");
				//system("PAUSE");

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {

						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
					}
					else {
						if (b_on_adaptive_local_refinement_mesh) {
							calculation_main_end_time = clock(); // момент окончания счёта.
							unsigned int calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;

							// Общее время вычисления до возникновения критического сообщения.
							int im = 0, is = 0, ims = 0;
							im = (int)(calculation_main_seach_time / 60000); // минуты
							is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // секунды
							ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

							printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

							//printf("Would you like to save the result on the ALICE grid ? y/n\n");
							//ch_EXPORT_ALICE_ONLY = getchar(); // Здесь именно getchar();
							//ch_EXPORT_ALICE_ONLY = 'y';
						}
						if (ch_EXPORT_ALICE_ONLY == 'y') {
							// Экспорт в программу tecplot температуры.
							//С АЛИС сетки.
						    ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
						}

					}
				}

			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
			// Вычисление значения максимальной температуры внутри расчётной области и на её границах:
			//for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
			// 23 декабря 2015
			// На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
			// физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

			FILE* fp = NULL;
			
#ifdef MINGW_COMPILLER
			int err1 = 0;
			fp = fopen64("report.txt", "w");
			if (fp == NULL) err1 = 1;
#else
			errno_t err1 = 0;
			err1 = fopen_s(&fp, "report.txt", "w");
#endif
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//system("PAUSE");
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			bool bMechanical = false;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
			// Печатает расход жидкости через выходную границу потока.
			// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
			// проходящий через выходную границу потока. 28.10.2019
			report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);

			
			if (ch_EXPORT_ALICE_ONLY != 'y') {

				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

				// На реальной модели большой размерности (11.4М ячеек на АЛИС Coarse )
				// преобразует неприемлемо долго. 

				// 25.11.2017
				// 1. Получили данные о температуре и сетке на АЛИС.
				// 2. Сохранили их в оперативной памяти.
				// 3. Освободили память.
				// 4. Построили обычную декартовую прямоугольную сетку.
				// 5. Перенесли данные о температуре с АЛИС на обычную декартовую прямоугольную сетку.
				// 6. Визуализировали температуру и построенные на структурированной сетке тепловые потоки.
				// 7.1 Tecplot идеально визуализирует то что есть на структурированной сетке и в сечении и в объёме. 
				// 7.2 Визуализация на АЛИС сетке в объеме даже идеально точно найденного поля не удовлетворительна (глюк tecplotа).


				if (b_on_adaptive_local_refinement_mesh) {
					// 1. Получение x,y,z,T,nvtx, m_sizeT, m_size_nvtx.
					doublereal* x_buf = nullptr;
					doublereal* y_buf = nullptr;
					doublereal* z_buf = nullptr;
					doublereal* t_buf = nullptr;
					integer** nvtx_buf = nullptr;
					integer m_sizeT = 0, m_size_nvtx = 0;

					ANES_tecplot360_export_temperature_preobrazovatel(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx, operatingtemperature);

					// 2. Освобождение памяти.
					// Освобождение оперативной памяти.
					if (my_global_temperature_struct.xpos_copy != nullptr) {
						delete[] my_global_temperature_struct.xpos_copy;
						my_global_temperature_struct.xpos_copy = nullptr;
					}
					if (my_global_temperature_struct.ypos_copy != nullptr) {
						delete[] my_global_temperature_struct.ypos_copy;
						my_global_temperature_struct.ypos_copy = nullptr;
					}
					if (my_global_temperature_struct.zpos_copy != nullptr) {
						delete[] my_global_temperature_struct.zpos_copy;
						my_global_temperature_struct.zpos_copy = nullptr;
					}


					if (bsource_term_radiation_for_relax != nullptr) {
						delete[] bsource_term_radiation_for_relax; // Релаксация источниковых членов радиационных потоков.
						bsource_term_radiation_for_relax = nullptr;
					}
					if (b_buffer_correct_source != nullptr) {
						delete[] b_buffer_correct_source;
						b_buffer_correct_source = nullptr;
					}

					if (rthdsd_no_radiosity_patch != nullptr) {
						delete[] rthdsd_no_radiosity_patch;
						rthdsd_no_radiosity_patch = nullptr;
					}
					

					// Быстрая обработка нелинейных граничных условий в Румба 0.14 решателе.
					if (qnbc != nullptr) {
						delete[] qnbc;
						qnbc = nullptr;
						iadd_qnbc_maxelm = 0;
					}

					// Нужно освободить оперативную память из под всех структур данных:
					free_level1_temp(my_global_temperature_struct);
					free_level2_temp(my_global_temperature_struct); // освобождение памяти из под матриц.
										 // Освобождает память для LR начало.
					if (my_global_temperature_struct.rootWE != nullptr) {
						free_root(my_global_temperature_struct.rootWE, my_global_temperature_struct.iWE);
					}
					if (my_global_temperature_struct.rootSN != nullptr) {
						free_root(my_global_temperature_struct.rootSN, my_global_temperature_struct.iSN);
					}
					if (my_global_temperature_struct.rootBT != nullptr) {
						free_root(my_global_temperature_struct.rootBT, my_global_temperature_struct.iBT);
					}
					if (my_global_temperature_struct.rootWE != nullptr) {
						delete[] my_global_temperature_struct.rootWE;
						my_global_temperature_struct.rootWE = nullptr;
					}
					if (my_global_temperature_struct.rootSN != nullptr) {
						delete[] my_global_temperature_struct.rootSN;
						my_global_temperature_struct.rootSN = nullptr;
					}
					if (my_global_temperature_struct.rootBT != nullptr) {
						delete[] my_global_temperature_struct.rootBT;
						my_global_temperature_struct.rootBT = nullptr;
					}

					if (bvery_big_memory) {
						if (my_global_temperature_struct.database.x != nullptr) {
							free(my_global_temperature_struct.database.x);
						}
						if (my_global_temperature_struct.database.y != nullptr) {
							free(my_global_temperature_struct.database.y);
						}
						if (my_global_temperature_struct.database.z != nullptr) {
							free(my_global_temperature_struct.database.z);
						}
						if (my_global_temperature_struct.database.nvtxcell != nullptr) {
							for (integer i = 0; i <= 7; i++) {
								delete[] my_global_temperature_struct.database.nvtxcell[i];
							}
							delete[] my_global_temperature_struct.database.nvtxcell;
						}
						if (my_global_temperature_struct.database.ptr != nullptr) {
							if (my_global_temperature_struct.database.ptr[0] != nullptr) {
								delete[] my_global_temperature_struct.database.ptr[0];
							}
							if (my_global_temperature_struct.database.ptr[1] != nullptr) {
								delete[] my_global_temperature_struct.database.ptr[1];
							}
							delete[] my_global_temperature_struct.database.ptr;
						}
					}

					// Освобождение памяти для LR конец.
					free_level1_flow(f, flow_interior);
					free_level2_flow(f, flow_interior); // освобождение памяти из под матриц.

					if (sourse2Dproblem != nullptr) {
						delete[] sourse2Dproblem;
						sourse2Dproblem = nullptr;
					}
					if (conductivity2Dinsource != nullptr) {
						delete[] conductivity2Dinsource;
						conductivity2Dinsource = nullptr;
					}

					if (x_jacoby_buffer != nullptr) {
						// 30 октября 2016. 
						// В seidelsor2 сделан переключатель на метод нижней релаксации К.Г. Якоби.
						// Освобождение памяти из под Jacobi buffer.
						delete[] x_jacoby_buffer;
					}


					
					// Освобождение общей памяти в ILU буффере.
					//if (milu_gl_buffer.alu_copy != nullptr) delete[] milu_gl_buffer.alu_copy;
					//if (milu_gl_buffer.jlu_copy != nullptr) delete[] milu_gl_buffer.jlu_copy;
					//if (milu_gl_buffer.ju_copy != nullptr) delete[] milu_gl_buffer.ju_copy;
					//milu_gl_buffer.alu_copy = nullptr;
					//milu_gl_buffer.jlu_copy = nullptr;
					//milu_gl_buffer.ju_copy = nullptr;
					
					flow_interior = 0;

					// 3. Построение обычной сетки.

					b_on_adaptive_local_refinement_mesh = false;
					iCabinetMarker = 0;
					load_TEMPER_and_FLOW(my_global_temperature_struct, f, inx, iny, inz, xpos, ypos, zpos, flow_interior,
						b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
						dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iCabinetMarker);

					my_global_temperature_struct.operatingtemperature = operatingtemperature;


					// Эти копии данных нужны для полного восстановления данных.
					my_global_temperature_struct.inx_copy = inx;
					my_global_temperature_struct.iny_copy = iny;
					my_global_temperature_struct.inz_copy = inz;
					my_global_temperature_struct.operatingtemperature_copy = operatingtemperature;
					my_global_temperature_struct.xpos_copy = new doublereal[inx + 1];
					my_global_temperature_struct.ypos_copy = new doublereal[iny + 1];
					my_global_temperature_struct.zpos_copy = new doublereal[inz + 1];
					// Данная информация нужна для экономии оперативной памяти,
					// некоторые данные будут выгружены из озу а потом восстановлены 
					// путём вычиcления.
					for (integer i_7 = 0; i_7 < inx + 1; i_7++) {
						my_global_temperature_struct.xpos_copy[i_7] = xpos[i_7];
					}
					for (integer i_7 = 0; i_7 < iny + 1; i_7++) {
						my_global_temperature_struct.ypos_copy[i_7] = ypos[i_7];
					}
					for (integer i_7 = 0; i_7 < inz + 1; i_7++) {
						my_global_temperature_struct.zpos_copy[i_7] = zpos[i_7];
					}

					my_global_temperature_struct.free_temper_level1 = false; // чистая теплопроводность освобождение памяти необходимой для сборки матрицы после успешной сборки.
					my_global_temperature_struct.free_temper_level2 = false; // освобождение памяти под хранение матрицы при перезаписи её в SIMPLESPARSE формат.	


					// 4. Интерполяция для температуры.
					ALICE_2_Structural(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx, my_global_temperature_struct.operatingtemperature_copy);


					if (x_buf != nullptr) {
						delete[] x_buf;
						x_buf = nullptr;
					}
					if (y_buf != nullptr) {
						delete[] y_buf;
						y_buf = nullptr;
					}
					if (z_buf != nullptr) {
						delete[] z_buf;
						z_buf = nullptr;
					}
					if (t_buf != nullptr) {
						delete[] t_buf;
						t_buf = nullptr;
					}
					if (nvtx_buf != nullptr) {
						for (integer i_1 = 0; i_1 < 8; i_1++) {
							if (nvtx_buf[i_1] != nullptr) {
								delete[] nvtx_buf[i_1];
								nvtx_buf[i_1] = nullptr;
							}
						}
						delete[] nvtx_buf;
						nvtx_buf = nullptr;
					}
					m_sizeT = 0, m_size_nvtx = 0;
					// 5. Обычный экспорт в tecplot.
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
			}
			else {
				// Объёмная визуализация на АЛИС очень плохая даже с точным полем, поэтому
				// для качественной объёмной визуализации нужен переход на структурированную сетку.
				// В сечении на АЛИС сетке все хорошо по tecplotу.
				// Точность на АЛИС очень плоха. Можно допиться хорошей картинки для поля температур, но
				// плотности тепловых потоков никогда не будут найдены (представлены точно) только в случае
				// патологического случая дико мелкой сетки когда АЛИС очень подробная (модели большой размерности).

				// Экспорт в программу tecplot температуры.
				// С АЛИС сетки.
				//ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent,t,0,b,lb);
			}
			


		}

		// steady Temperature Finite Element Method
		if (1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE) && (eqin.itemper == 2)) {

			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

#pragma omp parallel for
			for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) {
				my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // инициализация.
			}

			// Решаем только Static Structural.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка:
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
				// Мы инициализируем случайной величиной чтобы избавится от неопределённости при первой сборке СЛАУ.
				//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
			}

			// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
			bsolid_static_only = true;
			bool bcleantemp = false;
			// Мы тоже решаем задачу чистой теплопроводности только
			// вторым температурным солвером.
			//if (2 == eqin.itemper)
			{
				bcleantemp = true;
				for (integer i = 0; i < flow_interior; i++) {
					if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
				}
				// если bcleantemp==true то мы решаем задачу чистой теплопередачи без учёта конвекции.
			}
			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой теплопроводности в трёхмерной области. 
				printf("solution of pure Static Structural...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//system("PAUSE");
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
				}
				else {
					// Загрузка распределения начальной скорости.
					
					FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
					int err_inicialization_data = 0;
					fp_inicialization_data = fopen64("load.txt", "r");
					if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
					errno_t err_inicialization_data = 0;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

					if (err_inicialization_data == 0) {
						// Открытие удачно и файл присутствует.
						if (f[0].maxelm > 0) {
							bmyconvective = true;
						}
						fclose(fp_inicialization_data);
					}
				}

				// if (flow_interior>0) bmyconvective=true;
				// массив отладочной информации,
				// конкретно для проверки подхода Рхи-Чоу 1983
				doublereal** rhie_chow = nullptr;
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

				bPhysics_stop = false;

				// Температура 19.05.2018

				doublereal* lstub = nullptr;
				integer maxelm_global_ret = 0;
				doublereal* t_for_Mechanical = nullptr;
				solve_Thermal(my_global_temperature_struct, f, matlist, w, lw, lu, b, lb, ls, m, false, 
					operatingtemperature, false, 0.0, lstub, lstub,
					maxelm_global_ret, 1.0, bAVLrealesation, t_for_Mechanical);


				/*
				// если flow_interior == 0 то f[0] просто формальный параметр
				solve_nonlinear_temp(f[0], f, t,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				bmyconvective, nullptr, 0.001, 0.001,
				false,
				matlist, 0,
				bprintmessage,
				gtdps, ltdp, 1.0, m,
				nullptr, // скорость с предыдущего временного слоя.
				nullptr); // массовый поток через границу с предыдущего временного слоя.
				// последний параметр равный 1.0 означает что мощность подаётся.
				*/
				// Вычисление массы модели.
				massa_cabinet(my_global_temperature_struct, f,  flow_interior,
					b, lb, operatingtemperature,
					matlist);


				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
					}
					else {
						// Экспорт в программу tecplot температуры.
						//С АЛИС сетки.
						ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
					}
				}

			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
											 // Вычисление значения максимальной температуры внутри расчётной области и на её границах:
											 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
											 // 23 декабря 2015
											 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
											 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

			doublereal totaldeform_max = -1.0e+30;
			if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
			{
				for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
			}

			FILE* fp = NULL;
			
#ifdef MINGW_COMPILLER
			int err1 = 0;
			fp = fopen64("report.txt", "w");
			if (fp == NULL) err1 = 1;
#else
			errno_t err1 = 0;
			err1 = fopen_s(&fp, "report.txt", "w");
#endif
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//system("PAUSE");
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			bool bMechanical = false;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);


		}

		// steady Static Structural
		if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) {


			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

#pragma omp parallel for
			for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) {
				my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // инициализация.
			}

			// Решаем только Static Structural.
			bonly_solid_calculation = true;


			// Здесь предполагается что мы решаем стационарную задачу чистой механики по нахождению деформированного состояния.
			bsolid_static_only = true;
			bool bcleantemp = false;

			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой 
				// механики по нахождению деформированного состояния в трёхмерной области. 
				printf("solution of pure Static Structural...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//system("PAUSE");
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx * starting_speed_Vx +
					starting_speed_Vy * starting_speed_Vy +
					starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
				}
				else {
					// Загрузка распределения начальной скорости.
					
					FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
					int err_inicialization_data = 0;
					fp_inicialization_data = fopen64("load.txt", "r");
					if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
					errno_t err_inicialization_data = 0;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

					if (err_inicialization_data == 0) {
						// Открытие удачно и файл присутствует.
						if (f[0].maxelm > 0) {
							bmyconvective = true;
						}
						fclose(fp_inicialization_data);
					}
				}


				bPhysics_stop = false;
				// Вызов солвера Static Structural.
				// Погашено 19.05.2018
				doublereal* stub_Mechanical = nullptr;
				solve_Structural(my_global_temperature_struct, w, lw, false, operatingtemperature, b, lb, lu,
					false, 1.0e9, stub_Mechanical, stub_Mechanical, 1.0, matlist, stub_Mechanical);
				bPhysics_stop = true;


				// Вычисление массы модели.
				massa_cabinet(my_global_temperature_struct, f, flow_interior,
					b, lb, operatingtemperature,
					matlist);


				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
					}
					else {
						// Экспорт в программу tecplot температуры.
						//С АЛИС сетки.
						ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
					}
				}

			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
											 // Вычисление значения максимальной температуры внутри расчётной области и на её границах:
											 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
											 // 23 декабря 2015
											 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
											 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

			doublereal totaldeform_max = -1.0e+30;
			if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
			{
				for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
		    }

			FILE* fp = NULL;
			
#ifdef MINGW_COMPILLER
			int err1 = 0;
			fp = fopen64("report.txt", "w");
			if (fp == NULL) err1 = 1;
#else
			errno_t err1 = 0;
			err1 = fopen_s(&fp, "report.txt", "w");
#endif
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//system("PAUSE");
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			bool bMechanical = true;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);

		}

		// steady Static Structural and Temperature (Thermal Stress).
		if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) {

			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

			for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++)
				my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // инициализация.

			// Решаем теплопередачу, а затем Static Structural.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка:
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
				// Мы инициализируем случайной величиной чтобы избавится от неопределённости при первой сборке СЛАУ.
				//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
			}

			// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
			bsolid_static_only = true;
			bool bcleantemp = false;
			if (eqin.itemper == 1) {
				bcleantemp = true;
				for (integer i = 0; i < flow_interior; i++) {
					if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
				}
				// если bcleantemp==true то мы решаем задачу чистой теплопередачи без учёта конвекции.
			}
			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой теплопроводности в трёхмерной области. 
				printf("solution of pure Static Structural...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//system("PAUSE");
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
				}
				else {
					// Загрузка распределения начальной скорости.
					
					FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
					int err_inicialization_data = 0;
					fp_inicialization_data = fopen64("load.txt", "r");
					if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
					errno_t err_inicialization_data = 0;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

					if (err_inicialization_data == 0) {
						// Открытие удачно и файл присутствует.
						if (f[0].maxelm > 0) {
							bmyconvective = true;
						}
						fclose(fp_inicialization_data);
					}
				}

				// if (flow_interior>0) bmyconvective=true;
				// массив отладочной информации,
				// конкретно для проверки подхода Рхи-Чоу 1983
				doublereal** rhie_chow = nullptr;
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



				integer* color = nullptr;
				integer dist_max = 3;
				calculate_color_for_temperature(color, my_global_temperature_struct,inx,xpos);

				
				// если flow_interior == 0 то f[0] просто формальный параметр
				solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					bmyconvective, nullptr, 0.001, 0.001,
					false,
					matlist, 0,
					bprintmessage,
					gtdps, ltdp, 1.0, // последний параметр равный 1.0 означает что мощность подаётся.
					m, nullptr, // скорость с предыдущего временного слоя.
					nullptr,
					lu, my_union, color, dist_max); // массовый поток через границу с предыдущего временного слоя.
				  
				delete[] color;
				
				

				// Температура 19.05.2018

				doublereal* t_for_Mechanical = nullptr;
				//doublereal* lstub = nullptr;
				//integer maxelm_global_ret = 0;
				
				//solve_Thermal(my_global_temperature_struct, f, matlist, w, lw, lu, b, lb, ls, m, false,
					//operatingtemperature, false, 0.0, lstub, lstub,
					//maxelm_global_ret, 1.0, bAVLrealesation, t_for_Mechanical);
					
				t_for_Mechanical = new doublereal[my_global_temperature_struct.maxnod + 2];
				{
					// Метод линейного порядка.
					doublereal min_x = 1e60;
					doublereal min_y = 1e60;
					doublereal min_z = 1e60;
					doublereal max_x = -1e60;
					doublereal max_y = -1e60;
					doublereal max_z = -1e60;

					for (integer i = 0; i < my_global_temperature_struct.maxnod; i++) {
						if (my_global_temperature_struct.pa[i].x < min_x) {
							min_x = my_global_temperature_struct.pa[i].x;
						}
						if (my_global_temperature_struct.pa[i].y < min_y) {
							min_y = my_global_temperature_struct.pa[i].y;
						}
						if (my_global_temperature_struct.pa[i].z < min_z) {
							min_z = my_global_temperature_struct.pa[i].z;
						}
						if (my_global_temperature_struct.pa[i].x > max_x) {
							max_x = my_global_temperature_struct.pa[i].x;
						}
						if (my_global_temperature_struct.pa[i].y > max_y) {
							max_y = my_global_temperature_struct.pa[i].y;
						}
						if (my_global_temperature_struct.pa[i].z > max_z) {
							max_z = my_global_temperature_struct.pa[i].z;
						}
					}

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

					doublereal* vol = new doublereal[my_global_temperature_struct.maxnod];

					for (integer i = 0; i < my_global_temperature_struct.maxnod; i++) {
						vol[i] = 0.0;
					}

					// Преобразование температуры с сетки МКО на сетку МКЭ.
					SECOND_ORDER_QUADRATIC_RECONSTRUCTA(my_global_temperature_struct.maxnod,
						my_global_temperature_struct.maxelm, my_global_temperature_struct.pa,
						my_global_temperature_struct.nvtx, vol, t_for_Mechanical, min_x, min_y, min_z, my_global_temperature_struct.potent,
						my_global_temperature_struct, eps_mashine, false);

					delete[] vol;

				}
				//bPhysics_stop = false;
				// Вызов солвера Static Structural.
				doublereal* stub_Mechanical = nullptr;
				solve_Structural(my_global_temperature_struct, w, lw, true, operatingtemperature, b, lb,lu, 
					false, 1.0e9, stub_Mechanical, stub_Mechanical, 1.0,matlist, t_for_Mechanical);
				//bPhysics_stop = true;

				delete[] t_for_Mechanical;

				// Вычисление массы модели.
				massa_cabinet(my_global_temperature_struct, f,  flow_interior,
					b, lb, operatingtemperature,
					matlist);

				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
					}
					else {
						// Экспорт в программу tecplot температуры.
						//С АЛИС сетки.
						ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
					}
				}

			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
											 // Вычисление значения максимальной температуры внутри расчётной области и на её границах:
											 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
											 // 23 декабря 2015
											 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
											 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

			doublereal totaldeform_max = -1.0e+30;
			if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
			{
				for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
			}

			FILE* fp = NULL;
			
#ifdef MINGW_COMPILLER
			int err1 = 0;
			fp = fopen64("report.txt", "w");
			if (fp == NULL) err1 = 1;
#else
			errno_t err1 = 0;
			err1 = fopen_s(&fp, "report.txt", "w");
#endif

			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//system("PAUSE");
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			bool bMechanical = true;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
			// Печатает расход жидкости через выходную границу потока.
			// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
			// проходящий через выходную границу потока. 28,10,2019
			report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);


		}



		/*
		if (b_on_adaptive_local_refinement_mesh) {
		printf("t.maxbound=%d\n", t.maxbound);
		printf("v dvuch shagah ot ALICE sborki. \n");
		system("PAUSE");
		exit(1);  // Режим глобального тестирования сеточного генератора адаптивной локально измельчённой сетки.
		}

		if (b_on_adaptive_local_refinement_mesh) {
		printf("Solve temperature is compleate. \n");
		system("PAUSE");
		exit(1);  // Режим глобального тестирования сеточного генератора адаптивной локально измельчённой сетки.
		}
		*/

		//system("pause");

		// Нестационарная теплопроводность.
		if (1 && ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER)||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL)||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
			

			// Решаем только теплопередачу в твёрдом теле.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			// Это прерывание работает только на ПТБШ позволяя ускорить вычисления.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка:
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			bglobal_unsteady_temperature_determinant = true;
			// Здесь предполагается что мы решаем задачу чистой теплопроводности.
			// последний 13 параметр bconvective
			// По результатам тестовых расчётов, наиболее сильная зависимость результата,
			// наблюдается от подробности расчётной сетки. Влияние порядка аппроксимации на 
			// границе области сказывается весьма слабо.
			// для фиксированной сетки и фиксированных остальных 
			// параметров изменение порядка аппроксимации с первого на второй
			// даёт изменение температуры в нужную сторону на 5 градусов на фоне перегрева в 167 градусов в статике.
			// Для сравнения, перегрев в icepak равен 120 градусам, что близко к значениям паспортного Rt 
			// (мощность 6.875Вт, RT=16K/W) и экспериментальным
			// данным.
			doublereal dbeta = 1.3333333;//1.0; // если 1.0 то первый порядок аппроксимации на границе.
			dbeta = 1.0; // более стабильное значение.
			// массив отладочной информации,
			// конкретно для проверки подхода Рхи-Чоу 1983
			doublereal** rhie_chow = nullptr;

			// 29.06.2020
			if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
				calculate_Network_T_unsteady(my_global_temperature_struct, f,
					b, lb, w, lw, s, ls, matlist);
				//printf("unsteady temperature calculation is finished...\n");
			    //system("PAUSE");
			    //exit(1);
			}

			//std::cout << "steady_or_unsteady_global_determinant=" << steady_or_unsteady_global_determinant << std::endl;
			//getchar();

			//solve_nonlinear_temp(f[0], f, t, rhie_chow, b, lb, s, ls, w, lw, dbeta, flow_interior, false, nullptr, 0.001, false);
			bool bsecond_T_solver = false;
			if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) {
				// Температурный солвер на основе поячеечной сборки матрицы 10.11.2018.
				bsecond_T_solver = true;
			}
			
			if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
				(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)) {
				
				bool bMechanical = false; // true - расчёт механических деформаций нестационарных МКЭ.
				bool bTemperature = true; // true - расчёт нестационарной теплопередачи.
				// bMechanical   && bTemperature   - расчёт нестационарной механики и теплопередачи совместно.

				if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) {
					bMechanical = true;
					bTemperature = false;
				}
				if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) {
					bMechanical = true;
					bTemperature = true;
				}
				

				unsteady_temperature_calculation(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					matlist,
					operatingtemperature,
					gtdps, ltdp, lu, my_union, bsecond_T_solver, inx, xpos, bTemperature, bMechanical); // нестационарный температурный солвер
					
			}

			// Вычисляет среднюю температуру в жидких ячейках.
			doublereal Tavg_fluid = 0.0;
			doublereal ic_avg_Temp_fluid = 0.0;
			for (integer i_7 = 0; i_7 < f[0].maxelm; i_7++) {
				if ((f[0].ptr[i_7] >= 0) && (f[0].ptr[i_7] < my_global_temperature_struct.maxelm)) {
					Tavg_fluid += my_global_temperature_struct.potent[f[0].ptr[i_7]];
					ic_avg_Temp_fluid += 1.0;
				}
			}
			if (ic_avg_Temp_fluid > 1.0e-30) {
				Tavg_fluid /= ic_avg_Temp_fluid;
				printf("average fluid temperature is %e\n", Tavg_fluid);
			}
			else {
				//printf("no fluid cell\n");
			}

			// Вычисление массы модели.
			massa_cabinet(my_global_temperature_struct, f,  flow_interior,
				b, lb, operatingtemperature,
				matlist);

			// 10.10.2017
			// Построение двумерного графика температуры вдоль структуры
			// на заданный момент времени.
			xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);

			if (!bsecond_T_solver) {
				if (!b_on_adaptive_local_refinement_mesh) {
					// экспорт результата вычисления в программу tecplot360:
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
				else {
					// Экспорт в программу tecplot температуры.
					//С АЛИС сетки.
					ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
				}


				doublereal tmaxfinish = -273.15;
				// Вычисление значения максимальной температуры внутри расчётной области и на её границах:
				for (integer i = 0; i < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i++) {
					tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);
				}
				FILE* fp = NULL;
				
#ifdef MINGW_COMPILLER
				int err1 = 0;
				fp = fopen64("report.txt", "w");
				if (fp == NULL) err1 = 1;
#else
				errno_t err1 = 0;
				err1 = fopen_s(&fp, "report.txt", "w");
#endif
				// создание файла для записи.
				if ((err1) != 0) {
					printf("Create File report.txt Error\n");
					//system("PAUSE");
					system("pause");
				}
				else {
					// запись заголовка
					fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
					fclose(fp);
				}
				// 1 - solver/solid_static/
				bool bMechanical = true;
				report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
				// Печатает расход жидкости через выходную границу потока.
				// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
				// проходящий через выходную границу потока. 28,10,2019
				report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);
			}
			else {
				printf("THIS IS SECOND UNSTEADY TEMPERATURE SOLVER ON ALL MESHES.\n");
				printf("NO EXPOPRT TECPLOT.\n");
				printf("NO PRINT REPORT.\n");
			}
			printf("calculation complete...\n");
			// system("PAUSE");
		}

		

		// экспорт результата вычисления в программу tecplot360:
		// можно использовать как проверку построенной сетки.
		if (false) {
			exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, 
				my_global_temperature_struct.ncell, 
				f, my_global_temperature_struct,
				flow_interior, 0, bextendedprint,
				0, b, lb);
			printf("read values. OK.\n");
			if (bwait) {
				//system("PAUSE"); // debug avtosave
				system("pause");
			}
		}

		// Fluid dynamic стационарные установившиеся течения.
		if ((1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY)) {
			
#ifdef _OPENMP
			// Предупреждение о невозможности расчёта cfd на openMP 08.05.2019.
			printf("CFD not work in OPENMP ON and bparallelismold is true.\n");
			printf("uskorenie ot OPENMP otsutstvuet. Rabotaen odnopotochnaq versiq.\n");
			printf("variable bparallelismold must be equal false.\n");
			//system("PAUSE");
#endif

			told_temperature_global_for_HOrelax = new doublereal[my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound];
			bSIMPLErun_now_for_temperature = true;

			if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
				// надо также проверить включено ли для fluid материалов приближение Буссинеска.
				bool bbussinesk_7 = false;
#pragma omp parallel for 
				for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
					integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
					if (ib > -1) {
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							integer i_7 = b[ib].imatid;
							if (matlist[i_7].bBussineskApproach) {
#pragma omp critical
								{
									if (bbussinesk_7 == false) {
										bbussinesk_7 = true;
									}
								}
							}
						}
					}
				}
				if (bbussinesk_7) {
					bSIMPLErun_now_for_natural_convection = true;
				}
			}

			const integer ISIZEf = f[0].maxelm + f[0].maxbound;

			bHORF = true;
			bPamendment_source_old = new doublereal[ISIZEf];
#pragma omp parallel for
			for (integer i5 = 0; i5 < ISIZEf; i5++) {
				bPamendment_source_old[i5] = 0.0;
			}
			// exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint);
			//system("PAUSE");
			if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
				// надо также проверить включено ли для fluid материалов приближение Буссинеска.
				bool bbussinesk_7 = false;
#pragma omp parallel for 
				for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
					integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
					if (ib > -1) {
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							integer i_7 = b[ib].imatid;
							if (matlist[i_7].bBussineskApproach) {
#pragma omp critical 
								{
									if (bbussinesk_7 == false) {
										bbussinesk_7 = true;
									}
								}
							}
						}
					}
				}
				if (bbussinesk_7) {
					printf("Bussinesk approach Operating Temperature=%e\n", f[0].OpTemp); // Operating Temperature);
				}
			}

			//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
			//system("PAUSE");

			// стационарный гидродинамический решатель.
			steady_cfd_calculation(breadOk,
				eqin, dgx, dgy, dgz,
				continity_start,
				inumber_iteration_SIMPLE,
				flow_interior, f, my_global_temperature_struct, b, lb,
				s, ls, w, lw, matlist,
				gtdps, ltdp, bextendedprint, lu, my_union,inx,xpos);
			// xyplot( f, 0, t);
			// boundarylayer_info(f, t, flow_interior, w, lw);
			// 2 - solver/conjugate_heat_transfer_static/
			bool bMechanical = false;
			
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0/*2*/, matlist, bMechanical);
			// Печатает расход жидкости через выходную границу потока.
			// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
			// проходящий через выходную границу потока. 28,10,2019
			report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);

			// Вычисление массы модели.
			massa_cabinet(my_global_temperature_struct, f,  flow_interior,
				b, lb, operatingtemperature,
				matlist);

			// Построение двумерного графика гидродинамических величин.
			// Позиция отрезка вдоль которого строится график передаётся 
			// из графического интерфейса. AliceMesh v.0.45.
			xyplot(f, flow_interior, my_global_temperature_struct);

			// 10.10.2017
			// Построение двумерного графика температуры вдоль структуры
			// на заданный момент времени.
			xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);

			// экспорт результата вычисления в программу tecplot360:
			if (!b_on_adaptive_local_refinement_mesh) {
				exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
			}
			else {
				ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
			}
			
			save_velocity_for_init(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior);
			
			// exporttecplotxy360T_3D_part2_rev(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,b,lb);
			delete[] bPamendment_source_old;
			bPamendment_source_old = nullptr;
			delete[] told_temperature_global_for_HOrelax;
			told_temperature_global_for_HOrelax = nullptr;
		}


		// нестационарный гидродинамический решатель: (НЕ РАБОТАЕТ)
		if ((1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY))) {


#ifdef _OPENMP
			// Предупреждение о невозможности расчёта cfd на openMP 08.05.2019.
			printf("Unsteady CFD not work in OPENMP ON and bparallelismold is true.\n");
			printf("uskorenie ot OPENMP otsutstvuet. Rabotaen odnopotochnaq versiq.\n");
			printf("variable bparallelismold must be equal false.\n");
			//system("PAUSE");
#endif

			told_temperature_global_for_HOrelax = new doublereal[my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound];
			bSIMPLErun_now_for_temperature = true;


			if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
				// надо также проверить включено ли для fluid материалов приближение Буссинеска.
				bool bbussinesk_7 = false;
#pragma omp parallel for
				for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
					integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
					if (ib > -1) {
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							integer i_7 = b[ib].imatid;
							if (matlist[i_7].bBussineskApproach) {
#pragma omp critical
								{
									if (bbussinesk_7 == false) {
										bbussinesk_7 = true;
									}
								}
							}
						}
					}
				}
				if (bbussinesk_7) {
					bSIMPLErun_now_for_natural_convection = true;
				}
			}

			const integer ISIZEf = f[0].maxelm + f[0].maxbound;
			bHORF = true;
			bPamendment_source_old = new doublereal[ISIZEf];
#pragma omp parallel for
			for (integer i5 = 0; i5 < ISIZEf; i5++) {
				bPamendment_source_old[i5] = 0.0;
			}

			// exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint);
			//system("PAUSE");
			if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
				// надо также проверить включено ли для fluid материалов приближение Буссинеска.
				bool bbussinesk_7 = false;
#pragma omp parallel for 
				for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
					integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
					if (ib > -1) {
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							integer i_7 = b[ib].imatid;
							if (matlist[i_7].bBussineskApproach) {
								if (bbussinesk_7 == false)
								{
#pragma omp critical
									{
										bbussinesk_7 = true;
									}
								}
							}
						}
					}
				}
				if (bbussinesk_7) {
					printf("Bussinesk approach Operating Temperature=%e\n", f[0].OpTemp); // Operating Temperature);
				}
			}

			//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
			//system("PAUSE");

			if (1) {
				printf("Sorry unsteady cfd calcuation dont support... 21.07.2019\n");
				printf("Your may send your message to kirill7785@mail.ru.\n");
				system("pause");
			}
			else {
				// нестационарный гидродинамический решатель:
				usteady_cfd_calculation(breadOk, eqin,
					dgx, dgy, dgz,
					continity_start,
					inumber_iteration_SIMPLE,
					flow_interior,
					f, my_global_temperature_struct,
					b, lb, s, ls,
					w, lw, matlist, gtdps, ltdp, bextendedprint, lu, my_union,inx,xpos);
			}

			//xyplot( f, 0, t);
			// boundarylayer_info(f, t, flow_interior, w, lw);
			// 2 - solver/conjugate_heat_transfer_static/
			bool bMechanical = false;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0/*2*/, matlist, bMechanical);

			// Вычисление массы модели.
			massa_cabinet(my_global_temperature_struct, f,  flow_interior,
				b, lb, operatingtemperature,
				matlist);

			// Построение двумерного графика гидродинамических величин.
			// Позиция отрезка вдоль которого строится график передаётся 
			// из графического интерфейса. AliceMesh v.0.45.
			xyplot(f, flow_interior, my_global_temperature_struct);

			// 10.10.2017
			// Построение двумерного графика температуры вдоль структуры
			// на заданный момент времени.
			xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);

			// экспорт результата вычисления в программу tecplot360:
			if (!b_on_adaptive_local_refinement_mesh) {
				exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
			}
			else {
				ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
			}

			save_velocity_for_init(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior);
			// exporttecplotxy360T_3D_part2_rev(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,b,lb);

			delete[] bPamendment_source_old;
			bPamendment_source_old = nullptr;
			delete[] told_temperature_global_for_HOrelax;
			told_temperature_global_for_HOrelax = nullptr;

		}


		if (continity_start != nullptr) {
			delete[] continity_start;
			continity_start = nullptr;
		}

		if (inumber_iteration_SIMPLE != nullptr) {
			delete[] inumber_iteration_SIMPLE;
			inumber_iteration_SIMPLE = nullptr;
		}

	

	// Освобождение оперативной памяти.
	if (my_global_temperature_struct.xpos_copy != nullptr) {
		delete[] my_global_temperature_struct.xpos_copy;
		my_global_temperature_struct.xpos_copy = nullptr;
	}
	if (my_global_temperature_struct.ypos_copy != nullptr) {
		delete[] my_global_temperature_struct.ypos_copy;
		my_global_temperature_struct.ypos_copy = nullptr;
	}
	if (my_global_temperature_struct.zpos_copy != nullptr) {
		delete[] my_global_temperature_struct.zpos_copy;
		my_global_temperature_struct.zpos_copy = nullptr;
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


	// Освобождение оперативной памяти.
	if (xpos != nullptr) {
		delete[] xpos;
		xpos = nullptr;
	}
	if (ypos != nullptr) {
		delete[] ypos;
		ypos = nullptr;
	}
	if (zpos != nullptr) {
		delete[] zpos;
		zpos = nullptr;
	}

	if (bsource_term_radiation_for_relax != nullptr) {
		delete[] bsource_term_radiation_for_relax; // Релаксация источниковых членов радиационных потоков.
		bsource_term_radiation_for_relax = nullptr;
	}
	if (b_buffer_correct_source != nullptr) {
		delete[] b_buffer_correct_source;
		b_buffer_correct_source = nullptr;
	}

	printf("free memory begin...\n");
	if (bwait) {
		//system("PAUSE");
		system("pause");
	}

	if (rthdsd_no_radiosity_patch != nullptr) {
		delete[] rthdsd_no_radiosity_patch;
		rthdsd_no_radiosity_patch = nullptr;
	}


	// Быстрая обработка нелинейных граничных условий в Румба 0.14 решателе.
	if (qnbc != nullptr) {
		delete[] qnbc;
		qnbc = nullptr;
		iadd_qnbc_maxelm = 0;
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

	for (integer i_7 = 0; i_7 < lb; i_7++) {
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
	delete[] b; delete[] s; delete[] w; // освобождение памяти
	b = nullptr;
	s = nullptr;
	w = nullptr;
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
	}
	delete[] matlist;
	matlist = nullptr;
	delete[] gtdps;
	gtdps = nullptr;
	if (eqin.fluidinfo != nullptr) {
		delete[] eqin.fluidinfo;
		eqin.fluidinfo = nullptr;
	}
	// Нужно освободить оперативную память из под всех структур данных:
	free_level1_temp(my_global_temperature_struct);
	free_level2_temp(my_global_temperature_struct); // освобождение памяти из под матриц.
	// Освобождает память для LR начало.
	free_root(my_global_temperature_struct.rootWE, my_global_temperature_struct.iWE);
	free_root(my_global_temperature_struct.rootSN, my_global_temperature_struct.iSN);
	free_root(my_global_temperature_struct.rootBT, my_global_temperature_struct.iBT);
	if (my_global_temperature_struct.rootWE != nullptr) {
		delete[] my_global_temperature_struct.rootWE;
		my_global_temperature_struct.rootWE = nullptr;
	}
	if (my_global_temperature_struct.rootSN != nullptr) {
		delete[] my_global_temperature_struct.rootSN;
		my_global_temperature_struct.rootSN = nullptr;
	}
	if (my_global_temperature_struct.rootBT != nullptr) {
		delete[] my_global_temperature_struct.rootBT;
		my_global_temperature_struct.rootBT = nullptr;
	}
	// Освобождение памяти для LR конец.
	free_level1_flow(f, flow_interior);
	free_level2_flow(f, flow_interior); // освобождение памяти из под матриц.

	delete[] f;
	f = nullptr;

	if (sourse2Dproblem != nullptr) {
		delete[] sourse2Dproblem;
		sourse2Dproblem = nullptr;
	}
	if (conductivity2Dinsource != nullptr) {
		delete[] conductivity2Dinsource;
		conductivity2Dinsource = nullptr;
	}

	if (x_jacoby_buffer != nullptr) {
		// 30 октября 2016. 
		// В seidelsor2 сделан переключатель на метод нижней релаксации К.Г. Якоби.
		// Освобождение памяти из под Jacobi buffer.
		delete[] x_jacoby_buffer;
		x_jacoby_buffer = nullptr;
	}

	if (bvery_big_memory) {
		if (my_global_temperature_struct.database.x != nullptr) {
			free(my_global_temperature_struct.database.x);
			my_global_temperature_struct.database.x = nullptr;
		}
		if (my_global_temperature_struct.database.y != nullptr) {
			free(my_global_temperature_struct.database.y);
			my_global_temperature_struct.database.y = nullptr;
		}
		if (my_global_temperature_struct.database.z != nullptr) {
			free(my_global_temperature_struct.database.z);
			my_global_temperature_struct.database.z = nullptr;
		}
		if (my_global_temperature_struct.database.nvtxcell != nullptr) {
			for (integer i = 0; i <= 7; i++) {
				delete[] my_global_temperature_struct.database.nvtxcell[i];
			}
			delete[] my_global_temperature_struct.database.nvtxcell;
			my_global_temperature_struct.database.nvtxcell = nullptr;
		}
		if (my_global_temperature_struct.database.ptr != nullptr) {
			if (my_global_temperature_struct.database.ptr[0] != nullptr) {
				delete[] my_global_temperature_struct.database.ptr[0];
			}
			if (my_global_temperature_struct.database.ptr[1] != nullptr) {
				delete[] my_global_temperature_struct.database.ptr[1];
			}
			delete[] my_global_temperature_struct.database.ptr;
			my_global_temperature_struct.database.ptr = nullptr;
		}
	}

	for (integer i63 = 0; i63 < lu; i63++) {
		// Нужно освободить оперативную память из под всех структур данных:
		free_level1_temp(my_union[i63].t);
		free_level2_temp(my_union[i63].t); // освобождение памяти из под матриц.

		// Освобождает память для LR начало.
		free_root(my_union[i63].t.rootWE, my_union[i63].t.iWE);
		free_root(my_union[i63].t.rootSN, my_union[i63].t.iSN);
		free_root(my_union[i63].t.rootBT, my_union[i63].t.iBT);
		if (my_union[i63].t.rootWE != nullptr) {
			delete[] my_union[i63].t.rootWE;
			my_union[i63].t.rootWE = nullptr;
		}
		if (my_union[i63].t.rootSN != nullptr) {
			delete[] my_union[i63].t.rootSN;
			my_union[i63].t.rootSN = nullptr;
		}
		if (my_union[i63].t.rootBT != nullptr) {
			delete[] my_union[i63].t.rootBT;
			my_union[i63].t.rootBT = nullptr;
		}

		// Освобождение памяти для LR конец.
		free_level1_flow(my_union[i63].f, my_union[i63].flow_interior);
		free_level2_flow(my_union[i63].f, my_union[i63].flow_interior); // освобождение памяти из под матриц.

		delete[] my_union[i63].f;
		my_union[i63].f = nullptr;

		if (bvery_big_memory) {
			if (my_union[i63].t.database.x != nullptr) {
				free(my_union[i63].t.database.x);
				my_union[i63].t.database.x = nullptr;
			}
			if (my_union[i63].t.database.y != nullptr) {
				free(my_union[i63].t.database.y);
				my_union[i63].t.database.y = nullptr;
			}
			if (my_union[i63].t.database.z != nullptr) {
				free(my_union[i63].t.database.z);
				my_union[i63].t.database.z = nullptr;
			}
			if (my_union[i63].t.database.nvtxcell != nullptr) {
				for (integer i = 0; i <= 7; i++) {
					delete[] my_union[i63].t.database.nvtxcell[i];
				}
				delete[] my_union[i63].t.database.nvtxcell;
				my_union[i63].t.database.nvtxcell = nullptr;
			}
			if (my_union[i63].t.database.ptr != nullptr) {
				if (my_union[i63].t.database.ptr[0] != nullptr) {
					delete[] my_union[i63].t.database.ptr[0];
				}
				if (my_union[i63].t.database.ptr[1] != nullptr) {
					delete[] my_union[i63].t.database.ptr[1];
				}
				delete[] my_union[i63].t.database.ptr;
				my_union[i63].t.database.ptr = nullptr;
			}
		}
	}

	// Освобождение общей памяти в ILU буфере.
	if (milu_gl_buffer.alu_copy != nullptr) delete[] milu_gl_buffer.alu_copy;
	if (milu_gl_buffer.jlu_copy != nullptr) delete[] milu_gl_buffer.jlu_copy;
	if (milu_gl_buffer.ju_copy != nullptr) delete[] milu_gl_buffer.ju_copy;
	milu_gl_buffer.alu_copy = nullptr;
	milu_gl_buffer.jlu_copy = nullptr;
	milu_gl_buffer.ju_copy = nullptr;

	free_QSBid(); // Для ускоренной работы функции myisblock_id.

	flow_interior = 0;
	printf("free memory finish...\n");

	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY) {
		// Был только вызов сеточного генератора.
		printf("Mesh generation procedure is finish.\n");
	}
	else {
		printf("Calculation procedure is finish.\n");
	}

	if (bwait) {
		printf("Please, press any key to exit...\n");
		//system("PAUSE");
		system("pause");
	}

	// Освобождение оперативной памяти из под закона изменения
	// шага по времени.
	if (glTSL.table_law_piecewise_constant != nullptr) {
		delete[] glTSL.table_law_piecewise_constant;
		glTSL.table_law_piecewise_constant = nullptr;
	}

	// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
	// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
	// при повторном обращении идет обращение к ячейки хеш-таблицы.
	// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
	// Хеш-таблицы для automatic
	// Освобождение оперативной памяти из под хеш-таблицы.
	delete[] shorter_hash_X;
	delete[] shorter_hash_Y;
	delete[] shorter_hash_Z;
	delete[] bshorter_hash_X;
	delete[] bshorter_hash_Y;
	delete[] bshorter_hash_Z;

	// Освобождаем память из под массива строк.
	//delete[] StringList;

	calculation_main_end_time = clock();
	calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;


	/*printf("time=%d statistic vorst=%3.2f %% \n",calculation_main_seach_time,(float)(100.0*calculation_vorst_seach_time/calculation_main_seach_time));
	system("PAUSE");
	*/



	// Общее время вычисления.
	int im = 0, is = 0, ims = 0;
	im = (int)(calculation_main_seach_time / 60000); // минуты
	is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // секунды
	ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

	printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

	if (1 && (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT)) {
		//system("pause");
		//46.321
		//31.655
	}

	return calculation_main_seach_time;
} // main_body

int main(void)
{

	// Поддержка кириллицы в консоли Windows
	// Меняет разделитель целой и дробной части.
	//setlocale(LC_ALL, "");
	system("mode con cols=166 lines=12000");

	char ch_EXPORT_ALICE_ONLY = 'y';

	if (0) {
		// Оптимизация. Перебор на основе трёх переменных:
		// d_my_optimetric1_6_12_2019;
		// d_my_optimetric2_6_12_2019;
		// d_my_optimetric3_6_12_2019.
		// с глобальной областью видимости.
		// При оптимизации эти переменные требуется внедрить в код и 
		// и задать диапазоны их изменения.
		// Три переменных
		//int in_optimetric = 21;
		int* time_optimetric = new int[1000];
		doublereal* op_comp = new doublereal[1000];
		int iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.33; d_my_optimetric1_6_12_2019 < 0.7; d_my_optimetric1_6_12_2019 += 0.01) {
			//for (d_my_optimetric2_6_12_2019 = 0.09; d_my_optimetric2_6_12_2019 < 0.95; d_my_optimetric2_6_12_2019 += 0.005) {
				//for (d_my_optimetric3_6_12_2019 = 0.11; d_my_optimetric3_6_12_2019 < 0.12; d_my_optimetric3_6_12_2019 += 0.01) {
					//printf("%e %e %e %d\n", d_my_optimetric1_6_12_2019, d_my_optimetric2_6_12_2019, d_my_optimetric3_6_12_2019, time_optimetric[iscan_optimetric]);
					//printf("%e %d\n", d_my_optimetric1_6_12_2019, time_optimetric[iscan_optimetric]);
					time_optimetric[iscan_optimetric] = main_body(ch_EXPORT_ALICE_ONLY);
					op_comp[iscan_optimetric] = d_my_optimetric2_6_12_2019;
					iscan_optimetric++;
					
				//}
			//}
		}

		printf("\n\n\n");
		iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.33; d_my_optimetric1_6_12_2019 < 0.7; d_my_optimetric1_6_12_2019 += 0.01) {
			//for (d_my_optimetric2_6_12_2019 = 0.09; d_my_optimetric2_6_12_2019 < 0.95; d_my_optimetric2_6_12_2019 += 0.005) {
				//for (d_my_optimetric3_6_12_2019 = 0.11; d_my_optimetric3_6_12_2019 < 0.12; d_my_optimetric3_6_12_2019 += 0.01) {
					//printf("%e %e %e %d\n", d_my_optimetric1_6_12_2019, d_my_optimetric2_6_12_2019, d_my_optimetric3_6_12_2019, time_optimetric[iscan_optimetric]);
					printf("%e CA=%e %d\n", d_my_optimetric1_6_12_2019, op_comp[iscan_optimetric], time_optimetric[iscan_optimetric]);
					iscan_optimetric++;
				//}
			//}
		}
		delete[] time_optimetric;
		delete[] op_comp;
	}
	if (0) {
		// Оптимизация. Перебор на основе одной переменной.
		// d_my_optimetric1_6_12_2019.
		// с глобальной областью видимости.
		// При оптимизации эту переменную требуется внедрить в код и 
		// и задать диапазон её изменения.
		// Одна переменная
		
		int* time_optimetric = new int[300];
		int iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.65; d_my_optimetric1_6_12_2019 < 0.8; d_my_optimetric1_6_12_2019 += 0.01) {
			printf("%e %d\n", d_my_optimetric1_6_12_2019, time_optimetric[iscan_optimetric]);
			time_optimetric[iscan_optimetric] = main_body(ch_EXPORT_ALICE_ONLY);
			iscan_optimetric++;
		}

		printf("\n\n\n");
		iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.65; d_my_optimetric1_6_12_2019 < 0.8; d_my_optimetric1_6_12_2019 += 0.01) {
			printf("%e %d\n", d_my_optimetric1_6_12_2019, time_optimetric[iscan_optimetric]);
			iscan_optimetric++;
		}
		delete[] time_optimetric;
	}
	if (1) {
		main_body(ch_EXPORT_ALICE_ONLY);
	}

	// Освобождаем память из под массива строк.
	delete[] StringList;


	system("pause");
	return 0;
}
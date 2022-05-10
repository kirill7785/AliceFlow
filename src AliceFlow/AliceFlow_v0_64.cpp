// AliceFlow_v0_64.cpp
// 
// llvm-13.0.0-win64.exe
//
// 06.02.2022 Метод П.Л. Чебышева для давления.
// Ключ /Qlong_double 16 бит вещественное число.
// 
// 30.01.2022. Сглаживатель Пафнутия Львовича Чебышева.
// 1. n шаговый метод c оптимальным набором параметроов релаксации tau(i), i=1..n  нам каждом шаге. 
// tau(i),  i=1..n вычисляются на основе нулей многочлена Чебышева степени n (64, 128 задачи механики МКЭ) и границ спектра оператора. (Самарский, Николаев).
// 2. Переупорядочивание нулей многочлена Чебышева для устойчивости. Например для n равного степени двойки. (Николаев).
// 3. Верхняя граница спектра lambda_max вычисляется по теореме Гершгорина о кругах. (ddemidov, Жуков, Феодоритова).
// Верхняя граница спектра один раз вычисляется на каждом уровне и запоминается.
// 4. Итерационное уточнение нижней границы спектра о формулам (Жуков, Феодоритова). Начальное значение нижней границы lambda_min=lambda_max/30.0. 
// На каждом уровне многосеточной иерархии свои верхняя и нижняя границы спектра матричного оператора данного уровня.
// Уточнение нижней границы спектра 2 раза за V цикл (по разу на нисходящей и восходящей ветви V цикла).
// 5. Использование итерационного процесса Чебышева в качестве сглаживателя в алгебраическом многосеточном методе (Румба.v.0.14).
// 
// Время сборки проекта на 2*xeon2630v4  45s без библиотеки AMGCL.
// 2min 17s с библиотекой AMGCL.
// Нажмите кнопку Tools -> Options, а затем выберите Projects and Solutions -> Build and Run. 
// Измените MSBuild project build output verbosity на Normal . Таким образом, он будет отображать
//  время, прошедшее в каждом проекте решения, который он создает. Но, к сожалению, нет суммы
//  прошедшего времени по всему проекту. Вы также увидите начатую сборку Timestamp.
// 
// Tools -> Options -> Project and Solutions -> VC++ Project Settings -> Build Timing – Yes
// Сервис->Параметры->Проекты и решения->Параметры проекта VC++->Время сборки – Да
// 
// Еще важное замечение, многопоточная или она же многопроцессорная компиляция работает только, если выключен параметр минимального перестроения. Что бы его выключить идём в настройки проекта:
// Project Name->Properties->С / C++->Code Generation->Enable Minimal Rebuild – No
// Проект->Свойства->С / C++->Создание кода->Включить минимальное перестроение – Нет
//
// А теперь включаем саму многопоточность :
// Project Name->Properties->С / C++->General->Multi - processor Compilation – Yes
// Проект->Свойства->С / C++->Общие->Многопроцессорная компиляция – Да
// 
// 1.01.2022. Регионы адаптации АЛИС расчётной сетки. Продолжение.
// 
// 2021 год
// 
// 31.12.2021  Регионы адаптации АЛИС расчётной сетки. Начало.
// 17 декабря 2021 приступил к разработке и кодированию асемблесов произвольного уровня вложенности на базе работающих асемблесов глубины вложенности 1.
// 18 декабря 2021. Асемблесы произвольного уровня вложенностии для теплопередачи введены в строй.
// Период обращения спутника связи на высоте 600км равен 96мин=5760c (Время одного витка вокруг Земли).
// 19 - 20 ноября 2021 Второй температурный солвер с асемблесом правильное сшитие границ различных подобластей на внешней границе асемблеса.
// Сеточный генратор блочных сеток заработал ещё летом 2018 года Сутоморе Черногория (Sutomore Montenegro 2018).
// 22 ноября 2021 ускорение поиска соседей за счёт хеш таблицы при сшитиии границ различных подобластей на внешней границе асемблеса.
// 26 ноября 2021 вылавливание и устранение утечек памяти - memory leak's.
// 16-17 октября 2021;  откомпилировал компилятором gcc (g++ 9.3) как под Ubuntu Linux так и под Windows 10. GNU project.
// 28.09.2021 PMIS применённое к квадрату матрицы, дальнобойная интерполяция (5 F узлов подряд).
// Дало уменьшение операторной сложности до 1.8 - 1.9 (даже до значения 1.6) вместо 2.45-2.52 (примерно 25%).
// 01.08.2021 Теперь в Visual Studio community edition 2022 Preview 2.
// 19.07.2021. Для метода конечных элементов как для Механики так и для теплопередачи в твёрдом теле без
// конвекции можно использовать ICCG решатель (Метод сопряжённых градиентов усиленный неполным разложением Холецкого).
// 18.07.2021. В одном и том же исполняемом файле стало возможно использовать библиотеку AMGCL ddemidov как на GPU, 
// так и на CPU по выбору пользователя без перекомпиляции приложения. 
// 20.06.2021. Мой алгебраический многосеточный метод РУМБА.v.0.14 догнал и перегнал по скорости
// код немецкого качества amg1r5 1984 года Джона Руге и Клауса Штубена. Код РУМБА.v.0.14
// написан заново с нуля и не содержит заимствований сторонних кодовых баз. 
// Заработала стационарная механика 14.03.2021 как в случае приложения заданной силы
// так и при тепловом расширении. 
// 
// 2020 год
// 
// Заработал нестационарный cfd решатель 25.12.2020.
// 22.12.2020 - К нестационарному графовому методу решения подключен алгебраический
// многосеточный метод РУМБА v.0.14 как опция на выбор между ним и amg1r5.
// 12.12.2020 - OpenGL визуализация моделей на с++ FWGL. С Z буфером.
// 17.11.2020 - Начало внедрения в код CAD_stl obj для CAD конструкторской 
// геометрии из *.stl файлов. 
// 24.06.2020 - 16 июля 2020 графовый температурный солвер - network_T solver.
// 2 июня 2020 нелинейный температурный солвер на базе amg1r5.
// 27 февраля 2020 PMIS aggregator, FF1 amg interpolation. Acomp = 3.0;
// Как в магистерской диссертации Джефри Баттлера 2006 года. 
// Авторы агрегатора PMIS Ульрике М. Янг и Ханс де Стерк 2004 год.
// 12 января 2020 iluk сглаживатель для amg1r5 алгоритма.
// 
// 2019 год
// 
// 10 октября 2019 запрограммировал модель турбулентности K-Omega SST Ментера (RANS). 
// 2 октября 2019 запрограммировал модель турбулентности Спаларта Аллмареса (RANS).
// 4 августа 2019; 03 ноября 2019 откомпилировал компилятором gcc (g++ 9.1). GNU project. 
// 7-8 мая 2019 присоединил алгебраический многосеточный метод amgcl Дениса Демидова.
// 6 апреля 2019 откомпилирована в visual studio community edition 2019 (open source).
// 25.03.2019 Начал использовать PVS-Studio 6.0
// 19 марта(03) 2019 заработала гидродинамика на АЛИС сетках.
// 
// 2017-2018 год.
// 
// 6.05.2018 LINK: fatal error LNK1102: недостаточно памяти 2015 VS community.
// Выход: компиляция проекта с опцией /bigobj
// Подсветка синтаксиса для cuda:
// Tools->Options->Text Editor->File Extension Ввести cu и нажать кнопку add.
// 9 июля 2017 переход на 64 битные целые int64_t.
// 15 апреля 2017 откомпилирована в vs community edition 2017 (open source).
//
// 2016 год.
//
// 1 октября 2016 откомпилировал на nvidia cuda 8.0. 
// АЛИС сетки введены в строй для теплопередачи в твёрдом теле. 
// 11 января 2016 года добавил cl_agl_amg_v0_14.
//
//  2015 год.
//
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


//#include "stdafx.h"

#define NO_OPENGL_GLFW

#ifndef NO_OPENGL_GLFW

#include <GLFW/glfw3.h>

#endif

const bool bBlasiusIM = false; // Задача Блазиуса Игорь Михайлович (Конденсатор).

// Убирает печать логов в constr struct.
bool print_log_message = true;
bool print_log_message_T2Solver = true;

//#define VISUAL_TUDIO_2008_COMPILLER

#ifdef VISUAL_TUDIO_2008_COMPILLER
#define nullptr NULL // для Visual Studio 2008
#endif

// Раскомментировать в случае если сборка приложения 
//осуществляется компилятором g++ (g++ 9.1) от GNU.
// g++ AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt
// gcc_log.txt - файл с диагностическими сообщениями компилятора.
//#define MINGW_COMPILLER 1

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

// Поддерживается только компилятором gcc.

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

//#define NOMINMAX
#include <stdlib.h> 
#include <omp.h> // OpenMP
//-Xcompiler "/openmp" 

//#define NOMINMAX // раскомментировать
//#include <stdio.h>
//#include <stdlib.h>
//#include <windows.h> // раскомментировать

// Раскомментировать в случае если ошибки NOT NUMBER отслеживаются.
//#define MY_DEBUG_NOT_NUMBER

int idevice_Tesla = 0; // 1 Tesla K80; 0 - GeForce GTX 1080 Ti, Quadro K6000.
int inumiterSIMPLE371 = -1; // Номер итерации SIMPLE алгоритма.


int pause_suppression = 0; // 0 - пауза в конце, 1 - пропуск паузы в конце. (пауза - getchar().)

unsigned int calculation_main_start_time_global_Depend = 0; // начало счёта мс.

bool CAD_GEOMETRY_OCTREE_MESHGEN = false;
bool SPARSE_MESHER = false;// true;// дополнительное разрежение расчётной сетки в случае амплитрона.

						   // OpenGL global variables for visualization

bool bmyconvective7248 = false;// Конвекция или нет.

// Для построения асемблесов произвольной вложенности. 
bool* bconstruction_union_now = nullptr;

#ifndef NO_OPENGL_GLFW

const bool binverse_color_black_2_white = false;// Замена чёрного фона на белый.

GLfloat scale_all = 1000.0f;
int SCREEN_WIDTH = 1280; // 640; // Ширина экрна для вывода графики.
int SCREEN_HEIGHT = 960; // 480; // Высота экрана для вывода графики.

						 // Позиция центра экрана, где размещается отображаемый объект.
GLfloat halfScreenWidth = SCREEN_WIDTH / 2.0f;
GLfloat halfScreenHeight = SCREEN_HEIGHT / 2.0f;

// Глубина по оси z где рисуется отображаемый объект.
const int abbys = 500;

// n_render - размерность массива точек pa_opengl или -1 если память под pa_opengl не выделена.
int n_render = -1;

// Рекции отрисовщика на opengl на нажатие кнопок клавиатуры, пользовательский ввод.
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);


GLfloat rotationX = 0.0f; // Поворот вокруг оси Ox в градусах.
GLfloat rotationY = 0.0f; // Поворот вокруг оси Oy в градусах.

						  // Синусы и косинусы углов поворота вычисляются заранее один раз 
						  // и потом просто используются при рендере сцены.
GLfloat cosAlf, cosBet, sinAlf, sinBet;

// OpenGL end global declaration.
#endif

// Одна клетка по Oz
// Дает возможность решать двумерные задачи (2D) трехмерным солвером (3D solver).
// Направление по оси Oz вырождается в одну клетку расчётной области.
//2D чертеж модели нарисован в плоскости OXY.
bool b_one_cell_z = false;

// Свободный параметр используемый для отладки 
// программы AliceFlov_v.0.48 и передаваемый из
// программы интерфейса
// AliceMesh_v.0.45.
double free_debug_parametr1 = 0.0;


// Число потоков исполнения. 
// Передаётся из интерфейса пользователя.
// Шина памяти насыщается после 4 потоков на одном процессоре intel xeon 2630 v4 без использования сглаживателя Чебышева.
int number_processors_global_var = 1;

int number_cores() {
	return number_processors_global_var;

	// один processor xeon 2630 v4 => 4 thread.
	// Шина памяти насыщается после 4 потоков.

	// два процессора intel xeon 2630 v4 в двухсокетной системе 8 потоков исполнения.
	// Метод П.Л.Чебышева в качестве сглаживателя на основе полинома Чебышева высокой степени 64, 128 
	// позволяет, например, на задаче механики в твердом  теле получить реальное ускорение вычислений 
	// на 16 или 32 потоках на двух процессорах intel xeon 2630 v4 за счёт того что значительная часть
	// вычислительной нагрузки перекладывается на сглаживатель Чебышева.
}

//using namespace System;

// Если закоментировано - Проект собран как стандартное консольное приложение windows.
//#define _CUDA_IS_ACTIVE_ 1

// 11.12.2021
// Не использовать миксовую точность в задачах с излучением. Точность полностью утрачивается.
typedef double real_mix_precision;
// Неизвестных врядли будет больше чем уместится в типе unsigned int,
// однако для прохода по иерархии матриц по прежнему используется тип 
// int64_t.
// INTEGER_MIX_PRECISION_IS_INT  true - говорит о том что для индексации в матрицах используется тип int а не тип int64_t.
// INTEGER_MIX_PRECISION_IS_INT  true - вовсе не означает что в РУМБЕ у нас не может быть последовательности матриц с 
// суммарным числом ненулевых элементов большим чем размерность типа int. Просто здесь можно сэкономить память без 
// ущерба для работы алгоритма РУМБА на больших размерностях.
#define INTEGER_MIX_PRECISION_IS_INT  true
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



// Программа могла быть откомпилирована на старом оборудовании с 
// с компилятором Microsoft Visual Studio 2008.
#ifdef VISUAL_TUDIO_2008_COMPILLER

doublereal fmin(doublereal a, doublereal b) {
	doublereal r = a;
	if (b < r) r = b;
	return r;
}

doublereal fmax(doublereal a, doublereal b) {
	doublereal r = a;
	if (b > r) r = b;
	return r;
}

template <typename VType>
bool isfinite(VType x) {
	//return std::isfinite(x);

	if (x != x) {
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


// 1 - 64 битные целые.
// Если считается механика с цилиндрами то число ненулевых
// элементов в матрице превышает размер типа int  и требуется 
// использовать тип int64_t 64 битного целого. 
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


// Меньше этого уровня граничный узел дробится.
// -1 не используется
// 01.01.2022
integer BonLevelDrobim = -1;// zmeevik 7.
integer isizearrAdaptRegion = 0;
typedef struct TADAPTREGION {
	double xS, xE, yS, yE, zS, zE;
	integer ilevelMeshAdaptRegion; // -1 не используется
} ADAPTREGION;
ADAPTREGION* AdaptRegion = nullptr;

doublereal** HighOrderTermRelaxation = nullptr;
doublereal* RelaxationPAMTerm = nullptr;

// Дробление в заданном регионе.
// Синтаксис вызова функции:
// AdaptivRegionDrobim(ilevel, i, j, k, xpos, ypos, zpos);
//
bool AdaptivRegionDrobim(char ilevel, int i, int j, int k,
	doublereal*& xpos, doublereal*& ypos, doublereal*& zpos)
{
	bool b = true;

	// Дробление основной сетки.
	if (ilevel < BonLevelDrobim) {
		return false;// дробим.
	}

	for (integer i_63 = 0; i_63 < isizearrAdaptRegion; ++i_63) {
		if ((xpos[i] >= AdaptRegion[i_63].xS) && (xpos[i] <= AdaptRegion[i_63].xE) &&
			(ypos[j] >= AdaptRegion[i_63].yS) && (ypos[j] <= AdaptRegion[i_63].yE) &&
			(zpos[k] >= AdaptRegion[i_63].zS) && (zpos[k] <= AdaptRegion[i_63].zE))
		{
			if (ilevel < AdaptRegion[i_63].ilevelMeshAdaptRegion) {
				return false;// дробим.
			}
		}
	}

	return b;
}

#define bStableVersion  false

#if !bStableVersion

typedef struct T_PAIR_OMP
{
	int s, e;
} PAIR_OMP;

PAIR_OMP s_par[10] = { {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1} };

const int MAX_STRING_LENGTH_ELL_THERMAL = 27; // 27 диагоналей.
const int MAX_STRING_LENGTH_ELL_MECHANICAL = 81; // 81 диагональ.
int MAX_STRING_LENGTH_ELL = MAX_STRING_LENGTH_ELL_MECHANICAL;

int ell_global_number_of_strings_26_11_2021 = -1;// Храним глобально количество аллоцированых строк в матрице.
doublereal** data_ell = nullptr;
int** coll_ell = nullptr;
int* coll_size = nullptr; // Максимальный размер строки.
int nnz_ell = 0;


// При обнаружении сбоя включить stable версия данной функции.
// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
template <typename MY_IND_TYPE>
void ell_to_CRS(doublereal*& val, MY_IND_TYPE*& col_ind, MY_IND_TYPE*& row_ptr, integer nodes);

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

enum class DEFAULT_CABINET_BOUNDARY_CONDITION {
	ADIABATIC_WALL_BC = 0, // теплоизоляция, однородное условие Неймана
	NEWTON_RICHMAN_BC = 1, // Теплоотдача конвекцией в подвижную среду теплоносителя, для имитации конвекции
	STEFAN_BOLCMAN_BC = 2, // Теплооотдача излучением на удаленный объект с заданной температурой
	MIX_CONDITION_BC = 3 // Смешанное условие - одновременно теплоотдача конвекцией в подвижную среду теплоносителя и излучение на удаленный объект с заданной температурой
};

//const unsigned int DIRICHLET_FAMILY=1;
//const unsigned int NEIMAN_FAMILY=2;
//const unsigned int NEWTON_RICHMAN_FAMILY=3;
//const unsigned int STEFAN_BOLCMAN_FAMILY=4;

enum class WALL_BOUNDARY_CONDITION {
	DIRICHLET_FAMILY = 1, // Условие Дирихле, заданная температура на стенке.
	NEIMAN_FAMILY = 2, // Теплоизолированная стенка, на стенке стоит однородное условие Неймана.
	NEWTON_RICHMAN_FAMILY = 3, // Теплоотдача для ячеек сетки на стенке в подвижную среду теплоносителя, для имитации конвекции вблизи стенки. Граничное тепловое условия Ньютона - Рихмана.
	STEFAN_BOLCMAN_FAMILY = 4 // Теплоотдача излучением для ячеек сетки на стенке  на удаленный объект с заданной температурой. Граничное условие Стефана Больцмана.
};

enum ORDER_INTERPOLATION { parabola_MNK = 0, line_MNK = 1, cubic_parabola = 2 };
enum class ORDER_DERIVATIVE { FIRST_ORDER = 1, SECOND_ORDER = 2 };

enum class INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 { ZERO_INIT = 0, RANDOM_INIT = 1 };

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

enum class LINE_DIRECTIONAL {
	X_LINE_DIRECTIONAL = 0,  // YZ plane
	Y_LINE_DIRECTIONAL = 1,  // XZ plane
	Z_LINE_DIRECTIONAL = 2   // XY plane
};

// Параметры преобразователя xyplot графиков для отчетов.
// 5.01.2018
typedef struct Tpatcher_for_print_in_report {
	// Отображение графика от минимальной позиции на оси абсцисс до
	// максимальной позиции на оси абсцисс
	doublereal fminimum, fmaximum;
	// в заданном направлении одной из трёх возможных координатных осей.
	LINE_DIRECTIONAL idir; // 0 - X, 1 - Y, 2 - Z.
	Tpatcher_for_print_in_report() {
		fminimum = -1.0e+30;
		fmaximum = 1.0e+30;
		idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; // 0 - X, 1 - Y, 2 - Z.
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


enum class AMG1R5_OUT_ITERATOR { NONE_only_amg1r5, BiCGStab_plus_amg1r5, FGMRes_plus_amg1r5, Non_Linear_amg1r5 };

typedef struct TAMG1r5_Data {
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

	AMG1R5_OUT_ITERATOR stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5; // BiCGStab + amg1r5.
} AMG1r5_Data;

AMG1r5_Data AMG1r5Dat;

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

int iP_perefirier_start = 0; // индекс с которого стартует алгоритм Катхилла Маки.

bool b_setup_CathilMC_Temp = false;
int* new_number_CathilMC_Temp = nullptr;
int* new_number_internal_CathilMC_Temp = nullptr;
int* new_number_bound_CathilMC_Temp = nullptr;
int* rev_number_CathilMC_Temp = nullptr;

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
									   // WENO взвешенные не осцилирующие.
const unsigned int UNEVEN_WENO5 = 1036; // 14.01.2021

bool bglobal_first_start_radiation = true;

// Если мы решаем нестационарную задачу теплопередачи в твердом теле.
bool bglobal_unsteady_temperature_determinant = false;

// Выбор сеточного генератора:
// simplemeshgen == 0 или unevensimplemeshgen ==1.
// По умолчанию используется simplemeshgen == 0.
enum class CONFORMAL_MESH_GENERATOR_SELECTOR {
	SIMPLEMESHGEN_MESHER = 0, // Рекомендуется для очень простой геометрии типо один кубик. Неравномерная сетка со сгущением сеточных линий.
	UNEVENSIMPLEMESHGEN_MESHER = 1, // Рекомендуется для очень простой геометрии типо один кубик. Еще более неравномерная сетка с еще большим  сгущением сеточных линий.
	COARSEMESHGEN_MESHER = 2 // Рекомендуется для сложных комплексных геометрий.
};
CONFORMAL_MESH_GENERATOR_SELECTOR iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; // обычный сеточный генератор.

																												  // нестационарное или стационарное моделирование.
enum class PHYSICAL_MODEL_SWITCH {
	STEADY_TEMPERATURE = 0, UNSTEADY_TEMPERATURE = 1,
	MESHER_ONLY = 2, CFD_STEADY = 3,
	STEADY_STATIC_STRUCTURAL = 5, STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE = 6,
	SECOND_TEMPERATURE_SOLVER = 7,
	PREOBRAZOVATEL_FOR_REPORT = 8, CFD_UNSTEADY = 9,
	NETWORK_T = 10, NETWORK_T_UNSTEADY = 11, UNSTEADY_STATIC_STRUCTURAL = 12,
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

// Используется для досрочного прерывания внутреннего цикла for( ; ; ) при нестационарных расчётах
// графовым методом в случае если алгоритм забуксует. 27.06.2021.
bool GLOBAL_identity_situation_in_RUMBA_for_NEtworkT_solver = false;

// Использовать ли адаптивные локально измельчённые расчётные сетки.
bool b_on_adaptive_local_refinement_mesh = false;
enum class TYPE_ALICE_MESH { ONE_PASS_COARSE_ALICE_MESH = 0, MULTI_PASS_MEDIUM_ALICE_MESH = 1 };
TYPE_ALICE_MESH itype_ALICE_Mesh = TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH;// Тип АЛИС сетки.

																				 // Впервые применено для модели
																				 // блок вум без теплоотвода на массивной железке (один из возможных циклов).
typedef struct TPiecewiseConstantTimeStepLaw {
	doublereal time, timestep, m;
	TPiecewiseConstantTimeStepLaw() {
		time = 0.0; // Момент времени с которого начинается очередной шаг по времени, с
		timestep = 0.0; // Величина шага по времени, с
		m = 0.0; // Множитель на который домножается тепловая мощность в каждом элементе (кубике).
	}
} PiecewiseConstantTimeStepLawTimeStepLaw;

// 0 - Linear, 
// 1 - Square Wave,
// 2 - Square Wave 2, 
// 3 - Hot Cold (Нагрев до заданного момента времени, а потом остывание заданное время. Шаг по времени двойной логарифмический.)
// 4 - Piecewise Constant 20.12.2019 - Закон изменения заданный пользователем.
enum class TIME_STEP_lAW_SELECTOR {
	LINEAR = 0, // Шаг по времени изменяется по закону возрастающей геометрической прогресии от микросекунд до заданной пользователем величины времени.
	SQUARE_WAVE = 1, // Чередование включений тепловой мощности и её выключения. Многократно до заданного момента времени. Закон с заданной длительностью имплульса tau и скважностью Q.
	SQUARE_WAVE2 = 2, // Усовершествованная версия закона чередование включений тепловой мощности и её выключения. Используется для "Космических" задач.
	HOT_COLD = 3, // Тепловая мощность включена до заданного момента времени, а потом выключена заданное время. Шаг по времени двойной логарифмический, для описания наиболее резких изменений температуры.
	PIECEWISE_CONSTANT = 4 // Закон изменения шагов по времени и множителей для тепловой мощности заданный пользователем.
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

	doublereal Factor_a_for_Linear; // Знаменатель возрастающей геометрической прогрессии = 1.0 + Factor_a_for_Linear.
	doublereal tau; // длительность импульса для Square Wave
					// 06_03_2017 скважность может быть и дробной.
	doublereal Q; // Скважность для Square Wave.
				  // Импульсный режим для альтернативного Square Wave 2.
				  // off_multiplyer - Тепловая мощность теперь может подаваться в режиме пауза. 
				  // Множитель off_multiplyer домножает пиковую тепловую мощность на себя и задаёт её в паузе.
	doublereal m1, tau1, tau2, tau_pause, T_all, off_multiplyer;
	integer n_cycle; // 20 Циклов.
					 // hot cold reshime (double linear)
	doublereal on_time_double_linear; // 3c включено.

									  // 4-ый закон изменения шага по времени.
									  // Пользовательский закон изменения шагов по времени и тепловых мощностей.
	integer n_string_PiecewiseConst;
	PiecewiseConstantTimeStepLawTimeStepLaw* table_law_piecewise_constant;

	TTimeStepLaw() {
		// 0 - Linear, 
		// 1 - Square Wave,
		// 2 - Square Wave 2, 
		// 3 - Hot Cold (Евдокимова Н.Л.)
		// 4 - Piecewise Constant 20.12.2019 - Закон изменения заданный пользователем.
		id_law = TIME_STEP_lAW_SELECTOR::LINEAR;
		Factor_a_for_Linear = 0.2;
		tau = 60.0E-6; // длительность импульса для Square Wave
					   // 06_03_2017 скважность может быть и дробной.
		Q = 10.0; // Скважность для Square Wave.
				  // Импульсный режим для альтернативного Square Wave 2.
		m1 = 1.0; tau1 = 0.0; tau2 = 0.0; tau_pause = 0.0; T_all = 0.0; off_multiplyer = 0.0;
		n_cycle = 20; // 20 Циклов.
					  // hot cold reshime (double linear)
		on_time_double_linear = 3.0; // 3c включено.
									 // 4 закон изменения шага по времени.
		n_string_PiecewiseConst = 0;
		table_law_piecewise_constant = nullptr;
	}
} TimeStepLaw;

TimeStepLaw glTSL;

// 24 декабря 2016. 
// Для ускорения счёта нелинейных задач в РУМБА 0.14 решателе.
typedef struct TQuickNonlinearBoundaryCondition
{
	doublereal emissivity; // коэффициент излучательной способности.
	doublereal ViewFactor; // Фактор видимости.
	doublereal Tamb, dS; // Температура среды и площадь контакта со средой.
	doublereal film_coefficient; // Коэффициент теплоотдачи при конвективном охлаждении стенки.
	bool bactive;
	bool bStefanBolcman_q_on;
	bool bNewtonRichman_q_on;

	TQuickNonlinearBoundaryCondition() {
		emissivity = 0.8;
		ViewFactor = 1.0; // Фактор видимости.
		Tamb = 20.0, dS = 0.0;
		film_coefficient = 3.0;
		bactive = false;
		bStefanBolcman_q_on = false;
		bNewtonRichman_q_on = false;
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
bool* sourse2Dproblem = nullptr;
doublereal* conductivity2Dinsource = nullptr;

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
enum class WHAT_VISIBLE_OPTION { FLUID_AND_SOLID_BODY_VISIBLE = 0, ONLY_SOLID_BODY_VISIBLE = 1 };
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
enum class SECOND_T_SOLVER_ID_SWITCH {
	BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER = 0,
	DIRECT_SECOND_T_SOLVER = 1, CAMG_RUMBA_v0_14_SECOND_T_SOLVER = 2,
	AMG1R5_SECOND_T_SOLVER = 3, AMGCL_SECONT_T_SOLVER = 4,
	ICCG_SECOND_T_SOLVER = 5
};
SECOND_T_SOLVER_ID_SWITCH iswitchsolveramg_vs_BiCGstab_plus_ILU6 = SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; // BiCGStab + ILU6.

bool bwait = false; // если false, то мы игнорируем getchar().
					// Если задать нечто отличное от 1e-10 то программа уходит очень долгий цикл
const doublereal admission = 1.0e-30; //1.0e-10 // для определения совпадения двух вещественных чисел.

									  // 1.0e-30 слишком мало. При значении 1.0e-30 стенки не идентифицируются. MCB - не присваивается. 27.11.2020.
const doublereal admission_bon_con = 1.0e-10; // для идентификации стенок на которых задано граничное условие.

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
// Префиксная неявная схема.
int iprefix_Scheme_Flow = UDS;
unsigned int iFLOWScheme = UDS; // Противопоточная первого порядка
unsigned int iTEMPScheme = UDS; // Противопоточная первого порядка
unsigned int iTURBScheme = UDS; // Противопоточная первого порядка

								// включает более быстро сходящийся алгоритм SIMPLEC
								// SIMPLEC Van Doormal and Raithby, 1984 год.
								// SIMPLEC отличается от SIMPLE только в двух вещах:
								// 1. В SIMPLEC не используется нижняя релаксация для давления при коррекции давления, т.е. alphaP=1.0.
								// 2. В SIMPLEC псевдо время пропорционально tau ~ alpha/(1-alpha), а в SIMPLE tau ~ alpha. 
								// В остальном алгоритмы полностью совпадают.
enum class SIMPLE_CFD_ALGORITHM { SIMPLE_Carretto = 0, SIMPLEC_Van_Doormal_and_Raithby = 1 };
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

#ifdef _OPENMP
const integer inumcore = 2; // число ядер процессора
#endif

const bool bparallelizm_old = false;

// Структура границ деления используется при параллельных расчётах и
// декомпозиции области на две подобласти:
typedef struct TPARBOUND
{
	integer ileft_start, ileft_finish; // Пределы узлов для левой области разбиения
	integer iright_start, iright_finish; // Пределы узлов для правой области разбиения
	integer iseparate_start, iseparate_finish; // Пределы узлов для области разделителя левого множества узлов от прапвого.
	bool active; // активность декомпозиции.

	TPARBOUND() {
		ileft_start = -1; ileft_finish = -2;
		iright_start = -1; iright_finish = -2;
		iseparate_start = -1; iseparate_finish = -2;
		active = false; // активность декомпозиции.
	}
} PARBOUND;


// Структура данных используемая для распараллеливания.
typedef struct TPARDATA {
	integer ncore; // 1, 2, 4, 8.
	integer* inumerate;
	// это для ncore==2;
	PARBOUND b0;
	// это для ncore==4;
	PARBOUND b00, b01;
	// это для ncore==8;
	PARBOUND b000, b001, b010, b011;
	TPARDATA() {
		ncore = 2; // 1, 2, 4, 8.
		inumerate = nullptr;
	}
} PARDATA;

PARDATA nd; // nd - nested desection.

			// используется в Румбе для ускорения счёта.
doublereal* rthdsd_no_radiosity_patch = nullptr;

// При дроблении (если bdroblenie4=true)
// каждая из шести граней иожет граничить с четырьмя соседними ячейками.
/*
typedef struct TALICE_PARTITION {
bool bdroblenie4;
int iNODE1, iNODE2, iNODE3, iNODE4;
TALICE_PARTITION() {
bdroblenie4=false;
iNODE1=-1; iNODE2=-1; iNODE3=-1; iNODE4=-1;
}
void define_structural_mesh_neighbour(int id) {
bdroblenie4 = false;
iNODE1 = id; iNODE2 = -1; iNODE3 = -1; iNODE4 = -1;
}
} ALICE_PARTITION;
*/



// Глобальное сгущение АЛИС сетки к источникам тепла. 
// Это должно увеличить точность нахождения температуры.
const bool b_thermal_source_refinement = true;

#include "adaptive_local_refinement_mesh.cpp" // АЛИС

#ifndef NO_OPENGL_GLFW
// Для визуализации расчетной модели с использованием библиотеки OpenGL

// Координаты узлов pa использующиеся при рендере сцены с помощью алгоритма Z-буфера.
TOCHKA* pa_opengl = nullptr; // исходные сохраненные не преобразованные координаты.
TOCHKA* pa_render = nullptr; // Преобразованные координаты вычисляющиеся заново при изменении scale_all.
doublereal minimum_val_for_render_pic, maximum_val_for_render_pic; // Минимальное и максимальное значение полевой величины использующееся для отрисовки.
int iNUMBER_FUNC_openGL = 1;
int iCURENT_FUNC_openGL = 0;
doublereal** potent_array_opengl = nullptr;
// Результат расчёта нестационарным CFD решателем можно просматривать на экране монитора,
// смотреть анимацию.
int iCFD_animation = 0;
int iNUMBER_ANIMATION_FUNCTIONS = -1;
int iCURENT_ANIMATION_FUNCTION = 0;
int iNUMBER_ANIMATION_CADERS = -1;
int iCURENT_ANIMATION_CADER = 0;
doublereal*** animation_sequence_functions_openGL = nullptr;

// Каркасная визуализация с удалением невидимых линий. Проверено 19.12.2020.
// Данная версия подходит для механической задачи, не содержит информации о блоках.
int DrawZbufferColor(int maxelm, int maxnod, int**& nvtx);

// OpenGL end.
#endif

// 21.04.2021
// Если b_adhesion_Mesh то мы сокращаем количество сеточных линий путём слипания некоторых из них.
// Под слипанием понимается что границы некоторых body block тел редактируются чтобы специально слипнутся.
// Это всё делается в рамках экономиии числа ячеек сетки так чтобы ликвидировать щели между телами.
// Но иногда щели важны и тогда переменную b_adhesion_Mesh надо положить в false.
bool b_adhesion_Mesh = true;

#include "constr_struct.cpp" // заполнение структур данных TEMPER и FLOW
#include "uniformsimplemeshgen.cpp" // сеточный генератор

#ifndef NO_OPENGL_GLFW

// Ограниченное количество используемых цветов.
int* mask_color_for_openGL = nullptr;

// Устанавливает по целочисленному индексу от 0 до 1020 
// текущее значение цвета на цветовой шкале.
void set_color_for_render_icolor(int icol)
{

	if ((0 <= icol) && (icol <= 255)) {
		// Синий голубой
		glColor3f(0.0, (icol) / 255.0, 1.0);
	}
	else if ((256 <= icol) && (icol <= 510)) {
		//голубой - зелёный
		glColor3f(0, 1.0, (255.0 - (icol - 255.0)) / 255.0);
	}
	else if ((511 <= icol) && (icol <= 765)) {

		// зелёный - желтый
		glColor3f((icol - 510) / 255.0, 1.0, 0.0);
	}
	else if ((766 <= icol) && (icol <= 1020)) {

		// Жёлтый - красный
		glColor3f(1.0, (255 - (icol - 765)) / 255.0, 0.0);
	};
}

#endif

// Также используется в amg1r5.
integer myi_max(integer ia, integer ib)
{
	integer ir;
	if (ia < ib) ir = ib;
	else ir = ia;
	return ir;
} // max

integer myi_min(integer ia, integer ib)
{
	integer ir;
	if (ia < ib) ir = ia;
	else ir = ib;
	return ir;
} // min

#ifndef NO_OPENGL_GLFW

  // Устанавливает цвет вешины при отрисовке
  // по передаваемому идентификатору функции id для семейства хранимых функций,
  // а также хранимых для него значений максимума и минимума функции.
void set_color_for_render(int id) {

	int icol = 0;
	if (iCFD_animation == 1) {
		icol = myi_min(1020, myi_max(0, round(1020 * ((animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][id] - minimum_val_for_render_pic) / (maximum_val_for_render_pic - minimum_val_for_render_pic)))));
	}
	else {
		icol = myi_min(1020, myi_max(0, round(1020 * ((potent_array_opengl[iCURENT_FUNC_openGL][id] - minimum_val_for_render_pic) / (maximum_val_for_render_pic - minimum_val_for_render_pic)))));
	}
	//if ((0 <= icol) && (icol <= 1020)) {
	//icol = mask_color_for_openGL[icol];// небольшое число различных цветов.
	//}
	set_color_for_render_icolor(icol);
}

// Устанавливает цвет вешины при отрисовке
// по передаваемому вещественному значению, а также хранимых для него значений максимума и минимума.
void set_color_for_render(doublereal potent_val) {

	int icol = myi_min(1020, myi_max(0, round(1020 * ((potent_val - minimum_val_for_render_pic) / (maximum_val_for_render_pic - minimum_val_for_render_pic)))));
	//if ((0 <= icol) && (icol <= 1020)) {
	//icol = mask_color_for_openGL[icol];// небольшое число различных цветов.
	//}
	set_color_for_render_icolor(icol);
}

// Каркасная визуализация с удалением невидимых линий. Грани заливаются цветом фона.
int DrawZbuffer(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block);
// Каркасная визуализация с удалением невидимых линий. Грани заливаются цветом передаваемой функции potent.
int DrawZbufferColor(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent);
// Каркасная визуализация с удалением невидимых линий и с освещением. Грани заливаются цветом передаваемой функции potent.
int DrawZbufferColorLight(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent);

#endif

// Фомин, Фомина о применении полинейного метода (прогонки) в сочетании с алгоритмом BiCGStab.
// Не рекомендуется к использованию. Лучше мультигрид и BiCGStab.
#include "my_LR.cpp" // полилинейный метод

// 8 января 2016.
const bool bvery_big_memory = true; // true нет записи в файл всё храним в оперативной памяти. Это существенно быстрее по скорости.

									// Данные определения пределов раньше находились в файле MenterSST.cpp.
const doublereal K_limiter_min = 1.0e-20; //1.0e-14; // 1.0e-14 Fluent limits
										  // 0.1 тоже подходит и сходимость стала лучше для k по сравнению со случаем 1.0.
										  // при значении 1.0E-20 расходимость. Значение k на уровне 1e-14 до 1e-10. небольшой блок вум с оребрением и воздушным охлаждением 1 0.5м/с.
										  // Справляется только метод сглаженной аггрегации ddemidov AMGCL.
										  // Вточном решении (сошедшимся) omega не опускалась ниже значения 1.5.
const doublereal Omega_limiter_min = 0.1; // 1.0; иначе будет расходимость. // 1.0e-20 Fluent limits
const doublereal Epsilon_limiter_min = 1.0e-20; // 1.0e-14; значение требует уточнения... 1.0e-20 Fluent limits

UNION* my_union = nullptr; // для объединения.

						   // Глобальное объявление
TEMPER my_global_temperature_struct;
int flow_interior = 0; // Суммарное число FLUID зон
FLOW* f = nullptr;



#include "my_linalg.cpp" // самописные функции линейной алгебры

// экспорт картинки в программу tecplot360
#include "my_export_tecplot3.cpp"

#include "my_material_properties.cpp" // библиотека реальных свойств материалов


// Для функций: 
// eqsolve_simple_gauss - решает СЛАУ методом исключения Гаусса
// eqsolv_simple_holesskii - решает СЛАУ методом разложения Холесского

// аппроксимация обобщённого уравнения конвекции-диффузии
// на совмещённой сетке
#include "pamendment3.cpp"

#include "shortest_distance.cpp" // вычисление кратчайшего расстояния до стенки


// Информация о формуле вычисления невязок заимствована из 
// icepak user guide.
typedef struct TFLUENT_RESIDUAL {
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
	doublereal res_turb_gamma_Langtry_Mentor; // Перемежаемость для модели Лангтрии Ментора.
	doublereal res_turb_Re_Theta_Langtry_Mentor; // Число Рейнольдса для  модели Лангтрии Ментора.

	TFLUENT_RESIDUAL() {
		// Данные невязки приводятся на момент конца решения СЛАУ.
		// невязки согласованные с программой FLUENT
		// т.е. вычисляемые по формуле fluent.
		res_vx = 1.0; // невязка X скорости
		res_vy = 1.0; // невязка Y скорости
		res_vz = 1.0; // невязка Z скорости
		res_no_balance = 1.0; // несбалансированные источники массы.
		operating_value_b = 1.0; // значение несбалансированных источников массы с предыдущей итерации.
		res_nusha = 1.0; // Невязка модифицированной кинетической турбулентной вязкости.
		res_turb_kinetik_energy = 1.0; // Невязка кинетической энергии турбулентных пульсаций в модели SST K-Omega.
		res_turb_omega = 1.0; // Невязка удельной скорости диссипации кинетической энергии турбулентных пульсаций в модели SST K-Omega.
		res_turb_kinetik_energy_std_ke = 1.0; // Невязка кинетической энергии турбулентных пульсаций в двухслойной модели k-epsilon
		res_turb_epsilon = 1.0; // Невязка скорости диссипации кинетической энергии турбулентных пульсаций в двухслойной модели k-epsilon
		res_turb_gamma_Langtry_Mentor = 1.0; // Перемежаемость для модели Лангтрии Ментора.
		res_turb_Re_Theta_Langtry_Mentor = 1.0; // Число Рейнольдса для  модели Лангтрии Ментора.
	}
} FLUENT_RESIDUAL;


// соединение решателя
#include "mysolverv0_03.cpp"


// нестационарный солвер для температуры
// на основе стационарного солвера,
// а также нестационарный солвер для 
// гидродинамики на основе стационарного солвера.
#include "my_unsteady_temperature.cpp"

// Препроцессинг для параллельной обработки.
#include "my_nested_dissection.cpp"

#include <ctime> // для замера времени выполнения.

integer ltdp = 0; // количество таблично заданных мощностей от температуры и смещения стока.
TEMP_DEP_POWER* gtdps = nullptr; // Garber temperature depend power sequence. 

								 // база данных материалов:
integer lmatmax = 0; // максимальное число материалов
TPROP* matlist = nullptr; // хранилище базы данных материалов

void check_data(TEMPER t) {
	if (t.potent != nullptr) {
		const integer isize0 = static_cast<integer>(t.maxelm) + static_cast<integer>(t.maxbound);
		for (integer i = 0; i <isize0; ++i) {
			if (t.potent[i] != t.potent[i]) {
				std::cout << "t.potent[" << i << "] is " << t.potent[i] << "\n";
				system("pause");
			}
		}
	}
} // check_data

int main_body(BLOCK*& b, int& lb, char ch_EXPORT_ALICE_ONLY = 'y') {
	//printLOGO();

	//system("PAUSE");	

	// количество блоков, источников и стенок, юнионов.
	//integer lb = 0;
	//BLOCK* b = nullptr;// список блоков

	int ls = 0, lw = 0, lu = 0;
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
	int inx = 30, iny = 30, inz = 30;
	int inxadd = -1, inyadd = -1, inzadd = -1;
	doublereal dgx = 0.0, dgy = 0.0, dgz = 0.0; // сила тяжести
	doublereal operatingtemperature = 20.0; // Operating Temperature 20.0 Град. С.

	lu = 0;
	// lu, my_union
	loadFromFile();
	premeshin("premeshin.txt", lmatmax, lb, ls, lw, matlist, b, s, w,
		dgx, dgy, dgz, inx, iny, inz, operatingtemperature, ltdp, gtdps, lu, my_union);

	freeStringList();

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif

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

	bool bCAD = true;
	for (int i60 = 1; i60 < lb; i60++) {
		if (b[i60].g.itypegeom != CAD_STL) {
			if (b[i60].itype != PHYSICS_TYPE_IN_BODY::HOLLOW)
			{
				bCAD = false;
			}
		}
	}
	if (lb == 1) {
		// Есть только кабинет 
		bCAD = false;
	}

	integer iCabinetMarker = 0;
	if (bCAD) {
		// Для CAD геометрии просто задаём равномерную сетку в объеме расчетной области.

		int isize = inx;
		xpos = new doublereal[static_cast<integer>(isize) + 1];
		xpos[0] = b[0].g.xS;
		doublereal dstep = (b[0].g.xE - b[0].g.xS) / (1.0 * isize);
		for (int i60 = 1; i60 <= isize; i60++) {
			xpos[i60] = b[0].g.xS + dstep * i60;
		}
		xpos[isize] = b[0].g.xE;

		isize = iny;
		ypos = new doublereal[static_cast<integer>(isize) + 1];
		ypos[0] = b[0].g.yS;
		dstep = (b[0].g.yE - b[0].g.yS) / (1.0 * isize);
		for (int i60 = 1; i60 <= isize; i60++) {
			ypos[i60] = b[0].g.yS + dstep * i60;
		}
		ypos[isize] = b[0].g.yE;

		isize = inz;
		zpos = new doublereal[static_cast<integer>(isize) + 1];
		zpos[0] = b[0].g.zS;
		dstep = (b[0].g.zE - b[0].g.zS) / (1.0 * isize);
		for (int i60 = 1; i60 <= isize; i60++) {
			zpos[i60] = b[0].g.zS + dstep * i60;
		}
		zpos[isize] = b[0].g.zE;
	}
	else {


		bconstruction_union_now = new bool[static_cast<integer>(lu) + 1];
		for (integer iu = 0; iu <= lu; iu++) {
			bconstruction_union_now[iu] = false;
		}
		bconstruction_union_now[iCabinetMarker] = true;

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
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER" << "\n";
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER" << "\n";
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER" << "\n";
				break;
			default:
				std::cout << "error: your mesh generator is undefined " << "\n";
				break;
			}


			system("pause");
			exit(1);
		}

		//system("pause");
		//for (int i73 = 0; i73 <= inx; i73++) {

			//std::cout << "x[" << i73 << "]=" << xpos[i73] << std::endl;
		//}
		//system("pause");

		bool b_we_shold_be_continue_union = true;
		int i_level_union = 1;

		while (b_we_shold_be_continue_union) {

			// Строим расчётную сетку в объединениях.
			for (integer iu = 0; iu < lu; iu++) {
				if (my_union[iu].active)
				{
					my_union[iu].inxadd = -1;
					my_union[iu].inyadd = -1;
					my_union[iu].inzadd = -1;
					my_union[iu].xposadd = nullptr;
					my_union[iu].yposadd = nullptr;
					my_union[iu].zposadd = nullptr;
					my_union[iu].xpos = nullptr;
					my_union[iu].ypos = nullptr;
					my_union[iu].zpos = nullptr;
					//std::cout << "union.inx=" << my_union[iu].inx << "union.iny=" << my_union[iu].iny << "union.inz=" << my_union[iu].inz << std::endl;
					//system("pause"); проверка пройдена.
					my_union[iu].t.inx_copy = -1;
					my_union[iu].t.iny_copy = -1;
					my_union[iu].t.inz_copy = -1;
					if (my_union[iu].t.xpos_copy != nullptr) {
						delete[] my_union[iu].t.xpos_copy;
						my_union[iu].t.xpos_copy = nullptr;
					}
					if (my_union[iu].t.ypos_copy != nullptr) {
						delete[] my_union[iu].t.ypos_copy;
						my_union[iu].t.ypos_copy = nullptr;
					}
					if (my_union[iu].t.zpos_copy != nullptr) {
						delete[] my_union[iu].t.zpos_copy;
						my_union[iu].t.zpos_copy = nullptr;
					}
					my_union[iu].t.xpos_copy = nullptr;
					my_union[iu].t.ypos_copy = nullptr;
					my_union[iu].t.zpos_copy = nullptr;
					integer iup1 = iu + 1;
					switch (my_union[iu].iswitchMeshGenerator) {
					case CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER: simplemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
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
			}


			for (integer iu = 0; iu < lu; iu++) {
				if (bconstruction_union_now[my_union[iu].iunion_parent + 1]) {
					bconstruction_union_now[iu + 1] = true;
				}
			}

			b_we_shold_be_continue_union = false;

			for (integer iu = 0; iu <= lu; iu++) {
				if (!bconstruction_union_now[iu]) {
					// Необработанные уровни юнионов найдены. Продолжаем.
					b_we_shold_be_continue_union = true;
				}
			}

			//std::cout << "i level union " << i_level_union << "succsefully construct.\n";
			//system("pause");

			i_level_union++;
		}

		//std::cout << "all level union " << i_level_union << "succsefully construct.\n";
		//system("pause");

		// В файле add_line.txt при его наличии содержатся сеточные линии дополнительные 
		// заданные пользователем вручную для модели.
		FILE* fp_add_line = NULL;

		// создание файла для чтения.
#ifdef MINGW_COMPILLER
		fp_add_line = fopen64("add_line.txt", "r");
		int err_add_line = 0;
		if (fp_add_line != NULL) {
			err_add_line = 0;
		}
		else {
			err_add_line = 1; // ошибка открытия.
		}
#else
		errno_t err_add_line;
		err_add_line = fopen_s(&fp_add_line, "add_line.txt", "r");
#endif

		if (fp_add_line != NULL) {
			if (err_add_line != 0) {
				//printf("Open File add_line Error\n");
				//system("pause");
				//exit(0);
				// return bfound;
			}
			else
			{

				//printf("incomming");
				//system("pause");
				// Это избранные сеточные линиии добавленные пользователем вручную.

				int ix_add0 = 0;
				int iy_add0 = 0;
				int iz_add0 = 0;
#ifndef MINGW_COMPILLER
				fscanf_s(fp_add_line, "%d", &ix_add0);
#else
				fscanf(fp_add_line, "%d", &ix_add0);
#endif

				for (int i0c = 0; i0c < ix_add0; i0c++) {
					float fin0 = 0.0;
#ifndef MINGW_COMPILLER
					fscanf_s(fp_add_line, "%f", &fin0);
#else
					fscanf(fp_add_line, "%f", &fin0);
#endif

					// Добавляем глобальные сеточные линии кабинета.
					if ((fin0 >= b[0].g.xS) && (fin0 <= b[0].g.xE)) {
						addboundary_rudiment(xpos, inx, fin0, YZ_PLANE, b, lb, w, lw, s, ls);
					}
					Sort_method(xpos, inx);

				}
				Sort_method(xpos, inx);

#ifndef MINGW_COMPILLER
				fscanf_s(fp_add_line, "%d", &iy_add0);
#else
				fscanf(fp_add_line, "%d", &iy_add0);
#endif

				for (int i0c = 0; i0c < iy_add0; i0c++) {
					float fin0 = 0.0;
#ifndef MINGW_COMPILLER
					fscanf_s(fp_add_line, "%f", &fin0);
#else
					fscanf(fp_add_line, "%f", &fin0);
#endif

					// Добавляем глобальные сеточные линии кабинета.
					if ((fin0 >= b[0].g.yS) && (fin0 <= b[0].g.yE)) {
						addboundary_rudiment(ypos, iny, fin0, XZ_PLANE, b, lb, w, lw, s, ls);
					}
					Sort_method(ypos, iny);

				}
				Sort_method(ypos, iny);

#ifndef MINGW_COMPILLER
				fscanf_s(fp_add_line, "%d", &iz_add0);
#else
				fscanf(fp_add_line, "%d", &iz_add0);
#endif

				for (int i0c = 0; i0c < iz_add0; i0c++) {
					float fin0 = 0.0;
#ifndef MINGW_COMPILLER
					fscanf_s(fp_add_line, "%f", &fin0);
#else
					fscanf(fp_add_line, "%f", &fin0);
#endif


					// Добавляем глобальные сеточные линии кабинета.
					if ((fin0 >= b[0].g.zS) && (fin0 <= b[0].g.zE)) {
						addboundary_rudiment(zpos, inz, fin0, XY_PLANE, b, lb, w, lw, s, ls);
					}
					Sort_method(zpos, inz);

				}
				Sort_method(zpos, inz);

				fclose(fp_add_line);
			}

		}


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
					xpos_copy = new doublereal[2 * (static_cast<integer>(inx) + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < inx; i74++) {
						xpos_copy[2 * i74] = xpos[i74];
						xpos_copy[2 * i74 + 1] = 0.5 * (xpos[i74] + xpos[i74 + 1]);
					}
					xpos_copy[2 * (inx + 1) - 2] = xpos[inx];
					delete[] xpos;
					xpos = nullptr;
					xpos = new doublereal[2 * (static_cast<integer>(inx) + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < 2 * (static_cast<integer>(inx) + 1) - 1; i74++) {
						xpos[i74] = xpos_copy[i74];
					}
					delete[] xpos_copy;
					xpos_copy = nullptr;
					inx = 2 * (inx + 1) - 2;
				}

				for (integer jiter = 0; jiter < jiterM; jiter++) {
					xpos_copy = new doublereal[2 * (static_cast<integer>(iny) + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < iny; i74++) {
						xpos_copy[2 * i74] = ypos[i74];
						xpos_copy[2 * i74 + 1] = 0.5 * (ypos[i74] + ypos[i74 + 1]);
					}
					xpos_copy[2 * (iny + 1) - 2] = ypos[iny];
					delete[] ypos;
					ypos = nullptr;
					ypos = new doublereal[2 * (static_cast<integer>(iny) + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < 2 * (static_cast<integer>(iny) + 1) - 1; i74++) {
						ypos[i74] = xpos_copy[i74];
					}
					delete[] xpos_copy;
					xpos_copy = nullptr;
					iny = 2 * (iny + 1) - 2;
				}

				for (integer jiter = 0; jiter < jiterM; jiter++) {
					xpos_copy = new doublereal[2 * (static_cast<integer>(inz) + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < inz; i74++) {
						xpos_copy[2 * i74] = zpos[i74];
						xpos_copy[2 * i74 + 1] = 0.5 * (zpos[i74] + zpos[i74 + 1]);
					}
					xpos_copy[2 * (inz + 1) - 2] = zpos[inz];
					delete[] zpos;
					zpos = nullptr;
					zpos = new doublereal[2 * (static_cast<integer>(inz) + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < 2 * (static_cast<integer>(inz) + 1) - 1; i74++) {
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
		integer maxelm_loc = (static_cast<integer>(inx) + 1) * (static_cast<integer>(iny) + 1) * (static_cast<integer>(inz) + 1);
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
				doublereal t_1 = NormaV(xpos, static_cast<integer>(inx) + 1);
				doublereal t_2 = NormaV(ypos, static_cast<integer>(iny) + 1);
				doublereal t_3 = NormaV(zpos, static_cast<integer>(inz) + 1);
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



	//for (int i92 = 0; i92 < inz; i92++) {
		//std::cout << i92 << " " << zpos[i92] << " \n";
		//system("pause");
	//}

	iCabinetMarker = 0;
	load_TEMPER_and_FLOW(my_global_temperature_struct, f, inx, iny, inz, xpos, ypos, zpos, flow_interior,
		b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
		dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iCabinetMarker);

	my_global_temperature_struct.operatingtemperature = operatingtemperature;


	for (integer iu = 0; iu < lu; iu++) {
		if (my_union[iu].active) {
			// Только активные асемблесы с mesh assembles separatelly.

			integer iup1 = iu + 1;
			load_TEMPER_and_FLOW(my_union[iu].t, my_union[iu].f,
				my_union[iu].inx, my_union[iu].iny, my_union[iu].inz,
				my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
				my_union[iu].flow_interior,
				b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
				dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iup1);

			my_union[iu].t.operatingtemperature = operatingtemperature;
		}
	}




	// Эти копии данных нужны для полного восстановления данных.
	my_global_temperature_struct.inx_copy = inx;
	my_global_temperature_struct.iny_copy = iny;
	my_global_temperature_struct.inz_copy = inz;
	my_global_temperature_struct.operatingtemperature_copy = operatingtemperature;
	my_global_temperature_struct.xpos_copy = new doublereal[static_cast<integer>(inx) + 1];
	my_global_temperature_struct.ypos_copy = new doublereal[static_cast<integer>(iny) + 1];
	my_global_temperature_struct.zpos_copy = new doublereal[static_cast<integer>(inz) + 1];
	// Данная информация нужна для экономии оперативной памяти,
	// некоторые данные будут выгружены из озу а потом восстановлены 
	// путём вычисления.
#pragma omp parallel for
	for (integer i_7 = 0; i_7 < static_cast<integer>(inx) + 1; i_7++) {
		my_global_temperature_struct.xpos_copy[i_7] = xpos[i_7];
	}
#pragma omp parallel for
	for (integer i_7 = 0; i_7 < static_cast<integer>(iny) + 1; i_7++) {
		my_global_temperature_struct.ypos_copy[i_7] = ypos[i_7];
	}
#pragma omp parallel for
	for (integer i_7 = 0; i_7 < static_cast<integer>(inz) + 1; i_7++) {
		my_global_temperature_struct.zpos_copy[i_7] = zpos[i_7];
	}

	for (integer iu = 0; iu < lu; iu++) {
		if (my_union[iu].active) {
			// Эти копии данных нужны для полного восстановления данных.
			my_union[iu].t.inx_copy = my_union[iu].inx;
			my_union[iu].t.iny_copy = my_union[iu].iny;
			my_union[iu].t.inz_copy = my_union[iu].inz;
			my_union[iu].t.operatingtemperature_copy = operatingtemperature;
			if (my_union[iu].t.xpos_copy != nullptr) {
				delete[] my_union[iu].t.xpos_copy;
				my_union[iu].t.xpos_copy = nullptr;
			}
			if (my_union[iu].t.ypos_copy != nullptr) {
				delete[] my_union[iu].t.ypos_copy;
				my_union[iu].t.ypos_copy = nullptr;
			}
			if (my_union[iu].t.zpos_copy != nullptr) {
				delete[] my_union[iu].t.zpos_copy;
				my_union[iu].t.zpos_copy = nullptr;
			}
			my_union[iu].t.xpos_copy = new doublereal[static_cast<integer>(my_union[iu].inx) + 1];
			my_union[iu].t.ypos_copy = new doublereal[static_cast<integer>(my_union[iu].iny) + 1];
			my_union[iu].t.zpos_copy = new doublereal[static_cast<integer>(my_union[iu].inz) + 1];
			// Данная информация нужна для экономии оперативной памяти,
			// некоторые данные будут выгружены из озу а потом восстановлены 
			// путём вычисления.
#pragma omp parallel for
			for (integer i_7 = 0; i_7 < static_cast<integer>(my_union[iu].inx) + 1; i_7++) {
				my_union[iu].t.xpos_copy[i_7] = my_union[iu].xpos[i_7];
			}
#pragma omp parallel for
			for (integer i_7 = 0; i_7 < static_cast<integer>(my_union[iu].iny) + 1; i_7++) {
				my_union[iu].t.ypos_copy[i_7] = my_union[iu].ypos[i_7];
			}
#pragma omp parallel for
			for (integer i_7 = 0; i_7 < static_cast<integer>(my_union[iu].inz) + 1; i_7++) {
				my_union[iu].t.zpos_copy[i_7] = my_union[iu].zpos[i_7];
			}
		}
	}

	// Освобождение оперативной памяти из под octree дерева.
	if (b_on_adaptive_local_refinement_mesh) {
		printf("free octree start...\n");
		//system("PAUSE");
		//system("PAUSE");
		integer maxelm_loc = (static_cast<integer>(inx) + 1) * (static_cast<integer>(iny) + 1) * (static_cast<integer>(inz) + 1);
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
	for (integer i = 0; i < inx; ++i) {
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
	for (integer i = 0; i < iny; ++i) {
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
	for (integer i = 0; i < inz; ++i) {
		if ((zpos[i + 1] - zpos[i]) < minimum_gap) minimum_gap = (zpos[i + 1] - zpos[i]);
		if ((zpos[i + 1] - zpos[i]) < d_zalipanie) {
			//printf("error: zalipanie po Z: zpos[%lld]=%e zpos[%lld]=%e inz=%lld\n", i, zpos[i], i + 1, zpos[i + 1], inz);
			std::cout << "error: zalipanie po Z: zpos[" << i << "]=" << zpos[i] << " zpos[" << i + 1 << "]=" << zpos[i + 1] << " inz=" << inz << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap Z=" << minimum_gap << std::endl;

#pragma omp parallel for
	for (integer iP = 0; iP < my_global_temperature_struct.maxelm; ++iP) {
		if ((my_global_temperature_struct.nvtx[0][iP] == 0) ||
			(my_global_temperature_struct.nvtx[1][iP] == 0) ||
			(my_global_temperature_struct.nvtx[2][iP] == 0) ||
			(my_global_temperature_struct.nvtx[3][iP] == 0) ||
			(my_global_temperature_struct.nvtx[4][iP] == 0) ||
			(my_global_temperature_struct.nvtx[5][iP] == 0) ||
			(my_global_temperature_struct.nvtx[6][iP] == 0) ||
			(my_global_temperature_struct.nvtx[7][iP] == 0))
		{
			printf("nvtx[%lld]: %d %d %d %d %d %d %d %d \n", iP, my_global_temperature_struct.nvtx[0][iP] - 1, my_global_temperature_struct.nvtx[1][iP] - 1, my_global_temperature_struct.nvtx[2][iP] - 1, my_global_temperature_struct.nvtx[3][iP] - 1, my_global_temperature_struct.nvtx[4][iP] - 1, my_global_temperature_struct.nvtx[5][iP] - 1, my_global_temperature_struct.nvtx[6][iP] - 1, my_global_temperature_struct.nvtx[7][iP] - 1);
		}
	}
#else
	doublereal minimum_gap = 1.0e60;
	for (integer i = 0; i < inx; ++i) {
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
	for (integer i = 0; i < iny; ++i) {
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
	for (integer i = 0; i < inz; ++i) {
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
	for (integer i = 0; i < flow_interior; ++i) {
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

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
		sourse2Dproblem = new bool[my_global_temperature_struct.maxbound];
		conductivity2Dinsource = new doublereal[my_global_temperature_struct.maxbound];
		// Нижняя релаксация источникового члена при радиационных потоках.
		bsource_term_radiation_for_relax = new doublereal[my_global_temperature_struct.maxelm];
		for (integer i_init = 0; i_init < my_global_temperature_struct.maxelm; i_init++) bsource_term_radiation_for_relax[i_init] = 0.0;
		b_buffer_correct_source = new doublereal[my_global_temperature_struct.maxelm];
	}

	// невязка continity будет измеряться по отношению к уровню 1e0.
	doublereal* continity_start = nullptr;
	continity_start = new doublereal[flow_interior];
	if (continity_start == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for continity start in main...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	for (integer i = 0; i < flow_interior; ++i) continity_start[i] = 1.0;

	int* inumber_iteration_SIMPLE = nullptr;
	inumber_iteration_SIMPLE = new int[flow_interior];
	if (nullptr == inumber_iteration_SIMPLE) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for inumber_iteration_SIMPLE in main...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	for (integer i = 0; i < flow_interior; ++i) inumber_iteration_SIMPLE[i] = 0; // начальная итерация алгоритма SIMPLE для каждой FLUID зоны.

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

#ifndef NO_OPENGL_GLFW
	scale_all = fmin(SCREEN_WIDTH / (fabs(b[0].g.xE - b[0].g.xS)), SCREEN_HEIGHT / (fabs(b[0].g.yE - b[0].g.yS)));

	//halfScreenWidth -= scale_all * 0.5 * (b[0].g.xE + b[0].g.xS);
	//halfScreenHeight -= scale_all * 0.5 * (b[0].g.yE + b[0].g.yS);

	halfScreenWidth = SCREEN_WIDTH / 2.0;
	halfScreenHeight = SCREEN_HEIGHT / 2.0;

#endif

	for (integer i = 0; i < my_global_temperature_struct.maxbound; ++i) {
		if (my_global_temperature_struct.border_neighbor[i].iB > my_global_temperature_struct.maxbound + my_global_temperature_struct.maxelm) {
			printf("MCB=%lld  node error %lld maxelm=%d, maxbound=%d\n", my_global_temperature_struct.border_neighbor[i].MCB, i, my_global_temperature_struct.maxelm, my_global_temperature_struct.maxbound);
			system("pause");
		}
		if (my_global_temperature_struct.border_neighbor[i].iI > my_global_temperature_struct.maxelm) {
			printf("MCB=%lld node error %lld maxelm=%d, maxbound=%d\n", my_global_temperature_struct.border_neighbor[i].MCB, i, my_global_temperature_struct.maxelm, my_global_temperature_struct.maxbound);
			system("pause");
		}
	}

	/*for (integer i=0; i<lw; ++i) {
	printf("%e  \n",w[i].Tamb);
	}
	//exporttecplotxy360T_3D(t.maxelm, t.ncell, t.nvtx, t.nvtxcell, t.pa, t.potent);
	exporttecplotxy360T_3D_part2(t.maxelm, t.potent);
	system("PAUSE"); // debug
	*/

	int* inumerate = new int[my_global_temperature_struct.maxelm];
	for (int i = 0; i < my_global_temperature_struct.maxelm; ++i) inumerate[i] = -1;


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

	// true - разбиение намного более лучшего качества.
	// более сбалансированное, т.к. медиана вычисляется точно.
	const bool bMEDIAN = true;


	if (number_cores() == 2) {
		std::wcout << "inx=" << inx << "  iny=" << iny << "  inz=" << inz << std::endl;
		std::cout << "min_zazor_x = " << min_zazor_x << " min_zazor_y = " << min_zazor_y << " min_zazor_z = " << min_zazor_z << std::endl;
		int isepnum = 0;
		int iL = 0;
		int iR = 0;


		if ((inx >= iny) && (inx >= inz)) {
			doublereal dseparate = xpos[inx / 2];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].x;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparate = list[(my_global_temperature_struct.maxnod - 1) / 2];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfound = false;
				doublereal dcandidate = dseparate;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].x;
					if (fabs(dcandidate - dseparate) < min_zazor_x) bfound = true;

				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].x;
				}

				if (bfound) {
					inumerate[i] = 3;
					isepnum++;
				}
				else if (dcandidate < dseparate) {
					inumerate[i] = 1;
					iL++;
				}
				else {
					inumerate[i] = 2;
					iR++;
				}
			}


		}
		else if ((iny >= inx) && (iny >= inz)) {
			doublereal dseparate = ypos[iny / 2];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].y;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparate = list[(my_global_temperature_struct.maxnod - 1) / 2];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfound = false;
				doublereal dcandidate = dseparate;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].y;
					if (fabs(dcandidate - dseparate) < min_zazor_y) bfound = true;

				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].y;
				}

				if (bfound) {
					inumerate[i] = 3;
					isepnum++;
				}
				else if (dcandidate < dseparate) {
					inumerate[i] = 1;
					iL++;
				}
				else {
					inumerate[i] = 2;
					iR++;
				}
			}

		}
		else if ((inz >= inx) && (inz >= iny)) {
			doublereal dseparate = zpos[inz / 2];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].z;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparate = list[(my_global_temperature_struct.maxnod - 1) / 2];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfound = false;
				doublereal dcandidate = dseparate;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].z;
					if (fabs(dcandidate - dseparate) < min_zazor_z) bfound = true;

				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].z;
				}

				if (bfound) {
					inumerate[i] = 3;
					isepnum++;
				}
				else if (dcandidate < dseparate) {
					inumerate[i] = 1;
					iL++;
				}
				else {
					inumerate[i] = 2;
					iR++;
				}
			}

		}

		std::cout << "isepnum=" << isepnum << " L: " << iL << " R: " << iR << std::endl;
		for (int i = 0; i < my_global_temperature_struct.maxelm; ++i) if (inumerate[i] == -1) {
			std::cout << "nested desection error" << i << "\n"; system("pause");
		}
		//system("pause");
	}

	if (number_cores() == 3) {
		std::wcout << "inx=" << inx << "  iny=" << iny << "  inz=" << inz << std::endl;
		std::cout << "min_zazor_x = " << min_zazor_x << "  min_zazor_y = " << min_zazor_y << "  min_zazor_z = " << min_zazor_z << std::endl;
		int isepnum = 0;
		int iL = 0, iC = 0, iR = 0;


		if ((inx >= iny) && (inx >= inz)) {
			doublereal dseparateL = xpos[inx / 3];
			doublereal dseparateR = xpos[2 * inx / 3];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].x;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparateL = list[(my_global_temperature_struct.maxnod - 1) / 3];
				dseparateR = list[2 * (my_global_temperature_struct.maxnod - 1) / 3];
				delete[] list;
			}


			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfoundL = false;
				bool bfoundR = false;
				doublereal dcandidate = dseparateL;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].x;
					if (fabs(dcandidate - dseparateL) < min_zazor_x) bfoundL = true;

				}
				{
					dcandidate = dseparateR;
					for (int j = 0; j < 8; ++j) {
						int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

						dcandidate = my_global_temperature_struct.pa[i_1].x;
						if (fabs(dcandidate - dseparateR) < min_zazor_x) bfoundR = true;

					}
				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].x;
				}

				if (bfoundL) {
					inumerate[i] = 4;
					isepnum++;
				}
				else if (bfoundR) {
					inumerate[i] = 5;
					isepnum++;
				}
				else if (dcandidate < dseparateL) {
					inumerate[i] = 1;
					++iL;
				}
				else if (dcandidate > dseparateR) {
					inumerate[i] = 3;
					++iR;
				}
				else {
					inumerate[i] = 2;
					++iC;
				}
			}


		}
		else if ((iny >= inx) && (iny >= inz)) {
			doublereal dseparateL = ypos[iny / 3];
			doublereal dseparateR = ypos[2 * iny / 3];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].y;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparateL = list[(my_global_temperature_struct.maxnod - 1) / 3];
				dseparateR = list[2 * (my_global_temperature_struct.maxnod - 1) / 3];
				delete[] list;
			}


			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfoundL = false;
				bool bfoundR = false;
				doublereal dcandidate = dseparateL;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].y;
					if (fabs(dcandidate - dseparateL) < min_zazor_y) bfoundL = true;

				}
				{
					dcandidate = dseparateR;
					for (int j = 0; j < 8; ++j) {
						int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

						dcandidate = my_global_temperature_struct.pa[i_1].y;
						if (fabs(dcandidate - dseparateR) < min_zazor_y) bfoundR = true;

					}
				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].y;
				}


				if (bfoundL) {
					inumerate[i] = 4;
					isepnum++;
				}
				else if (bfoundR) {
					inumerate[i] = 5;
					isepnum++;
				}
				else if (dcandidate < dseparateL) {
					inumerate[i] = 1;
					++iL;
				}
				else if (dcandidate > dseparateR) {
					inumerate[i] = 3;
					++iR;
				}
				else {
					inumerate[i] = 2;
					++iC;
				}

			}

		}
		else if ((inz >= inx) && (inz >= iny)) {
			doublereal dseparateL = zpos[inz / 3];
			doublereal dseparateR = zpos[2 * inz / 3];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].z;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparateL = list[(my_global_temperature_struct.maxnod - 1) / 3];
				dseparateR = list[2 * (my_global_temperature_struct.maxnod - 1) / 3];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfoundL = false;
				bool bfoundR = false;
				doublereal dcandidate = dseparateL;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].z;
					if (fabs(dcandidate - dseparateL) < min_zazor_z) bfoundL = true;

				}
				{
					dcandidate = dseparateR;
					for (int j = 0; j < 8; ++j) {
						int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

						dcandidate = my_global_temperature_struct.pa[i_1].z;
						if (fabs(dcandidate - dseparateR) < min_zazor_z) bfoundR = true;

					}
				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].z;
				}



				if (bfoundL) {
					inumerate[i] = 4;
					isepnum++;
				}
				else if (bfoundR) {
					inumerate[i] = 5;
					isepnum++;
				}
				else if (dcandidate < dseparateL) {
					inumerate[i] = 1;
					++iL;
				}
				else if (dcandidate > dseparateR) {
					inumerate[i] = 3;
					++iR;
				}
				else {
					inumerate[i] = 2;
					++iC;
				}

			}

		}

		std::cout << "isepnum=" << isepnum << " L: " << iL << " C: " << iC << " R: " << iR << std::endl;
		for (int i = 0; i < my_global_temperature_struct.maxelm; ++i) if (inumerate[i] == -1) {
			std::cout << "nested desection error" << i << "\n"; system("pause");
		}
	}

	if ((number_cores() == 4) || (number_cores() == 5)) {
		std::wcout << "inx=" << inx << "  iny=" << iny << "  inz=" << inz << std::endl;
		std::cout << "min_zazor_x = " << min_zazor_x << " min_zazor_y = " << min_zazor_y << " min_zazor_z = " << min_zazor_z << std::endl;
		int isepnum = 0;


		if ((inx >= iny) && (inx >= inz)) {
			doublereal dseparate = xpos[inx / 2];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].x;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparate = list[(my_global_temperature_struct.maxnod - 1) / 2];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfound = false;
				doublereal dcandidate = dseparate;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].x;
					if (fabs(dcandidate - dseparate) < min_zazor_x) bfound = true;

				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].x;
				}

				if (bfound) {
					inumerate[i] = 3;
					isepnum++;
				}
				else if (dcandidate < dseparate) {
					inumerate[i] = 1;
				}
				else {
					inumerate[i] = 2;
				}





			}

			if (iny >= inz) {
				doublereal dseparate1 = ypos[iny / 2];

				if (bMEDIAN) {

					doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
					for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
					{
						list[i] = my_global_temperature_struct.pa[i].y;
					}
					std::sort(list, list + my_global_temperature_struct.maxnod - 1);
					dseparate1 = list[(my_global_temperature_struct.maxnod - 1) / 2];
					delete[] list;
				}


				for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
				{

					if ((inumerate[i] == 1) || (inumerate[i] == 2)) {

						bool bfound = false;
						doublereal dcandidate1 = dseparate1;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

							dcandidate1 = my_global_temperature_struct.pa[i_1].y;
							if (fabs(dcandidate1 - dseparate1) < min_zazor_y) bfound = true;

						}

						dcandidate1 = 0.0;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
							dcandidate1 += 0.125 * my_global_temperature_struct.pa[i_1].y;
						}

						if (bfound) {
							if (inumerate[i] == 1) {
								inumerate[i] = 5;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 8;
							}
							isepnum++;
						}
						else if (dcandidate1 < dseparate1) {
							if (inumerate[i] == 1) {
								inumerate[i] = 4;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 7;
							}
						}
						else {
							if (inumerate[i] == 1) {
								inumerate[i] = 6;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 9;
							}
						}
					}
				}

			}
			else if (inz >= iny) {
				doublereal dseparate1 = zpos[inz / 2];

				if (bMEDIAN) {

					doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
					for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
					{
						list[i] = my_global_temperature_struct.pa[i].z;
					}
					std::sort(list, list + my_global_temperature_struct.maxnod - 1);
					dseparate1 = list[(my_global_temperature_struct.maxnod - 1) / 2];
					delete[] list;
				}

				for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
				{

					if ((inumerate[i] == 1) || (inumerate[i] == 2)) {

						bool bfound = false;
						doublereal dcandidate1 = dseparate1;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

							dcandidate1 = my_global_temperature_struct.pa[i_1].z;
							if (fabs(dcandidate1 - dseparate1) < min_zazor_z) bfound = true;

						}

						dcandidate1 = 0.0;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
							dcandidate1 += 0.125 * my_global_temperature_struct.pa[i_1].z;
						}

						if (bfound) {
							if (inumerate[i] == 1) {
								inumerate[i] = 5;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 8;
							}
							isepnum++;
						}
						else if (dcandidate1 < dseparate1) {
							if (inumerate[i] == 1) {
								inumerate[i] = 4;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 7;
							}
						}
						else {
							if (inumerate[i] == 1) {
								inumerate[i] = 6;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 9;
							}
						}
					}
				}

			}


		}
		else if ((iny >= inx) && (iny >= inz)) {
			doublereal dseparate = ypos[iny / 2];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].y;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparate = list[(my_global_temperature_struct.maxnod - 1) / 2];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfound = false;
				doublereal dcandidate = dseparate;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].y;
					if (fabs(dcandidate - dseparate) < min_zazor_y) bfound = true;

				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].y;
				}

				if (bfound) {
					inumerate[i] = 3;
					isepnum++;
				}
				else if (dcandidate < dseparate) {
					inumerate[i] = 1;
				}
				else {
					inumerate[i] = 2;
				}
			}

			if (inx >= inz) {
				doublereal dseparate1 = xpos[inx / 2];

				if (bMEDIAN) {

					doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
					for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
					{
						list[i] = my_global_temperature_struct.pa[i].x;
					}
					std::sort(list, list + my_global_temperature_struct.maxnod - 1);
					dseparate1 = list[(my_global_temperature_struct.maxnod - 1) / 2];
					delete[] list;
				}

				for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
				{

					if ((inumerate[i] == 1) || (inumerate[i] == 2)) {

						bool bfound = false;
						doublereal dcandidate1 = dseparate1;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

							dcandidate1 = my_global_temperature_struct.pa[i_1].x;
							if (fabs(dcandidate1 - dseparate1) < min_zazor_x) bfound = true;

						}

						dcandidate1 = 0.0;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
							dcandidate1 += 0.125 * my_global_temperature_struct.pa[i_1].x;
						}

						if (bfound) {
							if (inumerate[i] == 1) {
								inumerate[i] = 5;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 8;
							}
							isepnum++;
						}
						else if (dcandidate1 < dseparate1) {
							if (inumerate[i] == 1) {
								inumerate[i] = 4;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 7;
							}
						}
						else {
							if (inumerate[i] == 1) {
								inumerate[i] = 6;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 9;
							}
						}
					}
				}

			}
			else if (inz >= inx) {
				doublereal dseparate1 = zpos[inz / 2];

				if (bMEDIAN) {

					doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
					for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
					{
						list[i] = my_global_temperature_struct.pa[i].z;
					}
					std::sort(list, list + my_global_temperature_struct.maxnod - 1);
					dseparate1 = list[(my_global_temperature_struct.maxnod - 1) / 2];
					delete[] list;
				}

				for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
				{

					if ((inumerate[i] == 1) || (inumerate[i] == 2)) {

						bool bfound = false;
						doublereal dcandidate1 = dseparate1;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

							dcandidate1 = my_global_temperature_struct.pa[i_1].z;
							if (fabs(dcandidate1 - dseparate1) < min_zazor_z) bfound = true;

						}

						dcandidate1 = 0.0;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
							dcandidate1 += 0.125 * my_global_temperature_struct.pa[i_1].z;
						}

						if (bfound) {
							if (inumerate[i] == 1) {
								inumerate[i] = 5;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 8;
							}
							isepnum++;
						}
						else if (dcandidate1 < dseparate1) {
							if (inumerate[i] == 1) {
								inumerate[i] = 4;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 7;
							}
						}
						else {
							if (inumerate[i] == 1) {
								inumerate[i] = 6;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 9;
							}
						}
					}
				}

			}

		}
		else if ((inz >= inx) && (inz >= iny)) {
			doublereal dseparate = zpos[inz / 2];

			if (bMEDIAN) {

				doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
				{
					list[i] = my_global_temperature_struct.pa[i].z;
				}
				std::sort(list, list + my_global_temperature_struct.maxnod - 1);
				dseparate = list[(my_global_temperature_struct.maxnod - 1) / 2];
				delete[] list;
			}

			for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
			{
				bool bfound = false;
				doublereal dcandidate = dseparate;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

					dcandidate = my_global_temperature_struct.pa[i_1].z;
					if (fabs(dcandidate - dseparate) < min_zazor_z) bfound = true;

				}

				dcandidate = 0.0;
				for (int j = 0; j < 8; ++j) {
					int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
					dcandidate += 0.125 * my_global_temperature_struct.pa[i_1].z;
				}

				if (bfound) {
					inumerate[i] = 3;
					isepnum++;
				}
				else if (dcandidate < dseparate) {
					inumerate[i] = 1;
				}
				else {
					inumerate[i] = 2;
				}
			}

			if (inx >= iny) {
				doublereal dseparate1 = xpos[inx / 2];

				if (bMEDIAN) {

					doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
					for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
					{
						list[i] = my_global_temperature_struct.pa[i].x;
					}
					std::sort(list, list + my_global_temperature_struct.maxnod - 1);
					dseparate1 = list[(my_global_temperature_struct.maxnod - 1) / 2];
					delete[] list;
				}

				for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
				{

					if ((inumerate[i] == 1) || (inumerate[i] == 2)) {

						bool bfound = false;
						doublereal dcandidate1 = dseparate1;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

							dcandidate1 = my_global_temperature_struct.pa[i_1].x;
							if (fabs(dcandidate1 - dseparate1) < min_zazor_x) bfound = true;

						}

						dcandidate1 = 0.0;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
							dcandidate1 += 0.125 * my_global_temperature_struct.pa[i_1].x;
						}

						if (bfound) {
							if (inumerate[i] == 1) {
								inumerate[i] = 5;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 8;
							}
							isepnum++;
						}
						else if (dcandidate1 < dseparate1) {
							if (inumerate[i] == 1) {
								inumerate[i] = 4;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 7;
							}
						}
						else {
							if (inumerate[i] == 1) {
								inumerate[i] = 6;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 9;
							}
						}
					}
				}

			}
			else if (iny >= inx) {
				doublereal dseparate1 = ypos[iny / 2];

				if (bMEDIAN) {

					doublereal* list = new doublereal[my_global_temperature_struct.maxnod];
					for (int i = 0; i < my_global_temperature_struct.maxnod; ++i)
					{
						list[i] = my_global_temperature_struct.pa[i].y;
					}
					std::sort(list, list + my_global_temperature_struct.maxnod - 1);
					dseparate1 = list[(my_global_temperature_struct.maxnod - 1) / 2];
					delete[] list;
				}


				for (int i = 0; i < my_global_temperature_struct.maxelm; ++i)
				{

					if ((inumerate[i] == 1) || (inumerate[i] == 2)) {

						bool bfound = false;
						doublereal dcandidate1 = dseparate1;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;

							dcandidate1 = my_global_temperature_struct.pa[i_1].y;
							if (fabs(dcandidate1 - dseparate1) < min_zazor_y) bfound = true;

						}

						dcandidate1 = 0.0;
						for (int j = 0; j < 8; ++j) {
							int i_1 = my_global_temperature_struct.nvtx[j][i] - 1;
							dcandidate1 += 0.125 * my_global_temperature_struct.pa[i_1].y;
						}

						if (bfound) {
							if (inumerate[i] == 1) {
								inumerate[i] = 5;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 8;
							}
							isepnum++;
						}
						else if (dcandidate1 < dseparate1) {
							if (inumerate[i] == 1) {
								inumerate[i] = 4;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 7;
							}
						}
						else {
							if (inumerate[i] == 1) {
								inumerate[i] = 6;
							}
							if (inumerate[i] == 2) {
								inumerate[i] = 9;
							}
						}
					}
				}

			}

		}

		std::cout << "isepnum=" << isepnum << std::endl;
		for (int i = 0; i < my_global_temperature_struct.maxelm; ++i) if (inumerate[i] == -1) {
			std::cout << "nested desection error" << i << "\n"; system("pause");
		}
	}


	// 29.01.2017
	// if (1 && steady_or_unsteady_global_determinant == MESHER_ONLY)  
	// То мы просто вызываем мешер не вызывая солвера.
	if (1 && ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T))) {
		// Замер времени.
		unsigned int calculation_start_time = 0; // начало счёта мс.
		unsigned int calculation_end_time = 0; // окончание счёта мс.
		unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

		calculation_start_time = clock(); // момент начала счёта.


#ifndef NO_OPENGL_GLFW
		pa_opengl = new TOCHKA[my_global_temperature_struct.database.maxelm];
		pa_render = new TOCHKA[my_global_temperature_struct.database.maxelm];
		n_render = my_global_temperature_struct.database.maxelm;
		for (int i = 0; i < n_render; ++i) {
			pa_opengl[i].x = my_global_temperature_struct.database.x[i];
			pa_opengl[i].y = my_global_temperature_struct.database.y[i];
			pa_opengl[i].z = my_global_temperature_struct.database.z[i];

			pa_render[i].x = halfScreenWidth + scale_all * my_global_temperature_struct.database.x[i];
			pa_render[i].y = halfScreenHeight + scale_all * my_global_temperature_struct.database.y[i];
			pa_render[i].z = -abbys + scale_all * my_global_temperature_struct.database.z[i];
		}

		DrawZbuffer(b, lb, my_global_temperature_struct.database.ncell, my_global_temperature_struct.database.maxelm, my_global_temperature_struct.database.nvtxcell, my_global_temperature_struct.whot_is_block);
		delete[] pa_opengl;
		delete[] pa_render;
		n_render = -1;

#endif

		// Экспортирует внешнюю поверхность геометрии пользователя в .stl формате.
		// 09.09.2019.
		// Закоментировал 27.05.2020
		//if (steady_or_unsteady_global_determinant == MESHER_ONLY) {
		//export_User_Geom_in_STL_format(my_global_temperature_struct);
		//}

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

																					  // 24.06.2020
																					  // Графовый метод решения уравнения теплопередачи.
		if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T) {
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
					// Экспорт в двумерную программу AliceFlow2D.
					exportAliceFlow2D(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, zpos, inz, ls, lw, w);
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


	/*doublereal x51 = Tochka_position_X0_for_XY_Plot;
	doublereal y51 = Tochka_position_Y0_for_XY_Plot;
	doublereal z51 = Tochka_position_Z0_for_XY_Plot;

	doublereal dist_max_CM = 1.0e30;

	for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm; ++i7) {

		TOCHKA pointP;
		center_cord3D(i7, my_global_temperature_struct.nvtx, my_global_temperature_struct.pa, pointP, 100);

		if ((pointP.x - x51) * (pointP.x - x51) + (pointP.y - y51) * (pointP.y - y51) + (pointP.z - z51) * (pointP.z - z51) < dist_max_CM) {

			dist_max_CM = (pointP.x - x51) * (pointP.x - x51) + (pointP.y - y51) * (pointP.y - y51) + (pointP.z - z51) * (pointP.z - z51);
			iP_perefirier_start = i7;

		}
	}*/


	// steady Temperature Finite Volume Method
	if (1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE) &&
		(1 == eqin.itemper)) {






		// Замер времени.
		unsigned int calculation_start_time = 0; // начало счёта мс.
		unsigned int calculation_end_time = 0; // окончание счёта мс.
		unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

		calculation_start_time = clock(); // момент начала счёта.

#pragma omp parallel for
		for (int i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; ++i7)
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
			//for (integer i7 = 0; i7<t.maxelm + t.maxbound; ++i7) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
		}

		// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
		bsolid_static_only = true;
		bool bcleantemp = false;
		// Условие уже проверялось выше.
		//if (1 == eqin.itemper) 
		{
			bcleantemp = true;
			for (integer i = 0; i < flow_interior; ++i) {
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

			if (bmyconvective) {
				bmyconvective7248 = true;
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
				// Есть подозрения что функция calculate_color_for_temperature() устарела, поэтому я отключил данный код 27.12.2021
				// if ((!b_on_adaptive_local_refinement_mesh) && (number_cores() == 2) && (my_amg_manager.lfil < 3)) {
				// Т.е. это работало для распараллеливания на два потока исполнения в случае структурированной сетки
				// и применения BiCGStab + ILU(2) алгоритма.
				//calculate_color_for_temperature(color, my_global_temperature_struct, inx, xpos);





				// если flow_interior == 0 то f[0] просто формальный параметр  
				solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					bmyconvective, nullptr, 0.001, 0.001,
					false,
					matlist, 0,
					bprintmessage,
					gtdps, ltdp, 1.0, 1.0, // последний параметр равный 1.0 означает что мощность подаётся.
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
					my_global_temperature_struct.potent = new doublereal[static_cast<integer>(my_global_temperature_struct.maxelm) + static_cast<integer>(my_global_temperature_struct.maxbound)];
					for (int i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; ++i7) my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // инициализация.
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
				calculate_color_for_temperature(color, my_global_temperature_struct, inx, xpos);




				// если flow_interior == 0 то f[0] просто формальный параметр  
				solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					bmyconvective, nullptr, 0.001, 0.001,
					false,
					matlist, 0,
					bprintmessage,
					gtdps, ltdp, 1.0, 1.0, m,
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

#ifndef NO_OPENGL_GLFW
					pa_opengl = new TOCHKA[my_global_temperature_struct.database.maxelm];
					pa_render = new TOCHKA[my_global_temperature_struct.database.maxelm];
					n_render = my_global_temperature_struct.database.maxelm;
					for (int i = 0; i < n_render; ++i) {
						pa_opengl[i].x = my_global_temperature_struct.database.x[i];
						pa_opengl[i].y = my_global_temperature_struct.database.y[i];
						pa_opengl[i].z = my_global_temperature_struct.database.z[i];

						pa_render[i].x = halfScreenWidth + scale_all * my_global_temperature_struct.database.x[i];
						pa_render[i].y = halfScreenHeight + scale_all * my_global_temperature_struct.database.y[i];
						pa_render[i].z = -abbys + scale_all * my_global_temperature_struct.database.z[i];
					}


					doublereal dmin = 1.0e30;
					doublereal dmax = -1.0e30;

					for (int i37 = 0; i37 < my_global_temperature_struct.maxelm; i37++) {
						if (my_global_temperature_struct.potent[i37] > dmax) {
							dmax = my_global_temperature_struct.potent[i37];
						}
						if (my_global_temperature_struct.potent[i37] < dmin) {
							dmin = my_global_temperature_struct.potent[i37];
						}
					}

					minimum_val_for_render_pic = dmin;
					maximum_val_for_render_pic = dmax;

					//DrawZbufferColor(b, lb, my_global_temperature_struct.database.ncell, my_global_temperature_struct.database.maxelm, my_global_temperature_struct.database.nvtxcell, my_global_temperature_struct.whot_is_block, my_global_temperature_struct.potent);
					DrawZbufferColorLight(b, lb, my_global_temperature_struct.database.ncell, my_global_temperature_struct.database.maxelm, my_global_temperature_struct.database.nvtxcell, my_global_temperature_struct.whot_is_block, my_global_temperature_struct.potent);
					delete[] pa_opengl;
					delete[] pa_render;
					n_render = -1;

#endif

					calculation_end_time = clock(); // момент окончания счёта.
					calculation_seach_time = calculation_end_time - calculation_start_time;
					unsigned int im = 0, is = 0, ims = 0;
					im = (unsigned int)(calculation_seach_time / 60000); // минуты
					is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
					ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

					printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);
					printf("export to tecplot start...\n");

					if (pause_suppression == 0)
					{
						// Только не вслучае оптимизации.
						// 

					    // экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
					}
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
						printf("export to tecplot start...\n");

						//printf("Would you like to save the result on the ALICE grid ? y/n\n");
						//ch_EXPORT_ALICE_ONLY = getchar(); // Здесь именно system("pause");
						//ch_EXPORT_ALICE_ONLY = 'y';
					}

					if (ch_EXPORT_ALICE_ONLY == 'y') {
						// Экспорт в программу tecplot температуры.
						//С АЛИС сетки.
						if (pause_suppression == 0)
						{
							// Только не вслучае оптимизации.
							// 

							ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
						}
					}

				}
			}

		}

		doublereal tmaxfinish = -273.15; // абсолютный ноль.
										 // Вычисление значения максимальной температуры внутри расчётной области и на её границах:
										 //for (integer i = 0; i < t.maxelm + t.maxbound; ++i) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 декабря 2015
										 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
										 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

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

				if (pause_suppression == 0)
				{
					// Только не вслучае оптимизации.
					// 

					ANES_tecplot360_export_temperature_preobrazovatel(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx, operatingtemperature);

				}

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
						for (integer i = 0; i <= 7; ++i) {
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
				my_global_temperature_struct.xpos_copy = new doublereal[static_cast<integer>(inx) + 1];
				my_global_temperature_struct.ypos_copy = new doublereal[static_cast<integer>(iny) + 1];
				my_global_temperature_struct.zpos_copy = new doublereal[static_cast<integer>(inz) + 1];
				// Данная информация нужна для экономии оперативной памяти,
				// некоторые данные будут выгружены из озу а потом восстановлены 
				// путём вычиcления.
				for (integer i_7 = 0; i_7 < static_cast<integer>(inx) + 1; i_7++) {
					my_global_temperature_struct.xpos_copy[i_7] = xpos[i_7];
				}
				for (integer i_7 = 0; i_7 < static_cast<integer>(iny) + 1; i_7++) {
					my_global_temperature_struct.ypos_copy[i_7] = ypos[i_7];
				}
				for (integer i_7 = 0; i_7 < static_cast<integer>(inz) + 1; i_7++) {
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
					for (integer i_1 = 0; i_1 < 8; ++i_1) {
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
		for (int i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; ++i7) {
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
			//for (integer i7 = 0; i7<t.maxelm + t.maxbound; ++i7) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
		}

		// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
		bsolid_static_only = true;
		bool bcleantemp = false;
		// Мы тоже решаем задачу чистой теплопроводности только
		// вторым температурным солвером.
		//if (2 == eqin.itemper)
		{
			bcleantemp = true;
			for (integer i = 0; i < flow_interior; ++i) {
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
			//bool bprintmessage = true; // печатать ли сообщения на консоль.

			//doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
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
			//doublereal** rhie_chow = nullptr;
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
			//02.11.2021
			// Эти массивы используются во втором температурном солвере.
			m.val = nullptr;
			m.row_ptr = nullptr;
			m.col_ind = nullptr;

			bPhysics_stop = false;

			// Температура 19.05.2018


			doublereal* lstub = nullptr;
			integer maxelm_global_ret = 0;
			doublereal* t_for_Mechanical = nullptr;
			solve_Thermal(my_global_temperature_struct, f, matlist, w, lw, lu, b, lb, ls, m, false,
				operatingtemperature, false, 0.0, lstub, lstub,
				maxelm_global_ret, 1.0, 1.0, bAVLrealesation, t_for_Mechanical, inumerate, true, -1);


			if (t_for_Mechanical != nullptr) {
				delete[] t_for_Mechanical;
			}
			t_for_Mechanical = nullptr;

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

			doublereal mass_s = 0.0;
			// Вычисление массы модели.
			mass_s += massa_cabinet(my_global_temperature_struct, f, flow_interior,
				b, lb, operatingtemperature,
				matlist);
			for (int i_37 = 0; i_37 < lu; i_37++) {
				if (my_union[i_37].active) {
					std::cout << "union " << i_37 << std::endl;
					mass_s += massa_cabinet(my_union[i_37].t, my_union[i_37].f, flow_interior,
						b, lb, operatingtemperature,
						matlist);
				}
			}
			std::cout << "massa=" << mass_s << " kg" << std::endl;

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
										 //for (integer i = 0; i < t.maxelm + t.maxbound; ++i) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 декабря 2015
										 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
										 // физичнее не смущать людей и приводить температуру только во внутренних КО. 

		// Неверно. Данная максимальная температура не уитывает асемблесы.
		for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		doublereal totaldeform_max = -1.0e+30;
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
		}

		//FILE* fp = NULL;

//#ifdef MINGW_COMPILLER
		//int err1 = 0;
		//fp = fopen64("report.txt", "w");
		//if (fp == NULL) err1 = 1;
//#else
		//errno_t err1 = 0;
		//err1 = fopen_s(&fp, "report.txt", "w");
//#endif
		// создание файла для записи.
		//if ((err1) != 0) {
			//printf("Create File report.txt Error\n");
			//system("PAUSE");
			//system("pause");
		//}
		//else {
			// запись заголовка
			// Не учитывает наличие асемблесов здесь. 
			// Печать перенесена внутрь функции solve_Thermal.
			//fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
			//fclose(fp);
		//}
		// 1 - solver/solid_static/
		bool bMechanical = false;
		if (eqin.itemper != 2) {
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
		}
		else {



			printf("THIS IS SECOND STEADY TEMPERATURE SOLVER ON ALL MESHES.\n");
			printf("NO EXPOPRT TECPLOT.\n");
			printf("NO PRINT REPORT.\n");
		}


	}

	// steady Static Structural
	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) {


		// Замер времени.
		unsigned int calculation_start_time = 0; // начало счёта мс.
		unsigned int calculation_end_time = 0; // окончание счёта мс.
		unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

		calculation_start_time = clock(); // момент начала счёта.







#pragma omp parallel for
		for (int i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; ++i7) {
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
			//bool bprintmessage = true; // печатать ли сообщения на консоль.

			//doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
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
				false, 1.0e9, 1.0e9, 1.0e9, 1.0e9, stub_Mechanical, stub_Mechanical, stub_Mechanical, stub_Mechanical, 1.0, matlist, stub_Mechanical, inumerate);
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
										 //for (integer i = 0; i < t.maxelm + t.maxbound; ++i) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 декабря 2015
										 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
										 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		doublereal totaldeform_max = -1.0e+30;
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
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

		for (int i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; ++i7)
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
			//for (integer i7 = 0; i7<t.maxelm + t.maxbound; ++i7) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
		}

		// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
		bsolid_static_only = true;
		bool bcleantemp = false;
		if (eqin.itemper == 1) {
			bcleantemp = true;
			for (integer i = 0; i < flow_interior; ++i) {
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
			calculate_color_for_temperature(color, my_global_temperature_struct, inx, xpos);


			// если flow_interior == 0 то f[0] просто формальный параметр
			solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				bmyconvective, nullptr, 0.001, 0.001,
				false,
				matlist, 0,
				bprintmessage,
				gtdps, ltdp, 1.0, 1.0, // последний параметр равный 1.0 означает что мощность подаётся.
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
			//maxelm_global_ret, 1.0, bAVLrealesation, t_for_Mechanical, -1);

			t_for_Mechanical = new doublereal[static_cast<integer>(my_global_temperature_struct.maxnod) + 2];
			{
				// Метод линейного порядка.
				doublereal min_x = 1e60;
				doublereal min_y = 1e60;
				doublereal min_z = 1e60;
				doublereal max_x = -1e60;
				doublereal max_y = -1e60;
				doublereal max_z = -1e60;

				for (int i = 0; i < my_global_temperature_struct.maxnod; ++i) {
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

				for (integer i = 0; i < my_global_temperature_struct.maxnod; ++i) {
					vol[i] = 0.0;
				}

				// Преобразование температуры с сетки МКО на сетку МКЭ.
				SECOND_ORDER_QUADRATIC_RECONSTRUCTA(my_global_temperature_struct.maxnod,
					my_global_temperature_struct.maxelm, my_global_temperature_struct.pa,
					my_global_temperature_struct.nvtx, vol, t_for_Mechanical, min_x, min_y, min_z, my_global_temperature_struct.potent,
					my_global_temperature_struct, eps_mashine, false, 1.0e-2);

				delete[] vol;

			}
			//bPhysics_stop = false;
			// Вызов солвера Static Structural.
			doublereal* stub_Mechanical = nullptr;
			solve_Structural(my_global_temperature_struct, w, lw, true, operatingtemperature, b, lb, lu,
				false, 1.0e9, 1.0e9, 1.0e9, 1.0e9, stub_Mechanical, stub_Mechanical, stub_Mechanical, stub_Mechanical, 1.0, matlist, t_for_Mechanical, inumerate);
			//bPhysics_stop = true;

			delete[] t_for_Mechanical;

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
										 //for (integer i = 0; i < t.maxelm + t.maxbound; ++i) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 декабря 2015
										 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
										 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		doublereal totaldeform_max = -1.0e+30;
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			for (integer i = 0; i < my_global_temperature_struct.maxelm; ++i) {
				totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
			}
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
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
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
		//system("pause");

		//solve_nonlinear_temp(f[0], f, t, rhie_chow, b, lb, s, ls, w, lw, dbeta, flow_interior, false, nullptr, 0.001, false);
		bool bsecond_T_solver = false;
		//if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) 
		if (eqin.itemper == 2)
		{
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

			if (f != nullptr) {
				unsteady_temperature_calculation(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					matlist,
					operatingtemperature,
					gtdps, ltdp, lu, my_union, bsecond_T_solver, inx, xpos, bTemperature, bMechanical, inumerate); // нестационарный температурный солвер
			}
			else {
				std::cout << " pointer f (fluid) is nullptr, can not start unsteady_temperature_calculation().\n";
				system("pause");
			}

		}

		// Вычисляет среднюю температуру в жидких ячейках.
		doublereal Tavg_fluid = 0.0;
		doublereal ic_avg_Temp_fluid = 0.0;
		if (f != nullptr) {
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
		}
		else {
			std::cout << "not found fluid cell node\n.";
		}

		// Вычисление массы модели.
		doublereal mass_s = 0.0;
		mass_s += massa_cabinet(my_global_temperature_struct, f, flow_interior,
			b, lb, operatingtemperature,
			matlist);

		for (int i_37 = 0; i_37 < lu; i_37++) {
			if (my_union[i_37].active) {
				std::cout << "union " << i_37 << std::endl;
				mass_s += massa_cabinet(my_union[i_37].t, my_union[i_37].f, flow_interior,
					b, lb, operatingtemperature,
					matlist);
			}
		}
		std::cout << "massa=" << mass_s << " kg" << std::endl;

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
			for (int i = 0; i < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; ++i) {
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

		// Печатает расход жидкости через выходную границу потока.
		// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
		// проходящий через выходную границу потока. 28,10,2019
		report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);
		//system("pause");

		told_temperature_global_for_HOrelax = new doublereal[static_cast<integer>(my_global_temperature_struct.maxelm) + static_cast<integer>(my_global_temperature_struct.maxbound)];
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

		const integer ISIZEf = static_cast<integer>(f[0].maxelm) + static_cast<integer>(f[0].maxbound);

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
			gtdps, ltdp, bextendedprint, lu, my_union, inx, xpos);
		// xyplot( f, 0, t);
		// boundarylayer_info(f, t, flow_interior, w, lw);
		// 2 - solver/conjugate_heat_transfer_static/
		bool bMechanical = false;

		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0/*2*/, matlist, bMechanical);
		// Печатает расход жидкости через выходную границу потока.
		// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
		// проходящий через выходную границу потока. 28,10,2019
		report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);

		{
			// Код занесения информации в Trials.txt который используется в параметрике setup.

		   doublereal power_diss_message_06_10_2018 = 0.0;
		   doublereal balancet = 0.0;
		   {

			bool bprintmessage = true;
			doublereal pdiss = 0.0;

#pragma omp parallel for reduction(+:pdiss)
			for (integer i = 0; i < ls; ++i) {
				if (s[i].power < 0.0) {
					//printf("warning source [%lld] is negative power = %e\n",i, s[i].power);
					std::cout << "warning source [" << i << "] is negative power = " << s[i].power << std::endl;
				}
				pdiss += s[i].power;
			}
			//for (integer i = 0; i < lb; ++i) {
				//pdiss += b[i].Sc*(fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS));
			//}
			// 19 november 2016.
			// Обновление мощности тепловыделения во всех внутренних узлах.

#pragma omp parallel for reduction(+:pdiss)
			for (integer i47 = 0; i47 < my_global_temperature_struct.maxelm; i47++) {
				// Скорость в том что значение не вычисляется как раньше а просто хранится.
				integer ib = my_global_temperature_struct.whot_is_block[i47];

				my_global_temperature_struct.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, my_global_temperature_struct.potent[i47]);
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
				volume3D(i47, my_global_temperature_struct.nvtx, my_global_temperature_struct.pa, dx, dy, dz);
				if (my_global_temperature_struct.Sc[i47] * dx * dy * dz < 0.0) {
					//printf("ERROR!!!  control volume [%lld] is negative power = %e\n", i47, my_global_temperature_struct.Sc[i47] * dx*dy*dz);
					std::cout << "ERROR!!!  control volume [" << i47 << "] is negative power =" << (my_global_temperature_struct.Sc[i47] * dx * dy * dz) << std::endl;
					//system("PAUSE");
				}
				pdiss += my_global_temperature_struct.Sc[i47] * dx * dy * dz;
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
			if (bprintmessage) {
				std::cout << "power generation is equal=" << pdiss << " W" << std::endl;
			}
			if (fabs(d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL - pdiss)
				> 0.01 * d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL) {

				// Корректируем значение тепловой мощности так чтобы она в 
				// точности была равна заданной пользователем.
				doublereal m_power = 1.0;
				if (fabs(pdiss) > 1.0e-30) {
					m_power = d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL / pdiss;
					if (bprintmessage) {
						printf("Correction power: m_power=%e\n", m_power);
						printf("Apriority power %e calculate real power %e\n", d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL, pdiss);
					}
					// m_power - множитель на который нужно домнодить ненулевую внутрипрограммную
					// тепловую мощность, чтобы в точности получить заданную пользователем мощность.
					pdiss = 0.0;

#pragma omp parallel for reduction(+:pdiss)
					for (integer i = 0; i < ls; ++i) {
						if (s[i].power < 0.0) {
							//printf("warning source [%lld] is negative power = %e\n",i, s[i].power);
							std::cout << "warning source [" << i << "] is negative power = " << s[i].power << std::endl;
						}
						s[i].power *= m_power;
						pdiss += s[i].power;
					}
					//for (integer i = 0; i < lb; ++i) {
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
								//system("pause");
							}

						}
					}

#pragma omp parallel for reduction(+:pdiss)
					for (integer i47 = 0; i47 < my_global_temperature_struct.maxelm; i47++) {
						// Скорость в том что значение не вычисляется как раньше а просто хранится.
						integer ib = my_global_temperature_struct.whot_is_block[i47];



						/*if (my_global_temperature_struct.Sc[i47] > 0.0) {
							printf("%e ", t.Sc[i47]);
							doublereal t12= get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, my_global_temperature_struct.potent[i47]);
							printf("%e ", t12);
							system("pause");
						}*/
						my_global_temperature_struct.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, my_global_temperature_struct.potent[i47]);
						// вычисление размеров текущего контрольного объёма:
						doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
						volume3D(i47, my_global_temperature_struct.nvtx, my_global_temperature_struct.pa, dx, dy, dz);
						if (my_global_temperature_struct.Sc[i47] * dx * dy * dz < 0.0) {
							//printf("ERROR!!!  control volume [%lld] is negative power = %e\n", i47, t.Sc[i47] * dx*dy*dz);
							std::cout << "ERROR!!!  control volume [" << i47 << "] is negative power =" << (my_global_temperature_struct.Sc[i47] * dx * dy * dz) << std::endl;
							//system("PAUSE");
						}
						pdiss += my_global_temperature_struct.Sc[i47] * dx * dy * dz;
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
					std::cout << "Real power dissipation = " << pdiss << std::endl;
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
			if (pdiss > 0.0) {
				doublereal square_bolc = 0.0;
				doublereal emissivity = 1.0;

#pragma omp parallel for reduction(+:square_bolc) 
				for (integer i = 0; i < lw; ++i) {
					if (w[i].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
						switch (w[i].iPlane) {
						case XY_PLANE: square_bolc += fabs(w[i].g.xE - w[i].g.xS) * fabs(w[i].g.yE - w[i].g.yS); break;
						case XZ_PLANE: square_bolc += fabs(w[i].g.xE - w[i].g.xS) * fabs(w[i].g.zE - w[i].g.zS); break;
						case YZ_PLANE: square_bolc += fabs(w[i].g.yE - w[i].g.yS) * fabs(w[i].g.zE - w[i].g.zS); break;
						}
						// Здесь мы предполагаем что на всех излучающих поверхностях излучающая способность одна и таже.
						// Если это не так то возникнет ошибка.
					}
				}
				for (integer i = 0; i < lw; ++i) {
					if (w[i].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
						emissivity = w[i].emissivity;
					}
				}

				if (fabs(square_bolc) > 1e-23) {
					//printf("Pdiss=%e, S=%e\n",pdiss, square_bolc);
					std::cout << "Pdiss=" << pdiss << ", S=" << square_bolc << std::endl;
					balancet = sqrt(sqrt((pdiss / (square_bolc * STEFAN_BOLCMAN_CONST * emissivity)))) - 273.15;
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
				for (integer i = 0; i < lw; ++i) {
					if (w[i].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
						// условие идеального теплоотвода обнаружено.

						bDirichlet = bDirichlet || true;

					}
				}
				if (!bDirichlet) {
					printf("negative power and the lack of Dirichlet conditions \n");
					//system("pause");
					balancet = -272.15;
				}
			}
			power_diss_message_06_10_2018 = pdiss;
		}

		   if (fabs(power_diss_message_06_10_2018) > 1.0e-20) {

			doublereal tmin = get_min_array_elm(my_global_temperature_struct.potent, my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound);
			doublereal tmax = get_max_array_elm(my_global_temperature_struct.potent, my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound);

			bool bprintmessage = true;

			if (bprintmessage) {
				printf("Intermediate maximum temperature in default interior\n");
				//printf("is equal %e  \u00B0C.\n", tmaxloc);
				std::cout << "is equal " << tmax << "  \u00B0C." << std::endl;
			}

			if (bprintmessage) {
				printf("Intermediate minimum temperature in default interior\n");
				//printf("is equal %e  \u00B0C.\n", tminloc);
				std::cout << "is equal " << tmin << "  \u00B0C." << std::endl;
			}

			std::cout << "Thermal resistance =" << (tmax - tmin) / power_diss_message_06_10_2018 << " \u00B0C/W." << std::endl;

			FILE* fp_trials = NULL;

			// создание файла для добавления.
#ifdef MINGW_COMPILLER
			fp_trials = fopen64("Trials.txt", "a");
			int err_trials = 0;
			if (fp_trials != NULL) {
				err_trials = 0;
			}
			else {
				err_trials = 1; // ошибка открытия.
			}
#else
			errno_t err_trials;
			err_trials = fopen_s(&fp_trials, "Trials.txt", "a");
#endif

			if (fp_trials != NULL) {
				if (err_trials != 0) {
					//printf("Open File add_line Error\n");
					//system("pause");
					//exit(0);
					// return bfound;
				}
				else
				{
					fprintf(fp_trials, "%e \n", (tmax - tmin) / power_diss_message_06_10_2018);
					fclose(fp_trials);
				}
			}
		}
	    }

		if (pause_suppression == 0)
		{
			// Только не вслучае оптимизации.
			// 

		    // Вычисление массы модели.
			massa_cabinet(my_global_temperature_struct, f, flow_interior,
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

		}

		save_velocity_for_init(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior);

		// exporttecplotxy360T_3D_part2_rev(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,b,lb);
		delete[] bPamendment_source_old;
		bPamendment_source_old = nullptr;
		delete[] told_temperature_global_for_HOrelax;
		told_temperature_global_for_HOrelax = nullptr;
	}


	// нестационарный гидродинамический решатель: 
	if ((1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY))) {


#ifdef _OPENMP
		// Предупреждение о невозможности расчёта cfd на openMP 08.05.2019.
		//printf("Unsteady CFD not work in OPENMP ON and bparallelismold is true.\n");
		//printf("uskorenie ot OPENMP otsutstvuet. Rabotaen odnopotochnaq versiq.\n");
		//printf("variable bparallelismold must be equal false.\n");
		//system("PAUSE");
#endif

		told_temperature_global_for_HOrelax = new doublereal[static_cast<integer>(my_global_temperature_struct.maxelm) + static_cast<integer>(my_global_temperature_struct.maxbound)];
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

		const integer ISIZEf = static_cast<integer>(f[0].maxelm) + static_cast<integer>(f[0].maxbound);
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

		if (0) {
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
				w, lw, matlist, gtdps, ltdp, bextendedprint, lu, my_union, inx, xpos);
		}

		//xyplot( f, 0, t);
		// boundarylayer_info(f, t, flow_interior, w, lw);
		// 2 - solver/conjugate_heat_transfer_static/
		bool bMechanical = false;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0/*2*/, matlist, bMechanical);

		// Вычисление массы модели.
		massa_cabinet(my_global_temperature_struct, f, flow_interior,
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


	delete[] inumerate;

	if (continity_start != nullptr) {
		delete[] continity_start;
		continity_start = nullptr;
	}

	if (inumber_iteration_SIMPLE != nullptr) {
		delete[] inumber_iteration_SIMPLE;
		inumber_iteration_SIMPLE = nullptr;
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


	// Быстрая обработка нелинейных граничных условий в Румба v. 0.14 решателе.
	if (qnbc != nullptr) {
		delete[] qnbc;
		qnbc = nullptr;
		iadd_qnbc_maxelm = 0;
	}



	/*
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
	delete[] b;
	b = nullptr;
	*/

	delete[] s; delete[] w; // освобождение памяти
	s = nullptr;
	w = nullptr;

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

	free_global_1(matlist, lmatmax, x_jacoby_buffer, bvery_big_memory, my_global_temperature_struct, amgGM, xposadd, yposadd, zposadd);



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
				for (integer i = 0; i <= 7; ++i) {
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

		// 26.11.2021
		if (my_union[i63].t.xpos_copy != nullptr) {
			delete[] my_union[i63].t.xpos_copy;
			my_union[i63].t.xpos_copy = nullptr;
		}
		if (my_union[i63].t.ypos_copy != nullptr) {
			delete[] my_union[i63].t.ypos_copy;
			my_union[i63].t.ypos_copy = nullptr;
		}
		if (my_union[i63].t.zpos_copy != nullptr) {
			delete[] my_union[i63].t.zpos_copy;
			my_union[i63].t.zpos_copy = nullptr;
		}

	}

	// Освобождение памяти из под UNIONов.
	delete[] my_union;
	my_union = nullptr;

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
	if (shorter_hash_X != nullptr) {
		delete[] shorter_hash_X;
		shorter_hash_X = nullptr;
	}
	if (shorter_hash_Y != nullptr) {
		delete[] shorter_hash_Y;
		shorter_hash_Y = nullptr;
	}
	if (shorter_hash_Z != nullptr) {
		delete[] shorter_hash_Z;
		shorter_hash_Z = nullptr;
	}
	if (bshorter_hash_X != nullptr) {
		delete[] bshorter_hash_X;
		bshorter_hash_X = nullptr;
	}
	if (bshorter_hash_Y != nullptr) {
		delete[] bshorter_hash_Y;
		bshorter_hash_Y = nullptr;
	}
	if (bshorter_hash_Z != nullptr) {
		delete[] bshorter_hash_Z;
		bshorter_hash_Z = nullptr;
	}

	// Освобождаем память из под массива строк.
	//delete[] StringList;

	calculation_main_end_time = clock();
	calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;


	/*printf("time=%d statistic vorst=%3.2f %% \n",calculation_main_seach_time,static_cast<float>(100.0*calculation_vorst_seach_time/calculation_main_seach_time));
	system("PAUSE");
	*/


	if (b_setup_CathilMC_Temp) {

		b_setup_CathilMC_Temp = false;

		delete[] new_number_CathilMC_Temp;
		delete[] new_number_internal_CathilMC_Temp;
		delete[] new_number_bound_CathilMC_Temp;
		delete[] rev_number_CathilMC_Temp;
	}


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

int main_solverv0_48(BLOCK*& b, int& lb)
{

	// Поддержка кириллицы в консоли Windows
	setlocale(LC_CTYPE, "rus"); // вызов функции настройки локали
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
			time_optimetric[iscan_optimetric] = main_body(b, lb, ch_EXPORT_ALICE_ONLY);
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
			time_optimetric[iscan_optimetric] = main_body(b, lb, ch_EXPORT_ALICE_ONLY);
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
		main_body(b, lb, ch_EXPORT_ALICE_ONLY);
	}

	// Освобождаем память из под массива строк.
	delete[] StringList;

	if (pause_suppression == 0)
	{
		// Не optimetric значит делаем паузу в конце расчёта.
		system("pause");
	}
	return 0;
}

#ifndef NO_OPENGL_GLFW

void DrawPrism(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawCylinder(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawPolygon(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawCADobj(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawKarkas(BLOCK*& b, integer lb, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {

	bool b_flag = true;
	for (integer i = 0; i < lb; ++i) {
		if (b_flag) {
			if (binverse_color_black_2_white) {
				glColor3f(0.0, 0.0, 0.0); // чёрный
			}
			else {
				glColor3f(1.0, 1.0, 1.0); // белый
			}
			b_flag = false;
		}
		if (b[i].g.itypegeom == PRISM) {
			DrawPrism(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
		}
		if (b[i].g.itypegeom == CYLINDER) {
			DrawCylinder(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
		}
		if (b[i].g.itypegeom == POLYGON) {
			DrawPolygon(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
		}
		if (b[i].g.itypegeom == CAD_STL)
		{
			glColor3f(0.502, 0.502, 0.502); // серый для CAD_STL obj
			DrawCADobj(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
			b_flag = true;
		}
	}

}

#endif

int main(void)
{




#ifndef NO_OPENGL_GLFW
	int icol_now = 0;
	mask_color_for_openGL = new int[1021];
	for (int i = 0; i < 1021; ++i) {
		mask_color_for_openGL[i] = icol_now;
		if (i % 72 == 0) icol_now = i;// Всего 14 различных цветов.
	}
#endif

	// количество блоков, источников и стенок, юнионов.
	int lb = 0;
	BLOCK* b = nullptr;// список блоков

					   // Вызываем солвер.
	main_solverv0_48(b, lb);

#ifndef NO_OPENGL_GLFW

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 

					  //GLfloat halfScreenWidth = SCREEN_WIDTH / 2;
					  //GLfloat halfScreenHeight = SCREEN_HEIGHT / 2;

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_LIGHTING); // рассеяный свет.
	glEnable(GL_LIGHT0); // источник света номера 0.
	glEnable(GL_COLOR_MATERIAL); // Учитываются цвета рисуемых объектов.

								 /* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{

		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
		}
		//glClearColor(0.7f, 1.0f, 0.7f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();

		float position[] = { 0,0,1,0 };
		glLightfv(GL_LIGHT0, GL_POSITION, position);

		glTranslatef(halfScreenWidth, halfScreenHeight, -500);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, 500);


		DrawKarkas(b, lb, halfScreenWidth, halfScreenHeight, -500.0, scale_all);


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

#endif

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
	delete[] b;
	b = nullptr;

#ifndef NO_OPENGL_GLFW
	delete[] mask_color_for_openGL;
#endif

	return 0;
}



#ifndef NO_OPENGL_GLFW

// 16.12.2020
// Соединять ли две точки сетки чёрной линией при прорисовке.
// Каркасная визуализация геометрии с помощью алгоритма Z буфера.
bool CheckLine1(int il1, int il2)
{
	if ((il1 != 4) &&
		(il1 != 8) &&
		(il2 != 4) &&
		(il2 != 8) &&
		(il1 == il2))
	{
		return true;
	}
	else
	{
		if (((il1 == 1) && ((il2 == 2) || (il2 == 3) || (il2 == 5)))
			|| ((il2 == 1) && ((il1 == 2) || (il1 == 3) || (il1 == 5))))
		{
			return true;
		}
		else if ((il1 == 2) && (il2 == 2)) {
			return true;
		}
		else if (((il1 == 3) && (il2 == 2)) ||
			((il2 == 3) && (il1 == 2))) {

			return true;
		}
		else if (((il1 == 7) && (il2 == 5)) ||
			((il2 == 7) && (il1 == 5)))
		{
			return true;
		}
		else if (((il1 == 7) && (il2 == 3)) ||
			((il2 == 7) && (il1 == 3)))
		{
			return  true;
		}
		//if (((il1==6)&&(il2==3))||
		//((il2==6)&&(il1==3))) 
		//{
		// return true;
		//}
		//else
		else if (((il1 == 6) && (il2 == 5)) ||
			((il2 == 6) && (il1 == 5)))
		{
			return  true;
		}
		else
		{
			if (((il1 == 6) && (il2 == 7)) ||
				((il2 == 6) && (il1 == 7)))
			{
				return true;
			}
			else
			{
				return  false;
			}
		}

	}
	return false;
} // CheckLine

  // Ускоренная проверка прорисовки линии с помощью хеш таблицы.
const int CHECK_LINE_SIZE = 9;
bool* checkline1 = nullptr, * checkline2 = nullptr;
bool CheckLine(int il1, int il2)
{
	int key = il1 * CHECK_LINE_SIZE + il2;
	// Для ускорения выполнения хеш таблица насчитывается заранее один раз.
	//if (checkline1[key]) {
	return checkline2[key];
	//}
	/*else {
	checkline1[key] = true;
	checkline2[key] = CheckLine1(il1, il2);
	return checkline2[key];
	}
	*/
}

//  оптимизировано быстродействие в ущерб общности и понятности 18.12.2020.
// Данная функция многократно вызывается при рендере, а рендер критичен к быстродействию.
bool binvisible_face_detect(doublereal nx, doublereal ny, doublereal nz)
{
	//if (0) {
	//return true;
	//}
	//else {


	//bool br;
	//GLfloat mx[3][3];
	//GLfloat my[3][3];
	//doublereal mz0[3][3];
	//doublereal mz[3][3];
	//GLfloat mr[3][3];
	//doublereal mr1[3][3];

	//int i88, j88, k88;
	GLfloat /*nx1, ny1,*/ nz1;

	// (nx,ny,nz) - нормаль.
	//
	//br = true;
	//

	//sinAlf, cosAlf, sinBet, cosBet - вычисляются один раз.

	// glRotatef(180.0*Gam0/3.141,0.0,0.0,1.0); // z apriory
	//glRotatef(180.0*Alf/3.141,1.0,0.0,0.0); // x
	//glRotatef(180.0*Bet/3.141,0.0,1.0,0.0); // y
	//glRotatef(180.0*Gam/3.141,0.0,0.0,1.0); // z
	// Матрица поворота вокруг оси Oz
	//mz0[0][0]: = cos(Gam0); mz0[0][1]: = -sin(Gam0); mz0[0][2]: = 0.0;
	//mz0[1][0]: = sin(Gam0); mz0[1][1]: = cos(Gam0); mz0[1][2]: = 0.0;
	//mz0[2][0]: = 0.0; mz0[2][1]: = 0.0; mz0[2][2]: = 1.0;
	/*
	mz0[0][0] = 1.0; mz0[0][1] = 0.0; mz0[0][2] = 0.0;
	mz0[1][0] = 0.0; mz0[1][1] = 1.0; mz0[1][2] = 0.0;
	mz0[2][0] = 0.0; mz0[2][1] = 0.0; mz0[2][2] = 1.0;
	*/
	// Матрица поворота вокруг оси Ox
	//mx[0][0] = 1.0; mx[0][1] = 0.0; mx[0][2] = 0.0;
	//mx[1][0] = 0.0; mx[1][1] = cosAlf; mx[1][2] = -sinAlf;
	//mx[2][0] = 0.0; mx[2][1] = sinAlf; mx[2][2] = cosAlf;
	// Матрица поворота вокруг оси Oy
	//my[0][0] = cosBet; my[0][1] = 0.0; my[0][2] = sinBet;
	//my[1][0] = 0.0; my[1][1] = 1.0; my[1][2] = 0.0;
	//my[2][0] = -sinBet; my[2][1] = 0.0; my[2][2] = cosBet;
	// Матрица поворота вокруг оси Oz
	//mz[0][0] = cos(Gam); mz[0][1] = -sin(Gam); mz[0][2] = 0.0;
	//mz[1][0] = sin(Gam); mz[1][1] = cos(Gam); mz[1][2] = 0.0;
	//mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;

	/*
	mz[0][0] = 1.0; mz[0][1] = 0.0; mz[0][2] = 0.0;
	mz[1][0] = 0.0; mz[1][1] = 1.0; mz[1][2] = 0.0;
	mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;
	*/

	//mr[0][0] = 0.0; mr[0][1] = 0.0; mr[0][2] = 0.0;
	//mr[1][0] = 0.0; mr[1][1] = 0.0; mr[1][2] = 0.0;
	/*mr[2][0] = 0.0; mr[2][1] = 0.0; mr[2][2] = 0.0;


	//for (i88 = 0; i88 < 3; i88++)
	i88 = 2; // Нужна только Z компонента.
	{
	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr[i88][j88] = mr[i88][j88] + mx[i88][k88] * my[k88][j88];
	}
	}
	}*/
	/*
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr[i88][j88] = mr[i88][j88] + mz0[i88][k88] * mx[k88][j88];
	}
	}
	}

	mr1[0][0] = 0.0; mr1[0][1] = 0.0; mr1[0][2] = 0.0;
	mr1[1][0] = 0.0; mr1[1][1] = 0.0; mr1[1][2] = 0.0;
	mr1[2][0] = 0.0; mr1[2][1] = 0.0; mr1[2][2] = 0.0;
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * my[k88][j88];
	}
	}
	}
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	mr[i88][j88] = mr1[i88][j88];
	mr1[i88][j88] = 0.0;
	}
	}
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * mz[k88][j88];
	}
	}
	}
	nx1 = 0.0;
	ny1 = 0.0;
	nz1 = 0.0;
	nx1 = mr1[0][0] * nx + mr1[0][1] * ny + mr1[0][2] * nz;
	ny1 = mr1[1][0] * nx + mr1[1][1] * ny + mr1[1][2] * nz;
	nz1 = mr1[2][0] * nx + mr1[2][1] * ny + mr1[2][2] * nz;
	*/
	//nx1 = 0.0;
	//ny1 = 0.0;
	//nz1 = 0.0;
	//nx1 = mr[0][0] * nx + mr[0][1] * ny + mr[0][2] * nz;
	//ny1 = mr[1][0] * nx + mr[1][1] * ny + mr[1][2] * nz;
	//nz1 = mr[2][0] * nx + mr[2][1] * ny + mr[2][2] * nz;

	// Всё вычисляется вручную, ни одного лишнего действия.
	// При этом синусы и косинусы посчитаны заранее и запомнены.
	// Данная функция вызывается колосальное число раз и быстродействие очень важно.
	// 18.12.2020.
	nz1 = -sinBet * cosAlf * nx + sinAlf * ny + cosAlf * cosBet * nz;

	//if (nx*nx1+ny*ny1+nz*nz1>-1.0e-3)
	if (nz1 > -1.0e-3)
	{
		return true;
	}
	else
	{
		return false;
	}


	//}
}

//  работает медленно.
bool binvisible_face_detectold(doublereal nx, doublereal ny, doublereal nz)
{
	if (0) {
		return true;
	}
	else {


		bool br;
		doublereal mx[3][3];
		doublereal my[3][3];
		doublereal mz0[3][3];
		doublereal mz[3][3];
		doublereal mr[3][3];
		doublereal mr1[3][3];

		int i88, j88, k88;
		doublereal nx1, ny1, nz1;

		// (nx,ny,nz) - нормаль.
		//
		br = true;
		//

		//sinAlf, cosAlf, sinBet, cosBet - вычисляются один раз.

		// glRotatef(180.0*Gam0/3.141,0.0,0.0,1.0); // z apriory
		//glRotatef(180.0*Alf/3.141,1.0,0.0,0.0); // x
		//glRotatef(180.0*Bet/3.141,0.0,1.0,0.0); // y
		//glRotatef(180.0*Gam/3.141,0.0,0.0,1.0); // z
		// Матрица поворота вокруг оси Oz
		//mz0[0][0]: = cos(Gam0); mz0[0][1]: = -sin(Gam0); mz0[0][2]: = 0.0;
		//mz0[1][0]: = sin(Gam0); mz0[1][1]: = cos(Gam0); mz0[1][2]: = 0.0;
		//mz0[2][0]: = 0.0; mz0[2][1]: = 0.0; mz0[2][2]: = 1.0;

		mz0[0][0] = 1.0; mz0[0][1] = 0.0; mz0[0][2] = 0.0;
		mz0[1][0] = 0.0; mz0[1][1] = 1.0; mz0[1][2] = 0.0;
		mz0[2][0] = 0.0; mz0[2][1] = 0.0; mz0[2][2] = 1.0;
		// Матрица поворота вокруг оси Ox
		mx[0][0] = 1.0; mx[0][1] = 0.0; mx[0][2] = 0.0;
		mx[1][0] = 0.0; mx[1][1] = cosAlf; mx[1][2] = -sinAlf;
		mx[2][0] = 0.0; mx[2][1] = sinAlf; mx[2][2] = cosAlf;
		// Матрица поворота вокруг оси Oy
		my[0][0] = cosBet; my[0][1] = 0.0; my[0][2] = sinBet;
		my[1][0] = 0.0; my[1][1] = 1.0; my[1][2] = 0.0;
		my[2][0] = -sinBet; my[2][1] = 0.0; my[2][2] = cosBet;
		// Матрица поворота вокруг оси Oz
		//mz[0][0] = cos(Gam); mz[0][1] = -sin(Gam); mz[0][2] = 0.0;
		//mz[1][0] = sin(Gam); mz[1][1] = cos(Gam); mz[1][2] = 0.0;
		//mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;

		mz[0][0] = 1.0; mz[0][1] = 0.0; mz[0][2] = 0.0;
		mz[1][0] = 0.0; mz[1][1] = 1.0; mz[1][2] = 0.0;
		mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;

		mr[0][0] = 0.0; mr[0][1] = 0.0; mr[0][2] = 0.0;
		mr[1][0] = 0.0; mr[1][1] = 0.0; mr[1][2] = 0.0;
		mr[2][0] = 0.0; mr[2][1] = 0.0; mr[2][2] = 0.0;
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				for (k88 = 0; k88 < 3; k88++) {

					mr[i88][j88] = mr[i88][j88] + mz0[i88][k88] * mx[k88][j88];
				}
			}
		}

		mr1[0][0] = 0.0; mr1[0][1] = 0.0; mr1[0][2] = 0.0;
		mr1[1][0] = 0.0; mr1[1][1] = 0.0; mr1[1][2] = 0.0;
		mr1[2][0] = 0.0; mr1[2][1] = 0.0; mr1[2][2] = 0.0;
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				for (k88 = 0; k88 < 3; k88++) {

					mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * my[k88][j88];
				}
			}
		}
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				mr[i88][j88] = mr1[i88][j88];
				mr1[i88][j88] = 0.0;
			}
		}
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				for (k88 = 0; k88 < 3; k88++) {

					mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * mz[k88][j88];
				}
			}
		}
		nx1 = 0.0;
		ny1 = 0.0;
		nz1 = 0.0;
		nx1 = mr1[0][0] * nx + mr1[0][1] * ny + mr1[0][2] * nz;
		ny1 = mr1[1][0] * nx + mr1[1][1] * ny + mr1[1][2] * nz;
		nz1 = mr1[2][0] * nx + mr1[2][1] * ny + mr1[2][2] * nz;

		//if (nx*nx1+ny*ny1+nz*nz1>-1.0e-3)
		if (nz1 > -1.0e-3)
		{
			br = true;
		}
		else
		{
			br = false;
		}

		if (br)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

void Drawnvtx(int**& nvtx, TOCHKA*& pa, int id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);


// Каркасная визуализация с удалением невидимых линий. Проверено 19.12.2020.
// Данная версия подходит для механической задачи, не содержит информации о блоках.
int DrawZbufferColor(int maxelm, int maxnod, int**& nvtx) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; ++i) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; ++i) {
		for (int j = 0; j < CHECK_LINE_SIZE; ++j) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool* [7];
	for (int i = 0; i < 7; ++i) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; ++i) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; ++i) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // Вся ячейка целиком!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++j) {
				//std::cout << nvtx[j][i] - 1 << " ";
				ipa_count[nvtx[j][i] - 1]++;
			}

			//system("pause");
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; ++j) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {


			bvisible_gran[6][j] = true;

			// XY bottom
			if ((ipa_count[inode1] <= bVisibleCount) &&
				(ipa_count[inode2] <= bVisibleCount) &&
				(ipa_count[inode3] <= bVisibleCount) &&
				(ipa_count[inode4] <= bVisibleCount)) {

				bvisible_gran[B_SIDE][j] = true;
			}

			// XY Top
			if ((ipa_count[inode5] <= bVisibleCount) &&
				(ipa_count[inode6] <= bVisibleCount) &&
				(ipa_count[inode7] <= bVisibleCount) &&
				(ipa_count[inode8] <= bVisibleCount)) {

				bvisible_gran[T_SIDE][j] = true;
			}

			// XZ SSIDE min Y
			if ((ipa_count[inode1] <= bVisibleCount) &&
				(ipa_count[inode2] <= bVisibleCount) &&
				(ipa_count[inode6] <= bVisibleCount) &&
				(ipa_count[inode5] <= bVisibleCount)) {

				bvisible_gran[S_SIDE][j] = true;
			}

			// XZ SSIDE max Y
			if ((ipa_count[inode4] <= bVisibleCount) &&
				(ipa_count[inode8] <= bVisibleCount) &&
				(ipa_count[inode7] <= bVisibleCount) &&
				(ipa_count[inode3] <= bVisibleCount)) {

				bvisible_gran[N_SIDE][j] = true;
			}

			// YZ SSIDE max X
			if ((ipa_count[inode2] <= bVisibleCount) &&
				(ipa_count[inode3] <= bVisibleCount) &&
				(ipa_count[inode7] <= bVisibleCount) &&
				(ipa_count[inode6] <= bVisibleCount)) {

				bvisible_gran[E_SIDE][j] = true;
			}

			// YZ SSIDE min X
			if ((ipa_count[inode1] <= bVisibleCount) &&
				(ipa_count[inode5] <= bVisibleCount) &&
				(ipa_count[inode8] <= bVisibleCount) &&
				(ipa_count[inode4] <= bVisibleCount))
			{
				bvisible_gran[W_SIDE][j] = true;
			}


		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
		}
		//glColor3f(0.0, 0.0, 0.0); // Чёрный
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);

		glLineWidth(1);

		// Рисовать тут.
		for (int j = 0; j < maxelm; ++j) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						set_color_for_render(inode4);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(inode3);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(inode2);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(inode1);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode8);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(inode7);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						set_color_for_render(inode6);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(inode5);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode6])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode5);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						set_color_for_render(inode1);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(inode2);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(inode6);
						//glColor3f(1.0, 1.0, 1.0);
						///glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode2],
							ipa_count[inode6])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode8);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(inode4);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(inode3);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(inode7);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode8],
							ipa_count[inode4])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode6);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(inode2);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(inode3);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(inode7);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}


						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode6])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode8);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(inode4);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(inode1);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(inode5);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode4],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; ++i) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbufferColor


  // Каркасная визуализация с удалением невидимых линий и с освещением. 
int DrawZbufferColorLight(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; ++i) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; ++i) {
		for (int j = 0; j < CHECK_LINE_SIZE; ++j) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool* [7];
	for (int i = 0; i < 7; ++i) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; ++i) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; ++i) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // Вся ячейка целиком!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// Использование быстродействующей хеш таблицы whot_is_block значительно
			// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//system("pause");
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//system("pause");
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//system("pause");
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//system("pause");
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//system("pause");
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//system("pause");
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//system("pause");
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//system("pause");
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {


				for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++j) {
					//std::cout << nvtx[j][i] - 1 << " ";
					ipa_count[nvtx[j][i] - 1]++;
				}
			}
			//system("pause");
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; ++j) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// Использование быстродействующей хеш таблицы whot_is_block значительно
			// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//system("pause");
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//system("pause");
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//system("pause");
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//system("pause");
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//system("pause");
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//system("pause");
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//system("pause");
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//system("pause");
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

				bvisible_gran[6][j] = true;

				// XY bottom
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount)) {

					bvisible_gran[B_SIDE][j] = true;
				}

				// XY Top
				if ((ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount)) {

					bvisible_gran[T_SIDE][j] = true;
				}

				// XZ SSIDE min Y
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount)) {

					bvisible_gran[S_SIDE][j] = true;
				}

				// XZ SSIDE max Y
				if ((ipa_count[inode4] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount)) {

					bvisible_gran[N_SIDE][j] = true;
				}

				// YZ SSIDE max X
				if ((ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount)) {

					bvisible_gran[E_SIDE][j] = true;
				}

				// YZ SSIDE min X
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount))
				{
					bvisible_gran[W_SIDE][j] = true;
				}

			}
		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_LIGHTING); // рассеяный свет.
	glEnable(GL_LIGHT0); // источник света номера 0.
	glEnable(GL_COLOR_MATERIAL); // Учитываются цвета рисуемых объектов.

								 /* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{


		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);// Белый фон
		}
		//glColor3f(0.0, 0.0, 0.0); // Чёрный
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();

		float position[] = { 0,0,1,0 };
		glLightfv(GL_LIGHT0, GL_POSITION, position);

		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);



		glLineWidth(1);

		// Рисовать тут.
		for (int j = 0; j < maxelm; ++j) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						/*
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode2])
						)
						{

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode2])
						)
						{

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode4])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode4])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode5],
						ipa_count[inode6])
						)
						{

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode6],
						ipa_count[inode7])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode7],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode5],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode2])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode2],
						ipa_count[inode6])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode6],
						ipa_count[inode5])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode5])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode4])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode7])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode7],
						ipa_count[inode8])
						) {
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode8],
						ipa_count[inode4])
						) {
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// белый
						}


						if (CheckLine(ipa_count[inode3],
						ipa_count[inode2])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode7])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode7],
						ipa_count[inode6])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode6],
						ipa_count[inode2])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode4])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode4],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode5],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode5])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; ++i) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbufferColorLight


  // Каркасная визуализация с удалением невидимых линий. Проверено 19.12.2020.
int DrawZbufferColor(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; ++i) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; ++i) {
		for (int j = 0; j < CHECK_LINE_SIZE; ++j) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool* [7];
	for (int i = 0; i < 7; ++i) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; ++i) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; ++i) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // Вся ячейка целиком!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// Использование быстродействующей хеш таблицы whot_is_block значительно
			// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//system("pause");
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//system("pause");
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//system("pause");
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//system("pause");
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//system("pause");
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//system("pause");
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//system("pause");
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//system("pause");
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {


				for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++j) {
					//std::cout << nvtx[j][i] - 1 << " ";
					ipa_count[nvtx[j][i] - 1]++;
				}
			}
			//system("pause");
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; ++j) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// Использование быстродействующей хеш таблицы whot_is_block значительно
			// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//system("pause");
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//system("pause");
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//system("pause");
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//system("pause");
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//system("pause");
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//system("pause");
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//system("pause");
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//system("pause");
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

				bvisible_gran[6][j] = true;

				// XY bottom
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount)) {

					bvisible_gran[B_SIDE][j] = true;
				}

				// XY Top
				if ((ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount)) {

					bvisible_gran[T_SIDE][j] = true;
				}

				// XZ SSIDE min Y
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount)) {

					bvisible_gran[S_SIDE][j] = true;
				}

				// XZ SSIDE max Y
				if ((ipa_count[inode4] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount)) {

					bvisible_gran[N_SIDE][j] = true;
				}

				// YZ SSIDE max X
				if ((ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount)) {

					bvisible_gran[E_SIDE][j] = true;
				}

				// YZ SSIDE min X
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount))
				{
					bvisible_gran[W_SIDE][j] = true;
				}

			}
		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{


		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);// Белый фон
		}
		//glColor3f(0.0, 0.0, 0.0); // Чёрный
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);

		glLineWidth(1);

		// Рисовать тут.
		for (int j = 0; j < maxelm; ++j) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode6])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						///glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode2],
							ipa_count[inode6])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode8],
							ipa_count[inode4])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}


						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode6])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode4],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; ++i) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbufferColor

  // Каркасная визуализация с удалением невидимых линий. Проверено 17.12.2020.
int DrawZbuffer(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; ++i) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; ++i) {
		for (int j = 0; j < CHECK_LINE_SIZE; ++j) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool* [7];
	for (int i = 0; i < 7; ++i) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; ++i) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; ++i) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // Вся ячейка целиком!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// Использование быстродействующей хеш таблицы whot_is_block значительно
			// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//system("pause");
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//system("pause");
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//system("pause");
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//system("pause");
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//system("pause");
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//system("pause");
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//system("pause");
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//system("pause");
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {


				for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++j) {
					//std::cout << nvtx[j][i] - 1 << " ";
					ipa_count[nvtx[j][i] - 1]++;
				}
			}
			//system("pause");
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; ++j) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// Использование быстродействующей хеш таблицы whot_is_block значительно
			// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//system("pause");
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//system("pause");
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//system("pause");
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//system("pause");
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//system("pause");
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//system("pause");
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//system("pause");
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//system("pause");
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

				bvisible_gran[6][j] = true;

				// XY bottom
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount)) {

					bvisible_gran[B_SIDE][j] = true;
				}

				// XY Top
				if ((ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount)) {

					bvisible_gran[T_SIDE][j] = true;
				}

				// XZ SSIDE min Y
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount)) {

					bvisible_gran[S_SIDE][j] = true;
				}

				// XZ SSIDE max Y
				if ((ipa_count[inode4] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount)) {

					bvisible_gran[N_SIDE][j] = true;
				}

				// YZ SSIDE max X
				if ((ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount)) {

					bvisible_gran[E_SIDE][j] = true;
				}

				// YZ SSIDE min X
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount))
				{
					bvisible_gran[W_SIDE][j] = true;
				}

			}
		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);// Белый фон.
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);

		glLineWidth(3);

		// Рисовать тут.
		for (int j = 0; j < maxelm; ++j) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // Чёрный
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// белый
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// белый
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode6])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// белый
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						//glColor3f(1.0, 1.0, 1.0);
						///glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode2],
							ipa_count[inode6])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// белый
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode8],
							ipa_count[inode4])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// белый
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}



						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode6])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// белый
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // Чёрный
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// белый
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode4],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; ++i) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbuffer

#endif

  // Рабочая сохраненная копия каркасного рендера.
int main_copy1(void)
{

	// количество блоков, источников и стенок, юнионов.
	int lb = 0;
	BLOCK* b = nullptr;// список блоков

					   // Вызываем солвер.
	main_solverv0_48(b, lb);

#ifndef NO_OPENGL_GLFW

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 

					  //GLfloat halfScreenWidth = SCREEN_WIDTH / 2;
					  //GLfloat halfScreenHeight = SCREEN_HEIGHT / 2;



					  /* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{

		//glClearColor(0.7f, 1.0f, 0.7f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -500);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, 500);

		DrawKarkas(b, lb, halfScreenWidth, halfScreenHeight, -500.0, scale_all);

		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

#endif


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
	delete[] b;
	b = nullptr;

	return 0;
} // Сохранённая копия рендера №1

#ifndef NO_OPENGL_GLFW

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// std::cout << key << std::endl;

	const GLfloat rotationSpeed = 10;
	doublereal dmin = 1.0e30;
	doublereal dmax = -1.0e30;

	// actions are GLFW_PRESS, GLFW_RELEASE or GLFW_REPEAT
	if (action == GLFW_PRESS || action == GLFW_REPEAT)
	{
		switch (key)
		{
		case GLFW_KEY_UP:
			rotationX -= rotationSpeed;
			rotationX = ((int)rotationX) % 360;
			break;
		case GLFW_KEY_DOWN:
			rotationX += rotationSpeed;
			rotationX = ((int)rotationX) % 360;
			break;
		case GLFW_KEY_RIGHT:
			rotationY += rotationSpeed;
			rotationY = ((int)rotationY) % 360;
			break;
		case GLFW_KEY_LEFT:
			rotationY -= rotationSpeed;
			rotationY = ((int)rotationY) % 360;
			break;
		case GLFW_KEY_R:
			// следующий кадр.
			if (animation_sequence_functions_openGL != nullptr)
			{
				if (iCURENT_ANIMATION_CADER < iNUMBER_ANIMATION_CADERS - 1) {
					iCURENT_ANIMATION_CADER++;
				}
				else {
					iCURENT_ANIMATION_CADER = 0;
				}

				dmin = 1.0e30;
				dmax = -1.0e30;


				for (int i37 = 0; i37 < n_render; i37++) {
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] > dmax) {
						dmax = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] < dmin) {
						dmin = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
				}

				minimum_val_for_render_pic = dmin;
				maximum_val_for_render_pic = dmax;
			}
			break;
		case GLFW_KEY_F:
			// Предыдущий кадр.
			if (animation_sequence_functions_openGL != nullptr)
			{
				if (iCURENT_ANIMATION_CADER > 0) {
					iCURENT_ANIMATION_CADER--;
				}
				else {
					iCURENT_ANIMATION_CADER = iNUMBER_ANIMATION_CADERS - 1;
				}

				dmin = 1.0e30;
				dmax = -1.0e30;


				for (int i37 = 0; i37 < n_render; i37++) {
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] > dmax) {
						dmax = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] < dmin) {
						dmin = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
				}

				minimum_val_for_render_pic = dmin;
				maximum_val_for_render_pic = dmax;
			}
			break;
			/*
			// Это так не работает. Там изометрия все сложно.
			case GLFW_KEY_O:
			halfScreenHeight += 0.05 * SCREEN_HEIGHT;
			break;
			case GLFW_KEY_L:
			halfScreenHeight -= 0.05 * SCREEN_HEIGHT;
			break;
			case GLFW_KEY_U:
			halfScreenWidth += 0.05 * SCREEN_WIDTH;
			break;
			case GLFW_KEY_Y:
			halfScreenWidth -= 0.05 * SCREEN_WIDTH;
			break;
			*/
		case GLFW_KEY_S: // удалить объект
			scale_all = scale_all / 1.2f;
			for (int i = 0; i < n_render; ++i) {
				pa_render[i].x = halfScreenWidth + scale_all * pa_opengl[i].x;
				pa_render[i].y = halfScreenHeight + scale_all * pa_opengl[i].y;
				pa_render[i].z = -abbys + scale_all * pa_opengl[i].z;
			}
			break;
		case GLFW_KEY_W: // приблизить обьект
			scale_all = scale_all * 1.2f;
			for (int i = 0; i < n_render; ++i) {
				pa_render[i].x = halfScreenWidth + scale_all * pa_opengl[i].x;
				pa_render[i].y = halfScreenHeight + scale_all * pa_opengl[i].y;
				pa_render[i].z = -abbys + scale_all * pa_opengl[i].z;
			}
			break;
		case GLFW_KEY_E: // смена анимируемой функциии.

			if (iCFD_animation == 1) {

				if (animation_sequence_functions_openGL != nullptr)
				{

					if (iCURENT_ANIMATION_FUNCTION < iNUMBER_ANIMATION_FUNCTIONS - 1) {
						iCURENT_ANIMATION_FUNCTION++;
					}
					else {
						iCURENT_ANIMATION_FUNCTION = 0;
					}

					dmin = 1.0e30;
					dmax = -1.0e30;


					for (int i37 = 0; i37 < n_render; i37++) {
						if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] > dmax) {
							dmax = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
						}
						if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] < dmin) {
							dmin = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
						}
					}

					minimum_val_for_render_pic = dmin;
					maximum_val_for_render_pic = dmax;
				}

			}
			else {
				if (potent_array_opengl != nullptr) {

					if (iCURENT_FUNC_openGL < iNUMBER_FUNC_openGL - 1) {
						iCURENT_FUNC_openGL++;
					}
					else {
						iCURENT_FUNC_openGL = 0;
					}
					dmin = 1.0e30;
					dmax = -1.0e30;


					for (int i37 = 0; i37 < n_render; i37++) {
						if (potent_array_opengl[iCURENT_FUNC_openGL][i37] > dmax) {
							dmax = potent_array_opengl[iCURENT_FUNC_openGL][i37];
						}
						if (potent_array_opengl[iCURENT_FUNC_openGL][i37] < dmin) {
							dmin = potent_array_opengl[iCURENT_FUNC_openGL][i37];
						}
					}

					minimum_val_for_render_pic = dmin;
					maximum_val_for_render_pic = dmax;
				}
			}
			break;

		}
		GLfloat Alf = 3.14159265 * rotationX / 180.0;// Угол в радианах.
		GLfloat Bet = 3.14159265 * rotationY / 180.0;// Угол в радианах.
		cosAlf = cos(Alf);
		sinAlf = sin(Alf);
		cosBet = cos(Bet);
		sinBet = sin(Bet);

	}

}


void DrawCADobj(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {
	//if (!b[id].g.bbigCADmodel) 
	{

		bool bvisible = true;

		if (id == 0) {
			bvisible = false;
		}
		if (b[id].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			bvisible = false;
		}


		if (bvisible) {

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

			gCAD_STL* tmp = b[id].g.root_CAD_STL;

			while (tmp != nullptr) {

				glBegin(GL_TRIANGLES); // Указываем что рисовать
				glNormal3f(tmp->n.x, tmp->n.y, tmp->n.z);
				glVertex3f(scale * tmp->pa.x + centerPosX, scale * tmp->pa.y + centerPosY, scale * tmp->pa.z + centerPosZ);
				glVertex3f(scale * tmp->pb.x + centerPosX, scale * tmp->pb.y + centerPosY, scale * tmp->pb.z + centerPosZ);
				glVertex3f(scale * tmp->pc.x + centerPosX, scale * tmp->pc.y + centerPosY, scale * tmp->pc.z + centerPosZ);
				glEnd(); // закончили рисовать

				tmp = tmp->next;
			}
		}
	}
}

void DrawPrism(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {

	bool bCAD = false;
	for (int i = 0; i < lb; ++i) {
		if (b[i].g.itypegeom == CAD_STL) {
			bCAD = true;
		}
	}

	bool bvisible = true;
	if (bCAD) {
		if (id == 0) {
			bvisible = false;
		}
		if (b[id].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			bvisible = false;
		}
	}

	if (bvisible) {

		GLfloat vertices[] =
		{
			// front face
			scale * b[id].g.xS + centerPosX, scale * b[id].g.yE + centerPosY, scale * b[id].g.zE + centerPosZ, // top left
			scale * b[id].g.xE + centerPosX, scale * b[id].g.yE + centerPosY, scale * b[id].g.zE + centerPosZ, // top right
			scale * b[id].g.xE + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zE + centerPosZ, // bottom right
			scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zE + centerPosZ, // bottom left

																											   // back face
																											   scale * b[id].g.xS + centerPosX, scale * b[id].g.yE + centerPosY, scale * b[id].g.zS + centerPosZ, // top left
																											   scale * b[id].g.xE + centerPosX, scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																											   scale * b[id].g.xE + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zS + centerPosZ, // bottom right
																											   scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zS + centerPosZ, // bottom left

																																																				  // left face
																																																				  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zE + centerPosZ, // top left
																																																				  scale * b[id].g.xS + centerPosX, scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																				  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yS + centerPosY, scale * b[id].g.zS + centerPosZ, // bottom right
																																																				  scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zE + centerPosZ, // bottom left

																																																																													 // right face
																																																																													 scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY, scale * b[id].g.zE + centerPosZ, // top left
																																																																													 scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																																													 scale * b[id].g.xE + centerPosX, scale * b[id].g.yS + centerPosY,  scale * b[id].g.zS + centerPosZ, // bottom right
																																																																													 scale * b[id].g.xE + centerPosX,  scale * b[id].g.yS + centerPosY,  scale * b[id].g.zE + centerPosZ, // bottom left

																																																																																																						  // top face
																																																																																																						  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zE + centerPosZ, // top left
																																																																																																						  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																																																																						  scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // bottom right
																																																																																																						  scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zE + centerPosZ, // bottom left

																																																																																																																															   // bottom face
																																																																																																																															   scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY,  scale * b[id].g.zE + centerPosZ, // top left
																																																																																																																															   scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																																																																																															   scale * b[id].g.xE + centerPosX,  scale * b[id].g.yS + centerPosY,  scale * b[id].g.zS + centerPosZ, // bottom right
																																																																																																																															   scale * b[id].g.xE + centerPosX,  scale * b[id].g.yS + centerPosY,  scale * b[id].g.zE + centerPosZ, // bottom left

		};

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, vertices);
		glDrawArrays(GL_QUADS, 0, 24);
		glDisableClientState(GL_VERTEX_ARRAY);

	}

}

void Drawnvtx(int**& nvtx, TOCHKA*& pa, int id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {



	GLfloat vertices[] =
	{
		// front face
		scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
		scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top right
		scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom right
		scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																			// back face
																																			scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top left
																																			scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																			scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																			scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom left

																																																																				// left face
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																																																																																					// right face
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																																																																																																																						  // top face
																																																																																																																																						  scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																																																																																						  scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																																																																																						  scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																																																																																						  scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																																																																																																																																																								// bottom face
																																																																																																																																																																								scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																																																																																																																								scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																																																																																																																								scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																																																																																																																								scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

	};

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertices);
	glDrawArrays(GL_QUADS, 0, 24);
	glDisableClientState(GL_VERTEX_ARRAY);

}

void DrawCylinder(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale)
{

	GLfloat angle33, dx33, dy33;
	GLfloat Hcyl = scale * b[id].g.Hcyl;
	GLfloat xC = scale * b[id].g.xC + centerPosX;
	GLfloat yC = scale * b[id].g.yC + centerPosY;
	GLfloat zC = scale * b[id].g.zC + centerPosZ;
	GLfloat R_out_cyl = scale * b[id].g.R_out_cyl;
	GLfloat R_out_cyl2 = 0.0;// scale* b[id].g.R_out_cyl;
	GLfloat R_in_cyl = scale * b[id].g.R_in_cyl;
	GLfloat R_in_cyl2 = 0.0; // scale* b[id].g.R_in_cyl;

							 // Cylinder
	switch (b[id].g.iPlane) {
	case XY_PLANE:
		// XY
		glBegin(GL_LINE_LOOP);
		for (int i33 = 0; i33 <= 29; i33++) {

			angle33 = 2.0 * 3.1415926 * i33 / 29.0;
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
		}
		glEnd();
		if (R_out_cyl2 <= 0.0) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl * cos(angle33);
				dy33 = R_out_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl2 * cos(angle33);
				dy33 = R_out_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
		}

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_in_cyl * cos(angle33);
				dy33 = R_in_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC);
			}
			glEnd();

			if (R_in_cyl2 <= 0.0) {

				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl * cos(angle33);
					dy33 = R_in_cyl * sin(angle33);
					glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
				}
				glEnd();
			}
			else
			{
				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl2 * cos(angle33);
					dy33 = R_in_cyl2 * sin(angle33);
					glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
				}
				glEnd();
			}
		}

		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.125;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.375;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.625;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.875;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.125;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.375;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.625;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.875;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
		}

		break;
	case XZ_PLANE:
		// XZ
		glBegin(GL_LINE_LOOP);
		for (int i33 = 0; i33 <= 29; i33++) {

			angle33 = 2.0 * 3.1415926 * i33 / 29.0;
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
		}
		glEnd();
		if (R_out_cyl2 <= 0.0) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl * cos(angle33);
				dy33 = R_out_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl2 * cos(angle33);
				dy33 = R_out_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
		}

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_in_cyl * cos(angle33);
				dy33 = R_in_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC, zC + dy33);
			}
			glEnd();
			if (R_in_cyl2 <= 0.0) {

				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl * cos(angle33);
					dy33 = R_in_cyl * sin(angle33);
					glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
				}
				glEnd();
			}
			else
			{
				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl2 * cos(angle33);
					dy33 = R_in_cyl2 * sin(angle33);
					glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
				}
				glEnd();
			}
		}

		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.125;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.375;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.625;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.875;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.125;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			};
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.375;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.625;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.875;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
		}


		break;
	case YZ_PLANE:
		// YZ
		glBegin(GL_LINE_LOOP);
		for (int i33 = 0; i33 <= 29; i33++) {

			angle33 = 2.0 * 3.1415926 * i33 / 29.0;
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
		}
		glEnd();
		if (R_out_cyl2 <= 0.0) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl * cos(angle33);
				dy33 = R_out_cyl * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl2 * cos(angle33);
				dy33 = R_out_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
		}
		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_in_cyl * cos(angle33);
				dy33 = R_in_cyl * sin(angle33);
				glVertex3f(xC, yC + dx33, zC + dy33);
			}
			glEnd();
			if (R_in_cyl2 <= 0.0) {

				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl * cos(angle33);
					dy33 = R_in_cyl * sin(angle33);
					glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
				}
				glEnd();
			}
			else
			{
				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl2 * cos(angle33);
					dy33 = R_in_cyl2 * sin(angle33);
					glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
				}
				glEnd();
			}
		}

		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.125;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.375;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.625;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.875;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.125;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.375;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.625;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.875;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
		}
		break;
	}

}

void DrawPolygon(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {



	// Polygon
	switch (b[id].g.iPlane_obj2) {
	case XY_PLANE:
		// XY
		for (int j3 = 0; j3 <= b[id].g.nsizei - 2; j3++)
		{
			// GL_LINE_LOOP
			glBegin(GL_LINE_LOOP);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1] + scale * b[id].g.hi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3] + scale * b[id].g.hi[j3]);
			glEnd();
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0] + scale * b[id].g.hi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1] + scale * b[id].g.hi[b[id].g.nsizei - 1]);
		glEnd();
		break;
	case XZ_PLANE:
		// XZ
		for (int j3 = 0; j3 <= b[id].g.nsizei - 2; j3++)
		{

			glBegin(GL_LINE_LOOP);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1] + scale * b[id].g.hi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3] + scale * b[id].g.hi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glEnd();
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0] + scale * b[id].g.hi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1] + scale * b[id].g.hi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glEnd();

		break;
	case YZ_PLANE:
		// YZ
		for (int j3 = 0; j3 <= b[id].g.nsizei - 2; j3++)
		{
			glBegin(GL_LINE_LOOP);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1] + scale * b[id].g.hi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3] + scale * b[id].g.hi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glEnd();
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0] + scale * b[id].g.hi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1] + scale * b[id].g.hi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glEnd();
		break;

	}
}

#endif



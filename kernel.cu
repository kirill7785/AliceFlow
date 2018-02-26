// AliceFlow_v0_24.cpp 
// 9 июля 2017 переход на 64 битные целые int64_t.
// 15 апреля 2017 откомпилирована в visual studio community edition 2017 (open source).
// 1 октября 2016 откомпилировал на nvidia cuda 8.0. 
// АЛИС сетки введены в строй для теплопередачи в твёрдом теле. 
// 11 января 2016 года добавил cl_agl_amg_v0_14.
// 15 августа 2015 года. Теперь в Visual Studio 2015.
// AliceFlow_v0_21.cpp
// 15 августа 2015. Теперь в Visual Studio 2013.
// AliceFlow_v0_20.cpp
// 14 августа 2015 Действительно правильное распараллеливание lusol и ilu2 decomposition на 2 потока.
// AliceFlow_v0_07.cpp: определяет точку входа для консольного приложения.
// AliceFlow_v0_07.cpp на основе  AliceFlow_v0_06.cpp, но теперь с LES моделью турбулентности.
// Программа не прошла тестирование и содержит ошибки.
// 17 апреля 2013 года. Правильное распараллеливание lusol_.
// 1 апреля 2013. Теперь в Visual Studio 2012.
//
// AliceFlow_v0_06.cpp :
// 3D программа AliceFlow_v0_06.cpp наследует свойства AliceFlowv0_05.cpp
// развивая их дальше, наращивая стабильность и функциональность.
// 
// Программа AliceFlowv0_05.cpp, 
// улучшенный вариант AliceFlowv0_03.cpp, преобразует 
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
//    данных необходимых для сборки матрцы СЛАУ.
// 3. Сборка матрицы СЛАУ.
// 4. Решатель СЛАУ.
// 5. Экспорт в визуализатор tecplot360.
// begin two 30 июня 2011 года.
// begin three 14 октября 2011 года. Теперь на Visual Studio 2010.
// begin four 12 марта 2012 года. (по мотивам книги Р. Лафоре - переход на ООП).
//
// Реальные размеры источника тепла в транзисторе типа TGF2023_*
// равны 0.2x120мкм толщиной 100 ангстрем (10-20нм).
// Для задачи Блазиуса без учёта краевых эфектов сначала 2-3 длины пластинки, в конце 3-5 длин пластинки. 

// Раскоментировать в случае если сборка приложения осуществляется компилятором gcc от GNU.
//#define MINGW_COMPILLER 1

#ifdef MINGW_COMPILLER
#include <stdio.h>
#endif

// для std::locale::global(std::locale("en_US.UTF-8"));
// Не работает.
//#include <locale.h>

//#include "stdafx.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <cinttypes> // для типа 64 битного целого числа int64_t


#include <stdlib.h> 
#include <omp.h> // OpenMP
//-Xcompiler "/openmp" 

//using namespace System;

// Вещественная арифметика.
#define doubleprecision 1
#if doubleprecision == 1
#define doublereal double //double // модель вещественного числа двойной точности
//#define doublereal Decimal // decimal
#else
#define doublereal float //float // модель вещественного числа одинарной точности
#endif

#define doubleintprecision 1



#if doubleintprecision == 1
#define integer int64_t
#else
#define integer int
#endif

// 9 september 2017.
// делать ли освобождения оперативной памяти и новые построения структур данных.
// Полигоны вызывают проблемы при перестроении сетки, т.к. сетка каждый раз строится по новому, а поле температур 
// остаётся старым. Другое решение данной проблемы состоит из перевыделения поля температур под новую сетку с 
// последующей переинтерполляцией значений температур на новую сетку.
integer ireconstruction_free_construct_alloc = 1; // 0 - off, 1 - on.
// Записывать ли ианимацию в текстовый файл
// по окончанию каждого нового шага по времени.
integer ianimation_write_on = 0; // 0 - off, 1 - on.

// Для достижения сходимости мы первые 300 итераций считаем со скоростью vel*rGradual_changes а потом 
// делаем перемасштабирование и досчитываем уже со скоростью vel.
doublereal rGradual_changes = 0.1; // 1.0 - не используется.

// инициализация компонент скорости константой.
// Единые значения для всей расчётной области.
// initialization value.
doublereal starting_speed_Vx = 0.0;
doublereal starting_speed_Vy = 0.0;
doublereal starting_speed_Vz = 0.0;

// Опорная линия для XY-Plot (variation Plot).
// Мы сохраняем опорную точку, через которую проходит линия, и направлеие 
// линии вдоль одной из осей декартовой прямоугольной системы координат.
doublereal Tochka_position_X0_for_XY_Plot = 0.0;
doublereal Tochka_position_Y0_for_XY_Plot = 0.0;
doublereal Tochka_position_Z0_for_XY_Plot = 0.0;
integer idirectional_for_XY_Plot = 0; // 0 - Ox axis. 

// При iVar==TEMP && lw==1 выход из солвера можно осуществлять когда максимальная температура между V циклами отличается менее 0.5K.
bool bPhysics_stop = false;
// Особое выделение количества оперативной памяти для ПТБШ.
bool bPhysics_PTBSH_memory = false;
// Решаем только теплопередачу в твёрдом теле :
bool bonly_solid_calculation = false;

// 3 августа 2015 схемы стали доступны через GUI пользователя
// в связи с чем объявление дентификаторов вынесено в самое начало кода.
// ограниченные схемы
#define UNEVEN_MUSCL 1017  // van Leer (1977)
#define UNEVEN_SOUCUP 1018 // MINMOD
#define UNEVEN_HLPA 1019
#define UNEVEN_SMART 1020 // Gaskell and Lau (1988)
#define UNEVEN_WACEB 1021
#define UNEVEN_SMARTER 1022
#define UNEVEN_STOIC 1023 // Darwish (1993)
#define UNEVEN_CLAM 1024
#define UNEVEN_OSHER 1025 // Chakravarthy and Osher (1983)
#define UNEVEN_VONOS 1026
#define UNEVEN_LPPA 1027
#define UNEVEN_EXPONENTIAL 1028
#define UNEVEN_SUPER_C 1029
#define UNEVEN_ISNAS 1030
#define UNEVEN_CUBISTA 1031
#define UNEVEN_GAMMA 1032 // схема с параметром beta_m
#define UNEVEN_COPLA 1033 // 1 08 2015
#define UNEVEN_SECBC 1034 // 2 08 2015 Yu et al., (2001b) Сингапур, Малазия.
#define UNEVEN_SGSD 1035 // 3 08 2015 Li and Tao (2002)

// Управление алгебраическим многосеточным методом из интерфейса.
typedef struct TMY_AMG_MANAGER {

	// 0 - не печатать портрет матрицы, 
	// 1 - печатать портрет матрицы.
	integer bTemperatureMatrixPortrait;
	integer bSpeedMatrixPortrait;
	integer bPressureMatrixPortrait;
	integer bStressMatrixPortrait;
	integer bMatrixPortrait;

	// fgmres(m_restart)
	integer m_restart;

	// lfil for BiCGStab+ILU2 and fgmres.
	integer lfil;

	// Temperature
	doublereal theta_Temperature;
	integer maximum_delete_levels_Temperature;
	integer nFinnest_Temperature, nu1_Temperature, nu2_Temperature;
	integer memory_size_Temperature;
	integer ilu2_smoother_Temperature; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Temperature;
	// Speed
	doublereal theta_Speed;
	integer maximum_delete_levels_Speed;
	integer nFinnest_Speed, nu1_Speed, nu2_Speed;
	integer memory_size_Speed;
	integer ilu2_smoother_Speed; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Speed;
	// Pressure
	doublereal theta_Pressure;
	integer maximum_delete_levels_Pressure;
	integer nFinnest_Pressure, nu1_Pressure, nu2_Pressure;
	integer memory_size_Pressure;
	integer ilu2_smoother_Pressure; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Pressure;
	// Stress
	doublereal theta_Stress;
	integer maximum_delete_levels_Stress;
	integer nFinnest_Stress, nu1_Stress, nu2_Stress;
	integer memory_size_Stress;
	integer ilu2_smoother_Stress; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Stress;
	// global
	bool bCFJacoby;
	integer iRunge_Kutta_smoother; // 3 - третьего порядка, 5 - пятого порядка, любое другое число не используется. 
	integer iFinnest_ilu; // 0 не используется, 1 - ilu0. Только на самой подробной сетке.
	// Использование iluk разложения на глубоких уровнях вложенности для которых
	// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
	bool b_ilu_smoothers_in_nnz_n_LE_6;
	doublereal theta; // strength threshold
	//integer maximum_levels; // максимальное количество уровней вложенности (уровни выше редуцируются).
	integer maximum_delete_levels; // Количество уровней отсекаемых снизу в области грубой сетки.
	integer nFinnest, nu1, nu2; // Количества сглаживаний.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure;
	integer memory_size; // В размерах матрицы А.
	// Для метода верхней релаксации в сглаживателе.
	integer number_interpolation_procedure; // идентификатор процедуры интерполляции.
	integer number_interpolation_procedure_Temperature;
	integer number_interpolation_procedure_Speed;
	integer number_interpolation_procedure_Pressure;
	integer number_interpolation_procedure_Stress;

	// 6 december 2016.
	integer itypemodifyinterpol; // номер модификации интерполляции.
	integer inumberadaptpass; // максимальное количество сканов-проходов с модификациями.
	doublereal gold_const, gold_const_Temperature, gold_const_Speed, gold_const_Pressure, gold_const_Stress;
	doublereal magic;
	doublereal F_to_F_Temperature, F_to_F_Speed, F_to_F_Pressure, F_to_F_Stress; // magic
	integer ilu2_smoother; // 0 - не использовать, 1 - использовать.
	// AMG Splitting (coarsening)
	// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
	// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
	integer icoarseningTemp, icoarseningSpeed, icoarseningPressure, icoarseningStress;
	integer icoarseningtype;
	// Stabilization BiCGStab.
	// 8.01.2017
	integer istabilizationTemp, istabilizationSpeed, istabilizationPressure, istabilizationStress;
	integer istabilization;
	// ipatch - номер патча.
	integer ipatch_number;

	integer iprint_log, iprint_log_Temperature, iprint_log_Speed, iprint_log_Pressure, iprint_log_Stress;

	// truncation for interpolation.
	integer itruncation_interpolation, itruncation_interpolation_Temperature, itruncation_interpolation_Speed, itruncation_interpolation_Pressure, itruncation_interpolation_Stress;
	double truncation_interpolation, truncation_interpolation_Temperature, truncation_interpolation_Speed, truncation_interpolation_Pressure, truncation_interpolation_Stress;

	// gmres smoother
	// Ю.Саад, Шульц [1986].
	bool b_gmresTemp, b_gmresSpeed, b_gmresPressure, b_gmresStress;
	bool b_gmres;

} MY_AMG_MANAGER;

MY_AMG_MANAGER my_amg_manager;

bool bglobal_first_start_radiation = true;

// Если мы решаем нестационарную задачу теплопередачи в твердом теле.
bool bglobal_unsteady_temperature_determinant = false;

// Выбор сеточного генератора :
// simplemeshgen == 0 или unevensimplemeshgen ==1.
// По умолчанию используется simplemeshgen == 0.
integer iswitchMeshGenerator = 0; // обычный сеточный генератор.
// нестационарное или стационарное моделирование.
integer steady_or_unsteady_global_determinant = 0; // 0 - steady, 1- unsteady.

// Использовать ли адаптивные локально измельчённые расчётные сетки.
bool b_on_adaptive_local_refinement_mesh = false;
integer itype_ALICE_Mesh = 1;// Тип АЛИС сетки.

typedef struct TTimeStepLaw
{
	integer id_law; // 0 - Linear, 1 - Square Wave.
	doublereal Factor_a_for_Linear;
	doublereal tau; // длительность импульса для Square Wave
	// 06_03_2017 скважность может быть и дробной.
	doublereal Q; // Скважность для Square Wave.
	// Импульсный режим для темы АППАРАТ.
	doublereal m1, tau1, tau2, tau_pause, T_all;
	integer n_cycle;
	// hot cold reshime (double linear)
	doublereal on_time_double_linear;

} TimeStepLaw;

TimeStepLaw glTSL;

// 24 декабря 2016. 
// Для ускорения счёта нелинейных задач в РУМБА 0.14 решателе.
typedef struct TQuickNonlinearBoundaryCondition {
	doublereal emissivity;
	doublereal Tamb, dS;
	doublereal film_coefficient;
	bool bactive;
	bool bStefanBolcman_q_on;
	bool bNewtonRichman_q_on;

} QuickNonlinearBoundaryCondition;

QuickNonlinearBoundaryCondition* qnbc = NULL;
integer iadd_qnbc_maxelm = 0; // для аддитивного сдвига
bool b_sign_on_nonlinear_bc = false;


// Считаем ли мы SIMPLE алгоритмом.
// Это нужно для более точной настройки невязки для уравнения теплопередачи.
// внутри BiCGStab_internal3 решателя.
bool bSIMPLErun_now_for_temperature = false;
// Это нужно для более точной настройки невязки для уравнения теплопередачи
// при расчёте amg1r5 алгоритмом задач с естественой конвекцией.
bool bSIMPLErun_now_for_natural_convection = false;
// Дополнительная нижняя релаксация для температуры.
doublereal* told_temperature_global_for_HOrelax = NULL;

/*
для внутренних плоских источников тепла организуется виртуальная грань.
эта грань общая для двух КО располагающихся по-бокам от источника.
принцип едиственности должен приводить к тому что теплопроводность грани должна быть единственна
при обработке обоих контрольных объемов примыкающих к данной грани. В качестве теплопроводности
берётся среднее геометрическое.
*/
bool *sourse2Dproblem=NULL;
doublereal *conductivity2Dinsource=NULL;

// дополнительная нижняя релаксация.
bool bHORF = false;
bool bdontstartsolver = false;
doublereal* bPamendment_source_old = NULL;
doublereal* bsource_term_radiation_for_relax = NULL;
doublereal* b_buffer_correct_source = NULL;
// Во избежании расходимости РУМБА_0.14 на первой же итерации.
bool bfirst_start_nonlinear_process = true;

// Условие Ньютона-Рихмана по дефолту для температуры.
integer adiabatic_vs_heat_transfer_coeff = 0; // 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
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

FILE* fp_radiation_log = NULL;
errno_t err_radiation_log;

// 1 - визуализация только твёрдого тела.
// 0 - визуализация и жидкости и твердого тела.
integer ionly_solid_visible = 0;

// переключение между алгебраическим многосеточным методом и алгоритмом Ван дер Ворста BiCGStab+ILU2.
// 0 - алгоритм BiCGStab + ILU2. 1 - алгоритм алгебраического многосеточного метода amg1r5.
// 2 - BiCGStab + ADI (Lr1sk).
integer iswitchsolveramg_vs_BiCGstab_plus_ILU2 = 0; // BiCGStab + ILU2.

bool bwait = false; // если false то мы игнорируем getchar().
#define admission 1.0e-10 // для определения совпадения двух вещественных чисел.

unsigned int calculation_vorst_seach_time = 0;

// Если температура канала превысит 
// температуру в 200 градусов Цельсия 
// то прибор Выйдет из строя (сгорит).
// В случае превышения температуры равной TEMPERATURE_FAILURE_DC
// на консоль печатается предупреждающее сообщение и после того
// как пользователь прочтёт сообщение осуществляется выход из программы.
#define TEMPERATURE_FAILURE_DC 5000.2


// Красавин Денис Андреевич Повышение порядка точности аппроксимации
// граничного условия. См. его диссертацию а также книжку С. Патанкара.
// BETA 1.0 4/3=1.333333333 6/5=1.2
#define BETA 1.0

// UDS  см. my_approx_convective.c
unsigned int iFLOWScheme = 2; // Противопоточная первого порядка
unsigned int iTEMPScheme = 2; // Противопоточная первого порядка

// включает более быстро сходящийся алгоритм SIMPLEC
// SIMPLEC Van Doormal and Raithby, 1984 год.
// SIMPLEC отличается от SIMPLE только в двух вещах:
// 1. В SIMPLEC не используется нижняя релаксация для давления при коррекции давления, т.е. alphaP=1.0.
// 2. В SIMPLEC псевдо время пропорционально tau ~ alpha/(1-alpha), а в SIMPLE tau ~ alpha. 
// В остальном алгоритмы полностью совпадают.
#define SIMPLE_Carretto 0 // алгоритм SIMPLE Carretto et al., 1973 используется по умолчанию.
#define SIMPLEC_Van_Doormal_and_Raithby 1 // алгоритм SIMPLEC Van Doormal and Raithby, 1984 год.
unsigned int iSIMPLE_alg = SIMPLE_Carretto;// SIMPLE_Carretto SIMPLEC_Van_Doormal_and_Raithby

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

// В этот файл будет записываться информация о работе
// линейных решателей СЛАУ.
FILE *fp_statistic_convergence=NULL;
// в этот файл будет осуществляться запись лога.
// лог нужен для анализа изменений внесённых в программу.
// он содержит информацию обо всех невязках получаемых в результате вычислительного процесса.
FILE *fp_log=NULL;

// используется для ускорения решения 
// задачи твёрдотельной теплопередачи.
bool bsolid_static_only = false;

const integer inumcore = 2; // число ядер процессора
const bool bparallelizm_old = false;

// Структура границ деления :
typedef struct TPARBOUND {
	integer ileft_start, ileft_finish, iright_start, iright_finish, iseparate_start, iseparate_finish;
	bool active; // активность декомпозиции.
} PARBOUND;


// Структура данных используемая для распараллеливания.
typedef struct TPARDATA {
	integer ncore; // 1, 2, 4, 8.
	integer *inumerate=NULL;
	// это для ncore==2;
	PARBOUND b0;
	// это для ncore==4;
	PARBOUND b00, b01;
	// это для ncore==8;
	PARBOUND b000, b001, b010, b011;
} PARDATA;

PARDATA nd;




// используется в Румбе для ускорения счёта.
doublereal* rthdsd_no_radiosity_patch = NULL;



#include "adaptive_local_refinement_mesh.cpp" // АЛИС
#include "constr_struct.cpp" // заполнение структур данных TEMPER и FLOW

// количество блоков, источников и стенок.
integer lb, ls, lw;
BLOCK* b = NULL; // блоков
SOURCE* s = NULL; // источников
WALL* w = NULL; // твёрдых стенок

doublereal *xpos = NULL, *ypos = NULL, *zpos = NULL;
doublereal *xposadd = NULL, *yposadd = NULL, *zposadd = NULL;

#include "my_LR.c" // полилинейный метод

#include "my_material_properties.c" // библиотека реальных свойств материалов

// аппроксимация обобщённого уравнения конвекции-диффузии
// на совмещённой сетке
#include "pamendment3.c"


#include "shortest_distance.cpp" // вычисление кратчайшего расстояния до стенки

// 8 января 2016.
const bool bvery_big_memory = true; // true нет записи в файл всё храним в оперативной памяти. Это существенно быстрее по скорости.

struct Tdatabase {
	doublereal *x=NULL, *y=NULL, *z=NULL; // координаты узлов.
	integer maxelm;
	integer** nvtxcell=NULL;
	integer ncell;
	// связь теплопередачи с гидродинамикой.
	integer **ptr=NULL;// для тестирования алгебраического многосеточного метода
};

Tdatabase database;

TEMPER t;
integer flow_interior; // Суммарное число FLUID зон
FLOW* f=NULL;

// экспорт картинки в программу tecplot360
#include "my_export_tecplot3.c"

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
} FLUENT_RESIDUAL;


// соединение решателя
#include "mysolverv0_03.c"




// нестационарный солвер для температуры
// на основе стационарного солвера,
// а также нестационарный солвер для 
// гидродинамики на основе стационарного солвера.
#include "my_unsteady_temperature.c"

// Препроцессинг для параллельной обработки.
#include "my_nested_dissection.cpp"

#include <ctime> // для замера времени выполнения.




integer ltdp; // количество таблично заданных мощностей от температуры и смещения стока.
TEMP_DEP_POWER* gtdps=NULL; // Garber temperature depend power sequence. 

// база данных материалов:
integer lmatmax; // максимальное число материалов
TPROP* matlist=NULL; // хранилище базы данных материалов





doublereal rterminate_residual_ICCG_Oh2(FLOW floc) {
	// точность дискретизации O(h!2)
	doublereal* resterm = new doublereal[floc.maxelm + floc.maxbound];
	for (integer i = 0; i<floc.maxelm + floc.maxbound; i++) {
		resterm[i] = 0.0; // инициализация.
	}

	for (integer iP = 0; iP<floc.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(iP, floc.nvtx, floc.pa, dx, dy, dz);
		doublereal dl = fmin(dx, fmin(dy, dz));
		resterm[iP] = 0.1*dl*dl; // O(h!2)
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
		iE = floc.sosedi[ESIDE][iP].iNODE1; iN = floc.sosedi[NSIDE][iP].iNODE1; iT = floc.sosedi[TSIDE][iP].iNODE1; iW = floc.sosedi[WSIDE][iP].iNODE1; iS = floc.sosedi[SSIDE][iP].iNODE1; iB = floc.sosedi[BSIDE][iP].iNODE1;
		// Если с одной из сторон стоит граница расчётной области
		// то соответствующая переменная равна true
		bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

		if (iE >= floc.maxelm) bE = true;
		if (iN >= floc.maxelm) bN = true;
		if (iT >= floc.maxelm) bT = true;
		if (iW >= floc.maxelm) bW = true;
		if (iS >= floc.maxelm) bS = true;
		if (iB >= floc.maxelm) bB = true;

		if ((bE) || (bW)) {
			dl = 0.5*dx;
			if (bE) resterm[iE] = 0.1*dl*dl; // O(h!2)
			if (bW) resterm[iW] = 0.1*dl*dl; // O(h!2)
		}
		if ((bN) || (bS)) {
			dl = 0.5*dy;
			if (bN) resterm[iN] = 0.1*dl*dl; // O(h!2)
			if (bS) resterm[iS] = 0.1*dl*dl; // O(h!2)
		}
		if ((bT) || (bB)) {
			dl = 0.5*dz;
			if (bT) resterm[iT] = 0.1*dl*dl; // O(h!2)
			if (bB) resterm[iB] = 0.1*dl*dl; // O(h!2)
		}
	}
	doublereal ret = Scal(resterm, resterm, floc.maxelm + floc.maxbound);
	delete[] resterm;
	resterm=NULL;
	return ret;
} // rterminate_residual_ICCG_Oh2

doublereal rterminate_residual_LR1sk_Oh3(FLOW floc) {
	// точность дискретизации O(h!2)
	doublereal* resterm = new doublereal[floc.maxelm + floc.maxbound];
	for (integer i = 0; i<floc.maxelm + floc.maxbound; i++) {
		resterm[i] = 0.0; // инициализация.
	}

	for (integer iP = 0; iP<floc.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(iP, floc.nvtx, floc.pa, dx, dy, dz);
		doublereal dl = fmin(dx, fmin(dy, dz));
		resterm[iP] = 0.1*dl*dl*dl; // O(h!3)
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
		iE = floc.sosedi[ESIDE][iP].iNODE1; iN = floc.sosedi[NSIDE][iP].iNODE1; iT = floc.sosedi[TSIDE][iP].iNODE1; iW = floc.sosedi[WSIDE][iP].iNODE1; iS = floc.sosedi[SSIDE][iP].iNODE1; iB = floc.sosedi[BSIDE][iP].iNODE1;
		// Если с одной из сторон стоит граница расчётной области
		// то соответствующая переменная равна true
		bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

		if (iE >= floc.maxelm) bE = true;
		if (iN >= floc.maxelm) bN = true;
		if (iT >= floc.maxelm) bT = true;
		if (iW >= floc.maxelm) bW = true;
		if (iS >= floc.maxelm) bS = true;
		if (iB >= floc.maxelm) bB = true;

		if ((bE) || (bW)) {
			dl = 0.5*dx;
			if (bE) resterm[iE] = 0.1*dl*dl*dl; // O(h!3)
			if (bW) resterm[iW] = 0.1*dl*dl*dl; // O(h!3)
		}
		if ((bN) || (bS)) {
			dl = 0.5*dy;
			if (bN) resterm[iN] = 0.1*dl*dl*dl; // O(h!3)
			if (bS) resterm[iS] = 0.1*dl*dl*dl; // O(h!3)
		}
		if ((bT) || (bB)) {
			dl = 0.5*dz;
			if (bT) resterm[iT] = 0.1*dl*dl*dl; // O(h!3)
			if (bB) resterm[iB] = 0.1*dl*dl*dl; // O(h!3)
		}
	}
	doublereal ret;
	//ret=Scal(resterm,resterm,floc.maxelm+floc.maxbound);
	ret = NormaV(resterm, floc.maxelm + floc.maxbound);
	// Освобождение оперативной памяти.
	if (resterm != NULL) {
		delete[] resterm;
		resterm = NULL;
	}
	return ret;
} // rterminate_residual_LR1sk_Oh3

doublereal rterminate_residual_LR1sk_temp_Oh3(TEMPER t) {
	// точность дискретизации O(h!2)
	// точность дискретизации теплопередачи в твёрдом теле O(h).
	doublereal* resterm = new doublereal[t.maxelm + t.maxbound];
	for (integer i = 0; i<t.maxelm + t.maxbound; i++) {
		resterm[i] = 0.0; // инициализация.
	}

	for (integer iP = 0; iP<t.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
		doublereal dl = fmin(dx, fmin(dy, dz));
		//resterm[iP]=0.1*dl*dl*dl; // O(h!3)
		resterm[iP] = dl; // O(h)
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
		iE = t.sosedi[ESIDE][iP].iNODE1; iN = t.sosedi[NSIDE][iP].iNODE1; iT = t.sosedi[TSIDE][iP].iNODE1; iW = t.sosedi[WSIDE][iP].iNODE1; iS = t.sosedi[SSIDE][iP].iNODE1; iB = t.sosedi[BSIDE][iP].iNODE1;
		// Если с одной из сторон стоит граница расчётной области
		// то соответствующая переменная равна true
		bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

		if (iE >= t.maxelm) bE = true;
		if (iN >= t.maxelm) bN = true;
		if (iT >= t.maxelm) bT = true;
		if (iW >= t.maxelm) bW = true;
		if (iS >= t.maxelm) bS = true;
		if (iB >= t.maxelm) bB = true;

		if ((bE) || (bW)) {
			dl = 0.5*dx;
			//if (bE) resterm[iE]=0.1*dl*dl*dl; // O(h!3)
			//if (bW) resterm[iW]=0.1*dl*dl*dl; // O(h!3)
			if (bE) resterm[iE] = dl; // O(h)
			if (bW) resterm[iW] = dl; // O(h)

		}
		if ((bN) || (bS)) {
			dl = 0.5*dy;
			//if (bN) resterm[iN]=0.1*dl*dl*dl; // O(h!3)
			//if (bS) resterm[iS]=0.1*dl*dl*dl; // O(h!3)
			if (bN) resterm[iN] = dl; // O(h)
			if (bS) resterm[iS] = dl; // O(h)

		}
		if ((bT) || (bB)) {
			dl = 0.5*dz;
			//if (bT) resterm[iT]=0.1*dl*dl*dl; // O(h!3)
			//if (bB) resterm[iB]=0.1*dl*dl*dl; // O(h!3)
			if (bT) resterm[iT] = dl; // O(h)
			if (bB) resterm[iB] = dl; // O(h)
		}
	}
	doublereal ret;
	//ret=Scal(resterm,resterm,f.maxelm+f.maxbound);
	ret = NormaV(resterm, t.maxelm + t.maxbound);
	// Освобождение оперативной памяти.
	if (resterm != NULL) {
		delete[] resterm;
		resterm = NULL;
	}
	return ret;
} // rterminate_residual_LR1sk_temp_Oh3



int main(void)
{
	getchar();



	

	// 22.01.2017 инициализация.
	eqin.fluidinfo = NULL;
	t.rootBT = NULL;
	t.rootSN = NULL;
	t.rootWE = NULL;

	// 29 10 2016.
	// Инициализация общей памяти в ILU буффере.
	milu_gl_buffer.alu_copy = NULL;
	milu_gl_buffer.jlu_copy = NULL;
	milu_gl_buffer.ju_copy = NULL;
	
	//fgmres(m_restart)
	my_amg_manager.m_restart = 20; // Количество итераций алгоритма fgmres перед перезапуском.

	// amg default settings:
	my_amg_manager.lfil = 2; // default value

	// Параметры собственного многосеточного метода о умолчанию.
	// Настройки решателя СЛАУ зависят от типа уравнения которое подаётся на вход:
	// симметричность <-> анизотропность, несимметричность, диффузионная задача <-> конвективная задача,
	// степень преобладания конвекции (число Рейнольдса).
	//my_amg_manager.maximum_levels = 20; // максимальное число уровней начиная с которого начинается усечение.
	my_amg_manager.maximum_delete_levels = 0; // Количество уровней отсекаемых в нижней части где грубая сетка.
	my_amg_manager.number_interpolation_procedure = 3; // номер интерполяционной процедуры.
	my_amg_manager.number_interpolation_procedure_Temperature = 3;
	my_amg_manager.number_interpolation_procedure_Speed = 3;
	my_amg_manager.number_interpolation_procedure_Pressure = 3;
	my_amg_manager.number_interpolation_procedure_Stress = 3;

	my_amg_manager.iCFalgorithm_and_data_structure=3; // 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Temperature=3;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Speed=3;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Pressure=3;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Stress=3;// 3-Treap.

	my_amg_manager.bTemperatureMatrixPortrait = 0; // NO_PRINT
	my_amg_manager.bSpeedMatrixPortrait = 0; // NO_PRINT
	my_amg_manager.bPressureMatrixPortrait = 0; // NO_PRINT
	my_amg_manager.bStressMatrixPortrait = 0; // NO_PRINT
	my_amg_manager.bMatrixPortrait = 0; // NO_PRINT


	my_amg_manager.nFinnest = 2; // число итераций на подробной сетке.
	my_amg_manager.nu1 = 1; // число предсглаживаний.
	my_amg_manager.nu2 = 2; // число пост сглаживаий.	
	my_amg_manager.memory_size = 9; // количество оперативной памяти в размерностях матрицы А.
	my_amg_manager.gold_const = 0.2; // Параметр верхней релаксации в сглаживателе.
	my_amg_manager.gold_const_Temperature = 0.2;
	my_amg_manager.gold_const_Speed = 0.2;
	my_amg_manager.gold_const_Pressure = 0.2;
	my_amg_manager.gold_const_Stress = 0.2;
	my_amg_manager.bCFJacoby = true; // CF-Jacobi smoothers 12% сокращение числа V циклов. 5.06.2017
	// Runge-Kutt smoother: 3 - третьего порядка, 5 - пятого порядка, любое другое число не используется.
	my_amg_manager.iRunge_Kutta_smoother = 0;
	my_amg_manager.iFinnest_ilu = 0; // 0 - не используется, 1 - используется, но только на самой подробной сетке.
	// Использование iluk разложения на глубоких уровнях вложенности для которых
	// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
	my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6 = false;
	my_amg_manager.theta = 0.24;
	my_amg_manager.magic = 0.4;
	my_amg_manager.F_to_F_Temperature = 0.4;
	my_amg_manager.F_to_F_Speed = 0.4; 
	my_amg_manager.F_to_F_Pressure = 0.4;
	my_amg_manager.F_to_F_Stress = 0.4;
	my_amg_manager.ilu2_smoother = 0; // 0 - не использовать, 1 - использовать.

	my_amg_manager.itypemodifyinterpol=0; // номер модификации интерполляции.
	my_amg_manager.inumberadaptpass=0; // максимальное количество сканов-проходов с модификациями.

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
	// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
	// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
	my_amg_manager.icoarseningTemp = 0; // standart
	my_amg_manager.icoarseningSpeed = 0; // standart
	my_amg_manager.icoarseningPressure=0; // standart
	my_amg_manager.icoarseningStress = 0; // standart
	my_amg_manager.icoarseningtype=0; // standart vs RS 2.
	// Stabilization BiCGStab.
	// 8.01.2017 Метод ван дер Ворста BiCGStab 
	// предобусловленный алгебраичесеким многосеточным методом.
	// 0 - используется просто алгебраический многосеточный метод без какого-либо привлечения алгоритмов подпространства Крылова,
	// 1 - Используется алгоритм Х. Ван дер Ворста BiCGStab [1992], предобусловленный алгебраическим многосеточным методом.
	// 2 - Используется алгоритм Саада и Шульца FGMRes [1986], предобусловленный алгебраическим многосеточным методом.
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
	my_amg_manager.iprint_log_Pressure=1;
	my_amg_manager.iprint_log_Stress = 1;

	// truncation for interpolation.
	// По умолчанию усечение интерполяции не используется.
	my_amg_manager.itruncation_interpolation = 0; // 0 - off
	my_amg_manager.itruncation_interpolation_Temperature = 0;
	my_amg_manager.itruncation_interpolation_Speed = 0;
	my_amg_manager.itruncation_interpolation_Pressure=0;
	my_amg_manager.itruncation_interpolation_Stress = 0;
	// 0.2 recomended Stuben.
	my_amg_manager.truncation_interpolation = 0.2; // 0.2 recomended default value.
	my_amg_manager.truncation_interpolation_Temperature = 0.2;
	my_amg_manager.truncation_interpolation_Speed = 0.2;
	my_amg_manager.truncation_interpolation_Pressure=0.2;
	my_amg_manager.truncation_interpolation_Stress = 0.2;

	// GMRES smoother.
	my_amg_manager.b_gmresTemp = false;
	my_amg_manager.b_gmresSpeed = false;
	my_amg_manager.b_gmresPressure=false;
	my_amg_manager.b_gmresStress = false;
	my_amg_manager.b_gmres=false;

	// Замер времени.
	unsigned int calculation_main_start_time = 0; // начало счёта мс.
	unsigned int calculation_main_end_time = 0; // окончание счёта мс.
	unsigned int calculation_main_seach_time = 0; // время выполнения участка кода в мс.

	calculation_main_start_time = clock(); // момент начала счёта.

	bool bextendedprint = false; // печать на граничных узлах расчитанных полей.

	// Диагностическая печать для двойного вакуумного промежутка.
	err_radiation_log = fopen_s(&fp_radiation_log, "log_radiation.txt", "a");
	if (err_radiation_log != 0) {
		printf("Error open file log.txt\n");
		printf("Please, press any key to continue...\n");
		//getchar();
		system("pause");
		exit(0);
	}

	//std::locale::global(std::locale("en_US.UTF-8"));
	system("mode con cols=126 lines=12000");
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
		

	printf("AliceFlow 3D x64 v0_07\n");
#ifdef _OPENMP 
	omp_set_num_threads(inumcore); // установка числа потоков
#endif

	errno_t err;
	err = fopen_s(&fp_log, "log.txt", "w");
	if (err != 0) {
		printf("Error open file log.txt\n");
		printf("Please, press any key to continue...\n");
		//getchar();
		system("pause");
		exit(0);
	}

	if (fp_log != NULL) {

		//ilu0_Saadtest();
		//printf("the end Saad ilu0 test\n");
		//getchar();

		// количество точек по каждой из осей.
		//integer inx=120, iny=64, inz=64;
		integer inx = 30, iny = 30, inz = 30;
		integer inxadd = -1, inyadd = -1, inzadd = -1;
		doublereal dgx = 0.0, dgy = 0.0, dgz = 0.0; // сила тяжести
		doublereal operatingtemperature = 20.0; // Operating Temperature 20.0 Град. С.

		

		premeshin("premeshin.txt", lmatmax, lb, ls, lw, matlist, b, s, w, dgx, dgy, dgz, inx, iny, inz, operatingtemperature,  ltdp, gtdps);
		if (iswitchMeshGenerator == 0) {
			simplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
		}
		else if (iswitchMeshGenerator == 1) {
			unevensimplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, matlist, dgx, dgy, dgz, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd); // генератор неравномерной сетки с автоматической балансировкой неравномерности.
		}
		else if (iswitchMeshGenerator == 2) {
			// Я стремился сделать coarse Mesh как в Icepak.
			// Реальные модели чрезвычайно многоэлементно-большеразмерные, а 
			// ресурсы персонального компьютера чрезвычайно слабы, т.к. cpu уперлись в 4ГГц а
			// правильное распараллеливание отдельная большая научная проблема.
			coarsemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
		}
		else {
#if doubleintprecision == 1
			printf("error : yuor mesh generator is undefined %lld\n", iswitchMeshGenerator);
#else
			printf("error : yuor mesh generator is undefined %d\n", iswitchMeshGenerator);
#endif
			
			system("pause");
			exit(1);
		}


		if (b_on_adaptive_local_refinement_mesh) {
			printf("starting ALICE\n");
			// Это слишком большая величина на многих промышленных задачах,
			// её нельзя задавать во избежании сбоя в операторе new или malloc.
			integer maxelm_loc = (inx + 1)*(iny + 1)*(inz + 1);
			bool bOkal=alice_mesh(xpos, ypos, zpos, inx, iny, inz, b, lb, lw, w, s, ls, maxelm_loc, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
			//system("PAUSE");
			
			if (itype_ALICE_Mesh == 1) {
				// Вызываем повторные генерации.

				while (!bOkal) {
					/* 3.09.2017
					АЛИС сетку лучше строить когда при дроблении ячейки соответствующая геометрическая длина делится
					ровно пополам. Раньше делился пополам целочисленный индекс в результате чего на нашей существенно неравномерной
					сетке АЛИС ячейки были очень сильно вытянутые. В результате точность аппроксимации чрезвычайно страдала. Теперь когда
					дробиться пополам именно геометрическая длина может нехватать ячеек сетки для балансировки сетки (т.е. будет невозможно
					без добавления новых сеточных линий сделать чтобы уровни соседних ячеек отличались не более чем на 1. Поэтому теперь
					такие невозможные для дробления ячейки ищутся в функции if_disbalnce(...) и для каждой такой ячейки принимается решение
					добавить соответствующую новую сеточную линию. Алгоритм освобождает память и возращается на исходные позиции и построение
					АЛИС сетки начинается заново только теперь базовая сетка уже содержит недостающие сеточные линии.
					*/


					// Нужно освободить память из под octree дерева и перестроить сетку.
					printf("free octree start...\n");
					//getchar();
					//system("PAUSE");
					free_octree(oc_global);
					delete[] my_ALICE_STACK;
					top_ALICE_STACK = 0;
					printf("free octree end...\n");
					// Новое построение расчётной сетки.
					delete[] xpos;
					xpos = NULL;
					inx = 0;
					delete[] ypos;
					ypos = NULL;
					iny = 0;
					delete[] zpos;
					zpos = NULL;
					inz = 0;

					printf("free xpos, ypos, zpos\n");
					//system("PAUSE");

					if (iswitchMeshGenerator == 0) {
						simplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
					}
					else if (iswitchMeshGenerator == 1) {
						unevensimplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, matlist, dgx, dgy, dgz, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd); // генератор неравномерной сетки с автоматической балансировкой неравномерности.
					}
					else if (iswitchMeshGenerator == 2) {
						// Я стремился сделать coarse Mesh как в Icepak.
						// Реальные модели чрезвычайно многоэлементно-большеразмерные, а 
						// ресурсы персонального компьютера чрезвычайно слабы, т.к. cpu уперлись в 4ГГц а
						// правильное распараллеливание отдельная большая научная проблема.
						coarsemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
					}
					else {
#if doubleintprecision == 1
						printf("error : yuor mesh generator is undefined %lld\n", iswitchMeshGenerator);
#else
						printf("error : yuor mesh generator is undefined %d\n", iswitchMeshGenerator);
#endif

						system("pause");
						exit(1);
					}
					// Новый запуск АЛИС меширования.

					printf("new construct xpos, ypos, zpos\n");
					//system("PAUSE");

					bOkal = alice_mesh(xpos, ypos, zpos, inx, iny, inz, b, lb, lw, w, s, ls, maxelm_loc, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);

					//getchar();
					//system("PAUSE");
				}
			}
			printf("end ALICE\n");
		}

		load_TEMPER_and_FLOW(t, f, inx, iny, inz, xpos, ypos, zpos, flow_interior, b, lb, lw, w, s, ls, operatingtemperature, matlist, bextendedprint, dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false);
		
		// Эти копии данных нужны для полного восстановления данных.
		t.inx_copy = inx;
		t.iny_copy = iny;
		t.inz_copy = inz;
		t.operatingtemperature_copy = operatingtemperature;
		t.xpos_copy = new doublereal[inx + 1];
		t.ypos_copy = new doublereal[iny + 1];
		t.zpos_copy = new doublereal[inz + 1];
		// Данная информармация нужна для экономии оперативной памяти,
		// некоторые данные будут выгружены из озу а потом восстановлены 
		// путём вычиления.
		for (integer i_7 = 0; i_7 < inx + 1; i_7++) {
			t.xpos_copy[i_7] = xpos[i_7];
		}
		for (integer i_7 = 0; i_7 < iny + 1; i_7++) {
			t.ypos_copy[i_7] = ypos[i_7];
		}
		for (integer i_7 = 0; i_7 < inz + 1; i_7++) {
			t.zpos_copy[i_7] = zpos[i_7];
		}

		// Освобождение оперативной памяти из под octree дерева.
		if (b_on_adaptive_local_refinement_mesh) {
			printf("free octree start...\n");
			//getchar();
			//system("PAUSE");
			free_octree(oc_global);
			delete[] my_ALICE_STACK;
			top_ALICE_STACK = 0;
			printf("free octree end...\n");
			//getchar();
			//system("PAUSE");
		}

		if (0) {
			xyplot(f, flow_interior, t);
			printf("after load temper and flow. OK.\n");
			//getchar(); // debug avtosave
			system("pause");
		}

		// На АЛИС сетке полинейный метод не является работоспособным.
		if (!b_on_adaptive_local_refinement_mesh) {
			// создаёт информацию о сеточных линиях для полилинейного метода LR:
			constr_line(f, flow_interior);  // для гидродинамики
			t.rootBT = NULL;
			t.rootSN = NULL;
			t.rootWE = NULL;
			constr_line_temp(t, b, lb); // для теплопроводности
			printf("LR preprocessing finish...\n");
		}

		// глобальная память для алгебраического многосеточного метода.

		amgGM.a = NULL;
		amgGM.f = NULL;
		amgGM.ia = NULL;
		amgGM.ig = NULL;
		amgGM.ja = NULL;
		amgGM.u = NULL;
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
		if (0 && (flow_interior == 1)) {
			calc_front(f, f[0], t, flow_interior, ls, lw, w, nd);
			// разделение выполнено !
			printf("separator compleate...\n");
			//getchar();
		}



		t.free_temper_level1 = false; // чистая теплопроводность освобождение памяти необходимой для сборки матрицы после успешной сборки.
		t.free_temper_level2 = false; // освобождение памяти под хранение матрицы при перезаписи её в SIMPLESPARSE формат.	

		printf("construction of all structures...\n");
		printf("mesh check start...\n");
#if doubleintprecision == 1
		for (integer i = 0; i < inx; i++) if (fabs(xpos[i + 1] - xpos[i]) < 1.0e-23)
			printf("error: zalipanie po X: xpos[%lld]=%e xpos[%lld]=%e inx=%lld\n", i, xpos[i], i + 1, xpos[i + 1], inx);
		for (integer i = 0; i < iny; i++) if (fabs(ypos[i + 1] - ypos[i]) < 1.0e-23)
			printf("error: zalipanie po X: ypos[%lld]=%e ypos[%lld]=%e iny=%lld\n", i, ypos[i], i + 1, ypos[i + 1], iny);
		for (integer i = 0; i < inz; i++) if (fabs(zpos[i + 1] - zpos[i]) < 1.0e-23)
			printf("error: zalipanie po X: zpos[%lld]=%e zpos[%lld]=%e inz=%lld\n", i, zpos[i], i + 1, zpos[i + 1], inz);
		for (integer iP = 0; iP < t.maxelm; iP++) {
			if ((t.nvtx[0][iP] == 0) || (t.nvtx[1][iP] == 0) || (t.nvtx[2][iP] == 0) || (t.nvtx[3][iP] == 0) || (t.nvtx[4][iP] == 0) || (t.nvtx[5][iP] == 0) || (t.nvtx[6][iP] == 0) || (t.nvtx[7][iP] == 0)) {
				printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, t.nvtx[0][iP] - 1, t.nvtx[1][iP] - 1, t.nvtx[2][iP] - 1, t.nvtx[3][iP] - 1, t.nvtx[4][iP] - 1, t.nvtx[5][iP] - 1, t.nvtx[6][iP] - 1, t.nvtx[7][iP] - 1);
			}
		}
#else
		for (integer i = 0; i < inx; i++) if (fabs(xpos[i + 1] - xpos[i]) < 1.0e-23)
			printf("error: zalipanie po X: xpos[%d]=%e xpos[%d]=%e inx=%d\n", i, xpos[i], i + 1, xpos[i + 1], inx);
		for (integer i = 0; i < iny; i++) if (fabs(ypos[i + 1] - ypos[i]) < 1.0e-23)
			printf("error: zalipanie po X: ypos[%d]=%e ypos[%d]=%e iny=%d\n", i, ypos[i], i + 1, ypos[i + 1], iny);
		for (integer i = 0; i < inz; i++) if (fabs(zpos[i + 1] - zpos[i]) < 1.0e-23)
			printf("error: zalipanie po X: zpos[%d]=%e zpos[%d]=%e inz=%d\n", i, zpos[i], i + 1, zpos[i + 1], inz);
		for (integer iP = 0; iP < t.maxelm; iP++) {
			if ((t.nvtx[0][iP] == 0) || (t.nvtx[1][iP] == 0) || (t.nvtx[2][iP] == 0) || (t.nvtx[3][iP] == 0) || (t.nvtx[4][iP] == 0) || (t.nvtx[5][iP] == 0) || (t.nvtx[6][iP] == 0) || (t.nvtx[7][iP] == 0)) {
				printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, t.nvtx[0][iP] - 1, t.nvtx[1][iP] - 1, t.nvtx[2][iP] - 1, t.nvtx[3][iP] - 1, t.nvtx[4][iP] - 1, t.nvtx[5][iP] - 1, t.nvtx[6][iP] - 1, t.nvtx[7][iP] - 1);
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
			printf("residual O(h!2) is equal=%e\n", f[i].resICCG);
			f[i].resLR1sk = rterminate_residual_LR1sk_Oh3(f[i]); // O(h!3)
			printf("residual O(h!3) is equal=%e\n", f[i].resLR1sk);
		}
		printf("TEMPERATURE\n");
		t.resLR1sk = rterminate_residual_LR1sk_temp_Oh3(t); // O(h!3)		
		printf("temp residual O(h!3) is equal=%e\n", t.resLR1sk);
		printf("mesh check.\n");
		if (bwait) {
			//getchar();
			system("pause");
		}
		
		sourse2Dproblem = new bool[t.maxbound];
		conductivity2Dinsource = new doublereal[t.maxbound];
		// Нижняя релаксация источникового члена при радиационых потоках.
		bsource_term_radiation_for_relax = new doublereal[t.maxelm];
		for (integer i_init = 0; i_init < t.maxelm; i_init++) bsource_term_radiation_for_relax[i_init] = 0.0;
		b_buffer_correct_source = new doublereal[t.maxelm];


		// невязка continity будет измеряться по отношению к уровню 1e0.
		doublereal* continity_start = NULL;
		continity_start = new doublereal[flow_interior];
		if (continity_start == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for continity start in main...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}
		for (integer i = 0; i < flow_interior; i++) continity_start[i] = 1.0;

		integer* inumber_iteration_SIMPLE = NULL;
		inumber_iteration_SIMPLE = new integer[flow_interior];
		if (inumber_iteration_SIMPLE == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for inumber_iteration_SIMPLE in main...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}
		for (integer i = 0; i < flow_interior; i++) inumber_iteration_SIMPLE[i] = 0; // начальная итерация алгоритма SIMPLE для каждой FLUID зоны.

		// считывание состояния расчёта из файла для возобновления расчёта
		bool breadOk = false;
		avtoreadvalue(f, t, flow_interior, inumber_iteration_SIMPLE, continity_start, breadOk);
		// Если считывание прошло неуспешно то breadOk==false и это значит что счёт начнётся заново со значений заданных при инициализации.

		if (b_on_adaptive_local_refinement_mesh) {
			// Инвариант корректности АЛИС сетки.
			printf("the invariant correctness...\n");
			ANES_ALICE_CORRECT(t.maxnod, t.pa, t.maxelm, t.nvtx);
		}

		// экспорт результата вычисления в программу tecplot360:
		// можно использовать как проверку построенной сетки.
		if (0) {
			exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,0);
			printf("read values. OK.\n");
			//getchar(); // debug avtosave
			system("pause");
		}



		/*for (integer i=0; i<lw; i++) {
		printf("%e  \n",w[i].Tamb);
		}
		//exporttecplotxy360T_3D(t.maxelm, t.ncell, t.nvtx, t.nvtxcell, t.pa, t.potent);
		exporttecplotxy360T_3D_part2(t.maxelm, t.potent);
		getchar(); // debug
		*/

		// 29.01.2017
		// if (1 && steady_or_unsteady_global_determinant == 2)  
		// То мы просто вызываем мешер не вызывая солвера.
		if (1 && steady_or_unsteady_global_determinant == 2) {
			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

			// Вычисление массы модели.
			massa_cabinet(t, f, inx, iny, inz,
				xpos, ypos, zpos, flow_interior,
				b, lb, operatingtemperature,
				matlist);

			calculation_end_time = clock(); // момент окончания счёта.
			calculation_seach_time = calculation_end_time - calculation_start_time;
			unsigned int im = 0, is = 0, ims = 0;
			im = (unsigned int)(calculation_seach_time / 60000); // минуты
			is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
			ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

			printf("time export to tecplot360 is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

			// Экспорт сетки в tecplot 360.
			if (1) {
				if (!b_on_adaptive_local_refinement_mesh) {
					// экспорт результата вычисления в программу tecplot360:
					exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
				}
				else {
					// Экспорт в программу техплот температуры.
					//С АЛИС сетки.
					ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t,0);
				}
			}

			
		}

		// steady
		if (1 && steady_or_unsteady_global_determinant == 0) {

			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

			for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = operating_temperature_for_film_coeff; // инициализация.

			// Решаем только теплопередачу в твёрдом теле.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка :
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			if (adiabatic_vs_heat_transfer_coeff == 1) {
				// Мы инициализируем случайной величиной чтобы избавится от неопределённости при первой сборке СЛАУ.
				//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
			}

			// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
			bsolid_static_only = true;
			bool bcleantemp = false;
			if (eqin.itemper == 1) {
				bcleantemp = true;
				integer i = 0; // счётчик цикла
				for (i = 0; i < flow_interior; i++) {
					if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
				}
				// если bcleantemp==true то мы решаем задачу чистой теплопередачи без учёта конвекции.
			}
			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой теплопроводности в трёхмерной области. 
				printf("solution of pure heat...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//getchar();
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx*starting_speed_Vx + starting_speed_Vy*starting_speed_Vy + starting_speed_Vz*starting_speed_Vz > 1.0e-30) {
					bmyconvective = true;
				}
				else {
					// Загрузка распределения начальной скорости.
					errno_t err_inicialization_data;
					FILE* fp_inicialization_data;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
					if (err_inicialization_data == 0) {
						// Открытие удачно и файл присутствует.
						bmyconvective = true;
						fclose(fp_inicialization_data);
					}
				}

				// if (flow_interior>0) bmyconvective=true;
				// массив отладочной информации,
				// конкретно для проверки подхода Рхи-Чоу
				doublereal **rhie_chow = NULL;
				QuickMemVorst m;
				m.ballocCRSt = false; // Выделять память
				m.bsignalfreeCRSt = true; // и сразу освобождать.
				// инициализация указателей.
				m.tval = NULL;
				m.tcol_ind = NULL;
				m.trow_ptr = NULL;
				m.tri = NULL;
				m.troc = NULL;
				m.ts = NULL;
				m.tt = NULL;
				m.tvi = NULL;
				m.tpi = NULL;
				m.tdx = NULL;
				m.tdax = NULL;
				m.ty = NULL;
				m.tz = NULL;
				m.ta = NULL;
				m.tja = NULL;
				m.tia = NULL;
				m.talu = NULL;
				m.tjlu = NULL;
				m.tju = NULL;
				m.tiw = NULL;
				m.tlevs = NULL;
				m.tw = NULL;
				m.tjw = NULL;
				m.icount_vel = 100000; // очень большое число.			
				
				
				
				// если flow_interior == 0 то f[0] просто формальный параметр  
				solve_nonlinear_temp(f[0], f, t,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					bmyconvective, NULL, 0.001, 0.001,
					false,
					matlist, 0,
					bprintmessage,
					gtdps, ltdp, 1.0, m,
					NULL, // скорость с предыдущего временного слоя. 
					NULL); // массовый поток через границу с предыдущего временного слоя.
				// последний параметр равный 1.0 означает что мощность подаётся.
				
				// Вычисление массы модели.
				massa_cabinet(t, f, inx, iny, inz,
					xpos, ypos, zpos, flow_interior,
					b, lb, operatingtemperature,
					matlist);

				// 10.10.2017
				// Построение двумерного графика вдоль структуры.
				xyplot_temp(t, t.potent);
				//printf("graphics writing sucseful\n");
				//getchar();

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,0);
					}
					else {
						// Экспорт в программу техплот температуры.
						//С АЛИС сетки.
						//ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent,t,0);
						
					}
				}
				//printf("marker stop\n");
				//getchar();
			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
			// Вычисление значения максимальной температуры внутри расчётной области и на её границах:
			//for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
			// 23 декабря 2015
			// На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
			// физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < t.maxelm; i++) tmaxfinish = fmax(tmaxfinish, t.potent[i]);
			FILE *fp;
			errno_t err1;
			err1 = fopen_s(&fp, "report.txt", "w");
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//getchar();
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			report_temperature(flow_interior, f, t, b, lb, s, ls, w, lw, 0);

			if (1) {

				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

			    // На реальной большразмерной модели (11.4М ячеек на АЛИС Coarse )
				// преобразует неприемлемо долго. 

				// 25.11.2017
				// 1. Получили данные о температуре и сетке на АЛИС.
				// 2. Сохранили их в оперативной памяти.
				// 3. Освободили память.
				// 4. Построили обычную декартовую прямоугольную сетку.
				// 5. Перенесли данные о температуре с АЛИС на обычную декартовую прямоугольную сетку.
				// 6. Визуализировали температуру и построенные на структурированной сетке тепловые потоки.
				// 7.1 Техплот идеально визуализирует то что есть на структурированной сетке и в сечении и в объёме. 
				// 7.2 Визуализация на АЛИС сетке в объеме даже идеально точно найденного поля не удовлетворитьельна (глюк техплота).


				if (b_on_adaptive_local_refinement_mesh) {
					// 1. Получение x,y,z,T,nvtx, m_sizeT, m_size_nvtx.
					doublereal *x_buf = NULL;
					doublereal *y_buf = NULL;
					doublereal *z_buf = NULL;
					doublereal *t_buf = NULL;
					integer **nvtx_buf = NULL;
					integer m_sizeT = 0, m_size_nvtx = 0;

					ANES_tecplot360_export_temperature_preobrazovatel(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx);

					// 2. Освобождение памяти.
					// Освобождение оперативной памяти.
					if (t.xpos_copy != NULL) {
						delete[] t.xpos_copy;
						t.xpos_copy = NULL;
					}
					if (t.ypos_copy != NULL) {
						delete[] t.ypos_copy;
						t.ypos_copy = NULL;
					}
					if (t.zpos_copy != NULL) {
						delete[] t.zpos_copy;
						t.zpos_copy = NULL;
					}


					if (bsource_term_radiation_for_relax != NULL) {
						delete[] bsource_term_radiation_for_relax; // Релаксация источниковых членов радиационных потоков.
						bsource_term_radiation_for_relax = NULL;
					}
					if (b_buffer_correct_source != NULL) {
						delete[] b_buffer_correct_source;
						b_buffer_correct_source = NULL;
					}

					if (rthdsd_no_radiosity_patch != NULL) {
						free(rthdsd_no_radiosity_patch);
					}
					rthdsd_no_radiosity_patch = NULL;

					// Быстрая обработка нелинейных граничных условий в Румба 0.14 решателе.
					if (qnbc != NULL) {
						delete[] qnbc;
						qnbc = NULL;
						iadd_qnbc_maxelm = 0;
					}

					// Нужно освободить оперативную память из под всех структур данных:
					free_level1_temp(t);
					free_level2_temp(t); // освобождение памяти из под матриц.
										 // Освобождает память для LR начало.
					if (t.rootWE != NULL) {
						free_root(t.rootWE, t.iWE);
					}
					if (t.rootSN != NULL) {
						free_root(t.rootSN, t.iSN);
					}
					if (t.rootBT != NULL) {
						free_root(t.rootBT, t.iBT);
					}
					if (t.rootWE != NULL) {
						delete[] t.rootWE;
						t.rootWE = NULL;
					}
					if (t.rootSN != NULL) {
						delete[] t.rootSN;
						t.rootSN = NULL;
					}
					if (t.rootBT != NULL) {
						delete[] t.rootBT;
						t.rootBT = NULL;
					}
					// Освобождение памяти для LR конец.
					free_level1_flow(f, flow_interior);
					free_level2_flow(f, flow_interior); // освобождение памяти из под матриц.

					if (sourse2Dproblem != NULL) {
						delete[] sourse2Dproblem;
						sourse2Dproblem = NULL;
					}
					if (conductivity2Dinsource != NULL) {
						delete[] conductivity2Dinsource;
						conductivity2Dinsource = NULL;
					}

					if (x_jacoby_buffer != NULL) {
						// 30 октября 2016. 
						// В seidelsor2 сделан переключатель на метод нижней релаксации К.Г. Якоби.
						// Освобождение памяти из под jacobi buffer.
						delete[] x_jacoby_buffer;
					}

					if (bvery_big_memory) {
						if (database.x != NULL) {
							free(database.x);
						}
						if (database.y != NULL) {
							free(database.y);
						}
						if (database.z != NULL) {
							free(database.z);
						}
						if (database.nvtxcell != NULL) {
							for (integer i = 0; i <= 7; i++) {
								delete[] database.nvtxcell[i];
							}
							delete[] database.nvtxcell;
						}
						if (database.ptr != NULL) {
							if (database.ptr[0] != NULL) {
								delete[] database.ptr[0];
							}
							if (database.ptr[1] != NULL) {
								delete[] database.ptr[1];
							}
							delete[] database.ptr;
						}
					}
					/*
					// Освобождение общей памяти в ILU буффере.
					if (milu_gl_buffer.alu_copy != NULL) delete[] milu_gl_buffer.alu_copy;
					if (milu_gl_buffer.jlu_copy != NULL) delete[] milu_gl_buffer.jlu_copy;
					if (milu_gl_buffer.ju_copy != NULL) delete[] milu_gl_buffer.ju_copy;
					milu_gl_buffer.alu_copy = NULL;
					milu_gl_buffer.jlu_copy = NULL;
					milu_gl_buffer.ju_copy = NULL;
					*/
					flow_interior = 0;

					// 3. Построение обычной сетки.

					b_on_adaptive_local_refinement_mesh = false;
					load_TEMPER_and_FLOW(t, f, inx, iny, inz, xpos, ypos, zpos, flow_interior, b, lb, lw, w, s, ls, operatingtemperature, matlist, bextendedprint, dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false);

					// Эти копии данных нужны для полного восстановления данных.
					t.inx_copy = inx;
					t.iny_copy = iny;
					t.inz_copy = inz;
					t.operatingtemperature_copy = operatingtemperature;
					t.xpos_copy = new doublereal[inx + 1];
					t.ypos_copy = new doublereal[iny + 1];
					t.zpos_copy = new doublereal[inz + 1];
					// Данная информармация нужна для экономии оперативной памяти,
					// некоторые данные будут выгружены из озу а потом восстановлены 
					// путём вычиления.
					for (integer i_7 = 0; i_7 < inx + 1; i_7++) {
						t.xpos_copy[i_7] = xpos[i_7];
					}
					for (integer i_7 = 0; i_7 < iny + 1; i_7++) {
						t.ypos_copy[i_7] = ypos[i_7];
					}
					for (integer i_7 = 0; i_7 < inz + 1; i_7++) {
						t.zpos_copy[i_7] = zpos[i_7];
					}

					t.free_temper_level1 = false; // чистая теплопроводность освобождение памяти необходимой для сборки матрицы после успешной сборки.
					t.free_temper_level2 = false; // освобождение памяти под хранение матрицы при перезаписи её в SIMPLESPARSE формат.	


					// 4. Интерполляция для температуры.
					ALICE_2_Structural(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx);


					if (x_buf != NULL) {
						delete[] x_buf;
						x_buf = NULL;
					}
					if (y_buf != NULL) {
						delete[] y_buf;
						y_buf = NULL;
					}
					if (z_buf != NULL) {
						delete[] z_buf;
						z_buf = NULL;
					}
					if (t_buf != NULL) {
						delete[] t_buf;
						t_buf = NULL;
					}
					if (nvtx_buf != NULL) {
						for (integer i_1 = 0; i_1 < 8; i_1++) {
							if (nvtx_buf[i_1] != NULL) {
								delete[] nvtx_buf[i_1];
								nvtx_buf[i_1] = NULL;
							}
						}
						delete[] nvtx_buf;
						nvtx_buf = NULL;
					}
					m_sizeT = 0, m_size_nvtx = 0;
					// 5. Обычный экспорт в техплот.
					exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
				}
			}
			else {
				// Объёмная визуализация на АЛИС очень плохая даже с точным полем, поэтому
				// для качественной объёмной визуализации нужен переход на структурированную сетку.
				// В сечении на АЛИС сетке все хорошо по техплоту.
				// Точность на АЛИС очень плоха. Можно допиться хорошей картинки для поля температур, но
				// плотности тепловых потоков никогда не будут найдены (представлены точно) только в случае
				// патологического случая дико мелкой сетки когда АЛИС очень подробная (большеэлементные модели).

				// Экспорт в программу техплот температуры.
				// С АЛИС сетки.
				//ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent,t,0);
			}

			

		}

		// steady Static Structural
		if (1 && steady_or_unsteady_global_determinant == 5) {

			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

			for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = operating_temperature_for_film_coeff; // инициализация.

		    // Решаем только Static Structural.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка :
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			if (adiabatic_vs_heat_transfer_coeff == 1) {
				// Мы инициализируем случайной величиной чтобы избавится от неопределённости при первой сборке СЛАУ.
				//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
			}

			// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
			bsolid_static_only = true;
			bool bcleantemp = false;
			if (eqin.itemper == 1) {
				bcleantemp = true;
				integer i = 0; // счётчик цикла
				for (i = 0; i < flow_interior; i++) {
					if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
				}
				// если bcleantemp==true то мы решаем задачу чистой теплопередачи без учёта конвекции.
			}
			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой теплопроводности в трёхмерной области. 
				printf("solution of pure Static Structural...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//getchar();
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx*starting_speed_Vx + starting_speed_Vy*starting_speed_Vy + starting_speed_Vz*starting_speed_Vz > 1.0e-30) {
					bmyconvective = true;
				}
				else {
					// Загрузка распределения начальной скорости.
					errno_t err_inicialization_data;
					FILE* fp_inicialization_data;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
					if (err_inicialization_data == 0) {
						// Открытие удачно и файл присутствует.
						bmyconvective = true;
						fclose(fp_inicialization_data);
					}
				}

				// if (flow_interior>0) bmyconvective=true;
				// массив отладочной информации,
				// конкретно для проверки подхода Рхи-Чоу
				doublereal **rhie_chow = NULL;
				QuickMemVorst m;
				m.ballocCRSt = false; // Выделять память
				m.bsignalfreeCRSt = true; // и сразу освобождать.
										  // инициализация указателей.
				m.tval = NULL;
				m.tcol_ind = NULL;
				m.trow_ptr = NULL;
				m.tri = NULL;
				m.troc = NULL;
				m.ts = NULL;
				m.tt = NULL;
				m.tvi = NULL;
				m.tpi = NULL;
				m.tdx = NULL;
				m.tdax = NULL;
				m.ty = NULL;
				m.tz = NULL;
				m.ta = NULL;
				m.tja = NULL;
				m.tia = NULL;
				m.talu = NULL;
				m.tjlu = NULL;
				m.tju = NULL;
				m.tiw = NULL;
				m.tlevs = NULL;
				m.tw = NULL;
				m.tjw = NULL;
				m.icount_vel = 100000; // очень большое число.

				bPhysics_stop = false;
				//bPhysics_stop = false;
				// Вызов солвера Static Structural.
				solve_Structural(t, w, lw, m, false, operatingtemperature);
				//bPhysics_stop = true;

				/*
				// если flow_interior == 0 то f[0] просто формальный параметр
				solve_nonlinear_temp(f[0], f, t,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				bmyconvective, NULL, 0.001, 0.001,
				false,
				matlist, 0,
				bprintmessage,
				gtdps, ltdp, 1.0, m,
				NULL, // скорость с предыдущего временного слоя.
				NULL); // массовый поток через границу с предыдущего временного слоя.
				// последний параметр равный 1.0 означает что мощность подаётся.
				*/
				// Вычисление массы модели.
				massa_cabinet(t, f, inx, iny, inz,
					xpos, ypos, zpos, flow_interior,
					b, lb, operatingtemperature,
					matlist);


				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
					}
					else {
						// Экспорт в программу техплот температуры.
						//С АЛИС сетки.
						ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t,0);
					}
				}

			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
											 // Вычисление значения максимальной температуры внутри расчётной области и на её границах:
											 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
											 // 23 декабря 2015
											 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
											 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < t.maxelm; i++) tmaxfinish = fmax(tmaxfinish, t.potent[i]);

			doublereal totaldeform_max = -1.0e+30;
			for (integer i = 0; i < t.maxelm; i++) totaldeform_max = fmax(totaldeform_max, t.total_deformation[TOTALDEFORMATION][i]);

			FILE *fp;
			errno_t err1;
			err1 = fopen_s(&fp, "report.txt", "w");
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//getchar();
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			report_temperature(flow_interior, f, t, b, lb, s, ls, w, lw, 0);

			

		}

		// steady Static Structural and Temperature (Thermal Stress).
		if (1 && steady_or_unsteady_global_determinant == 6) {

			// Замер времени.
			unsigned int calculation_start_time = 0; // начало счёта мс.
			unsigned int calculation_end_time = 0; // окончание счёта мс.
			unsigned int calculation_seach_time = 0; // время выполнения участка кода в мс.

			calculation_start_time = clock(); // момент начала счёта.

			for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = operating_temperature_for_film_coeff; // инициализация.

			// Решаем теплопередачу, а затем Static Structural.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка :
					// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
					bPhysics_PTBSH_memory = true;
				}
			}

			if (adiabatic_vs_heat_transfer_coeff == 1) {
				// Мы инициализируем случайной величиной чтобы избавится от неопределённости при первой сборке СЛАУ.
				//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
			}

			// Здесь предполагается что мы решаем стационарную задачу чистой теплопроводности.
			bsolid_static_only = true;
			bool bcleantemp = false;
			if (eqin.itemper == 1) {
				bcleantemp = true;
				integer i = 0; // счётчик цикла
				for (i = 0; i < flow_interior; i++) {
					if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
				}
				// если bcleantemp==true то мы решаем задачу чистой теплопередачи без учёта конвекции.
			}
			if (1 || bcleantemp) {
				// решение стационарной нелинейной (или линейной) задачи чистой теплопроводности в трёхмерной области. 
				printf("solution of pure Static Structural...\n");
				printf("please, press any key to continue...\n");
				if (bwait) {
					//getchar();
					system("pause");
				}

				// при тестировании рекомендуется обязательно печатать.
				bool bprintmessage = true; // печатать ли сообщения на консоль.

				doublereal dbeta = 1.0; // первый порядок аппроксимации на границе.
				bool bmyconvective = false;
				if (starting_speed_Vx*starting_speed_Vx + starting_speed_Vy*starting_speed_Vy + starting_speed_Vz*starting_speed_Vz > 1.0e-30) {
					bmyconvective = true;
				}
				else {
					// Загрузка распределения начальной скорости.
					errno_t err_inicialization_data;
					FILE* fp_inicialization_data;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
					if (err_inicialization_data == 0) {
						// Открытие удачно и файл присутствует.
						bmyconvective = true;
						fclose(fp_inicialization_data);
					}
				}

				// if (flow_interior>0) bmyconvective=true;
				// массив отладочной информации,
				// конкретно для проверки подхода Рхи-Чоу
				doublereal **rhie_chow = NULL;
				QuickMemVorst m;
				m.ballocCRSt = false; // Выделять память
				m.bsignalfreeCRSt = true; // и сразу освобождать.
										  // инициализация указателей.
				m.tval = NULL;
				m.tcol_ind = NULL;
				m.trow_ptr = NULL;
				m.tri = NULL;
				m.troc = NULL;
				m.ts = NULL;
				m.tt = NULL;
				m.tvi = NULL;
				m.tpi = NULL;
				m.tdx = NULL;
				m.tdax = NULL;
				m.ty = NULL;
				m.tz = NULL;
				m.ta = NULL;
				m.tja = NULL;
				m.tia = NULL;
				m.talu = NULL;
				m.tjlu = NULL;
				m.tju = NULL;
				m.tiw = NULL;
				m.tlevs = NULL;
				m.tw = NULL;
				m.tjw = NULL;
				m.icount_vel = 100000; // очень большое число.

				


				
				// если flow_interior == 0 то f[0] просто формальный параметр
				solve_nonlinear_temp(f[0], f, t,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				bmyconvective, NULL, 0.001, 0.001,
				false,
				matlist, 0,
				bprintmessage,
				gtdps, ltdp, 1.0, m,
				NULL, // скорость с предыдущего временного слоя.
				NULL); // массовый поток через границу с предыдущего временного слоя.
				// последний параметр равный 1.0 означает что мощность подаётся.
				
				//bPhysics_stop = false;
				// Вызов солвера Static Structural.
				solve_Structural(t, w, lw, m, true, operatingtemperature);
				//bPhysics_stop = true;

				// Вычисление массы модели.
				massa_cabinet(t, f, inx, iny, inz,
					xpos, ypos, zpos, flow_interior,
					b, lb, operatingtemperature,
					matlist);

				calculation_end_time = clock(); // момент окончания счёта.
				calculation_seach_time = calculation_end_time - calculation_start_time;
				unsigned int im = 0, is = 0, ims = 0;
				im = (unsigned int)(calculation_seach_time / 60000); // минуты
				is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // секунды
				ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

				printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

				if (1) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// экспорт результата вычисления в программу tecplot360:
						exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
					}
					else {
						// Экспорт в программу техплот температуры.
						//С АЛИС сетки.
						ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t,0);
					}
				}

			}

			doublereal tmaxfinish = -273.15; // абсолютный ноль.
											 // Вычисление значения максимальной температуры внутри расчётной области и на её границах:
											 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
											 // 23 декабря 2015
											 // На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
											 // физичнее не смущать людей и приводить температуру только во внутренних КО. 
			for (integer i = 0; i < t.maxelm; i++) tmaxfinish = fmax(tmaxfinish, t.potent[i]);

			doublereal totaldeform_max = -1.0e+30;
			for (integer i = 0; i < t.maxelm; i++) totaldeform_max = fmax(totaldeform_max, t.total_deformation[TOTALDEFORMATION][i]);

			FILE *fp;
			errno_t err1;
			err1 = fopen_s(&fp, "report.txt", "w");
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//getchar();
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			report_temperature(flow_interior, f, t, b, lb, s, ls, w, lw, 0);

			

		}



		/*
		if (b_on_adaptive_local_refinement_mesh) {
		printf("t.maxbound=%d\n", t.maxbound);
		printf("v dvuch shagah ot ALICE sborki. \n");
		getchar();
		exit(1);  // Режим глобального тестирования сеточного генератора адаптивной локально измельчённой сетки.
		}

		if (b_on_adaptive_local_refinement_mesh) {
		printf("Solve temperature is compleate. \n");
		getchar();
		exit(1);  // Режим глобального тестирования сеточного генератора адаптивной локально измельчённой сетки.
		}
		*/

		//system("pause");

		if (1 && steady_or_unsteady_global_determinant == 1) {
			// Нестационарная теплопроводность.
			
			// Решаем только теплопередачу в твёрдом теле.
			bonly_solid_calculation = true;

			// Включаем прекращение вычисления по физическому смыслу.
			// Это прерывание работает только на ПТБШ позволяя ускорить вычисления.
			if (lw == 1) {
				bPhysics_stop = true;
				if (lb < 11) {
					// Это стандартная подложка :
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
			// даёт изменение температуры в нужную сторону на 5 градусов на фоне перегрева в 167градусов в статике.
			// Для сравнения перегрев в icepak равен 120 градусам, что близко к значениям паспортного Rt 
			// (мощность 6.875Вт, RT=16K/W) и экспериментальным
			// данным.
			doublereal dbeta = 1.3333333;//1.0; // если 1.0 то первый порядок аппроксимации на границе.
			dbeta = 1.0; // более стабильное значение.
			// массив отладочной информации,
			// конкретно для проверки подхода Рхи-Чоу
			doublereal **rhie_chow = NULL;
			//solve_nonlinear_temp(f[0], f, t, rhie_chow, b, lb, s, ls, w, lw, dbeta, flow_interior, false, NULL, 0.001, false);
			unsteady_temperature_calculation(f[0], f, t,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				false, matlist,
				operatingtemperature,
				gtdps, ltdp); // нестационарный температурный солвер

			// Вычисление массы модели.
			massa_cabinet(t, f, inx, iny, inz,
				xpos, ypos, zpos, flow_interior,
				b, lb, operatingtemperature,
				matlist);

			if (!b_on_adaptive_local_refinement_mesh) {
				// экспорт результата вычисления в программу tecplot360:
				exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,0);
			}
			else {
				// Экспорт в программу техплот температуры.
				//С АЛИС сетки.
				ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent,t,0);
			}

			doublereal tmaxfinish = -273.15;
			// Вычисление значения максимальной температуры внутри расчётной области и на её границах:
			for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, t.potent[i]);
			FILE *fp;
			errno_t err1;
			err1 = fopen_s(&fp, "report.txt", "w");
			// создание файла для записи.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//getchar();
				system("pause");
			}
			else {
				// запись заголовка
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			report_temperature(flow_interior, f, t, b, lb, s, ls, w, lw, 0);

			printf("calculation complete...\n");
			// getchar();
		}

		fclose(fp_radiation_log);

		// экспорт результата вычисления в программу tecplot360:
		// можно использовать как проверку построенной сетки.
		if (false) {
			exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,0);
			printf("read values. OK.\n");
			if (bwait) {
				//getchar(); // debug avtosave
				system("pause");
			}
		}

		if ((1 && steady_or_unsteady_global_determinant == 3)) {
			// Fluid dynamic.

			told_temperature_global_for_HOrelax = new doublereal[t.maxelm + t.maxbound];
			bSIMPLErun_now_for_temperature = true;
			if (dgx*dgx + dgy*dgy + dgz*dgz > 1.0e-20) {
				// надо также проверить включено ли для fluid материалов приближение Буссинеска.
				bool bbussinesk_7 = false;
#pragma omp parallel for
				for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
					integer ib = t.whot_is_block[f[0].ptr[i_8]];
					if (ib > -1) {
						if (b[ib].itype == FLUID) {
							integer i_7 = b[ib].imatid;
							if (matlist[i_7].bBussineskApproach) bbussinesk_7 = true;
						}
					}
				}
				if (bbussinesk_7) {					
					bSIMPLErun_now_for_natural_convection = true;
				}
			}
			bHORF = true;
			bPamendment_source_old = new doublereal[f[0].maxelm + f[0].maxbound];
			for (integer i5 = 0; i5 < f[0].maxelm + f[0].maxbound; i5++) bPamendment_source_old[i5] = 0.0;
			// exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint);
			//getchar();
			if (dgx*dgx + dgy*dgy + dgz*dgz > 1.0e-20) {
				// надо также проверить включено ли для fluid материалов приближение Буссинеска.
				bool bbussinesk_7 = false;
#pragma omp parallel for
				for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
					integer ib = t.whot_is_block[f[0].ptr[i_8]];
					if (ib > -1) {
						if (b[ib].itype == FLUID) {
							integer i_7 = b[ib].imatid;
							if (matlist[i_7].bBussineskApproach) {
								bbussinesk_7 = true;
							}
						}
					}
				}
				if (bbussinesk_7) {
					printf("Bussinesk approach Operating Temperature=%e\n", f[0].OpTemp); // Operating Temperature);
				}
			}

			//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
			//getchar();

			// стационарный гидродинамический решатель.
			steady_cfd_calculation(breadOk,
				eqin, dgx, dgy, dgz,
				continity_start,
				inumber_iteration_SIMPLE,
				flow_interior, f, t, b, lb,
				s, ls, w, lw, matlist,
				gtdps, ltdp, bextendedprint);
			//xyplot( f, 0, t);
			// boundarylayer_info(f, t, flow_interior, w, lw);
			// 2 - solver/conjugate_heat_transfer_static/
			report_temperature(flow_interior, f, t, b, lb, s, ls, w, lw, 0/*2*/);

			// Вычисление массы модели.
			massa_cabinet(t, f, inx, iny, inz,
				xpos, ypos, zpos, flow_interior,
				b, lb, operatingtemperature,
				matlist);

			// экспорт результата вычисления в программу tecplot360:
			exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,0);
			save_velocity_for_init(t.maxelm, t.ncell, f, t, flow_interior);
			// exporttecplotxy360T_3D_part2_rev(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,b,lb);
			delete[] bPamendment_source_old;
			delete[] told_temperature_global_for_HOrelax;
		}
		if (0) {
			

			told_temperature_global_for_HOrelax = new doublereal[t.maxelm + t.maxbound];
			bSIMPLErun_now_for_temperature = true;
			bHORF = true;
			bPamendment_source_old = new doublereal[f[0].maxelm + f[0].maxbound];
			for (integer i5 = 0; i5 < f[0].maxelm + f[0].maxbound; i5++) bPamendment_source_old[i5] = 0.0;
			// нестационарный гидродинамический решатель :
			usteady_cfd_calculation(breadOk, eqin,
				dgx, dgy, dgz,
				continity_start,
				inumber_iteration_SIMPLE,
				flow_interior,
				f, t,
				b, lb, s, ls,
				w, lw, matlist, gtdps, ltdp, bextendedprint);
			delete[] bPamendment_source_old;
			delete[] told_temperature_global_for_HOrelax;

			// Вычисление массы модели.
			massa_cabinet(t, f, inx, iny, inz,
				xpos, ypos, zpos, flow_interior,
				b, lb, operatingtemperature,
				matlist);
		}

		fclose(fp_log); // закрытие файла лога.


		if (continity_start != NULL) {
			delete[] continity_start;
			continity_start = NULL;
		}

		if (inumber_iteration_SIMPLE != NULL) {
			delete[] inumber_iteration_SIMPLE;
			inumber_iteration_SIMPLE = NULL;
		}

	}

	// Освобождение оперативной памяти.
	if (t.xpos_copy != NULL) {
		delete[] t.xpos_copy;
		t.xpos_copy = NULL;
	}
	if (t.ypos_copy != NULL) {
		delete[] t.ypos_copy;
		t.ypos_copy = NULL;
	}
	if (t.zpos_copy != NULL) {
		delete[] t.zpos_copy;
		t.zpos_copy = NULL;
	}
	
	// Освобождение оперативной памяти.
	if (xposadd != NULL) {
		delete[] xposadd;
		xposadd = NULL;
	}
	if (yposadd != NULL) {
		delete[] yposadd;
		yposadd = NULL;
	}
	if (zposadd != NULL) {
		delete[] zposadd;
		zposadd = NULL;
	}

	
	// Освобождение оперативной памяти.
	if (xpos != NULL) {
		delete[] xpos;
		xpos = NULL;
	}
	if (ypos != NULL) {
		delete[] ypos;
		ypos = NULL;
	}
	if (zpos != NULL) {
		delete[] zpos;
		zpos = NULL;
	}

	if (bsource_term_radiation_for_relax != NULL) {
		delete[] bsource_term_radiation_for_relax; // Релаксация источниковых членов радиационных потоков.
		bsource_term_radiation_for_relax = NULL;
	}
	if (b_buffer_correct_source != NULL) {
		delete[] b_buffer_correct_source;
		b_buffer_correct_source = NULL;
	}

	printf("free memory begin...\n");
	if (bwait) {
		//getchar();
		system("pause");
	}

	if (rthdsd_no_radiosity_patch != NULL) {
		free(rthdsd_no_radiosity_patch);
	}
	rthdsd_no_radiosity_patch = NULL;

	// Быстрая обработка нелинейных граничных условий в Румба 0.14 решателе.
	if (qnbc != NULL) {
		delete[] qnbc;
		qnbc = NULL;
		iadd_qnbc_maxelm = 0;
	}

	// освобождение памяти из под amg1r5.
	if (amgGM.a != NULL) {
		delete amgGM.a;
		amgGM.a = NULL;
	}
	if (amgGM.ia != NULL) {
		delete amgGM.ia;
		amgGM.ia = NULL;
	}
	if (amgGM.ja != NULL) {
		delete amgGM.ja;
		amgGM.ja = NULL;
	}
	if (amgGM.u != NULL) {
		delete amgGM.u;
		amgGM.u = NULL;
	}
	if (amgGM.f != NULL) {
		delete amgGM.f;
		amgGM.f = NULL;
	}
	if (amgGM.ig != NULL) {
		delete amgGM.ig;
		amgGM.ig = NULL;
	}

	amgGM.nda = -1;
	amgGM.ndf = -1;
	amgGM.ndia = -1;
	amgGM.ndig = -1;
	amgGM.ndja = -1;
	amgGM.ndu = -1;

	for (integer i_7 = 0; i_7 < lb; i_7++) {
		if (b[i_7].temp_Sc != NULL) {
			delete[] b[i_7].temp_Sc;
			b[i_7].temp_Sc = NULL;
		}
		if (b[i_7].arr_Sc != NULL) {
			delete[] b[i_7].arr_Sc;
			b[i_7].arr_Sc = NULL;
		}
		if (b[i_7].g.hi != NULL) {
			delete[] b[i_7].g.hi;
			b[i_7].g.hi = NULL;
		}
		if (b[i_7].g.xi != NULL) {
			delete[] b[i_7].g.xi;
			b[i_7].g.xi = NULL;
		}
		if (b[i_7].g.yi != NULL) {
			delete[] b[i_7].g.yi;
			b[i_7].g.yi = NULL;
		}
		if (b[i_7].g.zi != NULL) {
			delete[] b[i_7].g.zi;
			b[i_7].g.zi = NULL;
		}
	}
	delete b; delete s; delete w; // освобождение памяти
	for (integer i_7 = 0; i_7 < lmatmax; i_7++) {
		if (matlist[i_7].arr_cp != NULL) {
			delete[] matlist[i_7].arr_cp;
			matlist[i_7].arr_cp = NULL;
     	}
		if (matlist[i_7].temp_cp != NULL) {
			delete[] matlist[i_7].temp_cp;
			matlist[i_7].temp_cp = NULL;
		}
		if (matlist[i_7].arr_lam != NULL) {
			delete[] matlist[i_7].arr_lam;
			matlist[i_7].arr_lam = NULL;
		}
		if (matlist[i_7].temp_lam != NULL) {
			delete[] matlist[i_7].temp_lam;
			matlist[i_7].temp_lam = NULL;
		}
	}
	delete[] matlist;
	delete[] gtdps;
	if (eqin.fluidinfo != NULL) {
		delete eqin.fluidinfo;
		eqin.fluidinfo = NULL;
	}
	// Нужно освободить оперативную память из под всех структур данных:
	free_level1_temp(t);
	free_level2_temp(t); // освобождение памяти из под матриц.
	// Освобождает память для LR начало.
	free_root(t.rootWE, t.iWE);
	free_root(t.rootSN, t.iSN);
	free_root(t.rootBT, t.iBT);
	if (t.rootWE != NULL) {
		delete[] t.rootWE;
		t.rootWE = NULL;
	}
	if (t.rootSN != NULL) {
		delete[] t.rootSN;
		t.rootSN = NULL;
	}
	if (t.rootBT != NULL) {
		delete[] t.rootBT;
		t.rootBT = NULL;
	}
	// Освобождение памяти для LR конец.
	free_level1_flow(f, flow_interior);
	free_level2_flow(f, flow_interior); // освобождение памяти из под матриц.

	delete[] f;
	f = NULL;

	if (sourse2Dproblem != NULL) {
		delete[] sourse2Dproblem;
		sourse2Dproblem = NULL;
	}
	if (conductivity2Dinsource != NULL) {
		delete[] conductivity2Dinsource;
		conductivity2Dinsource = NULL;
	}

	if (x_jacoby_buffer != NULL) {
		// 30 октября 2016. 
		// В seidelsor2 сделан переключатель на метод нижней релаксации К.Г. Якоби.
		// Освобождение памяти из под jacobi buffer.
		delete[] x_jacoby_buffer;
	}

	if (bvery_big_memory) {
		if (database.x != NULL) {
			free(database.x);
		}
		if (database.y != NULL) {
			free(database.y);
		}
		if (database.z != NULL) {
			free(database.z);
		}
		if (database.nvtxcell != NULL) {
			for (integer i = 0; i <= 7; i++) {
				delete[] database.nvtxcell[i];
			}
			delete[] database.nvtxcell;
		}
		if (database.ptr != NULL) {
			if (database.ptr[0] != NULL) {
				delete[] database.ptr[0];
			}
			if (database.ptr[1] != NULL) {
				delete[] database.ptr[1];
			}
			delete[] database.ptr;
		}
	}

	// Освобождение общей памяти в ILU буффере.
	if (milu_gl_buffer.alu_copy != NULL) delete[] milu_gl_buffer.alu_copy;
	if (milu_gl_buffer.jlu_copy != NULL) delete[] milu_gl_buffer.jlu_copy;
	if (milu_gl_buffer.ju_copy != NULL) delete[] milu_gl_buffer.ju_copy;
	milu_gl_buffer.alu_copy = NULL;
	milu_gl_buffer.jlu_copy = NULL;
	milu_gl_buffer.ju_copy = NULL;



	flow_interior = 0;
	printf("free memory finish...\n");

	if (1 && steady_or_unsteady_global_determinant == 2) {
		// Был только вызов сеточного генератора.
		printf("Mesh generation procedure is finish.\n");
	}
	else {
		printf("Calculation procedure is finish.\n");
	}
	printf("Please, press any key to exit...\n");
	if (bwait) {
		//getchar();
		system("pause");
	}

	

	calculation_main_end_time = clock();
	calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time;


	/*printf("time=%d statistic vorst=%3.2f %% \n",calculation_main_seach_time,(float)(100.0*calculation_vorst_seach_time/calculation_main_seach_time));
	getchar();
	*/

	
	
	// Общее время вычисления.
	int im=0, is=0, ims=0;
	im=(int)(calculation_main_seach_time/60000); // минуты
	is=(int)((calculation_main_seach_time-60000*im)/1000); // секунды
	ims=(int)((calculation_main_seach_time-60000*im-1000*is)/10); // миллисекунды делённые на 10
	
	printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10*ims);


	system("pause");
	return 0;
}
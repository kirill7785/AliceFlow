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
const doublereal rdivision_interval = 4.0;// 4.0 default 10.0
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
// при повторном обращении идет обращение к ячейки хеш таблицы.
// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
// Хеш таблицы для automatic
const integer isize_shorter_hash = 10000000;// Сделан значительный запас точности хеширования.
doublereal* shorter_hash_X = NULL;
doublereal* shorter_hash_Y = NULL;
doublereal* shorter_hash_Z = NULL;
bool* bshorter_hash_X = NULL;
bool* bshorter_hash_Y = NULL;
bool* bshorter_hash_Z = NULL;
// Параметры настройки 
const doublereal dcabinet_delta_multiplyer = 0.005; // 1/200
const doublereal dbody_dist_multiplyer = 12.5;
const doublereal dbody_R_out_multiplyer = 2.6;
const doublereal dbody_Hcyl_multiplyer = 10.0;
const doublereal dbody_LGATE_situation_multiplyer = 5.2;


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

// Мощность рассеиваемая в тепло зависит от
// температуры и смещения стока.
// Информация о зависимости хранится в таблично заданной 
// функция. Таблица взята из результатов расчёта в програме
// Г.З. Гарбера.
typedef struct TTEMP_DEP_POWER {
	char* sname = NULL; // имя текстового файла содержащего табличную функцию.

	integer intemp = 0; // количество различных дискретных значений температуры.
	integer inoffset_drain = 0; // количество различных дискретных хначений смещения стока.

	// таблица: левый столбец температуры,
	// верхняя строка без первого значения - смещения стока.
	// Остальные значения в таблице соответствуют рассеиваемой мощности.
	doublereal* rtemp = NULL; // значения температуры
	doublereal* roffset_drain = NULL; // значения смещения стока
	doublereal** rpower_table = NULL; // таблица мощностей.

} TEMP_DEP_POWER;

// источник тепла
typedef struct TSOURCE {
	// мощность источника тепла
	// в случае если мощность зависит от максимальной температуры и смещения стока то
	// это значение автоматически пересчитывается перед сборкой матрицы, 
	// так что в коде сборки граничного условия ничего менять не нужно.
	doublereal power = 0.0;
	// Мощность может зависеть от температуры и смещения стока.
	// Зависимость задаётся таблично на основе расчётов в программе Г.З. Гарбера.
	bool bgarber_depend = false; // задана ли зависимость рассеиваемой мощности от температуры и смещения стока.
	// В случае bgarber_depend == true, значение power является коэффициентом на который домножается 
	// истинное табличное значение мощности.
	integer igarber_depend = -1; // уникальный номер табличной зависимости. (нумерация начинается с нуля).
	doublereal roperation_offset_drain = 28.0; // рабочее значение смещения стока.
	doublereal power_multiplyer = 1.0; // на эту константу домножается реальная мощность, это удобно для корректирования
	// табличной мощности например когда нужно считать четверть транзистора.

	doublereal square = 0.0; // площадь источника тепла
	integer iPlane = XZ; // плоскость в которой лежит источник тепла
	// 1 - XY, 2 - XZ, 3 - YZ.
	GEOM g; // границы объекта.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id = 0;
} SOURCE;

// стенка (идеальный теплооотвод)
typedef struct TWALL {
	// 1 - Дирихле, 2 - однородное условие Неймана.
	integer ifamily = 2; // род краевого условия по температуре
	doublereal Tamb = 20.0, hf = 0.0; // значение температуры на идеальном теплоотводе и тепловой поток.
	doublereal emissivity = 0.8, film_coefficient = 3.0;
	bool bsymmetry = false, bpressure = false, bopening = false;
	doublereal Vx = 0.0, Vy = 0.0, Vz = 0.0, P = 0.0;
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
	integer ithermal_Stress_boundary_condition = 0;
	doublereal xForce = 0.0, yForce = 0.0, zForce = 0.0; // Три компоненты силы, приложенные к узлу в Ньютонах.
	integer iPlane = XZ; // плоскость в которой лежит стенка
	// 1 - XY, 2 - XZ, 3 - YZ.
	GEOM g; // границы объекта.
	// Фиксированная граница с нулевым смещением для
	// расчёта прочности конструкции.
	//bool bfixboundary; // true фиксированная, false свободная.

	// принадлежность объединению 
	// 0 - не принадлежит (принадлежит кабинету),
	// n > 0 принадлежит объединению с номером n.
	integer iunion_id = 0;
} WALL;



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

// Управление алгебраическим многосеточным методом из интерфейса.
typedef struct TMY_AMG_MANAGER {

	// Алгоритм сортировки используемый в многосеточном методе РУМБА.
	// 0 - Counting Sort, 1 - QUICKSORT, 3 - HEAPSORT.
	integer imySortAlgorithm = 0; // 0 - Counting Sort Default.

	// 0 - не печатать портрет матрицы, 
	// 1 - печатать портрет матрицы.
	integer bTemperatureMatrixPortrait = 0;
	integer bSpeedMatrixPortrait = 0;
	integer bPressureMatrixPortrait = 0;
	integer bStressMatrixPortrait = 0;
	integer bMatrixPortrait = 0;

	

	// Temperature
	doublereal theta_Temperature = 0.24;
	integer maximum_delete_levels_Temperature = 0;
	integer nFinnest_Temperature = 2, nu1_Temperature = 1, nu2_Temperature = 2;
	integer memory_size_Temperature = 13; //13*size(matrix A)
	integer ilu2_smoother_Temperature = 0; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Temperature = 2;
	// Speed
	doublereal theta_Speed = 0.24;
	integer maximum_delete_levels_Speed = 0;
	integer nFinnest_Speed = 2, nu1_Speed = 1, nu2_Speed = 2;
	integer memory_size_Speed = 13;
	integer ilu2_smoother_Speed = 0; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Speed = 2;
	// Pressure
	doublereal theta_Pressure = 0.24;
	integer maximum_delete_levels_Pressure = 0;
	integer nFinnest_Pressure = 2, nu1_Pressure = 1, nu2_Pressure = 2;
	integer memory_size_Pressure = 15;
	integer ilu2_smoother_Pressure = 0; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Pressure = 2;
	// Stress
	doublereal theta_Stress = 0.24;
	integer maximum_delete_levels_Stress = 0;
	integer nFinnest_Stress = 2, nu1_Stress = 1, nu2_Stress = 2;
	integer memory_size_Stress = 22;
	integer ilu2_smoother_Stress = 0; // 0 - не использовать, 1 - использовать.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 3.
	integer iCFalgorithm_and_data_structure_Stress = 2;
	// global
	bool bCFJacoby = true;
	integer iRunge_Kutta_smoother = 0; // 3 - третьего порядка, 5 - пятого порядка, любое другое число не используется. 
	integer iFinnest_ilu = 0; // 0 не используется, 1 - ilu0. Только на самой подробной сетке.
	// Использование iluk разложения на глубоких уровнях вложенности для которых
	// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
	bool b_ilu_smoothers_in_nnz_n_LE_6 = false;
	doublereal theta = 0.24; // strength threshold
	//integer maximum_levels; // максимальное количество уровней вложенности (уровни выше редуцируются).
	integer maximum_delete_levels = 0; // Количество уровней отсекаемых снизу в области грубой сетки.
	integer nFinnest = 2, nu1 = 1, nu2 = 2; // Количества сглаживаний.
	// 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap, 3 - Treap.
	// default - 2.
	integer iCFalgorithm_and_data_structure = 2;
	integer memory_size = 13; // В размерах матрицы А.
	// Для метода верхней релаксации в сглаживателе.
	const integer AMG1R5_IN_HOUSE = 1;// Собственная реализация интерполяции amg1r5.
	integer number_interpolation_procedure = AMG1R5_IN_HOUSE; // идентификатор процедуры интерполляции.
	integer number_interpolation_procedure_Temperature = AMG1R5_IN_HOUSE;
	integer number_interpolation_procedure_Speed = AMG1R5_IN_HOUSE;
	integer number_interpolation_procedure_Pressure = AMG1R5_IN_HOUSE;
	integer number_interpolation_procedure_Stress = AMG1R5_IN_HOUSE;

	// 6 december 2016.
	// Подлежит удалению Refactoring.
	//integer itypemodifyinterpol=0; // номер модификации интерполляции. // Подлежит удалению Refactoring.
	//integer inumberadaptpass=0; // максимальное количество сканов-проходов с модификациями. // Подлежит удалению Refactoring.
	//integer baglomeration_with_consistency_scaling = 0;
	// Принудительное усиление диагонали
	// в случае обнаружения внедиагональных
	// positive connections.
	integer bdiagonal_dominant = 1;


	// Параметры сглаживателя SOR Роуч.
	doublereal gold_const = 0.24;
	doublereal gold_const_Temperature = 0.24;
	doublereal gold_const_Speed = 0.24;
	doublereal gold_const_Pressure = 0.24;
	doublereal gold_const_Stress = 0.24;

	// Пороги отсечек
	// Отсечка же называется theta threshold. 
		
	doublereal magic = 0.4;
	doublereal F_to_F_Temperature = 0.4, F_to_F_Speed = 0.4, F_to_F_Pressure = 0.4, F_to_F_Stress = 0.4; // magic
	integer ilu2_smoother = 0; // 0 - не использовать, 1 - использовать.
	// AMG Splitting (coarsening)
	// Способ построения C-F разбиения : 0 - standart, 1 - RS2, 3 - ST classical standart, 4 - RS2 ST.
	// RS2 улучшенная версия построения C-F разбиения содержащая второй проход.
	// ST - на основе STRONG Transpose.
	integer icoarseningTemp = 3, icoarseningSpeed = 3, icoarseningPressure = 3, icoarseningStress = 3;
	integer icoarseningtype = 3;
	// Stabilization BiCGStab.
	// 8.01.2017
	// предобусловленный алгебраичесеким многосеточным методом.
	// 0 - используется просто алгебраический многосеточный метод без какого-либо привлечения алгоритмов подпространства Крылова,
	// 1 - Используется алгоритм Х. Ван дер Ворста BiCGStab [1992], предобусловленный алгебраическим многосеточным методом.
	// 2 - Используется алгоритм Саада и Шульца FGMRes [1986], предобусловленный алгебраическим многосеточным методом.
	// 3 - Нелинейный многосеточный метод (обновление правой части на каждом V цикле). Для нелинейных граничных условий.
	integer istabilizationTemp = 0, istabilizationSpeed = 0, istabilizationPressure = 0, istabilizationStress = 0; // 0 - none
	integer istabilization = 0; // 0 - none
	// ipatch - номер патча.
	integer ipatch_number = 7;

	integer iprint_log = 1, iprint_log_Temperature = 1, iprint_log_Speed = 1, iprint_log_Pressure = 1, iprint_log_Stress = 1;

	// truncation for interpolation.
	integer itruncation_interpolation = 0, itruncation_interpolation_Temperature = 0, itruncation_interpolation_Speed = 0, itruncation_interpolation_Pressure = 0, itruncation_interpolation_Stress = 0;
	doublereal truncation_interpolation = 0.2, truncation_interpolation_Temperature = 0.2, truncation_interpolation_Speed = 0.2, truncation_interpolation_Pressure = 0.2, truncation_interpolation_Stress = 0.2;

	// gmres smoother
	// Ю.Саад, Мартин Г. Шульц [1986].
	bool b_gmresTemp = false, b_gmresSpeed = false, b_gmresPressure = false, b_gmresStress = false;
	bool b_gmres = false;

	// fgmres(m_restart)
	integer m_restart = 20;

	// lfil for BiCGStab+ILU2 and fgmres.
	integer lfil = 2;

	// AMGCL parameters
	// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped_jacobi;
	// 4 - spai1; 5 - chebyshev; 6 - iluk, k=1; 7 - iluk, k=2;
	integer amgcl_smoother = 0; 
	integer amgcl_selector = 0; // 0 - Ruge-Stueben (amg1r5 analog); 1 - smoother aggregation.
	integer amgcl_iterator = 0; // 0 - BiCGStab; 1 - FGMRes;

} MY_AMG_MANAGER;

MY_AMG_MANAGER my_amg_manager;




void my_amg_manager_init() {
	
	// Алгоритм сортировки используемый в 
	// алгебраическом многосеточном методе РУМБА.
	// 0 - COUNTING SORT
	// 1 - QUICK SORT
	// 2 - HEAP SORT
	my_amg_manager.imySortAlgorithm = 0; // default value COUNTING SORT

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

	//my_amg_manager.baglomeration_with_consistency_scaling = 0;
	my_amg_manager.bdiagonal_dominant = 1;

	my_amg_manager.iCFalgorithm_and_data_structure = 3; // 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Temperature = 3;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Speed = 3;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Pressure = 3;// 3-Treap.
	my_amg_manager.iCFalgorithm_and_data_structure_Stress = 3;// 3-Treap.	

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

	// ВНИМАНИЕ !!! устаревшие параметры, они больше не используются. 24.09.2019
	// 0 - не используется, 1 используется ILU разложение в качестве предобуславливателя.
	my_amg_manager.iFinnest_ilu = 0; 
	// Использование iluk разложения на глубоких уровнях вложенности для которых
	// сеточный шаблон nnz/n имеет размер меньше либо равный 6 (шести).
	my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6 = false;
	// ВНИМАНИЕ !!! устаревшие параметры, они больше не используются. 24.09.2019


	my_amg_manager.theta = 0.24;
	my_amg_manager.magic = 0.4;
	my_amg_manager.F_to_F_Temperature = 0.4;
	my_amg_manager.F_to_F_Speed = 0.4;
	my_amg_manager.F_to_F_Pressure = 0.4;
	my_amg_manager.F_to_F_Stress = 0.4;
	my_amg_manager.ilu2_smoother = 0; // 0 - не использовать, 1 - использовать.

	// Устаревшие переменные, более не используются и подлежат удалению. 11.05.2019
	//my_amg_manager.itypemodifyinterpol=0; // номер модификации интерполляции.
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
	// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
	// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
	my_amg_manager.icoarseningTemp = 0; // standart
	my_amg_manager.icoarseningSpeed = 0; // standart
	my_amg_manager.icoarseningPressure = 0; // standart
	my_amg_manager.icoarseningStress = 0; // standart
	my_amg_manager.icoarseningtype = 0; // standart vs RS 2.
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

	//fgmres(m_restart)
	my_amg_manager.m_restart = 20; // Количество итераций алгоритма fgmres перед перезапуском.

	// amg default settings:
	my_amg_manager.lfil = 2; // default value


	// AMGCL parameters
	// 0 - spai0; 1 - ilu0; 2- gauss-seidel; 3 - damped-jacobi;
	// 4 - spai1; 5 - chebyshev; 6 - iluk, k=1; 7 - iluk, k=2;
	my_amg_manager.amgcl_smoother = 0; // 0 - spai0; default
	my_amg_manager.amgcl_selector = 1; // 0 - Ruge-Stueben (amg1r5 analog); 1 - smoother aggregation.
	my_amg_manager.amgcl_iterator = 0; // 0 - BiCGStab; 1 - FGMRes.

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
	case 5:// П.Роуч SOR
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
	printf("version v0.48 2009 - 2019\n");
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


integer ihash_key_shorter(doublereal g, doublereal min_g, doublereal max_g )
{
	return ((integer)(isize_shorter_hash*((g - min_g) / (1.05*(max_g - min_g)))));

} // ihash_key_73

//здесь будет управление допуском на местах.
doublereal shorter_length_for_simplificationX(doublereal g, BLOCK* &b, integer &lb, 
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {
	if (!bFULL_AUTOMATIC) {
		return shorter_length_for_simplificationX_BASIC;
	}
	else {

		doublereal db = 1.0e30;

		
			//integer key = (integer)(isize_shorter_hash*((g - b[0].g.xS) / (1.05*(b[0].g.xE - b[0].g.xS))));
			integer key = ihash_key_shorter(g, b[0].g.xS, b[0].g.xE);

			if (!bshorter_hash_X[key]) {

				for (integer i = 0; i < lb; i++) {
					if ((b[i].g.itypegeom == PRISM)||(b[i].g.itypegeom==POLYGON)) {
						doublereal dmult = dbody_R_out_multiplyer;
						if (((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.yE - b[i].g.yS)) && ((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.zE - b[i].g.zS))) {
							// Ситуация слой термопасты или клея.
							dmult = dbody_dist_multiplyer;
						}
						if (((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.yE - b[i].g.yS)) || ((5.0*(b[i].g.xE - b[i].g.xS)) < (b[i].g.zE - b[i].g.zS))) {
							// Ситуация палец затвора.
							dmult = dbody_LGATE_situation_multiplyer;
						}
						if ((g >= b[i].g.xS - dmult *(b[i].g.xE - b[i].g.xS)) && (g <= b[i].g.xE) + dmult *(b[i].g.xE - b[i].g.xS)) {
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.xE - b[i].g.xS)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.xE - b[i].g.xS));
						}
					}
					if (b[i].g.itypegeom == CYLINDER) {
						switch (b[i].g.iPlane) {
						case XY: case XZ:
							if (b[i].g.R_in_cyl < 1.0e-40) {
								if ((g >= b[i].g.xC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.xC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
									// g в зоне влияния блока сномером i.
									if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl));
								}
							}
							else {
								if ((g >= b[i].g.xC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.xC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
									// g в зоне влияния блока сномером i.
									if (((1.0 / rdivision_interval)*(b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_in_cyl));
									// g в зоне влияния блока сномером i.
									if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl));
								}
							}
							break;
						case YZ:
							if ((g >= b[i].g.xC - dbody_Hcyl_multiplyer *b[i].g.Hcyl) && (g <= b[i].g.xC + b[i].g.Hcyl + dbody_Hcyl_multiplyer *b[i].g.Hcyl)) {
								// g в зоне влияния блока сномером i.
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
				for (integer i = 0; i < lw; i++) {
					doublereal dmult = dbody_R_out_multiplyer;
					switch (w[i].g.iPlane) {						
					case XY: case XZ:
						if ((g >= w[i].g.xS - dmult * (w[i].g.xE - w[i].g.xS)) && (g <= w[i].g.xE) + dmult * (w[i].g.xE - w[i].g.xS)) {
							// g в зоне влияния стенки с номером i.
							if (((1.0 / rdivision_interval)*(w[i].g.xE - w[i].g.xS)) < db) db = ((1.0 / rdivision_interval)*(w[i].g.xE - w[i].g.xS));
						}
						break;
					}
				}

				for (integer i = 0; i < ls; i++) {
					doublereal dmult = dbody_R_out_multiplyer;
					switch (s[i].g.iPlane) {						
					case XY: case XZ:
						if ((g >= s[i].g.xS - dmult * (s[i].g.xE - s[i].g.xS)) && (g <= s[i].g.xE) + dmult * (s[i].g.xE - s[i].g.xS)) {
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

// Здесь будет управление допуском на местах.
doublereal shorter_length_for_simplificationY(doublereal g, BLOCK* &b, integer &lb,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {
	if (!bFULL_AUTOMATIC) {
		return shorter_length_for_simplificationY_BASIC;
	}
	else {
		doublereal db = 1.0e30;

		//integer key = (integer)(isize_shorter_hash*((g - b[0].g.yS) / (1.05*(b[0].g.yE - b[0].g.yS))));
		integer key = ihash_key_shorter(g, b[0].g.yS, b[0].g.yE);
		if (!bshorter_hash_Y[key]) {



		for (integer i = 0; i < lb; i++) {
			if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == POLYGON)) {
				doublereal dmult = dbody_R_out_multiplyer;
				if (((5.0*(b[i].g.yE - b[i].g.yS)) < (b[i].g.xE - b[i].g.xS)) && ((5.0*(b[i].g.yE - b[i].g.yS)) < (b[i].g.zE - b[i].g.zS))) {
					// Ситуация слой термопасты или клея.
					dmult = dbody_dist_multiplyer;
				}
				if (((5.0*(b[i].g.yE - b[i].g.yS)) < (b[i].g.xE - b[i].g.xS)) || ((5.0*(b[i].g.yE - b[i].g.yS)) < (b[i].g.zE - b[i].g.zS))) {
					// Ситуация палец затвора.
					dmult = dbody_LGATE_situation_multiplyer;
				}

				if ((g >= b[i].g.yS - dmult *(b[i].g.yE - b[i].g.yS)) && (g <= b[i].g.yE + dmult *(b[i].g.yE - b[i].g.yS))) {
					// g в зоне влияния блока сномером i.
					if (((1.0 / rdivision_interval)*(b[i].g.yE - b[i].g.yS)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.yE - b[i].g.yS));
				}
			}
			if (b[i].g.itypegeom == CYLINDER) {
				switch (b[i].g.iPlane) {
				case XY: case YZ:
					if (b[i].g.R_in_cyl < 1.0e-40) {
						if ((g >= b[i].g.yC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.yC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl));
						}
					}
					else {
						if ((g >= b[i].g.yC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.yC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_in_cyl));
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl));
						}
					}
					break;
				case XZ:
					if ((g >= b[i].g.yC- dbody_Hcyl_multiplyer *b[i].g.Hcyl) && (g <= b[i].g.yC + b[i].g.Hcyl+ dbody_Hcyl_multiplyer *b[i].g.Hcyl)) {
						// g в зоне влияния блока сномером i.
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
			if (db > dcabinet_delta_multiplyer *(fabs(b[0].g.yE - b[0].g.yS))) db = dcabinet_delta_multiplyer *(fabs(b[0].g.yE - b[0].g.yS));
		}

		for (integer i = 0; i < lw; i++) {
			doublereal dmult = dbody_R_out_multiplyer;
			switch (w[i].g.iPlane) {				
			case XY: case YZ:
				if ((g >= w[i].g.yS - dmult * (w[i].g.yE - w[i].g.yS)) && (g <= w[i].g.yE) + dmult * (w[i].g.yE - w[i].g.yS)) {
					// g в зоне влияния стенки с номером i.
					if (((1.0 / rdivision_interval)*(w[i].g.yE - w[i].g.yS)) < db) db = ((1.0 / rdivision_interval)*(w[i].g.yE - w[i].g.yS));
				}
				break;
			}
		}

		for (integer i = 0; i < ls; i++) {
			doublereal dmult = dbody_R_out_multiplyer;
			switch (s[i].g.iPlane) {				
			case XY: case YZ:
				if ((g >= s[i].g.yS - dmult * (s[i].g.yE - s[i].g.yS)) && (g <= s[i].g.yE) + dmult * (s[i].g.yE - s[i].g.yS)) {
					// g в зоне влияния источника с номером i.
					if (((1.0 / rdivision_interval)*(s[i].g.yE - s[i].g.yS)) < db) db = ((1.0 / rdivision_interval)*(s[i].g.yE - s[i].g.yS));
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

// Здесь будет управление допуском на местах.
doublereal shorter_length_for_simplificationZ(doublereal g, BLOCK* &b, integer &lb,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {

	if (!bFULL_AUTOMATIC) {
		return shorter_length_for_simplificationZ_BASIC;
	}
	else {
		doublereal db = 1.0e30;

		//integer key = (integer)(isize_shorter_hash*((g - b[0].g.zS) / (1.05*(b[0].g.zE - b[0].g.zS))));
		integer key = ihash_key_shorter(g, b[0].g.zS, b[0].g.zE);
		if (!bshorter_hash_Z[key]) {



		for (integer i = 0; i < lb; i++) {
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


				if ((g >= b[i].g.zS - dmult*(b[i].g.zE - b[i].g.zS)) && (g <= b[i].g.zE+ dmult*(b[i].g.zE - b[i].g.zS))) {
					// g в зоне влияния блока сномером i.
					if (((1.0 / rdivision_interval)*(b[i].g.zE - b[i].g.zS)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.zE - b[i].g.zS));
				}
			}
			if (b[i].g.itypegeom == CYLINDER) {
				switch (b[i].g.iPlane) {
				case XZ: case YZ:
					if (b[i].g.R_in_cyl < 1.0e-40) {
						if ((g >= b[i].g.zC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.zC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl));
						}
					}
					else {
						if ((g >= b[i].g.zC - dbody_R_out_multiplyer *b[i].g.R_out_cyl) && (g <= b[i].g.zC + dbody_R_out_multiplyer *b[i].g.R_out_cyl)) {
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_in_cyl));
							// g в зоне влияния блока сномером i.
							if (((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl)) < db) db = ((1.0 / rdivision_interval)*(b[i].g.R_out_cyl - b[i].g.R_in_cyl));
						}
					}
					break;
				case XY:
					if ((g >= b[i].g.zC- dbody_Hcyl_multiplyer *b[i].g.Hcyl) && (g <= b[i].g.zC + b[i].g.Hcyl+ dbody_Hcyl_multiplyer *b[i].g.Hcyl)) {
						// g в зоне влияния блока сномером i.
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
			if (db > dcabinet_delta_multiplyer *(fabs(b[0].g.zE - b[0].g.zS))) db = dcabinet_delta_multiplyer *(fabs(b[0].g.zE - b[0].g.zS));
		}

		for (integer i = 0; i < lw; i++) {
			doublereal dmult = dbody_R_out_multiplyer;
			switch (w[i].g.iPlane) {				
			case YZ: case XZ:
				if ((g >= w[i].g.zS - dmult * (w[i].g.zE - w[i].g.zS)) && (g <= w[i].g.zE) + dmult * (w[i].g.zE - w[i].g.zS)) {
					// g в зоне влияния стенки с номером i.
					if (((1.0 / rdivision_interval)*(w[i].g.zE - w[i].g.zS)) < db) db = ((1.0 / rdivision_interval)*(w[i].g.zE - w[i].g.zS));
				}
				break;
			}
		}

		for (integer i = 0; i < ls; i++) {
			doublereal dmult = dbody_R_out_multiplyer;
			switch (s[i].g.iPlane) {				
			case YZ: case XZ:
				if ((g >= s[i].g.zS - dmult * (s[i].g.zE - s[i].g.zS)) && (g <= s[i].g.zE) + dmult * (s[i].g.zE - s[i].g.zS)) {
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
			
			if (sqrt((x4 - x1)*(x4 - x1) + (y4 - y1)*(y4 - y1)) > 1.0e-40) {
				// Четвертая точка не должна совпадать с первой. Она может совпадать если у нас всего пять точек и
				// первый и последний отрезок лежат на одной прямой и имеют общую точку. Нам нужно корректно 
				// обработать такой случай.

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
// Переделал с printf("..."); на std::cout<<"..."<<std::endl; 14.08.2019.
void BODY_CHECK(BLOCK* &b, integer lb, WALL* &w, integer lw, SOURCE* &s, integer ls) {
	bool bOk = true;
	const doublereal dcheck_eps = 1e-10; // 1.0e-10 подобрано для tgf*

	for (integer i = lb - 1; i >= 0; i--) {
		// Проверка типа геометрии.
		CHECK_TYPE_GEOM(b[i].g.itypegeom);
		if (b[i].g.itypegeom == PRISM) {
			// Ошибка не число.
			if (b[i].g.xS != b[i].g.xS) {
				//printf("Error body[%lld] xS NOT VALUE=%e.\n", i, b[i].g.xS);
				std::cout<<"Error body["<<i<<"] xS NOT VALUE="<< b[i].g.xS <<"."<<std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.xE != b[i].g.xE) {
				//printf("Error body[%lld] xE NOT VALUE=%e.\n", i, b[i].g.xE);
				std::cout << "Error body[" << i << "] xE NOT VALUE=" << b[i].g.xE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.yS != b[i].g.yS) {
				//printf("Error body[%lld] yS NOT VALUE=%e.\n", i, b[i].g.yS);
				std::cout << "Error body[" << i << "] yS NOT VALUE=" << b[i].g.yS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.yE != b[i].g.yE) {
				//printf("Error body[%lld] yE NOT VALUE=%e.\n", i, b[i].g.yE);
				std::cout << "Error body[" << i << "] yE NOT VALUE=" << b[i].g.yE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.zS != b[i].g.zS) {
				//printf("Error body[%lld] zS NOT VALUE=%e.\n", i, b[i].g.zS);
				std::cout << "Error body[" << i << "] zS NOT VALUE=" << b[i].g.zS << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.zE != b[i].g.zE) {
				//printf("Error body[%lld] zE NOT VALUE=%e.\n", i, b[i].g.zE);
				std::cout << "Error body[" << i << "] zE NOT VALUE=" << b[i].g.zE << "." << std::endl;
				system("pause");
				bOk = false;
			}
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (b[i].g.xS >= b[i].g.xE) {
				//printf("body[%lld].xS=%e >= body[%lld].xE=%e\n", i, b[i].g.xS, i, b[i].g.xE);
				std::cout << "body["<<i<<"].xS="<< b[i].g.xS<<" >= body["<<i<<"].xE="<< b[i].g.xE<< std::endl;
				system("pause");
				doublereal temp = b[i].g.xE;
				b[i].g.xE = b[i].g.xS;
				b[i].g.xS = temp;
				bOk = false;
			}
			if (b[i].g.yS >= b[i].g.yE) {
				//printf("body[%lld].yS=%e >= body[%lld].yE=%e\n", i, b[i].g.yS, i, b[i].g.yE);
				std::cout << "body[" << i << "].yS=" << b[i].g.yS << " >= body[" << i << "].yE=" << b[i].g.yE << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				doublereal temp = b[i].g.yE;
				b[i].g.yE = b[i].g.yS;
				b[i].g.yS = temp;
				bOk = false;
			}
			if (b[i].g.zS >= b[i].g.zE) {
				//printf("body[%lld].zS=%e >= body[%lld].zE=%e\n", i, b[i].g.zS, i, b[i].g.zE);
				std::cout << "body[" << i << "].zS=" << b[i].g.zS << " >= body[" << i << "].zE=" << b[i].g.zE << std::endl;
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
					std::cout << "b["<<i<<"].g.xE="<< b[i].g.xE<<" > cabinet.xE="<< b[0].g.xE<<"."<< std::endl;
					system("pause");
					b[i].g.xE = b[0].g.xE; // исправление.
					bOk = false;
				}
				if (b[i].g.yE > b[0].g.yE) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.yE=%e > cabinet.yE=%e.\n", i, b[i].g.yE, b[0].g.yE);
					std::cout << "b[" << i << "].g.yE=" << b[i].g.yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
					system("pause");
					b[i].g.yE = b[0].g.yE; // исправление.
					bOk = false;
				}
				if (b[i].g.zE > b[0].g.zE) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.zE=%e > cabinet.zE=%e.\n", i, b[i].g.zE, b[0].g.zE);
					std::cout << "b[" << i << "].g.zE=" << b[i].g.zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
					system("pause");
					b[i].g.zE = b[0].g.zE; // исправление.
					bOk = false;
				}
				if (b[i].g.xS < b[0].g.xS) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.xS=%e < cabinet.xS=%e.\n", i, b[i].g.xS, b[0].g.xS);
					std::cout << "b[" << i << "].g.xS=" << b[i].g.xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
					system("pause");
					b[i].g.xS = b[0].g.xS; // исправление.
					bOk = false;
				}
				if (b[i].g.yS < b[0].g.yS) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.yS=%e < cabinet.yS=%e.\n", i, b[i].g.yS, b[0].g.yS);
					std::cout << "b[" << i << "].g.yS=" << b[i].g.yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
					std::cout << "" << std::endl;
					system("pause");
					b[i].g.yS = b[0].g.yS; // исправление.
					bOk = false;
				}
				if (b[i].g.zS < b[0].g.zS) {
					printf("ERROR. Your model is incorrect.\n");
					//printf("b[%lld].g.zS=%e < cabinet.zS=%e.\n", i, b[i].g.zS, b[0].g.zS);
					std::cout << "b[" << i << "].g.zS=" << b[i].g.zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
					system("pause");
					b[i].g.zS = b[0].g.zS; // исправление.
					bOk = false;
				}
			}
		}
		if (b[i].g.itypegeom == CYLINDER) {
			if (b[i].g.xC != b[i].g.xC) {
				//printf("Error body[%lld] xC NOT VALUE=%e.\n",i, b[i].g.xC);
				std::cout << "Error body[" << i << "] xC NOT VALUE=" << b[i].g.xC << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.yC != b[i].g.yC) {
				//printf("Error body[%lld] yC NOT VALUE=%e.\n", i, b[i].g.yC);
				std::cout << "Error body[" << i << "] yC NOT VALUE=" << b[i].g.yC << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.zC != b[i].g.zC) {
				//printf("Error body[%lld] zC NOT VALUE=%e.\n", i, b[i].g.zC);
				std::cout << "Error body[" << i << "] zC NOT VALUE=" << b[i].g.zC << "." << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.Hcyl != b[i].g.Hcyl) {
				//printf("Error body[%lld] Hcyl NOT VALUE=%e.\n", i, b[i].g.Hcyl);
				std::cout << "Error body[" << i << "] Hcyl NOT VALUE=" << b[i].g.Hcyl << "." << std::endl;
				std::cout << "" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_out_cyl != b[i].g.R_out_cyl) {
				//printf("Error body[%lld] R_out_cyl NOT VALUE=%e.\n", i, b[i].g.R_out_cyl);
				std::cout << "Error body[" << i << "] R_out_cyl NOT VALUE=" << b[i].g.R_out_cyl << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl != b[i].g.R_in_cyl) {
				//printf("Error body[%lld] R_in_cyl NOT VALUE=%e.\n", i, b[i].g.R_in_cyl);
				std::cout << "Error body[" << i << "] R_in_cyl NOT VALUE=" << b[i].g.R_in_cyl << "." << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.Hcyl <= 0.0) {
				//printf("Error: negative or zero Hcyl=%e of cylinder body[%lld]\n", b[i].g.Hcyl, i);
				std::cout << "Error: negative or zero Hcyl="<< b[i].g.Hcyl<<" of cylinder body["<<i<<"]" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_out_cyl <= 0.0) {
				//printf("Error: negative or zero R_out_cyl=%e of cylinder body[%lld]\n", b[i].g.R_out_cyl, i);
				std::cout << "Error: negative or zero R_out_cyl=" << b[i].g.R_out_cyl << " of cylinder body[" << i << "]" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl < 0.0) {
				printf("Error: negative R_in_cyl=%e of cylinder body[%lld]\n", b[i].g.R_in_cyl, i);
				std::cout << "Error: negative  R_in_cyl=" << b[i].g.R_in_cyl << " of cylinder body[" << i << "]" << std::endl;
				system("pause");
				bOk = false;
			}
			if (b[i].g.R_in_cyl > b[i].g.R_out_cyl) {
				//printf("Error: b[%lld].g.R_in_cyl=%e > b[%lld].g.R_out_cyl=%e\n", i, b[i].g.R_in_cyl, i, b[i].g.R_out_cyl);
				std::cout << "Error: b["<<i<<"].g.R_in_cyl="<< b[i].g.R_in_cyl<<" > b["<<i<<"].g.R_out_cyl="<< b[i].g.R_out_cyl << std::endl;
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
				
				if (xE_box_cyl > b[0].g.xE + dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] xE_box_cyl=%e > cabinet.xE=%e.\n", i, xE_box_cyl, b[0].g.xE);
					std::cout << "body["<<i<<"] xE_box_cyl="<< xE_box_cyl<<" > cabinet.xE="<< b[0].g.xE<<"." << std::endl;
					printf("tolerance=%e\n", xE_box_cyl - b[0].g.xE);
					system("pause");

					bOk = false;
				}
				if (yE_box_cyl > b[0].g.yE + dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] yE_box_cyl=%e > cabinet.yE=%e.\n", i, yE_box_cyl, b[0].g.yE);
					std::cout << "body[" << i << "] yE_box_cyl=" << yE_box_cyl << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
					printf("tolerance=%e\n", yE_box_cyl - b[0].g.yE);
					system("pause");

					bOk = false;
				}
				if (zE_box_cyl > b[0].g.zE+ dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] zE_box_cyl=%e > cabinet.zE=%e.\n", i, zE_box_cyl, b[0].g.zE);
					std::cout << "body[" << i << "] zE_box_cyl=" << zE_box_cyl << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
					printf("tolerance=%e\n", zE_box_cyl - b[0].g.zE);
					system("pause");

					bOk = false;
				}
				if (xS_box_cyl < b[0].g.xS- dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] xS_box_cyl=%e < cabinet.xS=%e.\n", i, xS_box_cyl, b[0].g.xS);
					std::cout << "body[" << i << "] xS_box_cyl=" << xS_box_cyl << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
					printf("tolerance=%e\n", b[0].g.xS - xS_box_cyl);
					system("pause");

					bOk = false;
				}
				if (yS_box_cyl < b[0].g.yS - dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] yS_box_cyl=%e < cabinet.yS=%e.\n", i, yS_box_cyl, b[0].g.yS);
					std::cout << "body[" << i << "] yS_box_cyl=" << yS_box_cyl << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
					printf("tolerance=%e\n", b[0].g.yS - yS_box_cyl);
					system("pause");

					bOk = false;
				}
				if (zS_box_cyl < b[0].g.zS - dcheck_eps) {
					printf("ERROR CYLINDER. Your model is incorrect.\n");
					//printf("body[%lld] zS_box_cyl=%e < cabinet.zS=%e.\n", i, zS_box_cyl, b[0].g.zS);
					std::cout << "body[" << i << "] zS_box_cyl=" << zS_box_cyl << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
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
					std::cout << "Error number points on polygon very small: b["<<i<<"].g.nsizei="<< b[i].g.nsizei<<"<2." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.xi == NULL) {
					//printf("Error: no memory allocate for b[%lld].g.xi array.\n",i);
					std::cout << "Error: no memory allocate for b["<<i<<"].g.xi array." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.yi == NULL) {
					//printf("Error: no memory allocate for b[%lld].g.yi array.\n", i);
					std::cout << "Error: no memory allocate for b["<<i<<"].g.yi array." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.zi == NULL) {
					//printf("Error: no memory allocate for b[%lld].g.zi array.\n", i);
					std::cout << "Error: no memory allocate for b["<<i<<"].g.zi array." << std::endl;
					system("pause");
					bOk = false;
				}

				if (b[i].g.hi == NULL) {
					//printf("Error: no memory allocate for b[%lld].g.hi array.\n", i);
					std::cout << "Error: no memory allocate for b["<<i<<"].g.hi array." << std::endl;
					system("pause");
					bOk = false;
				}

				for (integer ip = 0; ip < b[i].g.nsizei; ip++) {
					if (b[i].g.xi[ip] != b[i].g.xi[ip]) {
						//printf("Error body[%lld] xi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.xi[ip]);
						std::cout << "Error body["<<i<<"] xi["<<ip<<"] NOT VALUE="<< b[i].g.xi[ip]<<"." << std::endl;
						system("pause");
						bOk = false;
					}
					if (b[i].g.yi[ip] != b[i].g.yi[ip]) {
						//printf("Error body[%lld] yi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.yi[ip]);
						std::cout << "Error body[" << i << "] yi[" << ip << "] NOT VALUE=" << b[i].g.yi[ip] << "." << std::endl;
						system("pause");
						bOk = false;
					}
					if (b[i].g.zi[ip] != b[i].g.zi[ip]) {
						//printf("Error body[%lld] zi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.zi[ip]);
						std::cout << "Error body[" << i << "] zi[" << ip << "] NOT VALUE=" << b[i].g.zi[ip] << "." << std::endl;
						system("pause");
						bOk = false;
					}
					if (b[i].g.hi[ip] != b[i].g.hi[ip]) {
						//printf("Error body[%lld] hi[%lld] NOT VALUE=%e.\n", i, ip, b[i].g.hi[ip]);
						std::cout << "Error body[" << i << "] hi[" << ip << "] NOT VALUE=" << b[i].g.hi[ip] << "." << std::endl;
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
							//printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							std::cout << "Error body["<<i<<"].g.hi["<<ip<<"]="<< b[i].g.hi[ip]<<" is NEGATIV or ZERO." << std::endl;
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
							//printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							std::cout << "Error body[" << i << "].g.hi[" << ip << "]=" << b[i].g.hi[ip] << " is NEGATIV or ZERO." << std::endl;
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
							//printf("Error body[%lld].g.hi[%lld]=%e is NEGATIV or ZERO.\n", i, ip, b[i].g.hi[ip]);
							std::cout << "Error body[" << i << "].g.hi[" << ip << "]=" << b[i].g.hi[ip] << " is NEGATIV or ZERO." << std::endl;
							system("pause");
							bOk = false;
						}
						xS_box_polygon = b[i].g.xi[ip];
						xE_box_polygon = b[i].g.xi[ip] + b[i].g.hi[ip];
					}
					break;
				}

				if (xE_box_polygon > b[0].g.xE+ dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] xE_box_polygon=%e > cabinet.xE=%e.\n", i, xE_box_polygon, b[0].g.xE);
					std::cout << "body["<<i<<"] xE_box_polygon="<< xE_box_polygon<<" > cabinet.xE="<< b[0].g.xE<<"." << std::endl;
					system("pause");

					bOk = false;
				}
				if (yE_box_polygon > b[0].g.yE + dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] yE_box_polygon=%e > cabinet.yE=%e.\n", i, yE_box_polygon, b[0].g.yE);
					std::cout << "body[" << i << "] yE_box_polygon=" << yE_box_polygon << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (zE_box_polygon > b[0].g.zE + dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] zE_box_polygon=%e > cabinet.zE=%e.\n", i, zE_box_polygon, b[0].g.zE);
					std::cout << "body[" << i << "] zE_box_polygon=" << zE_box_polygon << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (xS_box_polygon < b[0].g.xS - dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] xS_box_polygon=%e < cabinet.xS=%e.\n", i, xS_box_polygon, b[0].g.xS);
					std::cout << "body[" << i << "] xS_box_polygon=" << xS_box_polygon << " > cabinet.xS=" << b[0].g.xS << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (yS_box_polygon < b[0].g.yS - dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] yS_box_polygon=%e < cabinet.yS=%e.\n", i, yS_box_polygon, b[0].g.yS);
					std::cout << "body[" << i << "] yS_box_polygon=" << yS_box_polygon << " > cabinet.yS=" << b[0].g.yS << "." << std::endl;
					system("pause");

					bOk = false;
				}
				if (zS_box_polygon < b[0].g.zS - dcheck_eps) {
					printf("ERROR POLYGON. Your model is incorrect.\n");
					//printf("body[%lld] zS_box_polygon=%e < cabinet.zS=%e.\n", i, zS_box_polygon, b[0].g.zS);
					std::cout << "body[" << i << "] zS_box_polygon=" << zS_box_polygon << " > cabinet.zS=" << b[0].g.zS << "." << std::endl;
					system("pause");

					bOk = false;
				}

			}

		}
    }

	for (integer i = 0; i < lw; i++) {

		// Ошибка не число.
		if (w[i].g.xS != w[i].g.xS) {
			//printf("Error wall[%lld] xS NOT VALUE=%e.\n", i, w[i].g.xS);
			std::cout << "Error wall["<<i<<"] xS NOT VALUE="<< w[i].g.xS<<"." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.xE != w[i].g.xE) {
			//printf("Error wall[%lld] xE NOT VALUE=%e.\n", i, w[i].g.xE);
			std::cout << "Error wall[" << i << "] xE NOT VALUE=" << w[i].g.xE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.yS != w[i].g.yS) {
			//printf("Error wall[%lld] yS NOT VALUE=%e.\n", i, w[i].g.yS);
			std::cout << "Error wall[" << i << "] yS NOT VALUE=" << w[i].g.yS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.yE != w[i].g.yE) {
			//printf("Error wall[%lld] yE NOT VALUE=%e.\n", i, w[i].g.yE);
			std::cout << "Error wall[" << i << "] yE NOT VALUE=" << w[i].g.yE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.zS != w[i].g.zS) {
			//printf("Error wall[%lld] zS NOT VALUE=%e.\n", i, w[i].g.zS);
			std::cout << "Error wall[" << i << "] zS NOT VALUE=" << w[i].g.zS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (w[i].g.zE != w[i].g.zE) {
			//printf("Error wall[%lld] zE NOT VALUE=%e.\n", i, w[i].g.zE);
			std::cout << "Error wall[" << i << "] zE NOT VALUE=" << w[i].g.zE << "." << std::endl;
			system("pause");
			bOk = false;
		}

		switch (w[i].iPlane) {
		case XY :
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.xS >= w[i].g.xE) {
				//printf("wall[%lld].xS=%e >= wall[%lld].xE=%e\n", i, w[i].g.xS, i, w[i].g.xE);
				std::cout << "wall["<<i<<"].xS="<< w[i].g.xS<<" >= wall["<<i<<"].xE="<< w[i].g.xE << std::endl;
				system("pause");
				doublereal temp = w[i].g.xE;
				w[i].g.xE = w[i].g.xS;
				w[i].g.xS = temp;
				bOk = false;
			}
			if (w[i].g.yS >= w[i].g.yE) {
				//printf("wall[%lld].yS=%e >= wall[%lld].yE=%e\n", i, w[i].g.yS, i, w[i].g.yE);
				std::cout << "wall[" << i << "].yS=" << w[i].g.yS << " >= wall[" << i << "].yE=" << w[i].g.yE << std::endl;
				system("pause");
				doublereal temp = w[i].g.yE;
				w[i].g.yE = w[i].g.yS;
				w[i].g.yS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.zS - w[i].g.zE) > 1.0e-40) {
				//printf("non union iso position on plane wall[%lld]: zS=%e zE=%e\n",i, w[i].g.zS, w[i].g.zE);
				std::cout << "non union iso position on plane wall[" << i << "]: zS=" << w[i].g.zS << " zE=" << w[i].g.zE << std::endl;
				system("pause");
				w[i].g.zS = w[i].g.zE;
				bOk = false;
			}
			break;
		case XZ :
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.xS >= w[i].g.xE) {
				//printf("wall[%lld].xS=%e >= wall[%lld].xE=%e\n", i, w[i].g.xS, i, w[i].g.xE);
				std::cout << "wall["<<i<<"].xS="<< w[i].g.xS<<" >= wall["<<i<<"].xE="<< w[i].g.xE  << std::endl;
				system("pause");
				doublereal temp = w[i].g.xE;
				w[i].g.xE = w[i].g.xS;
				w[i].g.xS = temp;
				bOk = false;
			}
			if (w[i].g.zS >= w[i].g.zE) {
				//printf("wall[%lld].zS=%e >= wall[%lld].zE=%e\n", i, w[i].g.zS, i, w[i].g.zE);
				std::cout << "wall[" << i << "].zS=" << w[i].g.zS << " >= wall[" << i << "].zE=" << w[i].g.zE << std::endl;
				system("pause");
				doublereal temp = w[i].g.zE;
				w[i].g.zE = w[i].g.zS;
				w[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.yS - w[i].g.yE) > 1.0e-40) {
				//printf("non union iso position on plane wall[%lld]: yS=%e yE=%e\n", i, w[i].g.yS, w[i].g.yE);
				std::cout << "non union iso position on plane wall[" << i << "]: yS=" << w[i].g.yS << " yE=" << w[i].g.yE << std::endl;
				system("pause");
				w[i].g.yS=w[i].g.yE;
				bOk = false;
			}
			break;
		case YZ :
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (w[i].g.yS >= w[i].g.yE) {
				//printf("wall[%lld].yS=%e >= wall[%lld].yE=%e\n", i, w[i].g.yS, i, w[i].g.yE);
				std::cout << "wall[" << i << "].yS=" << w[i].g.yS << " >= wall[" << i << "].yE=" << w[i].g.yE << std::endl;
				system("pause");
				doublereal temp = w[i].g.yE;
				w[i].g.yE = w[i].g.yS;
				w[i].g.yS = temp;
				bOk = false;
			}
			if (w[i].g.zS >= w[i].g.zE) {
				//printf("wall[%lld].zS=%e >= wall[%lld].zE=%e\n", i, w[i].g.zS, i, w[i].g.zE);
				std::cout << "wall[" << i << "].zS=" << w[i].g.zS << " >= wall[" << i << "].zE=" << w[i].g.zE << std::endl;
				system("pause");
				doublereal temp = w[i].g.zE;
				w[i].g.zE = w[i].g.zS;
				w[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(w[i].g.xS - w[i].g.xE) > 1.0e-40) {
				//printf("non union iso position on plane wall[%lld]: xS=%e xE=%e\n", i, w[i].g.xS, w[i].g.xE);
				std::cout << "non union iso position on plane wall[" << i << "]: xS=" << w[i].g.xS << " xE=" << w[i].g.xE << std::endl;
				system("pause");
				w[i].g.xS = w[i].g.xE;
				bOk = false;
			}
			break;
		}
		if (w[i].g.xE > b[0].g.xE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.xE=%e > cabinet.xE=%e.\n", i, w[i].g.xE, b[0].g.xE);
			std::cout << "wall["<<i<<"].g.xE="<< w[i].g.xE<<" > cabinet.xE="<< b[0].g.xE<<"." << std::endl;
			system("pause");
			w[i].g.xE = b[0].g.xE; // исправление.
			bOk = false;
		}
		if (w[i].g.yE > b[0].g.yE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.yE=%e > cabinet.yE=%e.\n", i, w[i].g.yE, b[0].g.yE);
			std::cout << "wall[" << i << "].g.yE=" << w[i].g.yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
			system("pause");
			w[i].g.yE = b[0].g.yE; // исправление.
			bOk = false;
		}
		if (w[i].g.zE > b[0].g.zE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.zE=%e > cabinet.zE=%e.\n", i, w[i].g.zE, b[0].g.zE);
			std::cout << "wall[" << i << "].g.zE=" << w[i].g.zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
			system("pause");
			w[i].g.zE = b[0].g.zE; // исправление.
			bOk = false;
		}
		if (w[i].g.xS < b[0].g.xS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.xS=%e < cabinet.xS=%e.\n", i, w[i].g.xS, b[0].g.xS);
			std::cout << "wall[" << i << "].g.xS=" << w[i].g.xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
			system("pause");
			w[i].g.xS = b[0].g.xS; // исправление.
			bOk = false;
		}
		if (w[i].g.yS < b[0].g.yS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.yS=%e < cabinet.yS=%e.\n", i, w[i].g.yS, b[0].g.yS);
			std::cout << "wall[" << i << "].g.yS=" << w[i].g.yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
			system("pause");
			w[i].g.yS = b[0].g.yS; // исправление.
			bOk = false;
		}
		if (w[i].g.zS < b[0].g.zS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("wall[%lld].g.zS=%e < cabinet.zS=%e.\n", i, w[i].g.zS, b[0].g.zS);
			std::cout << "wall[" << i << "].g.zS=" << w[i].g.zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
			system("pause");
			w[i].g.zS = b[0].g.zS; // исправление.
			bOk = false;
		}
	}

	for (integer i = 0; i < ls; i++) {

		// Ошибка не число.
		if (s[i].g.xS != s[i].g.xS) {
			//printf("Error source[%lld] xS NOT VALUE=%e.\n", i, s[i].g.xS);
			std::cout << "Error source[" << i << "] xS NOT VALUE=" << s[i].g.xS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.xE != s[i].g.xE) {
			//printf("Error source[%lld] xE NOT VALUE=%e.\n", i, s[i].g.xE);
			std::cout << "Error source[" << i << "] xE NOT VALUE=" << s[i].g.xE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.yS != s[i].g.yS) {
			//printf("Error source[%lld] yS NOT VALUE=%e.\n", i, s[i].g.yS);
			std::cout << "Error source[" << i << "] yS NOT VALUE=" << s[i].g.yS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.yE != s[i].g.yE) {
			//printf("Error source[%lld] yE NOT VALUE=%e.\n", i, s[i].g.yE);
			std::cout << "Error source[" << i << "] yE NOT VALUE=" << s[i].g.yE << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.zS != s[i].g.zS) {
			//printf("Error source[%lld] zS NOT VALUE=%e.\n", i, s[i].g.zS);
			std::cout << "Error source[" << i << "] zS NOT VALUE=" << s[i].g.zS << "." << std::endl;
			system("pause");
			bOk = false;
		}
		if (s[i].g.zE != s[i].g.zE) {
			//printf("Error source[%lld] zE NOT VALUE=%e.\n", i, s[i].g.zE);
			std::cout << "Error source[" << i << "] zE NOT VALUE=" << s[i].g.zE << "." << std::endl;
			system("pause");
			bOk = false;
		}

		switch (s[i].iPlane) {
		case XY:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.xS >= s[i].g.xE) {
				//printf("source[%lld].xS=%e >= source[%lld].xE=%e\n", i, s[i].g.xS, i, s[i].g.xE);
				std::cout << "source[" << i << "].xS=" << s[i].g.xS << " >= source[" << i << "].xE=" << s[i].g.xE << std::endl;
				system("pause");
				doublereal temp = s[i].g.xE;
				s[i].g.xE = s[i].g.xS;
				s[i].g.xS = temp;
				bOk = false;
			}
			if (s[i].g.yS >= s[i].g.yE) {
				//printf("source[%lld].yS=%e >= source[%lld].yE=%e\n", i, s[i].g.yS, i, s[i].g.yE);
				std::cout << "source[" << i << "].yS=" << s[i].g.yS << " >= source[" << i << "].yE=" << s[i].g.yE << std::endl;
				system("pause");
				doublereal temp = s[i].g.yE;
				s[i].g.yE = s[i].g.yS;
				s[i].g.yS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.zS - s[i].g.zE) > 1.0e-40) {
				//printf("non union iso position on plane source[%lld]: zS=%e zE=%e\n", i, s[i].g.zS, s[i].g.zE);
				std::cout << "non union iso position on plane source[" << i << "]: zS=" << s[i].g.zS << " zE=" << s[i].g.zE << std::endl;
				system("pause");
				s[i].g.zS = s[i].g.zE;
				bOk = false;
			}
			break;
		case XZ:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.xS >= s[i].g.xE) {
				//printf("source[%lld].xS=%e >= source[%lld].xE=%e\n", i, s[i].g.xS, i, s[i].g.xE);
				std::cout << "source[" << i << "].xS=" << s[i].g.xS << " >= source[" << i << "].xE=" << s[i].g.xE << std::endl;
				system("pause");
				doublereal temp = s[i].g.xE;
				s[i].g.xE = s[i].g.xS;
				s[i].g.xS = temp;
				bOk = false;
			}
			if (s[i].g.zS >= s[i].g.zE) {
				//printf("source[%lld].zS=%e >= source[%lld].zE=%e\n", i, s[i].g.zS, i, s[i].g.zE);
				std::cout << "source[" << i << "].zS=" << s[i].g.zS << " >= source[" << i << "].zE=" << s[i].g.zE << std::endl;
				system("pause");
				doublereal temp = s[i].g.zE;
				s[i].g.zE = s[i].g.zS;
				s[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.yS - s[i].g.yE) > 1.0e-40) {
				//printf("non union iso position on plane source[%lld]: yS=%e yE=%e\n", i, s[i].g.yS, s[i].g.yE);
				std::cout << "non union iso position on plane source[" << i << "]: yS=" << s[i].g.yS << " yE=" << s[i].g.yE << std::endl;
				system("pause");
				s[i].g.yS = s[i].g.yE;
				bOk = false;
			}
			break;
		case YZ:
			// Ошибка порядка следования Start должен быть строго меньше End.
			if (s[i].g.yS >= s[i].g.yE) {
				//printf("source[%lld].yS=%e >= source[%lld].yE=%e\n", i, s[i].g.yS, i, s[i].g.yE);
				std::cout << "source[" << i << "].yS=" << s[i].g.yS << " >= source[" << i << "].yE=" << s[i].g.yE << std::endl;
				system("pause");
				doublereal temp = s[i].g.yE;
				s[i].g.yE = s[i].g.yS;
				s[i].g.yS = temp;
				bOk = false;
			}
			if (s[i].g.zS >= s[i].g.zE) {
				//printf("source[%lld].zS=%e >= source[%lld].zE=%e\n", i, s[i].g.zS, i, s[i].g.zE);
				std::cout << "source[" << i << "].zS=" << s[i].g.zS << " >= source[" << i << "].zE=" << s[i].g.zE << std::endl;
				system("pause");
				doublereal temp = s[i].g.zE;
				s[i].g.zE = s[i].g.zS;
				s[i].g.zS = temp;
				bOk = false;
			}
			if (fabs(s[i].g.xS - s[i].g.xE) > 1.0e-40) {
				//printf("non union iso position on plane source[%lld]: xS=%e xE=%e\n", i, s[i].g.xS, s[i].g.xE);
				std::cout << "non union iso position on plane source[" << i << "]: xS=" << s[i].g.xS << " xE=" << s[i].g.xE << std::endl;
				system("pause");
				s[i].g.xS = s[i].g.xE;
				bOk = false;
			}
			break;
		}
		if (s[i].g.xE > b[0].g.xE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.xE=%e > cabinet.xE=%e.\n", i, s[i].g.xE, b[0].g.xE);
			std::cout << "source[" << i << "].g.xE=" << s[i].g.xE << " > cabinet.xE=" << b[0].g.xE << "." << std::endl;
			system("pause");
			s[i].g.xE = b[0].g.xE; // исправление.
			bOk = false;
		}
		if (s[i].g.yE > b[0].g.yE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.yE=%e > cabinet.yE=%e.\n", i, s[i].g.yE, b[0].g.yE);
			std::cout << "source[" << i << "].g.yE=" << s[i].g.yE << " > cabinet.yE=" << b[0].g.yE << "." << std::endl;
			system("pause");
			s[i].g.yE = b[0].g.yE; // исправление.
			bOk = false;
		}
		if (s[i].g.zE > b[0].g.zE) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.zE=%e > cabinet.zE=%e.\n", i, s[i].g.zE, b[0].g.zE);
			std::cout << "source[" << i << "].g.zE=" << s[i].g.zE << " > cabinet.zE=" << b[0].g.zE << "." << std::endl;
			system("pause");
			s[i].g.zE = b[0].g.zE; // исправление.
			bOk = false;
		}
		if (s[i].g.xS < b[0].g.xS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.xS=%e < cabinet.xS=%e.\n", i, s[i].g.xS, b[0].g.xS);
			std::cout << "source[" << i << "].g.xS=" << s[i].g.xS << " < cabinet.xS=" << b[0].g.xS << "." << std::endl;
			system("pause");
			s[i].g.xS = b[0].g.xS; // исправление.
			bOk = false;
		}
		if (s[i].g.yS < b[0].g.yS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.yS=%e < cabinet.yS=%e.\n", i, s[i].g.yS, b[0].g.yS);
			std::cout << "source[" << i << "].g.yS=" << s[i].g.yS << " < cabinet.yS=" << b[0].g.yS << "." << std::endl;
			system("pause");
			s[i].g.yS = b[0].g.yS; // исправление.
			bOk = false;
		}
		if (s[i].g.zS < b[0].g.zS) {
			printf("ERROR. Your model is incorrect.\n");
			//printf("source[%lld].g.zS=%e < cabinet.zS=%e.\n", i, s[i].g.zS, b[0].g.zS);
			std::cout << "source[" << i << "].g.zS=" << s[i].g.zS << " < cabinet.zS=" << b[0].g.zS << "." << std::endl;
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
	// 2 - RNG (LES).
	// 3 - Spalart-Allmares (RANS). 19.04.2019
	// 4 - K - Omega SST (RANS). 06.10.2019
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

	// Spalart-Allmares (RANS) [1992]. 19.04.2019
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

	// SST Ментер (RANS) [1993]. 03.10.2019
	// А.В. Гарбарук Моделирование турбулентности в расчётах сложных течений. Политех 2012.
	// Эмпирические константы, определяющие SST Ментера модель турбулентности.
	doublereal sigma_k1 = 0.85;
	doublereal sigma_omega1 = 0.5;
	doublereal beta1 = 0.075;
	doublereal sigma_k2 = 1.0;
	doublereal sigma_omega2 = 0.856;
	doublereal beta2 = 0.0828;
	doublereal beta_zvezda = 0.09;
	doublereal menter_a1 = 0.31;
	doublereal alpha1 = 0.555;
	doublereal alpha2 = 0.44;

} FLOWINFO;

// вся информация о полном наборе решаемых уравнений
typedef struct TEQUATIONINFO {
	// нужно ли решать уравнение теплопроводности
	// 0 - ненужно, 1 - нужно методом контрольного объема, 2  - методом конечных элементов.
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
	integer id=-1; // Уникальный номер объединения.

				// Из кабинета юнион видится как Hollow блок
				// в виде прмой прямоугольной призмы.
				// размеры передаются из интерфейса.
	doublereal xS=0.0, xE=0.0, yS=0.0, yE=0.0, zS=0.0, zE=0.0;

	// Внутренняя сетка union.
	// Union является кабинетом для своих внутренних блоков.
	// для внутреннего пользования.
	doublereal *xpos = NULL, *ypos = NULL, *zpos = NULL;
	doublereal *xposadd = NULL, *yposadd = NULL, *zposadd = NULL;

	// Для внутренней сетки (размерности).
	// Передаётся из интерфейса.
	integer inx=-1, iny=-1, inz=-1;
	// для внутреннего пользования.
	integer inxadd = -1, inyadd = -1, inzadd = -1;

	//  Тип сеточного генератора
	// передаётся из интерфейса.
	integer iswitchMeshGenerator = 2; // 2 - CoarseMeshGen

	// Локальные объявления.
	TEMPER t;
	integer flow_interior=-1; // Суммарное число FLUID зон
	FLOW* f = NULL;
} UNION;



// считывание параметров из 
// входного файла premeshin.txt
void premeshin_old(const char *fname, integer &lmatmax, integer &lb, integer &ls, integer &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w, 
	           doublereal &dgx, doublereal &dgy, doublereal &dgz, integer &inx, integer &iny, integer &inz, doublereal &operatingtemperature, 
			   integer &ltdp, TEMP_DEP_POWER* &gtdps, integer &lu, UNION* &my_union) {


	doublereal dmult = 1.0 / rdivision_interval;
	

	// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
	// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
	// при повторном обращении идет обращение к ячейки хеш таблицы.
	// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
	// Хеш таблицы для automatic
	// Инициализация хеш таблицы.
	for (integer i_1 = 0; i_1 < isize_shorter_hash; i_1++) {
		bshorter_hash_X[i_1] = false;
		bshorter_hash_Y[i_1] = false;
		bshorter_hash_Z[i_1] = false;
		shorter_hash_X[i_1] = 1.0e-10;
		shorter_hash_Y[i_1] = 1.0e-10;
		shorter_hash_Z[i_1] = 1.0e-10;
	}

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
			int din = 0;
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
			default : bsnap_TO_global = 1;  break;
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
			// Параметры импульсного режима для SquareWave 2 режима.
			fscanf(fp, "%f", &fin);
			if ((fin<=0.0)||(fin>=1.0)) {
				printf("error input parametr timestep law SquareWave2 multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
				printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;			
			fscanf(fp, "%f", &fin);
			if (fin<=0.0) {
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
			my_amg_manager.amgcl_smoother = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_selector = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_iterator = din;

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
				//printf("%e %e %e %e %e\n", matlist[i].rho, matlist[i].arr_cp[0], matlist[i].arr_lam[0], matlist[i].mu, matlist[i].beta_t);
				//std::cout << matlist[i].rho << " " << matlist[i].arr_cp[0] << " " << matlist[i].arr_lam[0] << " " << matlist[i].mu << " " << matlist[i].beta_t << std::endl;
				//printf("%d %d %d\n", matlist[i].blibmat, matlist[i].ilibident, matlist[i].ilawmu); // bBoussinesq не печатается
				//printf("%e %e %e %e %e %e\n", matlist[i].mumin, matlist[i].mumax, matlist[i].Amu, matlist[i].Bmu, matlist[i].Cmu, matlist[i].degreennmu);
				//std::cout << matlist[i].mumin << " " << matlist[i].mumax << " " << matlist[i].Amu << " " << matlist[i].Bmu << " " << matlist[i].Cmu << " " << matlist[i].degreennmu << std::endl;

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
				//printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				//std::cout << b[i].g.xS << " " << b[i].g.yS << " " << b[i].g.zS << " " << b[i].g.xE << " " << b[i].g.yE << " " << b[i].g.zE << std::endl;
				//printf("%d %d %d\n", b[i].imatid,  b[i].itype, b[i].ipower_time_depend);
				//printf("temperature depend power\n");
				//printf("t_C power_W\n");
				for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					//printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
					//std::cout << b[i].temp_Sc[i_54] << " " << b[i].arr_Sc[i_54] << std::endl;
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
				case XY: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				//printf("%e %d %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
				//std::cout << s[i].power << " " << s[i].iPlane << " " << s[i].g.xS << " " << s[i].g.yS << " " << s[i].g.zS << " " << s[i].g.xE << " " << s[i].g.yE << " " << s[i].g.zE << " " << s[i].square << std::endl;
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
				//printf("%d %e %e %d %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
				//std::cout << w[i].ifamily << " " << w[i].Tamb << " " << w[i].hf << " " << w[i].iPlane << " " << w[i].g.xS << " " << w[i].g.yS << " " << w[i].g.zS << " " << w[i].g.xE << " " << w[i].g.yE << " " << w[i].g.zE << std::endl;
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
					number_iteration_SIMPLE_algorithm = (integer)(din);

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					stabilization_amg1r5_algorithm = din;
				}
			}

			fclose(fp); // закрытие файла
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
				if (vol < 1.0e-40) {
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
			case XY:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-40) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
				}
				break;
			case XZ:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-40) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_in_cyl;
				}
				break;
			case YZ:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-40) {
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
			// Параметры импульсного режима для SquareWave 2 режима.
			fscanf_s(fp, "%f", &fin);
			if ((fin <= 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law SquareWave2 multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau2 = fin;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau_pause = fin;
			fscanf_s(fp, "%d", &din);
			glTSL.n_cycle = din;
			fscanf_s(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 Period must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.T_all = fin;
			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle * (2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl <= 0.0) {
				printf("error in parameters Square Wave SquareWave2 time step law.\n");
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
			my_amg_manager.amgcl_smoother = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_selector = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_iterator = din;

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
				// 2 - square wave 2 зависимость, 3 - hot cold режим.
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
					number_iteration_SIMPLE_algorithm = (integer)(din);

#if doubleintprecision == 1
					fscanf_s(fp, "%lld", &din);
#else
					fscanf_s(fp, "%d", &din);
#endif
					stabilization_amg1r5_algorithm = (integer)(din);

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
						// Параметры импульсного режима для SquareWave 2 режима.
		fscanf_s(fp, "%f", &fin);
		if ((fin <= 0.0) || (fin >= 1.0)) {
			printf("error input parametr timestep law SquareWave2 multiplyer\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.m1 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau1 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau2 = fin;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.tau_pause = fin;
		fscanf_s(fp, "%lld", &din);
		glTSL.n_cycle = din;
		fscanf_s(fp, "%f", &fin);
		if (fin <= 0.0) {
			printf("error input parametr timestep law SquareWave2 Period must be strongly positive\n");
			system("PAUSE");
			exit(1);
		}
		glTSL.T_all = fin;
		doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
		if (t_pause_gl <= 0.0) {
			printf("error in parameters Square Wave 2 time step law.\n");
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
		my_amg_manager.amgcl_smoother = din;

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.amgcl_selector = din;

		fscanf_s(fp, "%lld", &din);
		my_amg_manager.amgcl_iterator = din;

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
			// 2 - square wave 2 зависимость, 3 - hot cold режим.
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
				number_iteration_SIMPLE_algorithm = (integer)(din);

#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &din);
#else
				fscanf_s(fp, "%d", &din);
#endif
				stabilization_amg1r5_algorithm = (integer)(din);

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

				// На тот случай если геометрия состоит только из цилиндров.
				// 13.08.2019 Автомат настройки допусков сетки. 
				// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
				// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
				switch (b[i_1].g.iPlane) {
				case XY :
					if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.Hcyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-40) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_in_cyl;
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
					}
					break;
				case XZ :
					if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.Hcyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-40) {
						// Если внутренний радиус существует (задавался пользователем).
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.R_in_cyl;
						if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_in_cyl;
					}
					break;
				case YZ:
					if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult*b[i_1].g.Hcyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult*b[i_1].g.R_out_cyl;
					if (b[i_1].g.R_in_cyl > 1.0e-40) {
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
		// Запоминаем полное тепловыделение в модели.
		d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = dpower + dsoupow;
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
} // premeshin_old

// Код 2009 года. Из моей первой программы. Функции необхдимые для работы парсера.


// проверяет совпадают ли две строки:
// имя name и образец pattern
bool isname(char * name, char * pattern) {
	bool bret = false;
	int i; // счётчик
	int il1, il2; // длины строк

	i = 0;
	while (name[i] != '\0') i++;
	il1 = i;
	i = 0;
	while (pattern[i] != '\0') i++;
	il2 = i;

	if (il1 == il2) {
		//printf("%s \n", pattern);
		//getchar();
		bret = true;
		for (i = 0; i < il1; i++) {
			if (name[i] != pattern[i]) {
				bret = false;
				return bret;
			}
		}
	}
	//printf("%d \n", (int)(bret));
	return bret;
} // isname 2009


// анализирует строку символов buf длиной ilen.
bool ianalizestring(char *name0, char * buf, int ilen, int &iret)
{
	// состав анализируемой строки :
	// ключевое слово, символ равенства, число без знака.
	bool bfound = false;

	bool bstring = false;
	int inum=-1;
	int i; // счётчик цикла for
	for (i = 0; i < ilen; i++) {
		if (buf[i] == '=') {
			bstring = true;
			inum = i;
		}
	}
	if (bstring) {
		//
		//printf("%d\n",inum);

		// если строка нам подходит по структуре
		char * name = (char *)malloc(1024); // строка имя переменной
		for (i = 0; i < inum; i++) name[i] = buf[i];
		name[i] = '\n'; // конец имени
		name[i++] = '\0'; // конец строки
		//printf("%s \n", name); // debug
		char * value = (char *)malloc(1024); // строка значение переменной с имененм name
		i = inum + 1;
		int j = 0;
		while (buf[i] != '\n') {
			value[j++] = buf[i++];
		}
		value[j++] = '\n';
		value[j] = '\0'; // символ конца строки
		//printf("%s \n",value); // контроль считанного значения
		//getchar();

		if (isname(name, name0)) {
			iret = atoi(value);
			bfound = true;
		}

		free(name);
		free(value);
		return bfound;
		
	}
	
		return bfound;
} // ianalizestring по мотивам 2009

// анализирует строку символов buf длиной ilen.
bool fanalizestring(char *name0, char * buf, int ilen, float &fret)
{
	// состав анализируемой строки :
	// ключевое слово, символ равенства, число без знака.
	bool bfound = false;

	bool bstring = false;
	int inum=-1;
	int i; // счётчик цикла for
	for (i = 0; i < ilen; i++) {
		if (buf[i] == '=') {
			bstring = true;
			inum = i;
		}
	}
	if (bstring) {
		//
		//printf("%d\n",inum);

		// если строка нам подходит по структуре
		char * name = (char *)malloc(1024); // строка имя переменной
		for (i = 0; i < inum; i++) name[i] = buf[i];
		name[i] = '\n'; // конец имени
		name[i++] = '\0'; // конец строки
		//printf("%s \n", name); // debug
		char * value = (char *)malloc(1024); // строка значение переменной с имененм name
		i = inum + 1;
		int j = 0;
		while (buf[i] != '\n') {
			value[j++] = buf[i++];
		}
		value[j++] = '\n';
		value[j] = '\0'; // символ конца строки
		//printf("%s \n",value); // контроль считанного значения
		//getchar();

		if (isname(name, name0)) {
			fret = atof(value);
			bfound = true;
		}

		free(name);
		free(value);
		return bfound;

	}

	return bfound;
} // fanalizestring по мотивам 2009


// Массив строк для ускорения обработки файлов.
const integer ilimit_StringList = 1000000;
integer icurrentSize_StringList = 0;
integer icurrent_position__StringList = 0;
char ** StringList = new char*[ilimit_StringList];

// Загружает данные из входного файла в массив StringList.
// premeshin.txt
// реализовано 25 08 2019.
void loadFromFile()
{
	FILE *fp = NULL;
	errno_t err;

#ifdef MINGW_COMPILLER
	fp = fopen64("premeshin.txt", "r");
	if (fp != NULL) {
		err = 0;
	}
	else {
		err = 1; // ошибка открытия.
	}
#else
	err = fopen_s(&fp, "premeshin.txt", "r");
#endif

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

				StringList[icurrentSize_StringList] = new char[k+1];
				for (integer i_11 = 0; i_11 <= k; i_11++) {
					StringList[icurrentSize_StringList][i_11] = buf[i_11];
				}
				icurrentSize_StringList++;

				if (icurrentSize_StringList > ilimit_StringList) {
					printf("Error!!! buffer StringList is overflow...\n");
					printf("icurrentSize_StringList > ilimit_StringList==400K\n");
					system("pause");
					exit(1);
				}

				//printf("%s \n",buf);
				//getchar();
				b1 = false;
				k = 0;
			}
			if ((b1) && (bbeginstring == true)) {
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
		//getchar();	
	
		free(buf);
		fclose(fp);

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

// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool imakesource(char *name0, int &iret)
{
	bool bfound = false;
	
		
		
		char * buf = (char *)malloc(1024); // строка буфер
		
		bool bfound_loc = false;
		int iret_loc = -1;
		for (integer i_1 = icurrent_position__StringList; i_1 < icurrentSize_StringList; i_1++) {
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
						bfound_loc = false;
						bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
						if (bfound_loc) {
							bfound = true;
							iret = iret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						break;
						//printf("%s \n",buf);
						//getchar();	
					}
					if (c == '\n') {
						bbeginstring = false;
						buf[k++] = '\n';
						buf[k] = '\0';

						bfound_loc = false;
						bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
						if (bfound_loc) {
							bfound = true;
							iret = iret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						//printf("%s \n",buf);
						//getchar();
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
			for (integer i_1 = 0;  i_1 < icurrent_position__StringList; i_1++) {
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
							bfound_loc = false;
							bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
							if (bfound_loc) {
								bfound = true;
								iret = iret_loc;
								icurrent_position__StringList = i_1 + 1;
							}
							break;
							//printf("%s \n",buf);
							//getchar();	
						}
						if (c == '\n') {
							bbeginstring = false;
							buf[k++] = '\n';
							buf[k] = '\0';

							bfound_loc = false;
							bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
							if (bfound_loc) {
								bfound = true;
								iret = iret_loc;
								icurrent_position__StringList = i_1 + 1;
							}
							//printf("%s \n",buf);
							//getchar();
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
		
		return bfound;
	

} // imakesource по мотивам 2010


// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool fmakesource(char *name0, float &fret)
{
	bool bfound = false;



	char * buf = (char *)malloc(1024); // строка буфер

	bool bfound_loc = false;
	float fret_loc = 0.0;
	for (integer i_1 = icurrent_position__StringList; i_1 < icurrentSize_StringList; i_1++) {
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
					bfound_loc = false;
					bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
					if (bfound_loc) {
						bfound = true;
						fret = fret_loc;
						icurrent_position__StringList = i_1 + 1;
					}
					break;
					//printf("%s \n",buf);
					//getchar();	
				}
				if (c == '\n') {
					bbeginstring = false;
					buf[k++] = '\n';
					buf[k] = '\0';

					bfound_loc = false;
					bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
					if (bfound_loc) {
						bfound = true;
						fret = fret_loc;
						icurrent_position__StringList = i_1 + 1;
					}
					//printf("%s \n",buf);
					//getchar();
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

		// Повторяем всё снова.
		for (integer i_1 = 0;   i_1 < icurrent_position__StringList; i_1++) {
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
						bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = fret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						break;
						//printf("%s \n",buf);
						//getchar();	
					}
					if (c == '\n') {
						bbeginstring = false;
						buf[k++] = '\n';
						buf[k] = '\0';

						bfound_loc = false;
						bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
						if (bfound_loc) {
							bfound = true;
							fret = fret_loc;
							icurrent_position__StringList = i_1 + 1;
						}
						//printf("%s \n",buf);
						//getchar();
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

	return bfound;


} // fmakesource по мотивам 2010


// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool imakesource_old(char *name0, int &iret)
{
	bool bfound = false;
	FILE *fp=NULL;
	errno_t err;

#ifdef MINGW_COMPILLER
	fp = fopen64("premeshin.txt", "r");
	if (fp != NULL) {
		err = 0;
	}
	else {
		err = 1; // ошибка открытия.
	}
#else
	err = fopen_s(&fp, "premeshin.txt", "r");
#endif

	// создание файла для чтения.
	if (err  != 0) {
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
				bfound_loc = false;

				bfound_loc=ianalizestring(name0, buf, k - 1,iret_loc);
				if (bfound_loc) {
					bfound = true;
					iret = iret_loc;
					break;
				}
				//printf("%s \n",buf);
				//getchar();
				b1 = false;
				k = 0;
			}
			if ((b1) && (bbeginstring == true)) {
				buf[k++] = c;
				// printf("%c \n",c);
			}

		}// while
		if (!bfound) {
			buf[k++] = '\n';
			buf[k] = '\0';
			bfound_loc = false;
			bfound_loc = ianalizestring(name0, buf, k - 1, iret_loc);
			if (bfound_loc) {
				bfound = true;
				iret = iret_loc;
			}
			//printf("%s \n",buf);
			//getchar();	
		}
		free(buf);
		fclose(fp);

		return bfound;
	}

} // imakesource по мотивам 2010


// считывает данные из входного файла
// premeshin.txt
// реализовано 23 apr 2010.
// revised 21 августа 2019.
bool fmakesource_old(char *name0, float &fret)
{
	bool bfound = false;
	FILE *fp = NULL;
	errno_t err;
	// создание файла для чтения.
#ifdef MINGW_COMPILLER
	fp = fopen64("premeshin.txt", "r");
	if (fp != NULL) {
		err = 0;
	}
	else {
		err = 1; // ошибка открытия.
	}
#else
	err = fopen_s(&fp, "premeshin.txt", "r");
#endif

	if (err  != 0) {
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
		c = fgetc(fp); // читает один символ
		buf[k++] = c;
		bool bfound_loc = false;
		float fret_loc = 0.0;
		while (!feof(fp)) {
			b1 = true;
			c = fgetc(fp); // читает один символ
			if (c == '#') bbeginstring = false;
			if (c == '\n') {
				bbeginstring = true;
				buf[k++] = '\n';
				buf[k] = '\0';
				bfound_loc = false;
				bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
				if (bfound_loc) {
					bfound = true;
					fret = fret_loc;
					break;
				}
				//printf("%s \n",buf);
				//getchar();
				b1 = false;
				k = 0;
			}
			if ((b1) && (bbeginstring == true)) {
				buf[k++] = c;
				// printf("%c \n",c);
			}

		}// while
		if (!bfound) {
			buf[k++] = '\n';
			buf[k] = '\0';
			bfound_loc = false;
			bfound_loc = fanalizestring(name0, buf, k - 1, fret_loc);
			if (bfound_loc) {
				bfound = true;
				fret = fret_loc;
			}
			//printf("%s \n",buf);
			//getchar();	
		}
		free(buf);
		fclose(fp);

		return bfound;
	}

} // fmakesource по мотивам 2010


// считывание параметров из 
// входного файла premeshin.txt
void premeshin_new(const char *fname, integer &lmatmax, integer &lb, integer &ls, integer &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w,
	doublereal &dgx, doublereal &dgy, doublereal &dgz, integer &inx, integer &iny, integer &inz, doublereal &operatingtemperature,
	integer &ltdp, TEMP_DEP_POWER* &gtdps, integer &lu, UNION* &my_union) {


	doublereal dmult = 1.0 / rdivision_interval;


	// Так как в режиме bFULL_AUTOMATIC допуски определяются локально с
	// помощью тяжеловесной функции, то значения функции вычисляются лишь один раз, а
	// при повторном обращении идет обращение к ячейки хеш таблицы.
	// 20mm ПТБШ ускорился с 1мин 9с до 53с за счет режима bFULL_AUTOMATIC.
	// Хеш таблицы для automatic
	// Инициализация хеш таблицы.
	for (integer i_1 = 0; i_1 < isize_shorter_hash; i_1++) {
		bshorter_hash_X[i_1] = false;
		bshorter_hash_Y[i_1] = false;
		bshorter_hash_Z[i_1] = false;
		shorter_hash_X[i_1] = 1.0e-10;
		shorter_hash_Y[i_1] = 1.0e-10;
		shorter_hash_Z[i_1] = 1.0e-10;
	}

#ifdef MINGW_COMPILLER
	// eqin - информация о наборе решаемых уравнений.

	// dgx, dgy, dgz - вектор силы тяжести.
	// inx, iny, inz - количество точек по каждой из осей.

	FILE *fp;
	errno_t err1 = 0;
	fp = fopen64(fname, "r");
	if (err1 != 0) {
		printf("No input File premeshin.txt \n");
		//system("PAUSE");
		system("pause");

	}
	else
	{
		if (fp != NULL) {
			float fin = 0.0;
			int din = 0;
			doublereal scale = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания

			fscanf(fp, "%d", &din);
			ionly_solid_visible = din;
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
			operating_temperature_for_film_coeff = fin;


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
			Tochka_position_X0_for_XY_Plot = scale * fin;
			fscanf(fp, "%f", &fin);
			Tochka_position_Y0_for_XY_Plot = scale * fin;
			fscanf(fp, "%f", &fin);
			Tochka_position_Z0_for_XY_Plot = scale * fin;
			// Направление линии совпадает с направлением 
			// одной из осей декартовой прямоугольной системы координат:
			// 0 - Ox, 1 - Oy, 2 - Oz.
			fscanf(fp, "%d", &din);
			idirectional_for_XY_Plot = din;

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
			if ((fin <= 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law Factor a\n");
				system("PAUSE");
				exit(1);
			}
			fscanf(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law tau must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau = fin; // длительность импульса.
			fscanf(fp, "%f", &fin);
			glTSL.Q = fin;  // Скважность.
			// Параметры импульсного режима для SquareWave 2 режима.
			fscanf(fp, "%f", &fin);
			if ((fin <= 0.0) || (fin >= 1.0)) {
				printf("error input parametr timestep law SquareWave2 multiplyer\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.m1 = fin;
			fscanf(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau1 = fin;
			fscanf(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau2 = fin;
			fscanf(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.tau_pause = fin;
			fscanf(fp, "%d", &din);
			glTSL.n_cycle = din;
			fscanf(fp, "%f", &fin);
			if (fin <= 0.0) {
				printf("error input parametr timestep law SquareWave2 Period must be strongly positive\n");
				system("PAUSE");
				exit(1);
			}
			glTSL.T_all = fin;
			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl <= 0.0) {
				printf("error in parameters SquareWave2 time step law.\n");
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
			if (din == 0) {
				// Обычная структурированная сетка.
				b_on_adaptive_local_refinement_mesh = false;
			}
			else {
				// АЛИС
				b_on_adaptive_local_refinement_mesh = true;
			}
			fscanf(fp, "%d", &din);
			itype_ALICE_Mesh = din;
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
			my_amg_manager.amgcl_smoother = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_selector = din;

			fscanf(fp, "%d", &din);
			my_amg_manager.amgcl_iterator = din;

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
				matlist[i].orthotropy_multiplyer_x = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_y = fin;
				fscanf(fp, "%f", &fin);
				matlist[i].orthotropy_multiplyer_z = fin;
				// 5.08.2017.
				// Коэффициенты для задачи упругости.
				// Модуль Юнга и коэффициент Пуассона.
				doublereal Poissonratio = 0.154;
				doublereal Youngmodule = 217.5e9;
				fscanf(fp, "%f", &fin);
				Poissonratio = fin;
				fscanf(fp, "%f", &fin);
				Youngmodule = fin * 1e9;
				fscanf(fp, "%f", &fin);
				matlist[i].beta_t_solid = fin * 1E-6;
				// Коэффициенты Лямэ.
				doublereal E1_koef = Youngmodule / (1.0 - Poissonratio * Poissonratio);
				doublereal nu1_koef = Poissonratio / (1.0 - Poissonratio);
				matlist[i].mu_Lame = E1_koef / (2.0*(1.0 + nu1_koef));
				matlist[i].lambda_Lame = (E1_koef*nu1_koef) / (1.0 - nu1_koef * nu1_koef);
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
			for (i = 0; i < lb; i++) {

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
				b[i].g.xS = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.yS = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.zS = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.xE = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.yE = scale * fin;
				fscanf(fp, "%f", &fin);
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

				fscanf(fp, "%d", &din);
				b[i].g.iPlane = din;

				fscanf(fp, "%f", &fin);
				b[i].g.xC = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.yC = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.zC = scale * fin;
				fscanf(fp, "%f", &fin);
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


				fscanf(fp, "%f", &fin);
				b[i].g.R_out_cyl = scale * fin;
				fscanf(fp, "%f", &fin);
				b[i].g.R_in_cyl = scale * fin;

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
					b[i].g.hi[i73] = scale * fin;
					fscanf(fp, "%f", &fin);
					b[i].g.xi[i73] = scale * fin;
					fscanf(fp, "%f", &fin);
					b[i].g.yi[i73] = scale * fin;
					fscanf(fp, "%f", &fin);
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
				b[i].ipower_time_depend = din;
				// тип блока
				fscanf(fp, "%d", &din);
				b[i].itype = din;

				// печать считанных значений на консоль
				//printf("%e %e %e %e %e %e\n", b[i].g.xS, b[i].g.yS, b[i].g.zS, b[i].g.xE, b[i].g.yE, b[i].g.zE);
				//std::cout << b[i].g.xS << " " << b[i].g.yS << " " << b[i].g.zS << " " << b[i].g.xE << " " << b[i].g.yE << " " << b[i].g.zE << std::endl;
				//printf("%d %d %d\n", b[i].imatid,  b[i].itype, b[i].ipower_time_depend);
				//printf("temperature depend power\n");
				//printf("t_C power_W\n");
				for (integer i_54 = 0; i_54 < b[i].n_Sc; i_54++) {
					//printf("%e %e\n", b[i].temp_Sc[i_54], b[i].arr_Sc[i_54]);
					//std::cout << b[i].temp_Sc[i_54] << " " << b[i].arr_Sc[i_54] << std::endl;
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
				s[i].g.xS = scale * fin;
				fscanf(fp, "%f", &fin);
				s[i].g.yS = scale * fin;
				fscanf(fp, "%f", &fin);
				s[i].g.zS = scale * fin;
				fscanf(fp, "%f", &fin);
				s[i].g.xE = scale * fin;
				fscanf(fp, "%f", &fin);
				s[i].g.yE = scale * fin;
				fscanf(fp, "%f", &fin);
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
				case XY: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				//printf("%e %d %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
				//std::cout << s[i].power << " " << s[i].iPlane << " " << s[i].g.xS << " " << s[i].g.yS << " " << s[i].g.zS << " " << s[i].g.xE << " " << s[i].g.yE << " " << s[i].g.zE << " " << s[i].square << std::endl;
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
					w[i].emissivity = 0.0;
					w[i].film_coefficient = 0.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // первого рода
				case 2:  fscanf(fp, "%f", &fin);
					w[i].Tamb = 0.0;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// adiabatic wall
					w[i].emissivity = 0.0;
					w[i].film_coefficient = 0.0;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // однородное условие Неймана
				case 3:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// Newton-Richman condition, film coefficient.
					w[i].emissivity = 0.0;
					w[i].film_coefficient = fin;
					fscanf(fp, "%f", &fin);
					w[i].hf = 0.0;
					break; // Ньютон-Рихман.
				case 4:  fscanf(fp, "%f", &fin);
					w[i].Tamb = fin;
					fscanf(fp, "%f", &fin); // Stefan Bolcman
					// Stefan - Bolcman condition
					w[i].emissivity = fin;
					w[i].film_coefficient = 0.0;
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
				w[i].g.xS = scale * fin;
				fscanf(fp, "%f", &fin);
				w[i].g.yS = scale * fin;
				fscanf(fp, "%f", &fin);
				w[i].g.zS = scale * fin;
				fscanf(fp, "%f", &fin);
				w[i].g.xE = scale * fin;
				fscanf(fp, "%f", &fin);
				w[i].g.yE = scale * fin;
				fscanf(fp, "%f", &fin);
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
				//w[i].bfixboundary = false;// Свободная граница.
				//printf("%d %e %e %d %e %e %e %e %e %e\n", w[i].ifamily, w[i].Tamb, w[i].hf, w[i].iPlane, w[i].g.xS, w[i].g.yS, w[i].g.zS, w[i].g.xE, w[i].g.yE, w[i].g.zE);
				//std::cout << w[i].ifamily << " " << w[i].Tamb << " " << w[i].hf << " " << w[i].iPlane << " " << w[i].g.xS << " " << w[i].g.yS << " " << w[i].g.zS << " " << w[i].g.xE << " " << w[i].g.yE << " " << w[i].g.zE << std::endl;
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
				my_union[i].xS = scale * fin;
				fscanf(fp, "%f", &fin);
				my_union[i].xE = scale * fin;
				fscanf(fp, "%f", &fin);
				my_union[i].yS = scale * fin;
				fscanf(fp, "%f", &fin);
				my_union[i].yE = scale * fin;
				fscanf(fp, "%f", &fin);
				my_union[i].zS = scale * fin;
				fscanf(fp, "%f", &fin);
				my_union[i].zE = scale * fin;

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
					eqin.fluidinfo[i].xc = scale * fin;
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].yc = scale * fin;
					fscanf(fp, "%f", &fin);
					eqin.fluidinfo[i].zc = scale * fin;
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
					number_iteration_SIMPLE_algorithm = (integer)(din);

#if doubleintprecision == 1
					fscanf(fp, "%lld", &din);
#else
					fscanf(fp, "%d", &din);
#endif
					stabilization_amg1r5_algorithm = (integer)(din);

				}
			}

			fclose(fp); // закрытие файла
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
			if (b[i_1].n_Sc > 0) {
				doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
				doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS)*fabs(b[i_1].g.yE - b[i_1].g.yS)*fabs(b[i_1].g.zE - b[i_1].g.zS);
				if (vol < 1.0e-40) {
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
		if (b[i_1].g.itypegeom == CYLINDER) {
			// Cylinder

			// На тот случай если геометрия состоит только из цилиндров.
			// 13.08.2019 Автомат настройки допусков сетки. 
			// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
			// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
			switch (b[i_1].g.iPlane) {
			case XY:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-40) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
				}
				break;
			case XZ:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-40) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
				}
				break;
			case YZ:
				if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.Hcyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
				if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
				if (b[i_1].g.R_in_cyl > 1.0e-40) {
					// Если внутренний радиус существует (задавался пользователем).
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
					if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
				}
				break;
			}


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

					dpower += pdiss * vol;
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
	//printf("Blocks integral power =%e W\n", dpower);
	std::cout << "Blocks integral power =" << dpower << " W" << std::endl;
	printf("number of sources ls=%lld\n", ls);
	//printf("Sources integral power = %e W\n", dsoupow);
	std::cout << "Sources integral power =" << dsoupow << " W" << std::endl;
	//printf("Full total power = %e W\n", dpower + dsoupow);
	std::cout << "Full total power =" << (dpower + dsoupow) << " W" << std::endl;
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

	// eqin - информация о наборе решаемых уравнений.

	// dgx, dgy, dgz - вектор силы тяжести.
	// inx, iny, inz - количество точек по каждой из осей.

	FILE* fp=NULL;
	errno_t err1;
	if ((err1 = fopen_s(&fp, fname, "r")) != 0) {
		printf("No input File premeshin.txt \n");
		//system("PAUSE");
		system("pause");

	}
	else
	{
		if (fp != NULL) {

			fclose(fp);

			float fin = 0.0;
			integer din = 0;
			int idin = 0;
			doublereal scale = 1.0;
			doublereal dbuf; // для упорядочивания в порядке возрастания

			// Приостанавливает считывание если переменная не найдена.
			bool bSTOP_Reading = true; 


			if (imakesource("only_solid_visible", idin)) {
			    // Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					ionly_solid_visible = (integer)(idin);
					printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
				}
				else {
					printf("ionly_solid_visible must be equal 0 or 1. now value=%lld\n",idin);
					system("pause");
					ionly_solid_visible = 0;
				}
			}
			else {
				printf("WARNING!!! only_solid_visible not found in file premeshin.txt\n");
				ionly_solid_visible = 0; // Показываем всё
				printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
				if (bSTOP_Reading) system("pause");
			}
		

			if (fmakesource("mlength", fin)) {
				// Найдено успешно.
				if (scale <= 0.0) {
					printf("ERROR!!! scale must be >0.0\n");
					system("pause");
					exit(1);
				}
				scale = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! mlength not found in file premeshin.txt\n");
				scale = 1.0e-3; // mm
				printf("scale =%e\n", scale);
				if (bSTOP_Reading) system("pause");
			}


			
			if (imakesource("lmatmax", idin)) {
				// Найдено успешно.
				if (idin < 1) {
					printf("Error: lmatmax mast be > 0. now value=%lld\n",idin);
					system("pause");
					exit(1);
				}
				else {
					lmatmax = (integer)(idin);
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
					printf("Error: lb mast be > 0. now value=%lld\n", idin);
					system("pause");
					exit(1);
				}
				else {
					lb = (integer)(idin);
					//printf("lb =%lld\n", lb);
				}
			}
			else {
				printf("WARNING!!! lb not found in file premeshin.txt\n");
				lb = 0; // нет блоков
				printf("lb =%lld\n", lb);
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("ls", idin)) {
				// Найдено успешно.
				if (idin >= 0) {
					ls = (integer)(idin);
					//printf("ls =%lld\n", ls);
				}
				else {
					printf("Error: ls = %lld must be >=0\n",idin);
					system("pause");
					exit(1);
				}
			}
			else {
				printf("WARNING!!! ls not found in file premeshin.txt\n");
				ls = 0; // нет плоских бесконечно тонких источников тепла.
				printf("ls =%lld\n", ls);
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("lw", idin)) {
				// Найдено успешно.
				if (idin >= 0) {
					lw = (integer)(idin);
					//printf("lw =%lld\n", lw);
				}
				else {
					printf("Error: lw = %lld must be >=0\n", idin);
					system("pause");
					exit(1);
				}
			}
			else {
				printf("WARNING!!! lw not found in file premeshin.txt\n");
				lw = 0; // нет плоских бесконечно тонких источников тепла.
				printf("lw =%lld\n", lw);
				if (bSTOP_Reading) system("pause");
			}


			// количество уникальных данных с табличными данными 
			// по зависимости расеиваемой мощности от температуры.
			if (imakesource("iltdp", idin)) {
				// Найдено успешно.
				if (idin >= 0) {
				    ltdp = (integer)(idin);
				    //printf("ltdp =%lld\n", ltdp);
				}
				else {
					printf("Error: iltdp = %lld must be >=0\n", idin);
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
				dgx = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! gx not found in file premeshin.txt\n");
				dgx = 0.0; // Сила тяжести по направлению Ох отсутствует.
				printf("gx =%e\n", dgx);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("gy", fin)) {
				// Найдено успешно.
				dgy = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! gy not found in file premeshin.txt\n");
				dgy = 0.0; // Сила тяжести по направлению Оy отсутствует.
				printf("gy =%e\n", dgy);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("gz", fin)) {
				// Найдено успешно.
				dgz = (doublereal)(fin);
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
					inx = (integer)(idin);
					//printf("inx =%lld\n", inx);
				}
				else {
					printf("ERROR: inx=%lld\n",inx);
					system("pause");
				}
			}
			else {
				printf("WARNING!!! inx not found in file premeshin.txt\n");
				inx = 0; // нет информации о количестве узлов сетки в направлении оси Ох.
				printf("inx =%lld\n", inx);
				if (bSTOP_Reading) system("pause");
			}

		
			if (imakesource("iny", idin)) {
				// Найдено успешно.
				if (iny > 0) {
					iny = (integer)(idin);
					//printf("iny =%lld\n", iny);
				}
			    else {
				   printf("ERROR: inx=%lld\n", inx);
				   system("pause");
		     	}
			}
			else {
				printf("WARNING!!! iny not found in file premeshin.txt\n");
				iny = 0; // нет информации о количестве узлов сетки в направлении оси Оy.
				printf("iny =%lld\n", iny);
				if (bSTOP_Reading) system("pause");
			}


			if (imakesource("inz", idin)) {
				// Найдено успешно.
				if (inz>0) {
				   inz = (integer)(idin);
				   //printf("inz =%lld\n", inz);
				}
				else {
					printf("ERROR: inx=%lld\n", inx);
					system("pause");
				}
			}
			else {
				printf("WARNING!!! inz not found in file premeshin.txt\n");
				inz = 0; // нет информации о количестве узлов сетки в направлении оси Оz.
				printf("inz =%lld\n", inz);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("operatingtemperature", fin)) {
				// Найдено успешно.
				if (fin < -273.15) {
					fin = -273.15;
					printf("Error!!! operatingtemperature < -273.15 C\n");
					system("pause");
				}
				operatingtemperature = (doublereal)(fin);  // Operating Temperature
				operating_temperature_for_film_coeff = (doublereal)(fin);
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
				starting_speed_Vx = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! SpeedInitializationVx not found in file premeshin.txt\n");
				starting_speed_Vx = 0.0; // инициализируем нулевой скоростью в направлении оси Ох.
				printf("starting_speed_Vx =%e\n", starting_speed_Vx);
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("SpeedInitializationVy", fin)) {
				// Найдено успешно.
				starting_speed_Vy = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! SpeedInitializationVy not found in file premeshin.txt\n");
				starting_speed_Vy = 0.0; // инициализируем нулевой скоростью в направлении оси Оy.
				printf("starting_speed_Vy =%e\n", starting_speed_Vy);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("SpeedInitializationVz", fin)) {
				// Найдено успешно.
				starting_speed_Vz = (doublereal)(fin);
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
				Tochka_position_X0_for_XY_Plot = scale * ((doublereal)(fin));
			}
			else {
				printf("WARNING!!! XYPlotXo not found in file premeshin.txt\n");
				Tochka_position_X0_for_XY_Plot = 0.0; 
				printf("Tochka_position_X0_for_XY_Plot =%e\n", Tochka_position_X0_for_XY_Plot);
				if (bSTOP_Reading) system("pause");
			}

			
			if (fmakesource("XYPlotYo", fin)) {
				// Найдено успешно.
				Tochka_position_Y0_for_XY_Plot = scale * ((doublereal)(fin));
			}
			else {
				printf("WARNING!!! XYPlotYo not found in file premeshin.txt\n");
				Tochka_position_Y0_for_XY_Plot = 0.0;
				printf("Tochka_position_Y0_for_XY_Plot =%e\n", Tochka_position_Y0_for_XY_Plot);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("XYPlotZo", fin)) {
				// Найдено успешно.
				Tochka_position_Z0_for_XY_Plot = scale * ((doublereal)(fin));
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
				idirectional_for_XY_Plot  = (integer)(idin);
				//printf("iny =%lld\n", iny);
			}
			else {
				printf("WARNING!!! line_directional not found in file premeshin.txt\n");
				idirectional_for_XY_Plot = 0; // нет информации о количестве узлов сетки в направлении оси Оy.
				printf("idirectional_for_XY_Plot =%lld\n", idirectional_for_XY_Plot);
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
				etalon_max_size_ratio = (doublereal)(fin); // подробность расчётной сетки.
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
				etalon_max_size_ratio2 = (doublereal)(fin); 
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
					printf("Error!!! snap_to_grid must be equal 0, 1, 2 or 3. Current value = %lld\n",idin);
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

			// Выбор решающего устройства : либо amg1r5 либо BiCGStab+ILU2.
			if (imakesource("SolverSetting", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin < 14)) {
					iswitchsolveramg_vs_BiCGstab_plus_ILU2 = (integer)(idin);
					//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);
				}
				else {
					printf("Error !!! must be 0<= iswitchsolveramg_vs_BiCGstab_plus_ILU2 <14. current value iswitchsolveramg_vs_BiCGstab_plus_ILU2=%lld\n", idin);
					system("pause");
					idin = 0; // BiCGStab+ilu2 solver default.
					iswitchsolveramg_vs_BiCGstab_plus_ILU2 = (integer)(idin);
					//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);

				}
			}
			else {
				printf("WARNING!!! SolverSetting not found in file premeshin.txt\n");
				iswitchsolveramg_vs_BiCGstab_plus_ILU2 = 0; //BiCGStab +ILU* solver.
				printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);
				if (bSTOP_Reading) system("pause");
			}

			 
			// Выбор решающего устройства : либо РУМБА0.14 либо BiCGStab+ILU6.
			if (imakesource("StaticStructuralSolverSetting", idin)) {
				// Найдено успешно.
				iswitchsolveramg_vs_BiCGstab_plus_ILU6 = (integer)(idin);
				//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU6 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU6);
			}
			else {
				printf("WARNING!!! StaticStructuralSolverSetting not found in file premeshin.txt\n");
				iswitchsolveramg_vs_BiCGstab_plus_ILU6 = 0; //BiCGStab +ILU* solver.
				printf("iswitchsolveramg_vs_BiCGstab_plus_ILU6 =%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU6);
				if (bSTOP_Reading) system("pause");
			}
			

			if (imakesource("PressureVelocityCoupling", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					if (idin == 1) {
						// SIMPLEC algorithm.
						iSIMPLE_alg = SIMPLEC_Van_Doormal_and_Raithby;
					}
					else {
						// SIMPLE algorithm 1972.
						iSIMPLE_alg = SIMPLE_Carretto;
					}
					//printf("iSIMPLE_alg =%d\n", idin);
				}
				else {
					printf("ERROR!!! iSIMPLE_alg must be equal 0 or 1. Current value =%lld\n",idin);
					system("pause");
					// SIMPLE algorithm 1972.
					iSIMPLE_alg = SIMPLE_Carretto;
				}
			}
			else {
				printf("WARNING!!! PressureVelocityCoupling not found in file premeshin.txt\n");
				// SIMPLE algorithm 1972.
				iSIMPLE_alg = SIMPLE_Carretto;
				printf("iSIMPLE_alg = SIMPLE_Carretto\n");
				if (bSTOP_Reading) system("pause");
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
				default: iFLOWScheme = 2; break; // UDS самая стабильная схема.
				}
				//printf("iFLOWScheme =%lld\n", iFLOWScheme);
			}
			else {
				printf("WARNING!!! FlowScheme not found in file premeshin.txt\n");
				iFLOWScheme = 2; // UDS самая стабильная схема.
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
				default: iTEMPScheme = 2; break; // UDS самая стабильная схема.
				}
				//printf("iTEMPScheme  =%lld\n", iTEMPScheme);
			}
			else {
				printf("WARNING!!! SchemeTemperature not found in file premeshin.txt\n");
				iTEMPScheme = 2; // UDS самая стабильная схема.
				printf("iTEMPScheme = Upwind \n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("mesh_generator_algorithm", idin)) {
				// Найдено успешно.
				// Выбор сеточного генератора.
				iswitchMeshGenerator = (integer)(idin);
				//printf("iswitchMeshGenerator =%lld\n", iswitchMeshGenerator);
			}
			else {
				printf("WARNING!!! mesh_generator_algorithm not found in file premeshin.txt\n");
				iswitchMeshGenerator = 2; // Coarse Mesh.
				printf("iswitchMeshGenerator =CoarseMesh(2)\n");
				if (bSTOP_Reading) system("pause");
			}
			

			if (imakesource("Schedule", idin)) {
				// Найдено успешно.
				steady_or_unsteady_global_determinant = 2;
				din = (integer)(idin);
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
				//printf("steady_or_unsteady_global_determinant  =%lld\n", steady_or_unsteady_global_determinant );
			}
			else {
				printf("WARNING!!! Schedule not found in file premeshin.txt\n");
				steady_or_unsteady_global_determinant = 2; // Mesh Generator only
				printf("steady_or_unsteady_global_determinant =Mesh Generator only\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("TimeStep", idin)) {
				// Найдено успешно.
				// Закон изменения шага по времени.
				din = (integer)(idin);
				if ((din == 0) || (din == 1) || (din == 2) || (din == 3)) {
					glTSL.id_law = din;
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
				glTSL.id_law = 0; // default Linear.
				printf("glTSL.id_law =Linear\n");
				if (bSTOP_Reading) system("pause");
			}
			
			if (fmakesource("Factor_a_for_Linear", fin)) {
				// Найдено успешно.
				if ((fin <= 0.0) || (fin >= 1.0)) {
					printf("error input parametr timestep law Factor a=%e\n",fin);
					system("PAUSE");
					exit(1);
				}
				glTSL.Factor_a_for_Linear = (doublereal)(fin); // Factor_a
			}
			else {
				printf("WARNING!!! Factor_a_for_Linear not found in file premeshin.txt\n");
				glTSL.Factor_a_for_Linear = 0.3; // Знаменатель геометрической прогрессии.
				printf("glTSL.Factor_a_for_Linear = 0.3\n");
				if (bSTOP_Reading) system("pause");
			}


			
			
			if (fmakesource("tau", fin)) {
				// Найдено успешно.
				if (fin <= 0.0) {
					printf("error input parametr timestep law tau must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau = (doublereal)(fin); // длительность импульса.
				
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
				glTSL.Q = (doublereal)(fin); // Скважность.
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
				if ((fin <= 0.0) || (fin >= 1.0)) {
					printf("error input parametr timestep law SquareWave2 multiplyer\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.m1 = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! m1 not found in file premeshin.txt\n");
				glTSL.m1 = 0.33333; // SquareWave2 parameters.
				printf("glTSL.m1 =%e\n", glTSL.m1);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("tau1", fin)) {
				// Найдено успешно.
				if (fin <= 0.0) {
					printf("error input parametr timestep law SquareWave2 tau1 must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau1 = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! tau1 not found in file premeshin.txt\n");
				glTSL.tau1 = 180.0; // SquareWave2 parameters.
				printf("glTSL.tau1 =%e\n", glTSL.tau1);
				if (bSTOP_Reading) system("pause");
			}
					   			
			if (fmakesource("tau2", fin)) {
				// Найдено успешно.
				if (fin <= 0.0) {
					printf("error input parametr timestep law SquareWave2 tau2 must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau2 = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! tau2 not found in file premeshin.txt\n");
				glTSL.tau2 = 240.0; // SquareWave2 parameters.
				printf("glTSL.tau2 =%e\n", glTSL.tau2);
				if (bSTOP_Reading) system("pause");
			}

			if (fmakesource("tau_pause", fin)) {
				// Найдено успешно.
				if (fin <= 0.0) {
					printf("error input parametr timestep law SquareWave2 tau_pause must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.tau_pause = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! tau_pause not found in file premeshin.txt\n");
				glTSL.tau_pause = 5040.0; // SquareWave2 parameters.
				printf("glTSL.tau_pause =%e\n", glTSL.tau_pause);
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("n", idin)) {
				// Найдено успешно.
				glTSL.n_cycle = (integer)(idin);
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
				if (fin <= 0.0) {
					printf("error input parametr timestep law SquareWave2 T: Period must be strongly positive=%e\n",fin);
					system("PAUSE");
					exit(1);
				}
				glTSL.T_all = (doublereal)(fin);				
			}
			else {
				printf("WARNING!!! T not found in file premeshin.txt\n");
				glTSL.T_all = 86400.0; // SquareWave2 parameters.
				printf("glTSL.T_all =%e\n", glTSL.T_all);
				if (bSTOP_Reading) system("pause");
			}

			doublereal t_pause_gl = glTSL.T_all - glTSL.n_cycle*(2 * glTSL.tau1 + glTSL.tau2 + glTSL.tau_pause);
			if (t_pause_gl <= 0.0) {
				printf("error in parameters Square Wave SquareWave2 time step law.\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}

			
			if (fmakesource("on_time_double_linear", fin)) {
				// Найдено успешно.
				if (fin <= 0.0) {
					printf("error input parametr on_time_double_linear law hot cold reshime must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				glTSL.on_time_double_linear = (doublereal)(fin);
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
				if (fin <= 0.0) {
					printf("error input parametr globalEndTimeUnsteadyTemperatureCalculation must be strongly positive\n");
					system("PAUSE");
					exit(1);
				}
				globalEndTimeUnsteadyTemperatureCalculation = (doublereal)(fin);
			}
			else {
				printf("WARNING!!! EndTime not found in file premeshin.txt\n");
				globalEndTimeUnsteadyTemperatureCalculation = 2.0; // Продолжительность моделирования. 
				printf("globalEndTimeUnsteadyTemperatureCalculation =%e\n", globalEndTimeUnsteadyTemperatureCalculation);
				if (bSTOP_Reading) system("pause");
			}
			
			// Newton-Richman condition.
			if (imakesource("adiabatic_vs_heat_transfer_coeff", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) || (idin == 2) || (idin == 3)) {
					// 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
					adiabatic_vs_heat_transfer_coeff = (integer)(idin);
					//printf("adiabatic_vs_heat_transfer_coeff =%lld\n", adiabatic_vs_heat_transfer_coeff);
				}
				else {
					printf("Error!!! adiabatic_vs_heat_transfer_coeff must be equal 0 or 1 or 2 or 3. Current value %lld\n",idin);
					system("pause");
					adiabatic_vs_heat_transfer_coeff = 0; // adiabatic wall
				}
			}
			else {
				printf("WARNING!!! adiabatic_vs_heat_transfer_coeff not found in file premeshin.txt\n");
				adiabatic_vs_heat_transfer_coeff = 0; //  0 - adiabatic wall.
				printf("adiabatic_vs_heat_transfer_coeff =%lld\n", adiabatic_vs_heat_transfer_coeff);
				if (bSTOP_Reading) system("pause");
			}			
			
			if (fmakesource("filmcoefficient", fin)) {
				// Найдено успешно.
				film_coefficient = (doublereal)(fin);// Коэффициент теплоотдачи.
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
				itype_ALICE_Mesh = (integer)(idin);
				//printf("itype_ALICE_Mesh =%lld\n", itype_ALICE_Mesh);
			}
			else {
				printf("WARNING!!! version_ALICE_Mesh not found in file premeshin.txt\n");
				itype_ALICE_Mesh = 0; // Coarse Mesh.
				printf("itype_ALICE_Mesh =%lld\n", itype_ALICE_Mesh);
				if (bSTOP_Reading) system("pause");
			}

			
			if (imakesource("m_restart_gmres", idin)) {
				// Найдено успешно.
				my_amg_manager.m_restart = (integer)(idin);
				//printf("my_amg_manager.m_restart =%lld\n", my_amg_manager.m_restart);
			}
			else {
				printf("WARNING!!! m_restart_gmres not found in file premeshin.txt\n");
				my_amg_manager.m_restart = 20; // Количество рестартов агоритма Flexible GMRES.
				printf("my_amg_manager.m_restart =%lld\n", my_amg_manager.m_restart);
				if (bSTOP_Reading) system("pause");
			}


			// classical algebraic multigrid parameters:
			// only for my_agregat_amg.cu.
			if (imakesource("amg_manager_sorting_alg", idin)) {
				// Найдено успешно.
				my_amg_manager.imySortAlgorithm = (integer)(idin);
				//printf("my_amg_manager.imySortAlgorithm =%lld\n", my_amg_manager.imySortAlgorithm);
			}
			else {
				printf("WARNING!!! amg_manager_sorting_alg not found in file premeshin.txt\n");
				my_amg_manager.imySortAlgorithm = 0; // Counting Sort Саусвела.
				printf("my_amg_manager.imySortAlgorithm =%lld\n", my_amg_manager.imySortAlgorithm);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("inum_reduced_levels", idin)) {
				// Найдено успешно.
				my_amg_manager.maximum_delete_levels_Temperature = (integer)(idin);
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
				my_amg_manager.maximum_delete_levels_Speed  = (integer)(idin);
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
				my_amg_manager.maximum_delete_levels_Pressure  = (integer)(idin);
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
				my_amg_manager.maximum_delete_levels_Stress = (integer)(idin);
				//printf("my_amg_manager.maximum_delete_levels_Stress =%lld\n", my_amg_manager.maximum_delete_levels_Stress);
			}
			else {
				printf("WARNING!!! inum_reduced_levels_Stress not found in file premeshin.txt\n");
				my_amg_manager.maximum_delete_levels_Stress = 0; // нет удаляемых уровней.
				printf("my_amg_manager.maximum_delete_levels_Stress =%lld\n", my_amg_manager.maximum_delete_levels_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// type interpolation procedure :			
			if (imakesource("interpolation", idin)) {
				// Найдено успешно.
				my_amg_manager.number_interpolation_procedure_Temperature = (integer)(idin);
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
				my_amg_manager.number_interpolation_procedure_Speed = (integer)(idin);
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
				my_amg_manager.number_interpolation_procedure_Pressure = (integer)(idin);
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
				my_amg_manager.number_interpolation_procedure_Stress  = (integer)(idin);
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
				my_amg_manager.iCFalgorithm_and_data_structure_Temperature = (integer)(idin);
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Temperature =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Temperature);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Temperature not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Temperature = 3; //  3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Temperature =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Temperature);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("CFalgorithmandDataStruct_Speed", idin)) {
				// Найдено успешно.
				my_amg_manager.iCFalgorithm_and_data_structure_Speed = (integer)(idin);
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Speed =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Speed);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Speed not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Speed = 3; //  3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Speed =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("CFalgorithmandDataStruct_Pressure", idin)) {
				// Найдено успешно.
				my_amg_manager.iCFalgorithm_and_data_structure_Pressure = (integer)(idin);
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Pressure =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Pressure);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Pressure not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Pressure = 0; //  3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Pressure =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Pressure);
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("CFalgorithmandDataStruct_Stress", idin)) {
				// Найдено успешно.
				my_amg_manager.iCFalgorithm_and_data_structure_Stress  = (integer)(idin);
				//printf("my_amg_manager.iCFalgorithm_and_data_structure_Stress =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Stress);
			}
			else {
				printf("WARNING!!! CFalgorithmandDataStruct_Stress not found in file premeshin.txt\n");
				my_amg_manager.iCFalgorithm_and_data_structure_Stress = 3; // 3-Treap.
				printf("my_amg_manager.iCFalgorithm_and_data_structure_Stress =%lld\n", my_amg_manager.iCFalgorithm_and_data_structure_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// Устарело:
			//my_amg_manager.itypemodifyinterpol = din;
			//my_amg_manager.baglomeration_with_consistency_scaling = din;
			
			if (imakesource("modinterpol", idin)) {
				// Найдено успешно.
				my_amg_manager.bdiagonal_dominant = (integer)(idin);
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
					my_amg_manager.bTemperatureMatrixPortrait = (integer)(idin);
					//printf("my_amg_manager.bTemperatureMatrixPortrait =%lld\n", my_amg_manager.bTemperatureMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bTemperatureMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bTemperatureMatrixPortrait = 0; // 0 - NO_PRINT, 1 - PRINT.
				}
			}
			else {
				printf("WARNING!!! TemperatureMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bTemperatureMatrixPortrait = 0; // 0 - NO_PRINT, 1 - PRINT.
				printf("my_amg_manager.bTemperatureMatrixPortrait =%lld\n", my_amg_manager.bTemperatureMatrixPortrait);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("SpeedMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.bSpeedMatrixPortrait = (integer)(idin);
					//printf("my_amg_manager.bSpeedMatrixPortrait =%lld\n", my_amg_manager.bSpeedMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bSpeedMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bSpeedMatrixPortrait = 0; // 0 - NO_PRINT, 1 - PRINT.
				}
			}
			else {
				printf("WARNING!!! SpeedMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bSpeedMatrixPortrait = 0; //  0 - NO_PRINT, 1 - PRINT.
				printf("my_amg_manager.bSpeedMatrixPortrait =%lld\n", my_amg_manager.bSpeedMatrixPortrait);
				if (bSTOP_Reading) system("pause");
			}
						
			if (imakesource("PressureMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.bPressureMatrixPortrait = (integer)(idin);
					//printf("my_amg_manager.bPressureMatrixPortrait =%lld\n", my_amg_manager.bPressureMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bPressureMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bPressureMatrixPortrait = 0; // 0 - NO_PRINT, 1 - PRINT.
				}
			}
			else {
				printf("WARNING!!! PressureMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bPressureMatrixPortrait = 0; //  0 - NO_PRINT, 1 - PRINT.
				printf("my_amg_manager.bPressureMatrixPortrait =%lld\n", my_amg_manager.bPressureMatrixPortrait);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("StressMatrixPortrait", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.bStressMatrixPortrait = (integer)(idin);
					//printf("my_amg_manager.bStressMatrixPortrait =%lld\n", my_amg_manager.bStressMatrixPortrait);
				}
				else {
					printf("my_amg_manager.bStressMatrixPortrait must be equal 1 or 0. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.bStressMatrixPortrait = 0; // 0 - NO_PRINT, 1 - PRINT.
				}
			}
			else {
				printf("WARNING!!! StressMatrixPortrait not found in file premeshin.txt\n");
				my_amg_manager.bStressMatrixPortrait = 0; //  0 - NO_PRINT, 1 - PRINT.
				printf("my_amg_manager.bStressMatrixPortrait =%lld\n", my_amg_manager.bStressMatrixPortrait);
				if (bSTOP_Reading) system("pause");
			}

			// truncation of interpolation:
			if (imakesource("truncationT", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1)) {
					my_amg_manager.itruncation_interpolation_Temperature = (integer)(idin);
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
					my_amg_manager.itruncation_interpolation_Speed = (integer)(idin);
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
					my_amg_manager.itruncation_interpolation_Pressure = (integer)(idin);
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
					my_amg_manager.itruncation_interpolation_Stress = (integer)(idin);
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
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Temperature = (doublereal)(fin);
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
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Speed = (doublereal)(fin);
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
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Pressure = (doublereal)(fin);
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
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.truncation_interpolation_Stress = (doublereal)(fin);
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

			// number nFinnest sweeps :
			if (imakesource("nFinnest", idin)) {
				// Найдено успешно.
				my_amg_manager.nFinnest_Temperature = (integer)(idin);
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
				my_amg_manager.nFinnest_Speed = (integer)(idin);
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
				my_amg_manager.nFinnest_Pressure = (integer)(idin);
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
				my_amg_manager.nFinnest_Stress = (integer)(idin);
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
				my_amg_manager.nu1_Temperature = (integer)(idin);
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
				my_amg_manager.nu1_Speed = (integer)(idin);
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
				my_amg_manager.nu1_Pressure = (integer)(idin);
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
				my_amg_manager.nu1_Stress = (integer)(idin);
				//printf("my_amg_manager.nu1_Stress =%lld\n", my_amg_manager.nu1_Stress);
			}
			else {
				printf("WARNING!!! numberpresmoothersStress not found in file premeshin.txt\n");
				my_amg_manager.nu1_Stress = 1; // Количество предварительных сглаживаний.
				printf("my_amg_manager.nu1_Stress =%lld\n", my_amg_manager.nu1_Stress);
				if (bSTOP_Reading) system("pause");
			}

			// number postsweeps :
			if (imakesource("numberpostsweeps", idin)) {
				// Найдено успешно.
				my_amg_manager.nu2_Temperature = (integer)(idin);
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
				my_amg_manager.nu2_Speed = (integer)(idin);
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
				my_amg_manager.nu2_Pressure = (integer)(idin);
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
				my_amg_manager.nu2_Stress = (integer)(idin);
				//printf("my_amg_manager.nu2_Stress =%lld\n", my_amg_manager.nu2_Stress);
			}
			else {
				printf("WARNING!!! numberpostsweepsStress not found in file premeshin.txt\n");
				my_amg_manager.nu2_Stress = 2; // Количество пост сглаживаний.
				printf("my_amg_manager.nu2_Stress =%lld\n", my_amg_manager.nu2_Stress);
				if (bSTOP_Reading) system("pause");
			}


			// memory size :
		    if (imakesource("memorysize", idin)) {
				// Найдено успешно.
				my_amg_manager.memory_size_Temperature = (integer)(idin);
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
				my_amg_manager.memory_size_Speed = (integer)(idin);
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
				my_amg_manager.memory_size_Pressure = (integer)(idin);
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
				my_amg_manager.memory_size_Stress = (integer)(idin);
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
				din = (integer)(idin);
				my_amg_manager.ilu2_smoother_Temperature = din;
				my_amg_manager.b_gmresTemp = false;				

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
				printf("my_amg_manager.ilu2_smoother_Temperature =%lld\n", my_amg_manager.ilu2_smoother_Temperature);
				if (bSTOP_Reading) system("pause");
			}
					   		
			if (imakesource("smoothertypeSpeed", idin)) {
				// Найдено успешно.
				din = (integer)(idin);
				my_amg_manager.ilu2_smoother_Speed = din;
				my_amg_manager.b_gmresSpeed = false;				

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
				printf("my_amg_manager.ilu2_smoother_Speed =%lld\n", my_amg_manager.ilu2_smoother_Speed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("smoothertypePressure", idin)) {
				// Найдено успешно.
				din = (integer)(idin);
				my_amg_manager.ilu2_smoother_Pressure = din;
				my_amg_manager.b_gmresPressure = false;				

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
				printf("my_amg_manager.ilu2_smoother_Pressure =%lld\n", my_amg_manager.ilu2_smoother_Pressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("smoothertypeStress", idin)) {
				// Найдено успешно.
				din = (integer)(idin);
				my_amg_manager.ilu2_smoother_Stress = din;
				my_amg_manager.b_gmresStress = false;				

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
				printf("my_amg_manager.ilu2_smoother_Stress =%lld\n", my_amg_manager.ilu2_smoother_Stress);
				if (bSTOP_Reading) system("pause");
			}

			my_amg_manager.gold_const_Temperature = return_gold_const(my_amg_manager.ilu2_smoother_Temperature);
			my_amg_manager.gold_const_Speed = return_gold_const(my_amg_manager.ilu2_smoother_Speed);
			my_amg_manager.gold_const_Pressure = return_gold_const(my_amg_manager.ilu2_smoother_Pressure);
			my_amg_manager.gold_const_Stress = return_gold_const(my_amg_manager.ilu2_smoother_Stress);


			// strength threshold :
			
			if (fmakesource("threshold", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.theta_Temperature = (doublereal)(fin);
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
					my_amg_manager.theta_Speed = (doublereal)(fin);
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
					my_amg_manager.theta_Pressure = (doublereal)(fin);
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
					my_amg_manager.theta_Stress = (doublereal)(fin);
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

			// magic threshold :
			// magic <=> F_to_F
			if (fmakesource("magicT", fin)) {
				// Найдено успешно.
				if ((fin > 0.0) && (fin < 1.0)) {
					my_amg_manager.F_to_F_Temperature = (doublereal)(fin);
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
					my_amg_manager.F_to_F_Speed = (doublereal)(fin);
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
					my_amg_manager.F_to_F_Pressure = (doublereal)(fin);
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
					my_amg_manager.F_to_F_Stress = (doublereal)(fin);
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
			// Способ построения C-F разбиения : 0 - standart, 1 - RS 2.
			// RS 2 улучшенная версия построения C-F разбиения содержащая второй проход.
			if (imakesource("coarseningTemp", idin)) {
				// Найдено успешно.
				my_amg_manager.icoarseningTemp = (integer)(idin);
				//printf("my_amg_manager.icoarseningTemp =%lld\n", my_amg_manager.icoarseningTemp);
			}
			else {
				printf("WARNING!!! coarseningTemp not found in file premeshin.txt\n");
				my_amg_manager.icoarseningTemp = 0; // 0 - standart, 1 - RS 2.
				printf("my_amg_manager.icoarseningTemp =%lld\n", my_amg_manager.icoarseningTemp);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("coarseningSpeed", idin)) {
				// Найдено успешно.
				my_amg_manager.icoarseningSpeed = (integer)(idin);
				//printf("my_amg_manager.icoarseningSpeed =%lld\n", my_amg_manager.icoarseningSpeed);
			}
			else {
				printf("WARNING!!! coarseningSpeed not found in file premeshin.txt\n");
				my_amg_manager.icoarseningSpeed = 0; // 0 - standart, 1 - RS 2.
				printf("my_amg_manager.icoarseningSpeed =%lld\n", my_amg_manager.icoarseningSpeed);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("coarseningPressure", idin)) {
				// Найдено успешно.
				my_amg_manager.icoarseningPressure = (integer)(idin);
				//printf("my_amg_manager.icoarseningPressure =%lld\n", my_amg_manager.icoarseningPressure);
			}
			else {
				printf("WARNING!!! coarseningPressure not found in file premeshin.txt\n");
				my_amg_manager.icoarseningPressure = 0; // 0 - standart, 1 - RS 2.
				printf("my_amg_manager.icoarseningPressure =%lld\n", my_amg_manager.icoarseningPressure);
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("coarseningStress", idin)) {
				// Найдено успешно.
				my_amg_manager.icoarseningStress = (integer)(idin);
				//printf("my_amg_manager.icoarseningStress =%lld\n", my_amg_manager.icoarseningStress);
			}
			else {
				printf("WARNING!!! coarseningStress not found in file premeshin.txt\n");
				my_amg_manager.icoarseningStress = 0; // 0 - standart, 1 - RS 2.
				printf("my_amg_manager.icoarseningStress =%lld\n", my_amg_manager.icoarseningStress);
				if (bSTOP_Reading) system("pause");
			}

			// Если din==0 то просто алгебраический многосеточный метод без привлечения алгоритмов подпространства Крылова,
			// Если din==1, Stabilization BiCGStab.
			// 8.01.2017 Метод ван дер Ворста BiCGStab 
			// предобусловленный алгебраичесеким многосеточным методом.
			// 9.01.2018 Если din==2, FGMRes предобусловленный алгебраическим многосеточным методом.
			if (imakesource("StabilizationTemp", idin)) {
				// Найдено успешно.
				if ((idin == 0) || (idin == 1) || (idin == 2) || (idin == 3)) {
					my_amg_manager.istabilizationTemp = (integer)(idin);
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
					my_amg_manager.istabilizationSpeed = (integer)(idin);
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
				if ((idin == 0) || (idin == 1) || (idin == 2)) {
					my_amg_manager.istabilizationPressure = (integer)(idin);
					//printf("my_amg_manager.istabilizationPressure =%lld\n", my_amg_manager.istabilizationPressure);
				}
				else {
					printf("my_amg_manager.istabilizationPressure must be equal 2, 1 or 0. current_value=%d\n", idin);
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
				if ((idin == 0) || (idin == 1) || (idin == 2)) {
					my_amg_manager.istabilizationStress = (integer)(idin);
					//printf("my_amg_manager.istabilizationStress =%lld\n", my_amg_manager.istabilizationStress);
				}
				else {
					printf("my_amg_manager.istabilizationStress must be equal 2, 1 or 0. current_value=%d\n", idin);
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
					my_amg_manager.ipatch_number = (integer)(idin);
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
					my_amg_manager.iprint_log_Temperature = (integer)(idin);
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
					my_amg_manager.iprint_log_Speed = (integer)(idin);
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
					my_amg_manager.iprint_log_Pressure = (integer)(idin);
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
					my_amg_manager.iprint_log_Stress = (integer)(idin);
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

			if (imakesource("amgcl_smoother", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 7)) {
					my_amg_manager.amgcl_smoother = (integer)(idin);
					//printf("my_amg_manager.amgcl_smoother =%lld\n", my_amg_manager.amgcl_smoother);
				}
				else {
					printf("my_amg_manager.amgcl_smoother must be 0<= value <=7. current_value=%d\n", idin);
					system("pause");
					// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped-jacobi.
					// 4 - spai1; 5 - chebyshev; 6 - ilu1 (iluk,k=1); 7 - ilu2 (iluk,k=2);
					my_amg_manager.amgcl_smoother = 0; 
				}
			}
			else {
				printf("WARNING!!! amgcl_smoother not found in file premeshin.txt\n");
				// 0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped-jacobi.
                // 4 - spai1; 5 - chebyshev; 6 - ilu1 (iluk,k=1); 7 - ilu2 (iluk,k=2);
				my_amg_manager.amgcl_smoother = 0; 
				printf("my_amg_manager.amgcl_smoother =%lld\n", my_amg_manager.amgcl_smoother);
				printf("0 - spai0; 1 - ilu0; 2 - gauss-seidel; 3 - damped-jacobi;\n");
				printf("4 - spai1; 5 - chebyshev; 6 - ilu1 (iluk,k=1); 7 - ilu2 (iluk,k=2).\n");
				if (bSTOP_Reading) system("pause");
			}

			if (imakesource("amgcl_selector", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 1)) {
					my_amg_manager.amgcl_selector = (integer)(idin);
					//printf("my_amg_manager.amgcl_selector =%lld\n", my_amg_manager.amgcl_selector);
				}
				else {
					printf("my_amg_manager.amgcl_selector must be 0<= value <=1. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.amgcl_selector = 1; // 0 - Ruge-Stueben (analog amg1r5); 1 - smoother aggregation.
				}
			}
			else {
				printf("WARNING!!! amgcl_selector not found in file premeshin.txt\n");
				my_amg_manager.amgcl_selector = 1; // 0 - Ruge-Stueben (analog amg1r5); 1 - smoother aggregation.
				printf("my_amg_manager.amgcl_selector =%lld\n", my_amg_manager.amgcl_selector);
				printf(" 0 - Ruge-Stueben (analog amg1r5); 1 - smoother aggregation. \n");
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("amgcl_iterator", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 1)) {
					my_amg_manager.amgcl_iterator = (integer)(idin);
					//printf("my_amg_manager.amgcl_iterator =%lld\n", my_amg_manager.amgcl_iterator);
				}
				else {
					printf("my_amg_manager.amgcl_iterator must be 0<= value <=1. current_value=%d\n", idin);
					system("pause");
					my_amg_manager.amgcl_iterator = 0; // 0 - BiCGStab; 1 - FGMRes.
				}
			}
			else {
				printf("WARNING!!! amgcl_iterator not found in file premeshin.txt\n");
				my_amg_manager.amgcl_iterator = 0; // 0 - BiCGStab; 1 - FGMRes.
				printf("my_amg_manager.amgcl_iterator =%lld\n", my_amg_manager.amgcl_iterator);
				printf("0 - BiCGStab; 1 - FGMRes.\n");
				if (bSTOP_Reading) system("pause");
			}
			
			if (imakesource("lfil", idin)) {
				// Найдено успешно.
				if ((idin >= 0) && (idin <= 7)) {
					my_amg_manager.lfil = (integer)(idin);
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
					ireconstruction_free_construct_alloc = (integer)(idin);
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
					ianimation_write_on = (integer)(idin);
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

			char buffer[1000];

			// считывание базы материалов
			for (i = 0; i < lmatmax; i++) {
				// свойства материалов:
				// плотность
			
				char name0[1000] = "matherial";
				
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				
				strcat_s(name0,"_rho");
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].rho = (doublereal)(fin);
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
					matlist[i].n_cp = (integer)(idin);
					//printf("matlist[i].n_cp =%lld\n",matlist[i].n_cp );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					matlist[i].n_cp = 1; // .
					printf("matlist[i].n_cp =%lld\n", matlist[i].n_cp);
					if (bSTOP_Reading) system("pause");
				}

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
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "temp_cp");

					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].temp_cp[i_4] = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						matlist[i].temp_cp[i_4] = 30.0; // .
						printf(" matlist[%lld].temp_cp[%lld]=%e\n", i, i_4, matlist[i].temp_cp[i_4]);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "arr_cp");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].arr_cp[i_4] = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
						matlist[i].arr_cp[i_4] = 700.0; // .
						printf("matlist[%lld].arr_cp[%lld]  =%e\n",i, i_4, matlist[i].arr_cp[i_4]);
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
					matlist[i].n_lam = (integer)(idin);
					//printf("matlist[%lld].n_lam =%lld\n", i, matlist[i].n_lam);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].n_lam = 1; // Количество различных значений теплопропроводности при её зависимости от температуры.
					printf("matlist[%lld].n_lam =%lld\n",i, matlist[i].n_lam);
					if (bSTOP_Reading) system("pause");
				}

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
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "temp_lam");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);

					// Температура в C.
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].temp_lam[i_4] = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
						matlist[i].temp_lam[i_4] = 30.0; // Температура при которой задана теплопроводность.
						printf(" matlist[%lld].temp_lam[%lld]=%e\n", i, i_4, matlist[i].temp_lam[i_4]);
						if (bSTOP_Reading) system("pause");
					}

					
					name0[0] = '\0'; strcat_s(name0, "matherial");
					buffer[0] = '\0';
					_itoa_s(i, buffer, 10);
					strcat_s(name0, buffer);
					strcat_s(name0, "arr_lam");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);


					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						matlist[i].arr_lam[i_4] = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
						matlist[i].arr_lam[i_4] = 4.0; // ЭЧЭС default.
						printf("matlist[%lld].arr_lam[%lld] =%e\n", i, i_4, matlist[i].arr_lam[i_4]);
						if (bSTOP_Reading) system("pause");
					}
				}
				// ортотропность теплопроводности :
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_lam_x");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_x = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].orthotropy_multiplyer_x = 1.0; // no orthotropy .
					printf("matlist[%lld].orthotropy_multiplyer_x =%e\n",i, matlist[i].orthotropy_multiplyer_x);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_lam_y");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_y = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].orthotropy_multiplyer_y = 1.0; // no orthotropy .
					printf("matlist[%lld].orthotropy_multiplyer_y =%e\n", i, matlist[i].orthotropy_multiplyer_y);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mult_lam_z");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].orthotropy_multiplyer_z = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].orthotropy_multiplyer_z = 1.0; // no orthotropy .
					printf(" matlist[%lld].orthotropy_multiplyer_z=%e\n", i, matlist[i].orthotropy_multiplyer_z);
					if (bSTOP_Reading) system("pause");
				}
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

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "Poisson_ratio");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					Poissonratio  = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					Poissonratio  = 0.3; // default 0.3.
					printf(" Poissonratio =%e\n", Poissonratio);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "Young_Module");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					Youngmodule = (doublereal)(fin)* 1e9;
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					Youngmodule  = 200.0E9; // default 200.0e9.
					printf("Youngmodule =%e\n", Youngmodule);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "Linear_expansion_coefficient");
				
				// beta_t_solid*1E-6
				 
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].beta_t_solid = (doublereal)(fin)* 1E-6;
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].beta_t_solid = 1.0e-6; // default 1.0e-6.
					printf(" matlist[%lld].beta_t_solid =%e\n", i, matlist[i].beta_t_solid);
					if (bSTOP_Reading) system("pause");
				}
				// Коэффициенты Лямэ.
				//doublereal E1_koef = Youngmodule / (1.0 - Poissonratio*Poissonratio);
				//doublereal nu1_koef = Poissonratio / (1.0 - Poissonratio);
				//matlist[i].mu_Lame = E1_koef / (2.0*(1.0 + nu1_koef));
				//matlist[i].lambda_Lame = (E1_koef*nu1_koef) / (1.0 - nu1_koef*nu1_koef);
				// стр. 25 В.Н.Сидоров, В.В. Вершинин Метод конечных элементов в расчёте сооружений.
				//+ 19.10.2018 проверено.
				matlist[i].mu_Lame = Youngmodule / (2.0*(1.0 + Poissonratio));
				matlist[i].lambda_Lame = (Poissonratio*Youngmodule) / ((1.0 + Poissonratio)*(1.0 - 2.0*Poissonratio));
				//printf("E=%e N/m^2 mu=%e lambda=%e\n", Youngmodule, matlist[i].mu_Lame, matlist[i].lambda_Lame);
				//system("PAUSE");
				// коэффициент динамической вязкости
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "mu");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].mu = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].mu = 1.7894e-5; // air.
					printf("matlist[%lld].mu =%e\n", i, matlist[i].mu);
					if (bSTOP_Reading) system("pause");
				}

				// коэффициент линейного температурного расширения
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0';
				_itoa_s(i, buffer, 10);
				strcat_s(name0, buffer);
				strcat_s(name0, "beta_t");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].beta_t = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].beta_t = 0.003331; // air.
					printf(" matlist[%lld].beta_t =%e\n", i, matlist[i].beta_t);
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
					matlist[i].blibmat = (integer)(idin);
					//printf(" matlist[%lld].blibmat=%lld\n", i, matlist[i].blibmat);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].blibmat = 0; // не библиотечный материал.
					printf(" matlist[%lld].blibmat=%lld\n",i, matlist[i].blibmat);
					if (bSTOP_Reading) system("pause");
				}
				// номер материала в библиотеке
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ilibident");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].ilibident = (integer)(idin);
					//printf("matlist[%lld].ilibident =%lld\n", i, matlist[i].ilibident);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].ilibident = 0; // .
					printf(" matlist[%lld].ilibident=%lld\n", i, matlist[i].ilibident);
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
					din= (integer)(idin);
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
					printf("matlist[%lld].bBussineskApproach=true\n",i);
					if (bSTOP_Reading) system("pause");
				}
				// номер закона для зависимости динамической вязкости от напряжения сдвига
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ilawmu");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					matlist[i].ilawmu = (integer)(idin);
					//printf("matlist[%lld].ilawmu =%lld\n",i,  matlist[i].ilawmu);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].ilawmu = 0; // .
					printf(" matlist[%lld].ilawmu =%lld\n",i, matlist[i].ilawmu);
					if (bSTOP_Reading) system("pause");
				}
				// минимальное значение динамической вязкости
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "mumin");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].mumin = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].mumin = 0.0; // .
					printf("matlist[%lld].mumin =%e\n", i, matlist[i].mumin);
					if (bSTOP_Reading) system("pause");
				}

				// максимальное значение динамической вязкости
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "mumax");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].mumax = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].mumax = 1.0; // .
					printf("matlist[%lld].mumax =%e\n", i, matlist[i].mumax);
					if (bSTOP_Reading) system("pause");
				}
				// параметры модельных законов для зависимости вязкости от напряжения сдвига
				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Amu");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].Amu = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].Amu = 0.0; // .
					printf("matlist[%lld].Amu =%e\n", i, matlist[i].Amu);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Bmu");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].Bmu = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].Bmu = 0.0; // .
					printf(" matlist[%lld].Bmu =%e\n", i, matlist[i].Bmu);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Cmu");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].Cmu = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].Cmu = 0.0; // .
					printf("matlist[%lld].Cmu  =%e\n",i, matlist[i].Cmu );
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "matherial");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "degreennmu");
				// показатель степени
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					matlist[i].degreennmu = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n",name0);
					matlist[i].degreennmu = 1.0; // .
					printf("matlist[%lld].degreennmu =%e\n", i, matlist[i].degreennmu);
					if (bSTOP_Reading) system("pause");
				}

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
			for (i = 0; i < lb; i++) {

				char name0[1000] = "body";
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iunion");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].iunion_id  = (integer)(idin); // 0 == Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					//printf(" b[%lld].iunion_id=%lld\n",i, b[i].iunion_id);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].iunion_id = 0; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					printf(" b[%lld].iunion_id=%lld\n",i, b[i].iunion_id);
					if (bSTOP_Reading) system("pause");
				}


				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "igeometry_type");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din= (integer)(idin);
					CHECK_TYPE_GEOM(din);
					b[i].g.itypegeom = din;
					//printf(" b[%lld].g.itypegeom=%lld\n",i,b[i].g.itypegeom );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.itypegeom = 0; // 0 - Prism, 1 - Cylinder, 2 - Polygon.
					printf(" b[%lld].g.itypegeom=%lld\n", i, b[i].g.itypegeom);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bvisible");
								
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);
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
					printf(" b[%lld].bvisible=true\n",i );
					if (bSTOP_Reading) system("pause");
				}

				// геометрия
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.xS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.xS = scale*0.0; // .
					printf(" b[%lld].g.xS=%e\n",i, b[i].g.xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.yS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.yS = scale * 0.0; // .
					printf(" b[%lld].g.yS=%e\n",i, b[i].g.yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");
								
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.zS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.zS = scale * 0.0; // .
					printf(" b[%lld].g.zS =%e\n", i, b[i].g.zS);
					if (bSTOP_Reading) system("pause");
				}
	
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.xE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.xE = scale * 0.0; // .
					printf("b[%lld].g.xE =%e\n",i, b[i].g.xE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.yE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.yE = scale * 0.0; // .
					printf("b[%lld].g.yE =%e\n",i, b[i].g.yE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");
			
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.zE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.zE = scale * 0.0; // .
					printf("b[%lld].g.zE =%e\n",i, b[i].g.zE);
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
					b[i].g.iPlane = (integer)(idin);
					//printf("b[%lld].g.iPlane =%lld\n",i,b[i].g.iPlane );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.iPlane = 0; // .
					printf("b[%lld].g.iPlane =%lld\n", i, b[i].g.iPlane);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xC");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.xC = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.xC = scale * 0.0; // .
					printf("b[%lld].g.xC =%e\n",i, b[i].g.xC);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yC");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.yC = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.yC = scale * 0.0; // .
					printf("b[%lld].g.yC =%e\n", i, b[i].g.yC);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zC");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.zC = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.zC = scale * 0.0; // .
					printf("b[%lld].g.zC =%e\n", i, b[i].g.zC);
					if (bSTOP_Reading) system("pause");
				}
			
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Hcyl");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.Hcyl = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.Hcyl = scale * 0.0; // .
					printf("b[%lld].g.Hcyl =%e\n",i, b[i].g.Hcyl);
					if (bSTOP_Reading) system("pause");
				}
				
				if (b[i].g.Hcyl <= 0.0) {
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

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "R_out_cyl");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.R_out_cyl = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.R_out_cyl = scale * 0.0; // .
					printf("b[%lld].g.R_out_cyl =%e\n",i, b[i].g.R_out_cyl);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "R_in_cyl");
								 
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].g.R_in_cyl = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.R_in_cyl = scale * 0.0; // .
					printf("b[%lld].g.R_in_cyl =%e\n",i, b[i].g.R_in_cyl);
					if (bSTOP_Reading) system("pause");
				}

				// Polygon
				//printf("Polygon\n");
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iPlane_obj2");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].g.iPlane_obj2 = (integer)(idin);
					//printf("b[%lld].g.iPlane_obj2 =%lld\n", i, b[i].g.iPlane_obj2);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.iPlane_obj2 = 0; // .
					printf("b[%lld].g.iPlane_obj2 =%lld\n", i, b[i].g.iPlane_obj2);
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
					b[i].g.nsizei = (integer)(idin);
					//printf(" b[%lld].g.nsizei=%lld\n",i, b[i].g.nsizei );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].g.nsizei = 3; // .
					printf("b[%lld].g.nsizei =%lld\n", i, b[i].g.nsizei);
					if (bSTOP_Reading) system("pause");
				}
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
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "hi");

					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.hi[i73] = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.hi[i73] = scale * 0.0; // .
						printf(" b[%lld].g.hi[%lld]=%e\n", i,i73,  b[i].g.hi[i73]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "xi");
					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.xi[i73] = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.xi[i73] = scale * 0.0; // .
						printf("b[%lld].g.xi[%lld] =%e\n",i, i73, b[i].g.xi[i73]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "yi");
					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.yi[i73] = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.yi[i73] = scale * 0.0; // .
						printf("b[%lld].g.yi[%lld] =%e\n", i, i73, b[i].g.yi[i73]);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "zi");
					buffer[0] = '\0'; _itoa_s(i73,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].g.zi[i73] = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].g.zi[i73] = scale * 0.0; // .
						printf(" b[%lld].g.zi[%lld]=%e\n",i,i73, b[i].g.zi[i73]);
						if (bSTOP_Reading) system("pause");
					}

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
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissW");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissW = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissW = 0.0; // изоляция.
					printf("b[%lld].radiation.emissW =%e\n", i, b[i].radiation.emissW);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissE");				
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissE = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissE = 0.0; // изоляция.
					printf("b[%lld].radiation.emissE =%e\n", i, b[i].radiation.emissE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissS = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissS = 0.0; // изоляция.
					printf("b[%lld].radiation.emissS =%e\n", i, b[i].radiation.emissS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissN");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissN = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissN = 0.0; // изоляция.
					printf("b[%lld].radiation.emissN  =%e\n", i, b[i].radiation.emissN);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissB");			
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissB = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissB = 0.0; // изоляция.
					printf("b[%lld].radiation.emissB =%e\n", i, b[i].radiation.emissB);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "emissT");			
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					b[i].radiation.emissT = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].radiation.emissT = 0.0; // изоляция.
					printf("b[%lld].radiation.emissT =%e\n",i, b[i].radiation.emissT);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "binternalRadiation");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);
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
					printf("b[%lld].radiation.binternalRadiation =false\n", i);
					if (bSTOP_Reading) system("pause");
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
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "imatid");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].imatid = (integer)(idin);
					//printf(" b[%lld].imatid=%lld\n", i, b[i].imatid);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].imatid = 0; // идентификатор материала блока.
					printf(" b[%lld].imatid=%lld\n", i, b[i].imatid);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bCylinderFixed");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);
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
					printf("b[%lld].CylinderFixed = false\n", i);
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
					b[i].n_Sc = (integer)(idin);
					//printf("b[%lld].n_Sc =%lld\n",i,b[i].n_Sc );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].n_Sc = 1; // .
					printf(" b[%lld].n_Sc=%lld\n", i, b[i].n_Sc);
					if (bSTOP_Reading) system("pause");
				}

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
					name0[0] = '\0'; strcat_s(name0, "body");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "temp_power");
					buffer[0] = '\0'; _itoa_s(i_4,buffer,10); strcat_s(name0, buffer);

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						b[i].temp_Sc[i_4] = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].temp_Sc[i_4] = 30.0; // .
						printf("b[%lld].temp_Sc[%lld] =%e\n", i, i_4, b[i].temp_Sc[i_4]);
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
							b[i].arr_Sc[i_4] = (doublereal)(fin);
						}
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						b[i].arr_Sc[i_4] = 0.0; // .
						printf("b[%lld].arr_Sc[%lld] =%e\n",i,i_4, b[i].arr_Sc[i_4]);
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
				// 2 - square wave #2 зависимость, 3 - hot cold режим.
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ipower_time_depend");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].ipower_time_depend = (integer)(idin);
					//printf("b[%lld].ipower_time_depend =%lld\n", i, b[i].ipower_time_depend);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].ipower_time_depend = 0; // не зависит от времени.
					printf("b[%lld].ipower_time_depend =%lld\n", i, b[i].ipower_time_depend);
					if (bSTOP_Reading) system("pause");
				}
				// тип блока
				name0[0] = '\0'; strcat_s(name0, "body");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "itype");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					b[i].itype = (integer)(idin);
					//printf("b[%lld].itype  =%lld\n",i,b[i].itype  );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					b[i].itype = 0; // 0 - SOLID, 1 - HOLLOW, 2 - FLUID.
					printf(" b[%lld].itype =%lld\n", i, b[i].itype);
					if (bSTOP_Reading) system("pause");
				}

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

				char name0[1000] = "source";
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iunion");			
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					s[i].iunion_id = (integer)(idin);
					//printf(" s[%lld].iunion_id=%lld\n",i,s[i].iunion_id );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].iunion_id = 0; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					printf(" s[%lld].iunion_id=%lld\n",i, s[i].iunion_id);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Power");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].power = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].power = 0.0; // .
					printf(" s[%lld].power=%e\n",i, s[i].power);
					if (bSTOP_Reading) system("pause");
				}


				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "itempdep");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);
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
					printf("s[%lld].bgarber_depend = false\n",i );
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "id_table");

				// уникальный номер таблицы.
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					s[i].igarber_depend = (integer)(idin);
					//printf(" s[%lld].igarber_depend=%lld\n", i, s[i].igarber_depend);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].igarber_depend = 0; // .
					printf("s[%lld].igarber_depend =%lld\n", i, s[i].igarber_depend);
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
					s[i].roperation_offset_drain = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].roperation_offset_drain = 28.0; // .
					printf("s[%lld].roperation_offset_drain =%e\n",i, s[i].roperation_offset_drain);
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
					s[i].iPlane = (integer)(idin);
					//printf(" s[%lld].iPlane=%lld\n",i,s[i].iPlane );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].iPlane = 0; // .
					printf("s[%lld].iPlane =%lld\n",i, s[i].iPlane);
					if (bSTOP_Reading) system("pause");
				}
				// геометрия
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.xS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.xS = scale * 0.0; // .
					printf(" s[%lld].g.xS=%e\n",i, s[i].g.xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.yS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.yS = scale * 0.0; // .
					printf("s[%lld].g.yS  =%e\n", i, s[i].g.yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.zS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.zS = scale * 0.0; // .
					printf("s[%lld].g.zS =%e\n", i, s[i].g.zS);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.xE  = scale *  (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.xE = scale *  0.0; // .
					printf(" s[%lld].g.xE=%e\n",i, s[i].g.xE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.yE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.yE = scale *  0.0; // .
					printf(" s[%lld].g.yE =%e\n",i, s[i].g.yE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "source");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					s[i].g.zE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					s[i].g.zE = scale * 0.0; // .
					printf(" s[%lld].g.zE=%e\n",i, s[i].g.zE);
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
				case XY: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.yE - s[i].g.yS); break;
				case XZ: s[i].square = fabs(s[i].g.xE - s[i].g.xS)*fabs(s[i].g.zE - s[i].g.zS); break;
				case YZ: s[i].square = fabs(s[i].g.yE - s[i].g.yS)*fabs(s[i].g.zE - s[i].g.zS); break;
				default: break;
				}
				//printf("source %e %lld %e %e %e %e %e %e %e\n", s[i].power, s[i].iPlane, s[i].g.xS, s[i].g.yS, s[i].g.zS, s[i].g.xE, s[i].g.yE, s[i].g.zE, s[i].square);
			}

			// считывание твёрдых стенок

			for (i = 0; i < lw; i++) {

				char name0[1000] = "wall";
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iunion");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					w[i].iunion_id = (integer)(idin);
					//printf(" w[%lld].iunion_id=%lld\n", i, w[i].iunion_id);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].iunion_id = 0; // 0==Кабинет, номер АССЕМБЛЕСА которому принадлежит.
					printf(" w[%lld].iunion_id=%lld\n",i, w[i].iunion_id);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "family");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					w[i].ifamily = (integer)(idin);
					din= (integer)(idin); // for switch() см. далее.
					//printf("w[%lld].ifamily =%lld\n", i, w[i].ifamily);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].ifamily = 2; // Однородное условие Неймана.
					printf("w[%lld].ifamily =%lld\n", i, w[i].ifamily);
					if (bSTOP_Reading) system("pause");
				}

				switch (din) {
				case 1:  
					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "Tamb");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						w[i].Tamb = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%lld].Tamb =%e\n", i, w[i].Tamb);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "heat_transfer_coefficient_vs_emissivity");

					// Stefan Bolcman
					// termostability wall
					w[i].emissivity = 0.0;
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
						w[i].Tamb = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%lld].Tamb =%e\n", i, w[i].Tamb);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "wall");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "heat_transfer_coefficient_vs_emissivity");

					// Stefan Bolcman
					// adiabatic wall
					w[i].emissivity = 0.0;
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
						w[i].Tamb = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%lld].Tamb =%e\n", i, w[i].Tamb);
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
						w[i].film_coefficient = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].film_coefficient = 0.0; // .
						printf(" w[%lld].film_coefficient=%e\n",i, w[i].film_coefficient);
						if (bSTOP_Reading) system("pause");
					}

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
						w[i].Tamb = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].Tamb = 30.0; // .
						printf("w[%lld].Tamb =%e\n", i, w[i].Tamb);
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
						w[i].emissivity = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						w[i].emissivity = 0.0; // .
						printf("w[%lld].emissivity =%e\n",i, w[i].emissivity);
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
					din = (integer)(idin);
					if (din == 1) w[i].bsymmetry = true; else w[i].bsymmetry = false;
					//printf("w[%lld].bsymmetry =%lld\n",i, din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].bsymmetry = false; // .
					printf(" w[%lld].bsymmetry = false\n", i);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bpressure");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);
					if (din == 1) w[i].bpressure = true; else w[i].bpressure = false;
					//printf(" w[%lld].bpressure=%lld\n", i,din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].bpressure = false; // .
					printf("w[%lld].bpressure = false\n",i );
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "bopening");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);
					if (din == 1) w[i].bopening = true; else w[i].bopening = false;
					//printf("  w[i].bopening =%lld\n",i, din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].bopening = false; // .
					printf(" w[%lld].bopening = false\n",i );
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Vx");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].Vx = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].Vx = 0.0; // .
					printf(" w[%lld].Vx =%e\n", i, w[i].Vx);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Vy");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].Vy = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].Vy = 0.0; // .
					printf("w[%lld].Vy =%e\n", i, w[i].Vy);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "Vz");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].Vz = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].Vz = 0.0; // .
					printf("w[%lld].Vz =%e\n", i, w[i].Vz);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "P");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].P = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].P  = 0.0; // .
					printf("w[%lld].P =%e\n",i, w[i].P);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "ithermal_stress_boundary_condition");
				
				if (imakesource(name0, idin)) {
					// Найдено успешно.
					din = (integer)(idin);

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
					//printf("w[%lld].ithermal_Stress_boundary_condition =%lld\n",i,din );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].ithermal_Stress_boundary_condition = 0; // 0 - FREE.
					printf(" w[%lld].ithermal_Stress_boundary_condition=FREE BOUNDARY\n",i );
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xForce");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].xForce = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].xForce = 0.0; // .
					printf("w[%lld].xForce =%e\n",i, w[i].xForce);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yForce");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].yForce = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].yForce = 0.0; // .
					printf("w[%lld].yForce=%e\n", i, w[i].yForce);
					if (bSTOP_Reading) system("pause");
				}
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zForce");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].zForce = (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].zForce = 0.0; // .
					printf(" w[%lld].zForce=%e\n",i, w[i].zForce);
					if (bSTOP_Reading) system("pause");
				}				
				//printf("Force Fx=%e Fy=%e Fz=%e\n", w[i].xForce, w[i].yForce, w[i].zForce);
				//system("PAUSE");

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iPlane");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					w[i].iPlane = (integer)(idin);
					//printf("w[%lld].iPlane =%lld\n", i, w[i].iPlane);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].iPlane = 0; // .
					printf(" w[%lld].iPlane=%lld\n", i, w[i].iPlane);
					if (bSTOP_Reading) system("pause");
				}
				// геометрия
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.xS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.xS = scale * 0.0; // .
					printf("w[%lld].g.xS =%e\n", i, w[i].g.xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.yS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.yS = scale * 0.0; // .
					printf(" w[%lld].g.yS=%e\n", i, w[i].g.yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.zS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.zS = scale * 0.0; // .
					printf(" w[%lld].g.zS  =%e\n", i, w[i].g.zS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.xE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.xE = scale * 0.0; // .
					printf("w[%lld].g.xE  =%e\n", i, w[i].g.xE);
					if (bSTOP_Reading) system("pause");
				}
				
				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.yE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.yE = scale * 0.0; // .
					printf("w[%lld].g.yE =%e\n", i, w[i].g.yE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "wall");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");
				
				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					w[i].g.zE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					w[i].g.zE = scale * 0.0; // .
					printf("w[%lld].g.zE =%e\n", i, w[i].g.zE);
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

			if (imakesource("iu_number", idin)) {
				// Найдено успешно.
				lu = (integer)(idin);
				//printf(" =%lld\n", );
			}
			else {
				printf("WARNING!!! iu_number not found in file premeshin.txt\n");
				lu = 0; // .
				printf(" lu=%lld\n",lu );
				if (bSTOP_Reading) system("pause");
			}
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

				char name0[1000] = "assembles";
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].xS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].xS = scale * 0.0; // .
					printf(" my_union[%lld].xS=%e\n", i, my_union[i].xS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "xE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].xE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].xE = scale * 0.0; // .
					printf(" my_union[%lld].xE=%e\n",i, my_union[i].xE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].yS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].yS = scale * 0.0; // .
					printf(" my_union[%lld].yS=%e\n", i, my_union[i].yS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "yE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].yE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].yE = scale * 0.0; // .
					printf("my_union[%lld].yE =%e\n", i, my_union[i].yE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zS");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].zS = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].zS = scale * 0.0; // .
					printf(" my_union[%lld].zS=%e\n",i, my_union[i].zS);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "zE");

				if (fmakesource(name0, fin)) {
					// Найдено успешно.
					my_union[i].zE = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].zE = scale * 0.0; // .
					printf(" my_union[%lld].zE=%e\n", i, my_union[i].zE);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "identifire");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].id  = (integer)(idin);
					//printf(" my_union[%lld].id=%lld\n", i, my_union[i].id);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].id = i+1; // Уникальный идентификатор АССЕМБЛЕСА.
					printf(" my_union[%lld].id=%lld\n",i, my_union[i].id);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "inx");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].inx = (integer)(idin);
					//printf(" my_union[%lld].inx=%lld\n", i, my_union[i].inx);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].inx = 23; // .
					printf(" my_union[%lld].inx=%lld\n", i, my_union[i].inx);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "iny");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].iny = (integer)(idin);
					//printf(" my_union[%lld].iny=%lld\n", i, my_union[i].iny);
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].iny = 23; // .
					printf("my_union[%lld].iny =%lld\n", i, my_union[i].iny);
					if (bSTOP_Reading) system("pause");
				}

				name0[0] = '\0'; strcat_s(name0, "assembles");
				buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
				strcat_s(name0, "inz");

				if (imakesource(name0, idin)) {
					// Найдено успешно.
					my_union[i].inz = (integer)(idin);
					//printf(" my_union[%lld].inz =%lld\n", i, my_union[i].inz );
				}
				else {
					printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
					my_union[i].inz = 23; // .
					printf(" my_union[%lld].inz =%lld\n", i, my_union[i].inz);
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
					eqin.itemper = (integer)(idin);
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
				eqin.imaxflD = (integer)(idin);
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
					char name0[1000] = "egddata";
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "xc");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].xc = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].xc = scale * 0.0; // .
						printf(" eqin.fluidinfo[%lld].xc=%e\n",i, eqin.fluidinfo[i].xc);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "yc");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].yc = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].yc = scale * 0.0; // .
						printf(" eqin.fluidinfo[%lld].yc=%e\n", i, eqin.fluidinfo[i].yc);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "zc");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].zc = scale * (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].zc = scale * 0.0; // .
						printf(" eqin.fluidinfo[%lld].zc=%e\n", i, eqin.fluidinfo[i].zc);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "iflow");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].iflow = (integer)(idin);
						//printf("eqin.fluidinfo[%lld].iflow =%lld\n", i, eqin.fluidinfo[i].iflow);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].iflow = 0; // .
						printf(" eqin.fluidinfo[%lld].iflow=%lld\n", i, eqin.fluidinfo[i].iflow);
						if (bSTOP_Reading) system("pause");
					}

					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "iflowregime");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].iflowregime = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].iflowregime=%lld\n",i,eqin.fluidinfo[i].iflowregime );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].iflowregime = 0; // Laminar.
						printf(" eqin.fluidinfo[%lld].iflowregime=%lld\n", i, eqin.fluidinfo[i].iflowregime);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "iturbmodel");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].iturbmodel = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].iturbmodel=%lld\n", i, eqin.fluidinfo[i].iturbmodel);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].iturbmodel = 0; // Zero Equation Turbulence Model.
						printf(" eqin.fluidinfo[%lld].iturbmodel=%lld\n", i, eqin.fluidinfo[i].iturbmodel);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "SmagConst");

					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].Cs = (doublereal)(fin);// постоянная Смагоринского.
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].Cs = 0.151; // default 0.151.
						printf(" eqin.fluidinfo[%lld].Cs=%e\n", i, eqin.fluidinfo[i].Cs);
						if (bSTOP_Reading) system("pause");
					}
					
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "iDynamicStressGermano");
					
					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din = (integer)(idin);
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
						printf(" eqin.fluidinfo[%lld].bDynamic_Stress = false\n",i );
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "iLimitersCs");
					
					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din = (integer)(idin);
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
						printf(" eqin.fluidinfo[%lld].bLimiters_Cs = false\n",i );
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "minCs");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].rminCs = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rminCs = -1.0e20; // минимальное возможное значение постоянной Смагоринского.
						printf("eqin.fluidinfo[%lld].rminCs =%e\n", i, eqin.fluidinfo[i].rminCs);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "maxCs");
				
				
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].rmaxCs = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rmaxCs = 1.0e23; // максимальное возможное значение постоянной Смагоринского.
						printf(" eqin.fluidinfo[%lld].rmaxCs=%e\n",i, eqin.fluidinfo[i].rmaxCs);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "itypeFiltrGermano");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].itypeFiltrGermano = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].itypeFiltrGermano=%lld\n",i,eqin.fluidinfo[i].itypeFiltrGermano );
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].itypeFiltrGermano = 2; // тип фильтра в модели Германо 1991 года. 
						printf(" eqin.fluidinfo[%lld].itypeFiltrGermano=%lld\n", i, eqin.fluidinfo[i].itypeFiltrGermano);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "roughness");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].roughness = 1.0e-6*(doublereal)(fin); // шероховатость стенки в м.
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].roughness = 1.0e-6*10.0; // шероховатость стенки в м.
						printf(" eqin.fluidinfo[%lld].roughness=%e\n",i, eqin.fluidinfo[i].roughness);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "rRimult");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						// множитель корректирующий турбулентное число Ричардсона.
						eqin.fluidinfo[i].rRi_mult = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rRi_mult = 1.0; // множитель корректирующий турбулентное число Ричардсона.
						printf(" eqin.fluidinfo[%lld].rRi_mult=%e\n", i, eqin.fluidinfo[i].rRi_mult);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "rSelectiveAngle");
					
					if (fmakesource(name0, fin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].rSelectiveAngle = (doublereal)(fin);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].rSelectiveAngle = 15.0; // пороговое значение угла в модели Selective Smagorinsky.
						printf(" eqin.fluidinfo[%lld].rSelectiveAngle=%e\n", i, eqin.fluidinfo[i].rSelectiveAngle);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "ipowerroughness");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].ipowerroughness = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].ipowerroughness=%lld\n", i, eqin.fluidinfo[i].ipowerroughness);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].ipowerroughness = 2; // показатель степени в модели учёта шероховатости.
						printf(" eqin.fluidinfo[%lld].ipowerroughness=%lld\n", i, eqin.fluidinfo[i].ipowerroughness);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "itypefiltr");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].itypeSelectiveSmagorinsky_filtr=%lld\n", i, eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr = 2; // тип фильтра в модели Selective Smagorinsky.
						printf(" eqin.fluidinfo[%lld].itypeSelectiveSmagorinsky_filtr=%lld\n", i, eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "bfdelta");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].bfdelta = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].bfdelta=%lld\n", i, eqin.fluidinfo[i].bfdelta);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bfdelta = 1; // учёт неравномерности сетки.
						printf(" eqin.fluidinfo[%lld].bfdelta=%lld\n", i, eqin.fluidinfo[i].bfdelta);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "bSmagorinsky_Lilly");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].bSmagorinskyLilly = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].bSmagorinskyLilly=%lld\n", i, eqin.fluidinfo[i].bSmagorinskyLilly);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bSmagorinskyLilly = 1; // модель Смагоринского-Лиллу.
						printf(" eqin.fluidinfo[%lld].bSmagorinskyLilly=%lld\n", i, eqin.fluidinfo[i].bSmagorinskyLilly);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "bsurface_roughness");

					if (imakesource(name0, idin)) {
						// Найдено успешно.
						eqin.fluidinfo[i].bsurface_roughness = (integer)(idin);
						//printf(" eqin.fluidinfo[%lld].bsurface_roughness=%lld\n", i, eqin.fluidinfo[i].bsurface_roughness);
					}
					else {
						printf("WARNING!!! %s not found in file premeshin.txt\n", name0);
						eqin.fluidinfo[i].bsurface_roughness = 0; // учёт шероховатости стенки.
						printf(" eqin.fluidinfo[%lld].bsurface_roughness=%lld\n", i, eqin.fluidinfo[i].bsurface_roughness);
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "bSwirlamendment");
					
					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din= (integer)(idin);
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
						printf(" eqin.fluidinfo[%lld].bSwirlAmendment = true\n",i );
						if (bSTOP_Reading) system("pause");
					}
					name0[0] = '\0'; strcat_s(name0, "egddata");
					buffer[0] = '\0'; _itoa_s(i,buffer,10); strcat_s(name0, buffer);
					strcat_s(name0, "bSelectiveSmagorinsky");
					
					if (imakesource(name0, idin)) {
						// Найдено успешно.
						din= (integer)(idin);
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
						printf(" eqin.fluidinfo[%lld].bSelectiveSmagorinsky = false\n", i);
						if (bSTOP_Reading) system("pause");
					}

				}
				// Параметры преобразователя картинок для отчетов.
				// 5.01.2018
			
				if (fmakesource("report_min", fin)) {
					// Найдено успешно.
					pfpir.fminimum = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! report_min not found in file premeshin.txt\n");
					pfpir.fminimum = scale * 0.0; // .
					printf(" pfpir.fminimum =%e\n", pfpir.fminimum);
					if (bSTOP_Reading) system("pause");
				}
					
				if (fmakesource("report_max", fin)) {
					// Найдено успешно.
					pfpir.fmaximum = scale * (doublereal)(fin);
				}
				else {
					printf("WARNING!!! report_max not found in file premeshin.txt\n");
					pfpir.fmaximum = scale * 0.0; // .
					printf("pfpir.fmaximum =%e\n", pfpir.fmaximum);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("report_directional", idin)) {
					// Найдено успешно.
					pfpir.idir = (integer)(idin);
					//printf(" pfpir.idir=%lld\n", pfpir.idir);
				}
				else {
					printf("WARNING!!! report_directional not found in file premeshin.txt\n");
					pfpir.idir = 0; // XY.
					printf(" pfpir.idir=%lld\n", pfpir.idir);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("amg1r6_checker", idin)) {
					// Найдено успешно.
					AMG1R6_LABEL = (integer)(idin);
					//printf(" AMG1R6_LABEL=%lld\n",AMG1R6_LABEL );
				}
				else {
					printf("WARNING!!! amg1r6_checker not found in file premeshin.txt\n");
					AMG1R6_LABEL = 0; // default amg1r5.
					printf(" AMG1R6_LABEL=%lld\n", AMG1R6_LABEL);
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
					printf(" number_processors=%lld\n", number_processors_global_var);
					if (bSTOP_Reading) system("pause");
				}

				if (imakesource("number_iterations_SIMPLE_algorithm", idin)) {
					// Найдено успешно.
					number_iteration_SIMPLE_algorithm  = (int)(idin);
					//printf(" number_iterations_SIMPLE_algorithm=%lld\n",number_iteration_SIMPLE_algorithm );
				}
				else {
					printf("WARNING!!! number_iterations_SIMPLE_algorithm not found in file premeshin.txt\n");
					number_iteration_SIMPLE_algorithm = 0; // default 1000 iterations.
					printf(" number_iterations_SIMPLE_algorithm=%lld\n", number_iteration_SIMPLE_algorithm);
					if (bSTOP_Reading) system("pause");
				}
				
				if (imakesource("stabilization_amg1r5_algorithm", idin)) {
					// Найдено успешно.
					stabilization_amg1r5_algorithm = (integer)(idin);
					//printf(" stabilization_amg1r5_algorithm=%lld\n",stabilization_amg1r5_algorithm );
				}
				else {
					printf("WARNING!!! stabilization_amg1r5_algorithm not found in file premeshin.txt\n");
					stabilization_amg1r5_algorithm = 1; // default BiCGStab+amg1r5.
					printf(" stabilization_amg1r5_algorithm=%lld\n", stabilization_amg1r5_algorithm);
					if (bSTOP_Reading) system("pause");
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
					if (b[i_1].n_Sc > 0) {
						doublereal pdiss = get_power(b[i_1].n_Sc, b[i_1].temp_Sc, b[i_1].arr_Sc, 20.0);
						doublereal vol = fabs(b[i_1].g.xE - b[i_1].g.xS)*fabs(b[i_1].g.yE - b[i_1].g.yS)*fabs(b[i_1].g.zE - b[i_1].g.zS);
						if (vol < 1.0e-40) {
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
				if (b[i_1].g.itypegeom == CYLINDER) {
					// Cylinder

					// На тот случай если геометрия состоит только из цилиндров.
					// 13.08.2019 Автомат настройки допусков сетки. 
					// Некое разумное уточнение принебрежимо малой длины для упрощения (shorter_length_for_simplification*).
					// Она не может быть больше чем 10% от характерной длины объекта заданной пользователем.
					switch (b[i_1].g.iPlane) {
					case XY:
						if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.Hcyl;
						if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
						if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
						if (b[i_1].g.R_in_cyl > 1.0e-40) {
							// Если внутренний радиус существует (задавался пользователем).
							if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
							if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
						}
						break;
					case XZ:
						if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.Hcyl;
						if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_out_cyl;
						if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
						if (b[i_1].g.R_in_cyl > 1.0e-40) {
							// Если внутренний радиус существует (задавался пользователем).
							if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.R_in_cyl;
							if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
						}
						break;
					case YZ:
						if (0.1*b[i_1].g.Hcyl < shorter_length_for_simplificationX_BASIC) shorter_length_for_simplificationX_BASIC = dmult * b[i_1].g.Hcyl;
						if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_out_cyl;
						if (0.1*b[i_1].g.R_out_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_out_cyl;
						if (b[i_1].g.R_in_cyl > 1.0e-40) {
							// Если внутренний радиус существует (задавался пользователем).
							if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationY_BASIC) shorter_length_for_simplificationY_BASIC = dmult * b[i_1].g.R_in_cyl;
							if (0.1*b[i_1].g.R_in_cyl < shorter_length_for_simplificationZ_BASIC) shorter_length_for_simplificationZ_BASIC = dmult * b[i_1].g.R_in_cyl;
						}
						break;
					}


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

							dpower += pdiss * vol;
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
			printf("number of sources ls=%lld\n", ls);
			printf("Sources integral power = %e W\n", dsoupow);
			printf("Full total power = %e W\n", dpower + dsoupow);
			// Запоминаем полное тепловыделение в модели.
			d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = dpower + dsoupow;
			printf("number of blocks lb=%lld\n", lb);
			printf("PRISMS = %lld, CYLINDERS = %lld, POLYGONS = %lld\n", iprism, icyl, ipoly);
			printf("SOLID: %lld\n", isol);
			printf("HOLLOW: %lld\n", ihol);
			printf("FLUID: %lld\n", iflui);
			printf("number of walls lw=%lld\n", lw);
			printf("number of units lu=%lld\n", lu);

			fclose(fp); // закрытие файла
		}
	}

#endif
	printf("OK. \n");
} // premeshin_new


// считывание параметров из 
// входного файла premeshin.txt
void premeshin(const char *fname, integer &lmatmax, integer &lb, integer &ls, integer &lw, TPROP* &matlist, BLOCK* &b, SOURCE* &s, WALL* &w,
	doublereal &dgx, doublereal &dgy, doublereal &dgz, integer &inx, integer &iny, integer &inz, doublereal &operatingtemperature,
	integer &ltdp, TEMP_DEP_POWER* &gtdps, integer &lu, UNION* &my_union)
{

	FILE* fp = NULL;
	errno_t err1;

#ifdef MINGW_COMPILLER
	fp = fopen64(fname, "r");
	if (fp != NULL) {
		err1 = 0;
	}
	else {
		err1 = 1;
	}
#else
	err1 = fopen_s(&fp, fname, "r");
#endif

	if (err1 != 0) {
		printf("No input File premeshin.txt \n");
		system("pause");

	}
	else
	{
		if (fp != NULL) {

			fclose(fp);

			int idin = 0;

			if (imakesource("lb", idin)) {
				// Найдено успешно.
				lb = (integer)(idin);
				//printf("lb =%lld\n", lb);

				if (lb > 0) {
					// начало 21.08.2019 - окончание 23.08.2019
			        // Сдержит простейший парсер.
					premeshin_new(fname, lmatmax, lb, ls, lw, matlist, b, s, w, dgx, dgy, dgz, inx, iny, inz,
						operatingtemperature, ltdp, gtdps, lu, my_union);
				}
				else {
					printf("ERROR!!! your model incorrect...\n");
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
// Файл my_long_range_amg_interpolation.cpp
// содержит реализацию дальнобойной интерполяции расстояния 3.
// Её планируется использовать для PMIS огрубления применённого
// квадрату матрицы (соседи соседей). В теории это приводит к 
// дополнительному уменьшению операторной сложности. 
// Начало разработки 18.09.2020.
// Первый работоспособный вариант 20.09.2020

#pragma once
#ifndef MY_LONG_RANGE_AMG_INTERPOLATION_CPP 
#define MY_LONG_RANGE_AMG_INTERPOLATION_CPP 1

// Дальнобойная интерполяция на основе нтерполяционной процедуры №3.
// Протестированный вариант освобожденный от лишнего кода.
// За счёт освобождения от лишнего кода улучшена ясность кода.
// 12.05.2018 - основа положенная в реализацию интерполяции.
// 18.09.2020
template <typename doublerealT>
void my_long_range3_interpolation_procedure_number3A(integer& the_number_of_neighbors_that_are_not_C_nodes,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative, bool bsuffix_work)
{

	// ослабление добавления новых С узлов 01.10.2020
	bool weakening_the_addition_of_new_C_nodes = true;
	if (!b_on_adaptive_local_refinement_mesh) {
		// Данная модификация предназначена только для АЛИС сетки и
		// на структурированной сетке не работает.
		weakening_the_addition_of_new_C_nodes = false;
	}

	//printf("incomming long range%e\n", free_debug_parametr1);
	//getchar();

	// do_not_let_one_neighbor = false; допустим один С сосед а также полное отсутствие С соседей.
			// do_not_let_one_neighbor =true; добавляется новый С узел и делается рестарт.
	bool do_not_let_one_neighbor = true;// bsuffix_work; // true; // Новые С узлы обязательно нужно добавлять иначе сходимость деградирует.

	integer icount1_memo = icount1;

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // оптимальный выбор.
	doublerealT theta_strong_F = theta83; // 3 июля 2016
	

	// четвёртая версия интерполяции.

	// Здесь узел F не имеющий Strong С соседей сам становится С узлом.
	// Узел F имеющий одного Strong  С соседа обрабатывается с помощью сильных С соседей 
	// сильных F узлов.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	//integer icount_93 = 0;

	// the begining

	bool byes_add = false;
	if (do_not_let_one_neighbor) {
		// Быстрое добавление недостающих С узлов.
		if (1) {
			// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
			//doublerealTT magic = 0.4; // 0.4 optimum


			the_number_of_neighbors_that_are_not_C_nodes = 0;
			number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]) {


				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer icsosF = 0;

				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увеличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;

				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;

				if (bpositive_connections) {


#if doubleintprecision == 1
					//printf("i8=%lld n=%lld\n", i8, n_a[ilevel - 1]);
#else
					//printf("i8=%d n=%d\n", i8, n_a[ilevel - 1]);
#endif

				//getchar();


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполяции а значит и оператора Галёркина.
					doublerealT maxelem_threshold = -1.0;

					//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
								//if (this_is_C_node[Amat.j[is0]] ) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold) {
									maxelem_threshold = fabs(Amat.aij[is0]);
								}
								//}
							}
						}
					}
					else {
						maxelem_threshold = threshold_quick_all[Amat.i[ii1]];
					}
					// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.


					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]] ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
										SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// Подсчитываем количество соседей которые не являются С узлами.
								the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
							}
						}
					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
						// Следующий случай назван "безнадёжным".
						// Безнадёжный случай оставляется без изменений.
						// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
						// интерполяцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
					}
				}
				else {
					// only negative connections


					// Нам нужен разреженный оператор Галёркина.
					// 5 декабря 2015 года мы попробуем увеличить разреженность
					// оператора интерполяции а значит и оператора Галёркина.
					doublerealT maxelem_threshold = -1.0;

					//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
								//if (this_is_C_node[Amat.j[is0]] ) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold)) {
									maxelem_threshold = fabs(Amat.aij[is0]);
								}
								//}
							}
						}
					}
					else {
						maxelem_threshold = threshold_quick_only_negative[Amat.i[ii1]];
					}
					// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.


					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (this_is_C_node[Amat.j[is0]] ) {
								//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
									sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]] ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
										SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// Подсчитываем количество соседей которые не являются С узлами.
								the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
							}
						}
					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
						// Следующий случай назван "безнадёжным".
						// Безнадёжный случай оставляется без изменений.
						// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
						// интерполяцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
					}

				}

				/*
				9.05.2018
				1. Определили thresholdы.
				2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
				3. icsos - количество strong C. icsosF - количество strong F.
				*/

				// 1 января 2016 Один сосед это недостаточно.
				// Поэтому в случае одного соседа делаем такой узел С узлом.

				bool baddC = true;
				for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
					if (this_is_C_node[Amat.j[is0]]) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (weakening_the_addition_of_new_C_nodes) {
								// Если нет ни одного сильного С соседа
								// то для хорошей сходимости за минимальное
								// число итераций такой F узел надо делать С узлом.
								// Здесь в случае наличия у F узла без сильных С соседей
								// хотябы одного слабого С соседа (самого сильного из слабых)
								// то мы интерполируемся из этого
								// слабого С соседа и не добавляем новых С узлов.
								// Это действие ухудшает сходимость в три раза.
								baddC = false;
							}
						}
					}
				}

				if (baddC) {


					if (do_not_let_one_neighbor) {
						// Мы делаем из F узла C узел только в том случае если у нас нету С связей.
						//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (this_is_F_node[Amat.j[is0]]) {
								if (Amat.j[is0] != Amat.i[ii1]) {


									// 20 jan 2016.
									// Можно не добавлять узел в случае если есть не менее трёх сильных F соседей.


									if (fabs(sumP) < RealZERO) {
										// Это случай когда вообще нет сильных С соседей.

										// if (icsosF==0) полная потеря сходимости.
										//if (icsosF<2)// наблюдается замедление скорости сходимости и времени расчёта.
										{
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//printf("i8 is Dirichlet node\n");

											// В статье написано что операторную сложность PMIS нельзя уменьшить меньше
											// чем у меня в программе. Уменьшение операторной сложности или пропуск 
											// узлов в данном месте кода приводит к потери скорости сходимости.
											// Скорость сходимости деградирует очень быстро.
											//if ((icount_93 % ((integer)(free_debug_parametr1)) == 0))
											{
												// Добавление узлов нужно делать обязательно, иначе сходимость 
												// деградирует. 21.09.2020
												if (this_is_C_node[i8] == false) iadditionalCstatistic++;
												this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
												this_is_C_node[i8] = true;
												bweSholdbeContinue = true;
												byes_add = true; // Были добавления узлов.
												//exit(1);
												// здесь нужна непрямая интерполяция.
											}
											//icount_93++;
										}
									}


								}
							}
						}//end

					}
				}
			}

		}

	}



	/*
	9.05.2018
	1. Определили thresholdы.
	2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
	3. icsos - количество strong C. icsosF - количество strong F.
	4. Если SumP ==0 то узел становится С узлом. iadditionalCstatistic++;
	*/

	// Если do_not_let_one_neighbor =true, то byes_add==false.

	if (1&&(!byes_add)) {

		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 июля 2016 для конвективных задач
		// линейный прогноз не дает улучшений
		// ни по быстродействию ни по числу V циклов.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82; // 0.4 is recomended.

		doublerealT magic_FF = magic82; // 18.09.2020
		//magic_FF =  0.45; // 20.09.2020 //free_debug_parametr1;//

		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		bool* permission_to_use_F = new bool[n_a[ilevel - 1] + 1];
		for (integer i8 = 0; i8 <= n_a[ilevel - 1]; i8++) {
			permission_to_use_F[i8] = true; // инициализация, всё разрешено по умолчанию.
		}
		bool* bactiveF = new bool[n_a[ilevel - 1] + 1];
		for (integer i8 = 0; i8 <= n_a[ilevel - 1]; i8++) {
			bactiveF[i8] = true; // инициализация, всё разрешено по умолчанию.
		}

		if (bpositive_connections) {

			// Построение пролонгации для узлов которые составляют F-nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] ) {

				permission_to_use_F[i8] = false;

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполяции а значит и оператора Галёркина.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] ) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold) {
								maxelem_threshold = fabs(Amat.aij[is0]);
							}
							//}
						}
					}
				}
				else {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
					}
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увеличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


									        // Сумма SumPall формируется только если у зла F есть сильные С соседи через magic.

											// Рассматриваем Strong F связи.

											// Смысл:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											
											
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
													// Важно что threshold должен быть определен только по С узлам.
													// Иначе потеря сходимости.
													if (this_is_C_node[Amat.j[is0_loc]]) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
														}
													}
												}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;
											integer icsos_locall = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublerealT sumP_loc = 0.0;
											doublerealT sumP_locall = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat.j[is0_loc]] ) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else {

														//if (Amat.j[is0_loc] != istr_etalon) {
														// Подсчитываем количество соседей которые не являются С узлами.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
														//}

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_locall++;
														}

													}

												}
											}

											//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											sumP_locall += sumP_loc;

											if (fabs(sumP_loc) > 1.0e-33) {
												// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.
												bactiveF[Amat.j[is0]] = true;

												// Если эту строку закоментировать то сходимость
												// резко ухудшится. 247 итераций вместо 20 итераций.
												permission_to_use_F[Amat.j[is0]] = false;

												SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
												icount_StronglyF++;
												icsosF++;
											}
											else {
												bactiveF[Amat.j[is0]] = false;
											}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}


				if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
					bool baddC = true;
					doublerealT Cmax = -1.0; // Наиболее сильная С связь из слабых С связей.
					for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
						if (this_is_C_node[Amat.j[is0]]) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								baddC = false;
								if (Amat.abs_aij[is0] > Cmax) {
									Cmax = Amat.abs_aij[is0];
								}
							}
						}
					}
					if (!baddC) {
						maxelem_threshold_theta = 0.98*Cmax;
						for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
							if (this_is_C_node[Amat.j[is0]]) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									doublerealT aij_abs = Amat.abs_aij[is0];
									if (aij_abs > 0.98*Cmax) {
										sumP += aij_abs; // сумма модулей внедиагональных элементов которые принадлежат Strongly C узлам.
										icsos++;
									}
								}
							}
						}
					}

				}



				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
					// Следующий случай назван "безнадёжным".
					// Безнадёжный случай оставляется без изменений.
					// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
					// интерполяцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом
					// и к тому-же не имеющих сильных F соседов.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть по крайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// здесь нужна непрямая интерполяция.
										}
										else {
											// Это как раз единственная Strong C связь. 
											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = fabs(Amat.aij[is0]) / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n, P);
												}
											}

										}
									}

								}
								else
									if ((bactiveF[Amat.j[is0]])&&(this_is_F_node[Amat.j[is0]])) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
											// Рассматриваем Strong F связи.

											// Смысл:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											if (fabs(SumPall) < 1.0e-33) {
												printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
												//getchar();
												system("PAUSE");
												exit(1);
											}
											doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
													// Важно что threshold должен быть определен только по С узлам.
													// Иначе потеря сходимости.
													if (this_is_C_node[Amat.j[is0_loc]] ) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
														}
													}
												}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;
											integer icsos_locall = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublerealT sumP_loc = 0.0;
											doublerealT sumP_locall = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat.j[is0_loc]] ) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else if ((permission_to_use_F[Amat.j[is0_loc]])&&(this_is_F_node[Amat.j[is0_loc]] )) {

														//if (Amat.j[is0_loc] != istr_etalon) {
														// Подсчитываем количество соседей которые не являются С узлами.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
														//}

														if (Amat.j[is0_loc] != istr_etalon) {


															
																// F node

																integer iFpoint_loc = Amat.j[is0_loc];
																
																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// сканируем всех соседей iFpointeger_loc 
																// чтобы найти среди них С узлы.

																// Идентифицируем локальный порог.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2) {
																		// Важно что threshold должен быть определен только по С узлам.
																		// Иначе потеря сходимости.
																		if (this_is_C_node[Amat.j[is0_loc2]]  ) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {
																				maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																			}
																		}
																	}
																}

																

																if (maxelem_threshold_loc2  < 0) {
																	//printf("Error !!! not found C nodes");
																	//system("PAUSE");
																	//exit(1);
																}
																else {
																	// Узел F без С узлов значит он не учитывается.

																	sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																	icsos_locall++;
																}
														}

													}

												}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											sumP_locall += sumP_loc;

											if (fabs(sumP_loc) > 1.0e-33) {
												// Учитывается только в том случае если есть сильные С соседи.

												// В общем это двойная интерполяционная сумма 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


													// Внедиагональный элемент из множества С узлов.

													// Данная вставка должна существенно сохранять 
													// разреженность оператора Галёркина на глубоких 
													// сеточных уровнях.
													// Модификация 5 декабря 2015.
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {
															if (Amat.j[is0_loc] != istr_etalon) {

																P[icount1].j = i8;
																P[icount1].i = C_numerate[Amat.j[is0_loc]];
																//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																if (fabs(sumP_loc) < 1.0e-33) {
																	printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																P[icount1].aij = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;
																icount1++;
																if (icount1 >= nsizePR * n) {
																	printf("memory error!!!\n");
																	printf("not enough memory for the interpolation operator.\n");
																	//system("PAUSE");
																	//exit(1);
																	deallocate_prolongation(nsizePR, n, P);
																}

															}
														}
														else if ((permission_to_use_F[Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]]  )) {
															// F node

															integer iFpoint_loc = Amat.j[is0_loc];
															if (fabs(sumP_locall) < 1.0e-33) {
																printf("error 2.a ! division by zero. SumP_locall =%e\n", sumP_locall);
																//getchar();
																system("PAUSE");
																exit(1);
															}
															doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

															// сканируем всех соседей iFpointeger_loc 
															// чтобы найти среди них С узлы.

															// Идентифицируем локальный порог.
															doublerealT maxelem_threshold_loc2 = -1.0;
															//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
															integer ii1_loc2 = row_startA[iFpoint_loc];
															integer istr_etalon2 = Amat.i[ii1_loc2];															
															integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
															for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2) {
																	// Важно что threshold должен быть определен только по С узлам.
																	// Иначе потеря сходимости.
																	if (this_is_C_node[Amat.j[is0_loc2]]  ) {
																		if (Amat.j[is0_loc2] != istr_etalon2) {
																			maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																		}
																	}
																}
															}

															doublerealT maxelem_threshold_loc_magic2 = maxelem_threshold_loc2 * magic_FF;

															if (maxelem_threshold_loc2 < 0) {
																//printf("Error !!! not found C nodes");
																//system("PAUSE");
																//exit(1);
															}
															else {
																// Здесь maxelem_threshold_loc2 это модуль максимального внедиагонального 
																// элемента в строке среди С соседей локально.

																// Найти соседей данного F-node которые C-node.
																integer icsos_loc2 = 0;

																// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
																// Для каждого такого члена суммы увеличиваем счётчик iscos_loc2. 
																// По идее iscos_loc2 должно быть 2 и более.
																doublerealT sumP_loc2 = 0.0;
																//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																	// Будем отсекать сразу наибольшее число вариантов.
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2) {

																		if (this_is_C_node[Amat.j[is0_loc2]]  ) {

																			if (Amat.j[is0_loc2] != istr_etalon2) {

																				sumP_loc2 += fabs(Amat.aij[is0_loc2]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																				icsos_loc2++;
																			}

																		}
																		else {

																			//if (Amat.j[is0_loc2] != istr_etalon2) {

																		}

																	}
																}

																doublerealT maxelem_threshold_loc_magic_minus2 = -maxelem_threshold_loc_magic2;

																// В общем это двойная интерполяционная сумма 
																//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																	// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


																	// Внедиагональный элемент из множества С узлов.

																	// Данная вставка должна существенно сохранять 
																	// разреженность оператора Галёркина на глубоких 
																	// сеточных уровнях.
																	// Модификация 5 декабря 2015.
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2) {

																		if (this_is_C_node[Amat.j[is0_loc2]]  ) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {

																				P[icount1].j = i8;
																				P[icount1].i = C_numerate[Amat.j[is0_loc2]];
																				//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																				if (fabs(sumP_loc2) < 1.0e-33) {
																					printf("error 3.a ! division by zero. sumP_loc2 =%e\n", sumP_loc2);
																					//getchar();
																					system("PAUSE");
																					exit(1);
																				}
																				P[icount1].aij = multiplyer_nu_loc * fabs(Amat.aij[is0_loc2]) / sumP_loc2;
																				icount1++;
																				if (icount1 >= nsizePR * n) {
																					printf("memory error!!!\n");
																					printf("not enough memory for the interpolation operator.\n");
																					//system("PAUSE");
																					//exit(1);
																					deallocate_prolongation(nsizePR, n, P);
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
									}
							}

						}
					}

				}

				// Снимаем лок.
				permission_to_use_F[i8] = true;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						{
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


									// Сумма SumPall формируется только если у зла F есть сильные С соседи через magic.

									// Рассматриваем Strong F связи.

									// Смысл:
									// 



									//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
									// Это лишнее условие, ранее уже было проверено что мы имеем дело
									// с сильной F связью.
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

									integer iFpoint = Amat.j[is0];


									// сканируем всех соседей iFpointeger 
									// чтобы найти среди них С узлы.

									// Идентифицируем локальный порог.
									doublerealT maxelem_threshold_loc = -1.0;
									//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
									integer ii1_loc = row_startA[iFpoint];
									integer istr_etalon = Amat.i[ii1_loc];
									integer iend_for = -1;
									integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
									for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
										iend_for = is0_loc;
										if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
											// Важно что threshold должен быть определен только по С узлам.
											// Иначе потеря сходимости.
											if (this_is_C_node[Amat.j[is0_loc]]  ) {
												if (Amat.j[is0_loc] != istr_etalon) {
													maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
												}
											}
										}
									}

									doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
									// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

									// Найти соседей данного F-node которые C-node.
									integer icsos_loc = 0;
									integer icsos_locall = 0;

									// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
									// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
									doublerealT sumP_loc = 0.0;
									doublerealT sumP_locall = 0.0;
									//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
									for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

										// Будем отсекать сразу наибольшее число вариантов.
										if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

											if (this_is_C_node[Amat.j[is0_loc]]  ) {

												if (Amat.j[is0_loc] != istr_etalon) {

													sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
													icsos_loc++;
												}

											}
											else {

												//if (Amat.j[is0_loc] != istr_etalon) {
												// Подсчитываем количество соседей которые не являются С узлами.
												//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
												//}

												if (Amat.j[is0_loc] != istr_etalon) {

													sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
													icsos_locall++;
												}

											}

										}
									}

									//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

									sumP_locall += sumP_loc;
									if (fabs(sumP_loc) > 1.0e-33) {
										// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.
										permission_to_use_F[Amat.j[is0]] = true;
									}
									bactiveF[Amat.j[is0]] = true;// инициализация.
								}
							}
							 
						}
					}
				}

			}
			////
		}
		else {
			// only negative connections.

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				permission_to_use_F[i8] = false;

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполяции а значит и оператора Галёркина.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]]  ) {
							if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold)) {
								maxelem_threshold = fabs(Amat.aij[is0]);
							}
							//}
						}
					}
				}
				else {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
					}
					maxelem_threshold = threshold_quick_only_negative[istr_etalon1];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увеличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {


									// Рассматриваем Strong F связи.

											// Смысл:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

									integer iFpoint = Amat.j[is0];


									// сканируем всех соседей iFpointeger 
									// чтобы найти среди них С узлы.

									// Идентифицируем локальный порог.
									doublerealT maxelem_threshold_loc = -1.0;
									//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
									integer ii1_loc = row_startA[iFpoint];
									integer istr_etalon = Amat.i[ii1_loc];
									integer iend_for = -1;
									integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
									for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
										iend_for = is0_loc;
										if ((Amat.aij[is0_loc]<0.0)&&(fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
											// Важно что threshold должен быть определен только по С узлам.
											// Иначе потеря сходимости.
											if (this_is_C_node[Amat.j[is0_loc]]  ) {
												if (Amat.j[is0_loc] != istr_etalon) {
													maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
												}
											}
										}
									}

									doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
									// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

									// Найти соседей данного F-node которые C-node.
									integer icsos_loc = 0;
									integer icsos_locall = 0;

									// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
									// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
									doublerealT sumP_loc = 0.0;
									doublerealT sumP_locall = 0.0;
									//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
									for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

										// Будем отсекать сразу наибольшее число вариантов.
										if ((Amat.aij[is0_loc]<0.0)&&(fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

											if (this_is_C_node[Amat.j[is0_loc]]  ) {

												if (Amat.j[is0_loc] != istr_etalon) {

													sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
													icsos_loc++;
												}

											}
											else {

												//if (Amat.j[is0_loc] != istr_etalon) {
												// Подсчитываем количество соседей которые не являются С узлами.
												//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
												//}

												if (Amat.j[is0_loc] != istr_etalon) {

													sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
													icsos_locall++;
												}

											}

										}
									}

									//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

									sumP_locall += sumP_loc;

									if (fabs(sumP_loc) > 1.0e-33) {
										// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.

										bactiveF[Amat.j[is0]] = true;

										permission_to_use_F[Amat.j[is0]] = false;

										SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
									else {
										bactiveF[Amat.j[is0]] = false;
									}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
					// Следующий случай назван "безнадёжным".
					// Безнадёжный случай оставляется без изменений.
					// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
					// интерполяцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть по крайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// здесь нужна непрямая интерполяция.
										}
										else {
											// Это как раз единственная Strong C связь. 
											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < 1.0e-23) {
													printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = fabs(Amat.aij[is0]) / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n, P);
												}
											}

										}
									}

								}
								else
									if ((bactiveF[Amat.j[is0]]) && (this_is_F_node[Amat.j[is0]]  )) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
											// Рассматриваем Strong F связи.

											// Смысл:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
														}
													}
												}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;
											integer icsos_locall = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublerealT sumP_loc = 0.0;
											doublerealT sumP_locall = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else if ((permission_to_use_F[Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]]  ))  {

														if (Amat.j[is0_loc] != istr_etalon) {
															
															// F node

															integer iFpoint_loc = Amat.j[is0_loc];

															doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

															// сканируем всех соседей iFpointeger_loc 
															// чтобы найти среди них С узлы.

															// Идентифицируем локальный порог.
															doublerealT maxelem_threshold_loc2 = -1.0;
															//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
															integer ii1_loc2 = row_startA[iFpoint_loc];
															integer istr_etalon2 = Amat.i[ii1_loc2];															
															integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
															for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																if ((Amat.aij[is0_loc2]<0.0)&&(fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2)) {
																	// Важно что threshold должен быть определен только по С узлам.
																	// Иначе потеря сходимости.
																	if (this_is_C_node[Amat.j[is0_loc2]]  ) {
																		if (Amat.j[is0_loc2] != istr_etalon2) {
																			maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																		}
																	}
																}
															}



															if (maxelem_threshold_loc2 < 0) {
																//printf("Error !!! not found C nodes");
																//system("PAUSE");
																//exit(1);
															}
															else {
																// Узел F без С узлов значит он не учитывается.

																sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																icsos_locall++;
															}
															
														}
													}

												}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											sumP_locall += sumP_loc;

											if (fabs(sumP_loc) > 1.0e-33) {

												// Учитывается только в том случае если есть сильные С соседи.


												// В общем это двойная интерполяционная сумма 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


													// Внедиагональный элемент из множества С узлов.

													// Данная вставка должна существенно сохранять 
													// разреженность оператора Галёркина на глубоких 
													// сеточных уровнях.
													// Модификация 5 декабря 2015.
													if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {
															if (Amat.j[is0_loc] != istr_etalon) {

																P[icount1].j = i8;
																P[icount1].i = C_numerate[Amat.j[is0_loc]];
																//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																if (fabs(sumP_loc) < 1.0e-23) {
																	printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																P[icount1].aij = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_loc;
																icount1++;
																if (icount1 >= nsizePR * n) {
																	printf("memory error!!!\n");
																	printf("not enough memory for the interpolation operator.\n");
																	//system("PAUSE");
																	//exit(1);
																	deallocate_prolongation(nsizePR, n, P);
																}

															}
														}
														else if ((permission_to_use_F[Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]]  )) {
															// F node

															integer iFpoint_loc = Amat.j[is0_loc];
															if (fabs(sumP_locall) < 1.0e-33) {
																printf("error 2.a ! division by zero. SumP_locall =%e\n", sumP_locall);
																//getchar();
																system("PAUSE");
																exit(1);
															}
															doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

															// сканируем всех соседей iFpointeger_loc 
															// чтобы найти среди них С узлы.

															// Идентифицируем локальный порог.
															doublerealT maxelem_threshold_loc2 = -1.0;
															//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
															integer ii1_loc2 = row_startA[iFpoint_loc];
															integer istr_etalon2 = Amat.i[ii1_loc2];
															integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
															for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2)) {
																	// Важно что threshold должен быть определен только по С узлам.
																	// Иначе потеря сходимости.
																	if (this_is_C_node[Amat.j[is0_loc2]]  ) {
																		if (Amat.j[is0_loc2] != istr_etalon2) {
																			maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																		}
																	}
																}
															}

															doublerealT maxelem_threshold_loc_magic2 = maxelem_threshold_loc2 * magic_FF;


															if (maxelem_threshold_loc2 < 0) {
																//printf("Error !!! not found C nodes");
																//system("PAUSE");
																//exit(1);
															}
															else {
																// Здесь maxelem_threshold_loc2 это модуль максимального внедиагонального 
																// элемента в строке среди С соседей локально.


															// Найти соседей данного F-node которые C-node.
																integer icsos_loc2 = 0;

																// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
																// Для каждого такого члена суммы увеличиваем счётчик iscos_loc2. 
																// По идее iscos_loc2 должно быть 2 и более.
																doublerealT sumP_loc2 = 0.0;
																//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																	// Будем отсекать сразу наибольшее число вариантов.
																	if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2)) {

																		if (this_is_C_node[Amat.j[is0_loc2]]  ) {

																			if (Amat.j[is0_loc2] != istr_etalon2) {

																				sumP_loc2 += fabs(Amat.aij[is0_loc2]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																				icsos_loc2++;
																			}

																		}
																		else {

																			//if (Amat.j[is0_loc2] != istr_etalon2) {

																		}

																	}
																}

																doublerealT maxelem_threshold_loc_magic_minus2 = -maxelem_threshold_loc_magic2;

																// В общем это двойная интерполяционная сумма 
																//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																	// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


																	// Внедиагональный элемент из множества С узлов.

																	// Данная вставка должна существенно сохранять 
																	// разреженность оператора Галёркина на глубоких 
																	// сеточных уровнях.
																	// Модификация 5 декабря 2015.
																	if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2)) {

																		if (this_is_C_node[Amat.j[is0_loc2]]  ) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {

																				P[icount1].j = i8;
																				P[icount1].i = C_numerate[Amat.j[is0_loc2]];
																				//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																				if (fabs(sumP_loc2) < 1.0e-33) {
																					printf("error 3.a ! division by zero. sumP_loc2 =%e\n", sumP_loc2);
																					//getchar();
																					system("PAUSE");
																					exit(1);
																				}
																				P[icount1].aij = multiplyer_nu_loc * fabs(Amat.aij[is0_loc2]) / sumP_loc2;
																				icount1++;
																				if (icount1 >= nsizePR * n) {
																					printf("memory error!!!\n");
																					printf("not enough memory for the interpolation operator.\n");
																					//system("PAUSE");
																					//exit(1);
																					deallocate_prolongation(nsizePR, n, P);
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
									}
							}

						}
					}

				}

				// Снимаем лок.
				permission_to_use_F[i8] = true;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						{
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


									// Сумма SumPall формируется только если у зла F есть сильные С соседи через magic.

									// Рассматриваем Strong F связи.

									// Смысл:
									// 



									//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
									// Это лишнее условие, ранее уже было проверено что мы имеем дело
									// с сильной F связью.
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

									integer iFpoint = Amat.j[is0];


									// сканируем всех соседей iFpointeger 
									// чтобы найти среди них С узлы.

									// Идентифицируем локальный порог.
									doublerealT maxelem_threshold_loc = -1.0;
									//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
									integer ii1_loc = row_startA[iFpoint];
									integer istr_etalon = Amat.i[ii1_loc];
									integer iend_for = -1;
									integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
									for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
										iend_for = is0_loc;
										if ((Amat.aij[is0_loc]<0.0)&&(fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
											// Важно что threshold должен быть определен только по С узлам.
											// Иначе потеря сходимости.
											if (this_is_C_node[Amat.j[is0_loc]]  ) {
												if (Amat.j[is0_loc] != istr_etalon) {
													maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
												}
											}
										}
									}

									doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
									// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

									// Найти соседей данного F-node которые C-node.
									integer icsos_loc = 0;
									integer icsos_locall = 0;

									// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
									// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
									doublerealT sumP_loc = 0.0;
									doublerealT sumP_locall = 0.0;
									//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
									for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

										// Будем отсекать сразу наибольшее число вариантов.
										if ((Amat.aij[is0_loc]<0.0)&&(fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

											if (this_is_C_node[Amat.j[is0_loc]]) {

												if (Amat.j[is0_loc] != istr_etalon) {

													sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
													icsos_loc++;
												}

											}
											else {

												//if (Amat.j[is0_loc] != istr_etalon) {
												// Подсчитываем количество соседей которые не являются С узлами.
												//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
												//}

												if (Amat.j[is0_loc] != istr_etalon) {

													sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
													icsos_locall++;
												}

											}

										}
									}

									//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

									sumP_locall += sumP_loc;
									if (fabs(sumP_loc) > 1.0e-33) {
										// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.
										permission_to_use_F[Amat.j[is0]] = true;
									}
									bactiveF[Amat.j[is0]] = true;// инициализация.
								}
							}

						}
					}
				}

			} // end only negative connections

		}

		delete[] permission_to_use_F;
		delete[] bactiveF;


		

	}


	/*
	if (0) {
		// Надо читать теорию.
		// 26.09.2020 Сглаженная аггрегация.

		//Ak1* P_copy = new Ak1[nsizePR * n_a[0] + 1];

		//icount1
		//for (integer i_1 = 1; i_1 < icount1; i_1++) {
			//P_copy[i_1] = P[i_1];// Копия.
		//}

		bool* bvisible = new bool[n_a[ilevel - 1] + 1];
		for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
			bvisible[i_1] = false;
		}

		doublerealT* xP = new doublerealT[n_a[ilevel - 1] + 1];
		doublerealT* xP2 = new doublerealT[n_a[ilevel - 1] + 1];
		for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
			xP[i_1] = 0.0;
			xP2[i_1] = 0.0;
		}

		integer k_1 = icount1_memo;

		//for (integer i_1 = icount1_memo; i_1 < icount1; i_1++)
		for (integer i_1 = 1; i_1 < n_a[ilevel - 1]; i_1++) {
			if (this_is_F_node[i_1])
			{
				//if (!bvisible[P_copy[i_1].i]) 
				{

					// Инициализация.
					//for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++) {
						//xP[P_copy[i_1].j] = P_copy[i_1].aij;
					//}


					// smoother
					doublerealT omega = 0.6667;

					// Одна итерация.
					// Однопоточный вариант программы.
					//for (integer ii = 1; ii <= n_a[ilevel-1]; ii++)
					for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
					{
						//integer ii = P_copy[i_1].i;
						for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
						{
							//integer	ii = P_copy[i_1].j;
							integer ii = Amat.j[s_1];

							integer istr = ii;
							doublerealT rold = xP[istr];

							// нулевая правая часть.

							doublerealT rsum = 0.0;
							doublerealT rsum_ap = 0.0;
							integer is1 = row_startA[ii];
							integer is2 = row_startA[ii + 1] - 1;
							doublerealT ap = 0.0;
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ii1++)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];

								if (ipos != Amat.i[ii1]) {
									if (this_is_C_node[ipos]) {
										if (Amat.abs_aij[ii1] > theta * threshold_quick_all[Amat.i[ii1]]) {
											rsum += Amat.abs_aij[ii1] * xP[ipos];
											rsum_ap += Amat.abs_aij[ii1];
										}
									}
								}
								//else {
									//ap = Amat.aij[ii1];
								//}
							}
							xP2[istr] += rsum;
							//xP2[istr] *= ap;// делим на диагональный элемент, т.к. в
							//диагонали хранится обратная величина для диагонального коэффициента.
							xP2[istr] /= rsum_ap;

							// Запускается только после первого раза сразу запускать нельзя
							// т.к. начальное значение не является адекватным.
							xP2[istr] = omega * xP2[istr] + (1.0 - omega) * rold; // this is SOR
						}
					}

					// find_max
					doublerealT xPmax = -1.0;
					doublerealT sumxP = 0.0;
					for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
					{
						for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
						{
							//integer	ii = P_copy[i_1].j;
							integer ii = Amat.j[s_1];

							{
								xP[ii] = xP2[ii];
								if (xP[ii] > xPmax) {
									xPmax = xP[ii];
								}
								sumxP += xP[ii];
							}
						}
					}
					// Просеивание.
					doublerealT sumxPnew = 0.0;
					for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
					{
						for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
						{
							//integer	ii = P_copy[i_1].j;
							integer ii = Amat.j[s_1];


							{
								if ((this_is_F_node[ii]) || (xP[ii] <= 0.2 * xPmax)) {
									xP[ii] = 0.0;
								}
								else {
									sumxPnew += xP[ii];
								}
							}
						}
					}
					if (sumxPnew > 1.0e-30) {

						for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
						{
							for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
							{
								//integer	ii = P_copy[i_1].j;
								integer ii = Amat.j[s_1];

								{
									xP[ii] *= sumxP / sumxPnew;
									if ((this_is_C_node[ii]) && (xP[ii] > 1.0e-30)) {
										P[k_1].aij = xP[ii];
										P[k_1].i = P_copy[i_1].i;
										P[k_1].j = ii;
										k_1++;
										if (k_1 >= nsizePR * n) {
											printf("memory error!!!\n");
											printf("not enough memory for the interpolation operator.\n");
											//system("PAUSE");
											//exit(1);
											deallocate_prolongation(nsizePR, n, P);
										}
									}
								}
							}
						}

					}
					else {
						printf("ERROR !!! zero string in matrix P\n");
						printf("string number=%lld\n", P_copy[i_1].i);
						system("PAUSE");
						//this_is_C_node[P_copy[i_1].i] = true;
						//this_is_F_node[P_copy[i_1].i] = false;
						//P[k_1].aij = 1.0;
						//P[k_1].i = P_copy[i_1].i;
						//P[k_1].j = P_copy[i_1].i;
						//k_1++;
					}

					bvisible[P_copy[i_1].i] = true;

					for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
					{
						for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
						{
							//integer	ii = P_copy[i_1].j;
							integer ii = Amat.j[s_1];
							{
								xP[ii] = 0.0;
								xP2[ii] = 0.0;
							}
						}
					}
				}
			}
		}
		delete[] P_copy;
		delete[] bvisible;
		delete[] xP;
		delete[] xP2;

		icount1 = k_1;

	}
	*/

	//system("pause");


} // my_interpolation_procedure_number3A


  // Дальнобойная интерполяция на основе нтерполяционной процедуры №3.
  // Протестированный вариант освобожденный от лишнего кода.
  // За счёт освобождения от лишнего кода улучшена ясность кода.
  // 12.05.2018 - основа положенная в реализацию интерполяции.
  // 18.09.2020. 
  // 12.10.2020 Начало написания параллельной версии интерполяции с расстоянием 3.
template <typename doublerealT>
void my_long_range3_interpolation_procedure_number3A_parallelx8(integer& the_number_of_neighbors_that_are_not_C_nodes_gl,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor_gl,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic_gl,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF_gl,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative, bool bsuffix_work)
{

	// ослабление добавления новых С узлов 01.10.2020
	bool weakening_the_addition_of_new_C_nodes = true;
	if (!b_on_adaptive_local_refinement_mesh) {
		// Данная модификация предназначена только для АЛИС сетки и
		// на структурированной сетке не работает.
		weakening_the_addition_of_new_C_nodes = false;
	}


	//bool bprint_mesage_diagnostic = true;
	//if (my_amg_manager.iprint_log == 0) {
		//bprint_mesage_diagnostic = false;
	//}


#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#endif

	//printf("incomming long range%e\n", free_debug_parametr1);
	//getchar();

	// do_not_let_one_neighbor = false; допустим один С сосед а также полное отсутствие С соседей.
	// do_not_let_one_neighbor =true; добавляется новый С узел и делается рестарт.
	bool do_not_let_one_neighbor = true;// bsuffix_work; // true; // Новые С узлы обязательно нужно добавлять иначе сходимость деградирует.

	integer icount1_memo = icount1;

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // оптимальный выбор.
	doublerealT theta_strong_F = theta83; // 3 июля 2016
	doublerealT theta = theta_gl;


										  // четвёртая версия интерполяции.

										  // Здесь узел F не имеющий Strong С соседей сам становится С узлом.
										  // Узел F имеющий одного Strong  С соседа обрабатывается с помощью сильных С соседей 
										  // сильных F узлов.

										  //6interpolation 0.4 6.77 11 26 28.355
										  //6interpolation 0.45 6.6 10 27 28.151
										  //6interpolation 0.5 6.42 12 32 28.735
										  //4interpolation 0.4 3.7  52 24.736 // best
										  //4interpolation 0.3 3.78 13 59 27.525
										  //4interpolation 0.5 3.61 12 55 25.533
										  //4interpolation 0.45 3.65 10 63 30.24

										  //integer icount_93 = 0;

										  // the begining

	integer iadditionalCstatistic = 0;
	integer number_of_F_nodes_with_one_single_strong_C_neighborF = 0;

	bool byes_add = false;
	if (do_not_let_one_neighbor) {
		// Быстрое добавление недостающих С узлов.
		if (1) {
			// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
			//doublerealTT magic = 0.4; // 0.4 optimum


			integer the_number_of_neighbors_that_are_not_C_nodes = 0;
			integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
#pragma omp parallel for reduction(+: number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes, iadditionalCstatistic, number_of_F_nodes_with_one_single_strong_C_neighborF)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++)
			{
				if (this_is_F_node[i8]) {


					// Найти соседей данного F-node которые C-node.
					integer icsos = 0;
					integer icsosF = 0;

					// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
					// Для каждого такого члена суммы увеличиваем счётчик iscos. По идее iscos должно быть 2 и более.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;

					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;

					if (bpositive_connections) {


#if doubleintprecision == 1
						//printf("i8=%lld n=%lld\n", i8, n_a[ilevel - 1]);
#else
						//printf("i8=%d n=%d\n", i8, n_a[ilevel - 1]);
#endif

					//getchar();


					// Нам нужен разреженный оператор Галёркина.
					// 5 декабря 2015 года мы попробуем увеличить разреженность
					// оператора интерполяции а значит и оператора Галёркина.
						doublerealT maxelem_threshold = -1.0;

						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

						if (!btreshold_on_new_vetv) {
							for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
									//if (this_is_C_node[Amat.j[is0]] ) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold) {
										maxelem_threshold = fabs(Amat.aij[is0]);
									}
									//}
								}
							}
						}
						else {
							maxelem_threshold = threshold_quick_all[Amat.i[ii1]];
						}
						// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.


						doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
						doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (this_is_C_node[Amat.j[is0]]) {
									//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
										sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
										icsos++;
									}
								}
								else {
									if (this_is_F_node[Amat.j[is0]]) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
											SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
											icsosF++;
										}
									}
									// Подсчитываем количество соседей которые не являются С узлами.
									the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
								}
							}
						}
						if (icsos == 1) {
							number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																				   // Следующий случай назван "безнадёжным".
																				   // Безнадёжный случай оставляется без изменений.
																				   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																				   // интерполяцию на этих уровнях можно улучшить.
							if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
						}
					}
					else {
						// only negative connections


						// Нам нужен разреженный оператор Галёркина.
						// 5 декабря 2015 года мы попробуем увеличить разреженность
						// оператора интерполяции а значит и оператора Галёркина.
						doublerealT maxelem_threshold = -1.0;

						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

						if (!btreshold_on_new_vetv) {
							for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
									//if (this_is_C_node[Amat.j[is0]] ) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold)) {
										maxelem_threshold = fabs(Amat.aij[is0]);
									}
									//}
								}
							}
						}
						else {
							maxelem_threshold = threshold_quick_only_negative[Amat.i[ii1]];
						}
						// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.


						doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
						doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (this_is_C_node[Amat.j[is0]]) {
									//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
										sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
										icsos++;
									}
								}
								else {
									if (this_is_F_node[Amat.j[is0]]) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
											SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
											icsosF++;
										}
									}
									// Подсчитываем количество соседей которые не являются С узлами.
									the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
								}
							}
						}
						if (icsos == 1) {
							number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																				   // Следующий случай назван "безнадёжным".
																				   // Безнадёжный случай оставляется без изменений.
																				   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																				   // интерполяцию на этих уровнях можно улучшить.
							if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
						}

					}

					/*
					9.05.2018
					1. Определили thresholdы.
					2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
					3. icsos - количество strong C. icsosF - количество strong F.
					*/

					// 1 января 2016 Один сосед это недостаточно.
					// Поэтому в случае одного соседа делаем такой узел С узлом.

					bool baddC = true;
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (this_is_C_node[Amat.j[is0]]) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (weakening_the_addition_of_new_C_nodes) {
									// Если нет ни одного сильного С соседа
									// то для хорошей сходимости за минимальное
									// число итераций такой F узел надо делать С узлом.
									// Здесь в случае наличия у F узла без сильных С соседей
									// хотябы одного слабого С соседа (самого сильного из слабых)
									// то мы интерполируемся из этого
									// слабого С соседа и не добавляем новых С узлов.
									// Это действие ухудшает сходимость в три раза.
									baddC = false;
								}
							}
						}
					}

					if (baddC) {


						if (do_not_let_one_neighbor) {
							// Мы делаем из F узла C узел только в том случае если у нас нету С связей.
							//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
							for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								if (this_is_F_node[Amat.j[is0]]) {
									if (Amat.j[is0] != Amat.i[ii1]) {


										// 20 jan 2016.
										// Можно не добавлять узел в случае если есть не менее трёх сильных F соседей.


										if (fabs(sumP) < RealZERO) {
											// Это случай когда вообще нет сильных С соседей.

											// if (icsosF==0) полная потеря сходимости.
											//if (icsosF<2)// наблюдается замедление скорости сходимости и времени расчёта.
											{
												//printf("error interpolation zero diagonal sumP.\n");
												//printf("Fnode all sosed is F");
												//system("pause");
												//printf("i8 is Dirichlet node\n");

												// В статье написано что операторную сложность PMIS нельзя уменьшить меньше
												// чем у меня в программе. Уменьшение операторной сложности или пропуск 
												// узлов в данном месте кода приводит к потери скорости сходимости.
												// Скорость сходимости деградирует очень быстро.
												//if ((icount_93 % ((integer)(free_debug_parametr1)) == 0))
												{
													// Добавление узлов нужно делать обязательно, иначе сходимость 
													// деградирует. 21.09.2020
													if (this_is_C_node[i8] == false) iadditionalCstatistic++;
													this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
													this_is_C_node[i8] = true;
													//bweSholdbeContinue = true;

#pragma omp critical
													{// Важно делать именно так, это даёт низкую Acomp
														if (!byes_add) {
															byes_add = true;  // Были добавления узлов.
														}
													}

													
																	 //exit(1);
																	 // здесь нужна непрямая интерполяция.
												}
												//icount_93++;
											}
										}


									}
								}
							}//end

						}
					}
				}
			
			}
		
			number_of_F_nodes_with_one_single_strong_C_neighbor_gl = number_of_F_nodes_with_one_single_strong_C_neighbor;
			the_number_of_neighbors_that_are_not_C_nodes_gl = the_number_of_neighbors_that_are_not_C_nodes;

		}

	}


	//if (iadditionalCstatistic > 0) byes_add = true;
	bweSholdbeContinue = bweSholdbeContinue || byes_add;

	iadditionalCstatistic_gl += iadditionalCstatistic;
	iadditionalCstatistic = 0;

	number_of_F_nodes_with_one_single_strong_C_neighborF_gl += number_of_F_nodes_with_one_single_strong_C_neighborF;
	number_of_F_nodes_with_one_single_strong_C_neighborF = 0;


	/*
	9.05.2018
	1. Определили thresholdы.
	2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
	3. icsos - количество strong C. icsosF - количество strong F.
	4. Если SumP ==0 то узел становится С узлом. iadditionalCstatistic++;
	*/

	// Если do_not_let_one_neighbor =true, то byes_add==false.

	if (1 && (!byes_add)) {

		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 июля 2016 для конвективных задач
		// линейный прогноз не дает улучшений
		// ни по быстродействию ни по числу V циклов.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82; // 0.4 is recomended.

		const doublerealT magic_FF = magic82; // 18.09.2020
										//magic_FF =  0.45; // 20.09.2020 //free_debug_parametr1;//

		integer LOC_P_LIM_SIZE = -1;
		double dsize = 2.2*0.125*nsizePR*n_a[0];
		LOC_P_LIM_SIZE = (integer)(dsize);

		Ak1** P_loc = new Ak1*[8];
		integer* isize_P_loc = new integer[8];
		for (int itid = 0; itid < 8; itid++) {
			isize_P_loc[itid] = 0;

			P_loc[itid] = new Ak1[LOC_P_LIM_SIZE];
		}


		bool** permission_to_use_F = new bool*[8];
		bool** bactiveF = new bool*[8];


		for (int itid = 0; itid < 8; itid++) {

			permission_to_use_F[itid] = new bool[n_a[ilevel - 1] + 1];
#pragma omp parallel for
			for (integer i8 = 0; i8 <= n_a[ilevel - 1]; i8++) {
				permission_to_use_F[itid][i8] = true; // инициализация, всё разрешено по умолчанию.
			}
			bactiveF[itid] = new bool[n_a[ilevel - 1] + 1];
#pragma omp parallel for
			for (integer i8 = 0; i8 <= n_a[ilevel - 1]; i8++) {
				bactiveF[itid][i8] = true; // инициализация, всё разрешено по умолчанию.
			}
		}

		integer the_number_of_neighbors_that_are_not_C_nodes = 0;
		integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// Построение пролонгации для узлов которые составляют F-nodes.
			// Каждый F-nodes окружён C-nodes.
#pragma omp parallel for reduction(+ : number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes,  number_of_F_nodes_with_one_single_strong_C_neighborF, iadditionalCstatistic) schedule (guided)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) {

#ifdef _OPENMP 
				int itid = omp_get_thread_num();
#else
				int itid = 0;
#endif		

				if (this_is_F_node[i8]) {

					permission_to_use_F[itid][i8] = false;

					// Нам нужен разреженный оператор Галёркина.
					// 5 декабря 2015 года мы попробуем увеличить разреженность
					// оператора интерполяции а значит и оператора Галёркина.
					doublerealT maxelem_threshold = -1.0;
					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer istr_etalon1 = Amat.i[ii1];
					integer iend_for1 = -1;
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
							iend_for1 = is0;
							if (Amat.j[is0] != istr_etalon1) {
								// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
								//if (this_is_C_node[Amat.j[is0]] ) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold) {
									maxelem_threshold = fabs(Amat.aij[is0]);
								}
								//}
							}
						}
					}
					else {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
							iend_for1 = is0;
						}
						maxelem_threshold = threshold_quick_all[istr_etalon1];
					}
					// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

					// Найти соседей данного F-node которые C-node.
					integer icsos = 0;
					integer icsosF = 0;

					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


					// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
					// Для каждого такого члена суммы увеличиваем счётчик iscos. По идее iscos должно быть 2 и более.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;

					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
					for (integer is0 = ii1; is0 <= iend_for1; is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


										// Сумма SumPall формируется только если у зла F есть сильные С соседи через magic.

										// Рассматриваем Strong F связи.

										// Смысл:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// Это лишнее условие, ранее уже было проверено что мы имеем дело
										// с сильной F связью.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// сканируем всех соседей iFpointeger 
										// чтобы найти среди них С узлы.

										// Идентифицируем локальный порог.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
											iend_for = is0_loc;
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
												// Важно что threshold должен быть определен только по С узлам.
												// Иначе потеря сходимости.
												if (this_is_C_node[Amat.j[is0_loc]]) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

										// Найти соседей данного F-node которые C-node.
										integer icsos_loc = 0;
										integer icsos_locall = 0;

										// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
										// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// Будем отсекать сразу наибольшее число вариантов.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// Подсчитываем количество соседей которые не являются С узлами.
													//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
													//}

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;

										if (fabs(sumP_loc) > 1.0e-33) {
											// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.
											bactiveF[itid][Amat.j[is0]] = true;

											// Если эту строку закоментировать то сходимость
											// резко ухудшится. 247 итераций вместо 20 итераций.
											permission_to_use_F[itid][Amat.j[is0]] = false;

											SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
											icsosF++;
										}
										else {
											bactiveF[itid][Amat.j[is0]] = false;
										}
									}
								}
								// Подсчитываем количество соседей которые не являются С узлами.
								the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
							}
						}
					}


					if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
						bool baddC = true;
						doublerealT Cmax = -1.0; // Наиболее сильная С связь из слабых С связей.
						for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
							if (this_is_C_node[Amat.j[is0]]) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									baddC = false;
									if (Amat.abs_aij[is0] > Cmax) {
										Cmax = Amat.abs_aij[is0];
									}
								}
							}
						}
						if (!baddC) {
							maxelem_threshold_theta = 0.98*Cmax;
							for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
								if (this_is_C_node[Amat.j[is0]]) {
									if (Amat.j[is0] != Amat.i[ii1]) {
										doublerealT aij_abs = Amat.abs_aij[is0];
										if (aij_abs > 0.98*Cmax) {
											sumP += aij_abs; // сумма модулей внедиагональных элементов которые принадлежат Strongly C узлам.
											icsos++;
										}
									}
								}
							}
						}

					}



					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																			   // Следующий случай назван "безнадёжным".
																			   // Безнадёжный случай оставляется без изменений.
																			   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																			   // интерполяцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом
																								 // и к тому-же не имеющих сильных F соседов.
					}



					{

						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						if (1) {
							//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
							// Случай одного Strong C соседа в добавок к которому есть по крайней мере один Strong F сосед.
							//
							SumPall += sumP;

							for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// Нас интересуют только Strong связи.


									if (this_is_C_node[Amat.j[is0]]) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
										if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
											if (fabs(sumP) < RealZERO) {
												//printf("error interpolation zero diagonal sumP.\n");
												//printf("Fnode all sosed is F");
												//system("pause");
												//	printf("i8 is Dirichlet node\n");
												if (this_is_C_node[i8] == false) iadditionalCstatistic++;
												this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
												this_is_C_node[i8] = true;
#pragma omp critical
												{
													if (!bweSholdbeContinue) {
														bweSholdbeContinue = true;
													}
												}
												//exit(1);
												// здесь нужна непрямая интерполяция.
											}
											else {
												// Это как раз единственная Strong C связь. 
												// Внедиагональный элемент из множества С узлов.

												// Данная вставка должна существенно сохранять 
												// разреженность оператора Галёркина на глубоких 
												// сеточных уровнях.
												// Модификация 5 декабря 2015.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
												if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {

													Ak1 Ptmp;

													Ptmp.j = i8;
													Ptmp.i = C_numerate[Amat.j[is0]];
													//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
													if (fabs(SumPall) < 1.0e-33) {
														printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
														//getchar();
														system("PAUSE");
														exit(1);
													}
													Ptmp.aij = fabs(Amat.aij[is0]) / SumPall;

													if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

														P_loc[itid][isize_P_loc[itid]] = Ptmp;
														isize_P_loc[itid]++;
													}
													else {
														std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
														system("PAUSE");
													}

													
												}

											}
										}

									}
									else
										if ((bactiveF[itid][Amat.j[is0]]) && (this_is_F_node[Amat.j[is0]])) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
											if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
												// Рассматриваем Strong F связи.

												// Смысл:
												// 



												//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
												// Это лишнее условие, ранее уже было проверено что мы имеем дело
												// с сильной F связью.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat.j[is0];
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
												// сканируем всех соседей iFpointeger 
												// чтобы найти среди них С узлы.

												// Идентифицируем локальный порог.
												doublerealT maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												integer ii1_loc = row_startA[iFpoint];
												integer istr_etalon = Amat.i[ii1_loc];
												integer iend_for = -1;
												integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
													iend_for = is0_loc;
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
														// Важно что threshold должен быть определен только по С узлам.
														// Иначе потеря сходимости.
														if (this_is_C_node[Amat.j[is0_loc]]) {
															if (Amat.j[is0_loc] != istr_etalon) {
																maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
															}
														}
													}
												}

												doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

												// Найти соседей данного F-node которые C-node.
												integer icsos_loc = 0;
												integer icsos_locall = 0;

												// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
												// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
												doublerealT sumP_loc = 0.0;
												doublerealT sumP_locall = 0.0;
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// Будем отсекать сразу наибольшее число вариантов.
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]) {

															if (Amat.j[is0_loc] != istr_etalon) {

																sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																icsos_loc++;
															}

														}
														else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {

															//if (Amat.j[is0_loc] != istr_etalon) {
															// Подсчитываем количество соседей которые не являются С узлами.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
															//}

															if (Amat.j[is0_loc] != istr_etalon) {



																// F node

																integer iFpoint_loc = Amat.j[is0_loc];

																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// сканируем всех соседей iFpointeger_loc 
																// чтобы найти среди них С узлы.

																// Идентифицируем локальный порог.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2) {
																		// Важно что threshold должен быть определен только по С узлам.
																		// Иначе потеря сходимости.
																		if (this_is_C_node[Amat.j[is0_loc2]]) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {
																				maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																			}
																		}
																	}
																}



																if (maxelem_threshold_loc2 < 0) {
																	//printf("Error !!! not found C nodes");
																	//system("PAUSE");
																	//exit(1);
																}
																else {
																	// Узел F без С узлов значит он не учитывается.

																	sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																	icsos_locall++;
																}
															}

														}

													}
												}

												doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												sumP_locall += sumP_loc;

												if (fabs(sumP_loc) > 1.0e-33) {
													// Учитывается только в том случае если есть сильные С соседи.

													// В общем это двойная интерполяционная сумма 
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


														// Внедиагональный элемент из множества С узлов.

														// Данная вставка должна существенно сохранять 
														// разреженность оператора Галёркина на глубоких 
														// сеточных уровнях.
														// Модификация 5 декабря 2015.
														if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {

																	Ak1 Ptmp;

																	Ptmp.j = i8;
																	Ptmp.i = C_numerate[Amat.j[is0_loc]];
																	//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																	if (fabs(sumP_loc) < 1.0e-33) {
																		printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//getchar();
																		system("PAUSE");
																		exit(1);
																	}
																	Ptmp.aij = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																	if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																		P_loc[itid][isize_P_loc[itid]] = Ptmp;
																		isize_P_loc[itid]++;
																	}
																	else {
																		std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																		system("PAUSE");
																	}

																}
															}
															else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {
																// F node

																integer iFpoint_loc = Amat.j[is0_loc];
																if (fabs(sumP_locall) < 1.0e-33) {
																	printf("error 2.a ! division by zero. SumP_locall =%e\n", sumP_locall);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// сканируем всех соседей iFpointeger_loc 
																// чтобы найти среди них С узлы.

																// Идентифицируем локальный порог.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2) {
																		// Важно что threshold должен быть определен только по С узлам.
																		// Иначе потеря сходимости.
																		if (this_is_C_node[Amat.j[is0_loc2]]) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {
																				maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																			}
																		}
																	}
																}

																doublerealT maxelem_threshold_loc_magic2 = maxelem_threshold_loc2 * magic_FF;

																if (maxelem_threshold_loc2 < 0) {
																	//printf("Error !!! not found C nodes");
																	//system("PAUSE");
																	//exit(1);
																}
																else {
																	// Здесь maxelem_threshold_loc2 это модуль максимального внедиагонального 
																	// элемента в строке среди С соседей локально.

																	// Найти соседей данного F-node которые C-node.
																	integer icsos_loc2 = 0;

																	// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
																	// Для каждого такого члена суммы увеличиваем счётчик iscos_loc2. 
																	// По идее iscos_loc2 должно быть 2 и более.
																	doublerealT sumP_loc2 = 0.0;
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																		// Будем отсекать сразу наибольшее число вариантов.
																		if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {

																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					sumP_loc2 += fabs(Amat.aij[is0_loc2]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																					icsos_loc2++;
																				}

																			}
																			else {

																				//if (Amat.j[is0_loc2] != istr_etalon2) {

																			}

																		}
																	}

																	doublerealT maxelem_threshold_loc_magic_minus2 = -maxelem_threshold_loc_magic2;

																	// В общем это двойная интерполяционная сумма 
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																		// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


																		// Внедиагональный элемент из множества С узлов.

																		// Данная вставка должна существенно сохранять 
																		// разреженность оператора Галёркина на глубоких 
																		// сеточных уровнях.
																		// Модификация 5 декабря 2015.
																		if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {
																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					Ak1 Ptmp;

																					Ptmp.j = i8;
																					Ptmp.i = C_numerate[Amat.j[is0_loc2]];
																					//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																					if (fabs(sumP_loc2) < 1.0e-33) {
																						printf("error 3.a ! division by zero. sumP_loc2 =%e\n", sumP_loc2);
																						//getchar();
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = multiplyer_nu_loc * fabs(Amat.aij[is0_loc2]) / sumP_loc2;

																					if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																						P_loc[itid][isize_P_loc[itid]] = Ptmp;
																						isize_P_loc[itid]++;
																					}
																					else {
																						std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																						system("PAUSE");
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
										}
								}

							}
						}

					}

					// Снимаем лок.
					permission_to_use_F[itid][i8] = true;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
					for (integer is0 = ii1; is0 <= iend_for1; is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							{
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


										// Сумма SumPall формируется только если у зла F есть сильные С соседи через magic.

										// Рассматриваем Strong F связи.

										// Смысл:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// Это лишнее условие, ранее уже было проверено что мы имеем дело
										// с сильной F связью.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// сканируем всех соседей iFpointeger 
										// чтобы найти среди них С узлы.

										// Идентифицируем локальный порог.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
											iend_for = is0_loc;
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
												// Важно что threshold должен быть определен только по С узлам.
												// Иначе потеря сходимости.
												if (this_is_C_node[Amat.j[is0_loc]]) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

										// Найти соседей данного F-node которые C-node.
										integer icsos_loc = 0;
										integer icsos_locall = 0;

										// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
										// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// Будем отсекать сразу наибольшее число вариантов.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// Подсчитываем количество соседей которые не являются С узлами.
													//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
													//}

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;
										if (fabs(sumP_loc) > 1.0e-33) {
											// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.
											permission_to_use_F[itid][Amat.j[is0]] = true;
										}
										bactiveF[itid][Amat.j[is0]] = true;// инициализация.
									}
								}

							}
						}
					}

				}
			
			}
			////
		}
		else {
			// only negative connections.

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
#pragma omp parallel for reduction(+ : number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes,  number_of_F_nodes_with_one_single_strong_C_neighborF, iadditionalCstatistic)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) {

#ifdef _OPENMP 
				int itid = omp_get_thread_num();
#else
				int itid = 0;
#endif

				if (this_is_F_node[i8]) {

					permission_to_use_F[itid][i8] = false;

					// Нам нужен разреженный оператор Галёркина.
					// 5 декабря 2015 года мы попробуем увеличить разреженность
					// оператора интерполяции а значит и оператора Галёркина.
					doublerealT maxelem_threshold = -1.0;
					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer istr_etalon1 = Amat.i[ii1];
					integer iend_for1 = -1;
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
							iend_for1 = is0;
							if (Amat.j[is0] != istr_etalon1) {
								// Если закомментировано то смотрится максимальный внедиагональный элемент в строке.
								//if (this_is_C_node[Amat.j[is0]]  ) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold)) {
									maxelem_threshold = fabs(Amat.aij[is0]);
								}
								//}
							}
						}
					}
					else {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
							iend_for1 = is0;
						}
						maxelem_threshold = threshold_quick_only_negative[istr_etalon1];
					}
					// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

					// Найти соседей данного F-node которые C-node.
					integer icsos = 0;
					integer icsosF = 0;

					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


					// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
					// Для каждого такого члена суммы увеличиваем счётчик iscos. По идее iscos должно быть 2 и более.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
					for (integer is0 = ii1; is0 <= iend_for1; is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
									sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {


										// Рассматриваем Strong F связи.

										// Смысл:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// Это лишнее условие, ранее уже было проверено что мы имеем дело
										// с сильной F связью.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// сканируем всех соседей iFpointeger 
										// чтобы найти среди них С узлы.

										// Идентифицируем локальный порог.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
											iend_for = is0_loc;
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												// Важно что threshold должен быть определен только по С узлам.
												// Иначе потеря сходимости.
												if (this_is_C_node[Amat.j[is0_loc]]) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

										// Найти соседей данного F-node которые C-node.
										integer icsos_loc = 0;
										integer icsos_locall = 0;

										// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
										// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// Будем отсекать сразу наибольшее число вариантов.
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// Подсчитываем количество соседей которые не являются С узлами.
													//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
													//}

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;

										if (fabs(sumP_loc) > 1.0e-33) {
											// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.

											bactiveF[itid][Amat.j[is0]] = true;

											permission_to_use_F[itid][Amat.j[is0]] = false;

											SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
											icsosF++;
										}
										else {
											bactiveF[itid][Amat.j[is0]] = false;
										}
									}
								}
								// Подсчитываем количество соседей которые не являются С узлами.
								the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
							}
						}
					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
																			   // Следующий случай назван "безнадёжным".
																			   // Безнадёжный случай оставляется без изменений.
																			   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																			   // интерполяцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
					}



					{

						if (1) {
							//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
							// Случай одного Strong C соседа в добавок к которому есть по крайней мере один Strong F сосед.
							//
							SumPall += sumP;

							for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// Нас интересуют только Strong связи.


									if (this_is_C_node[Amat.j[is0]]) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
										if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
											if (fabs(sumP) < RealZERO) {
												//printf("error interpolation zero diagonal sumP.\n");
												//printf("Fnode all sosed is F");
												//system("pause");
												//	printf("i8 is Dirichlet node\n");
												if (this_is_C_node[i8] == false) iadditionalCstatistic++;
												this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
												this_is_C_node[i8] = true;

#pragma omp critical
												{
													if (!bweSholdbeContinue) {
														bweSholdbeContinue = true;
													}
												}
												
												//exit(1);
												// здесь нужна непрямая интерполяция.
											}
											else {
												// Это как раз единственная Strong C связь. 
												// Внедиагональный элемент из множества С узлов.

												// Данная вставка должна существенно сохранять 
												// разреженность оператора Галёркина на глубоких 
												// сеточных уровнях.
												// Модификация 5 декабря 2015.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
												if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {

													Ak1 Ptmp;

													Ptmp.j = i8;
													Ptmp.i = C_numerate[Amat.j[is0]];
													//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
													if (fabs(SumPall) < 1.0e-23) {
														printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
														//getchar();
														system("PAUSE");
														exit(1);
													}
													Ptmp.aij = fabs(Amat.aij[is0]) / SumPall;

													if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

														P_loc[itid][isize_P_loc[itid]] = Ptmp;
														isize_P_loc[itid]++;
													}
													else {
														std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
														system("PAUSE");
													}													
												}

											}
										}

									}
									else
										if ((bactiveF[itid][Amat.j[is0]]) && (this_is_F_node[Amat.j[is0]])) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
											if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
												// Рассматриваем Strong F связи.

												// Смысл:
												// 



												//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
												// Это лишнее условие, ранее уже было проверено что мы имеем дело
												// с сильной F связью.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat.j[is0];
												doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
												// сканируем всех соседей iFpointeger 
												// чтобы найти среди них С узлы.

												// Идентифицируем локальный порог.
												doublerealT maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												integer ii1_loc = row_startA[iFpoint];
												integer istr_etalon = Amat.i[ii1_loc];
												integer iend_for = -1;
												integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
													iend_for = is0_loc;
													if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
														if (this_is_C_node[Amat.j[is0_loc]]) {
															if (Amat.j[is0_loc] != istr_etalon) {
																maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
															}
														}
													}
												}

												doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

												// Найти соседей данного F-node которые C-node.
												integer icsos_loc = 0;
												integer icsos_locall = 0;

												// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
												// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
												doublerealT sumP_loc = 0.0;
												doublerealT sumP_locall = 0.0;
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// Будем отсекать сразу наибольшее число вариантов.
													if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]) {

															if (Amat.j[is0_loc] != istr_etalon) {

																sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																icsos_loc++;
															}

														}
														else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {

															if (Amat.j[is0_loc] != istr_etalon) {

																// F node

																integer iFpoint_loc = Amat.j[is0_loc];

																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// сканируем всех соседей iFpointeger_loc 
																// чтобы найти среди них С узлы.

																// Идентифицируем локальный порог.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																	if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2)) {
																		// Важно что threshold должен быть определен только по С узлам.
																		// Иначе потеря сходимости.
																		if (this_is_C_node[Amat.j[is0_loc2]]) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {
																				maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																			}
																		}
																	}
																}



																if (maxelem_threshold_loc2 < 0) {
																	//printf("Error !!! not found C nodes");
																	//system("PAUSE");
																	//exit(1);
																}
																else {
																	// Узел F без С узлов значит он не учитывается.

																	sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																	icsos_locall++;
																}

															}
														}

													}
												}

												doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												sumP_locall += sumP_loc;

												if (fabs(sumP_loc) > 1.0e-33) {

													// Учитывается только в том случае если есть сильные С соседи.


													// В общем это двойная интерполяционная сумма 
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


														// Внедиагональный элемент из множества С узлов.

														// Данная вставка должна существенно сохранять 
														// разреженность оператора Галёркина на глубоких 
														// сеточных уровнях.
														// Модификация 5 декабря 2015.
														if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {

																	Ak1 Ptmp;

																	Ptmp.j = i8;
																	Ptmp.i = C_numerate[Amat.j[is0_loc]];
																	//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																	if (fabs(sumP_loc) < 1.0e-23) {
																		printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//getchar();
																		system("PAUSE");
																		exit(1);
																	}
																	Ptmp.aij = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_loc;
																	

																	if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																		P_loc[itid][isize_P_loc[itid]] = Ptmp;
																		isize_P_loc[itid]++;
																	}
																	else {
																		std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																		system("PAUSE");
																	}

																}
															}
															else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {
																// F node

																integer iFpoint_loc = Amat.j[is0_loc];
																if (fabs(sumP_locall) < 1.0e-33) {
																	printf("error 2.a ! division by zero. SumP_locall =%e\n", sumP_locall);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// сканируем всех соседей iFpointeger_loc 
																// чтобы найти среди них С узлы.

																// Идентифицируем локальный порог.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); is0_loc2++) {
																	if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2)) {
																		// Важно что threshold должен быть определен только по С узлам.
																		// Иначе потеря сходимости.
																		if (this_is_C_node[Amat.j[is0_loc2]]) {
																			if (Amat.j[is0_loc2] != istr_etalon2) {
																				maxelem_threshold_loc2 = fabs(Amat.aij[is0_loc2]);
																			}
																		}
																	}
																}

																doublerealT maxelem_threshold_loc_magic2 = maxelem_threshold_loc2 * magic_FF;


																if (maxelem_threshold_loc2 < 0) {
																	//printf("Error !!! not found C nodes");
																	//system("PAUSE");
																	//exit(1);
																}
																else {
																	// Здесь maxelem_threshold_loc2 это модуль максимального внедиагонального 
																	// элемента в строке среди С соседей локально.


																	// Найти соседей данного F-node которые C-node.
																	integer icsos_loc2 = 0;

																	// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
																	// Для каждого такого члена суммы увеличиваем счётчик iscos_loc2. 
																	// По идее iscos_loc2 должно быть 2 и более.
																	doublerealT sumP_loc2 = 0.0;
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																		// Будем отсекать сразу наибольшее число вариантов.
																		if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2)) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {

																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					sumP_loc2 += fabs(Amat.aij[is0_loc2]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																					icsos_loc2++;
																				}

																			}
																			else {

																				//if (Amat.j[is0_loc2] != istr_etalon2) {

																			}

																		}
																	}

																	doublerealT maxelem_threshold_loc_magic_minus2 = -maxelem_threshold_loc_magic2;

																	// В общем это двойная интерполяционная сумма 
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); is0_loc2++) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; is0_loc2++) {

																		// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


																		// Внедиагональный элемент из множества С узлов.

																		// Данная вставка должна существенно сохранять 
																		// разреженность оператора Галёркина на глубоких 
																		// сеточных уровнях.
																		// Модификация 5 декабря 2015.
																		if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2)) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {
																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					Ak1 Ptmp;

																					Ptmp.j = i8;
																					Ptmp.i = C_numerate[Amat.j[is0_loc2]];
																					//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																					if (fabs(sumP_loc2) < 1.0e-33) {
																						printf("error 3.a ! division by zero. sumP_loc2 =%e\n", sumP_loc2);
																						//getchar();
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = multiplyer_nu_loc * fabs(Amat.aij[is0_loc2]) / sumP_loc2;
																					
																					if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																						P_loc[itid][isize_P_loc[itid]] = Ptmp;
																						isize_P_loc[itid]++;
																					}
																					else {
																						std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																						system("PAUSE");
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
										}
								}

							}
						}

					}

					// Снимаем лок.
					permission_to_use_F[itid][i8] = true;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
					for (integer is0 = ii1; is0 <= iend_for1; is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							{
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


										// Сумма SumPall формируется только если у зла F есть сильные С соседи через magic.

										// Рассматриваем Strong F связи.

										// Смысл:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// Это лишнее условие, ранее уже было проверено что мы имеем дело
										// с сильной F связью.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// сканируем всех соседей iFpointeger 
										// чтобы найти среди них С узлы.

										// Идентифицируем локальный порог.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
											iend_for = is0_loc;
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												// Важно что threshold должен быть определен только по С узлам.
												// Иначе потеря сходимости.
												if (this_is_C_node[Amat.j[is0_loc]]) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

										// Найти соседей данного F-node которые C-node.
										integer icsos_loc = 0;
										integer icsos_locall = 0;

										// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
										// Для каждого такого члена суммы увеличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// Будем отсекать сразу наибольшее число вариантов.
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// Подсчитываем количество соседей которые не являются С узлами.
													//the_number_of_neighbors_that_are_not_C_nodes_loc++; // подсчитываем проблемы интерполяции 
													//}

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;
										if (fabs(sumP_loc) > 1.0e-33) {
											// Учитываем этот дополнительный F узел только в том случае если у него есть сильные C соседи.
											permission_to_use_F[itid][Amat.j[is0]] = true;
										}
										bactiveF[itid][Amat.j[is0]] = true;// инициализация.
									}
								}

							}
						}
					}

				} // end only negative connections
			
			}
		}

		iadditionalCstatistic_gl += iadditionalCstatistic;

		// Обратное копирование результата.
		for (int itid = 0; itid < 8; itid++) {

			if (icount1 + isize_P_loc[itid] >= nsizePR * n) {
				printf("memory error!!!\n");
				printf("not enough memory for the interpolation operator.\n");
				//system("PAUSE");
				//exit(1);
				deallocate_prolongation(nsizePR, n, P);
			}

#pragma omp parallel for
			for (integer iscan_loc = 0; iscan_loc < isize_P_loc[itid]; iscan_loc++) {
				P[icount1 + iscan_loc] = P_loc[itid][iscan_loc];
			}
			icount1 += isize_P_loc[itid];

		}


		for (int itid = 0; itid < 8; itid++) {
			delete[] permission_to_use_F[itid];
			delete[] bactiveF[itid];
		}

		delete[] permission_to_use_F;
		delete[] bactiveF;


		number_of_F_nodes_with_one_single_strong_C_neighbor_gl = number_of_F_nodes_with_one_single_strong_C_neighbor;
		the_number_of_neighbors_that_are_not_C_nodes_gl = the_number_of_neighbors_that_are_not_C_nodes;

	}


	number_of_F_nodes_with_one_single_strong_C_neighborF_gl += number_of_F_nodes_with_one_single_strong_C_neighborF;

	/*
	if (0) {
	// Надо читать теорию.
	// 26.09.2020 Сглаженная аггрегация.

	//Ak1* P_copy = new Ak1[nsizePR * n_a[0] + 1];

	//icount1
	//for (integer i_1 = 1; i_1 < icount1; i_1++) {
	//P_copy[i_1] = P[i_1];// Копия.
	//}

	bool* bvisible = new bool[n_a[ilevel - 1] + 1];
	for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
	bvisible[i_1] = false;
	}

	doublerealT* xP = new doublerealT[n_a[ilevel - 1] + 1];
	doublerealT* xP2 = new doublerealT[n_a[ilevel - 1] + 1];
	for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
	xP[i_1] = 0.0;
	xP2[i_1] = 0.0;
	}

	integer k_1 = icount1_memo;

	//for (integer i_1 = icount1_memo; i_1 < icount1; i_1++)
	for (integer i_1 = 1; i_1 < n_a[ilevel - 1]; i_1++) {
	if (this_is_F_node[i_1])
	{
	//if (!bvisible[P_copy[i_1].i])
	{

	// Инициализация.
	//for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++) {
	//xP[P_copy[i_1].j] = P_copy[i_1].aij;
	//}


	// smoother
	doublerealT omega = 0.6667;

	// Одна итерация.
	// Однопоточный вариант программы.
	//for (integer ii = 1; ii <= n_a[ilevel-1]; ii++)
	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
	{
	//integer ii = P_copy[i_1].i;
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
	{
	//integer	ii = P_copy[i_1].j;
	integer ii = Amat.j[s_1];

	integer istr = ii;
	doublerealT rold = xP[istr];

	// нулевая правая часть.

	doublerealT rsum = 0.0;
	doublerealT rsum_ap = 0.0;
	integer is1 = row_startA[ii];
	integer is2 = row_startA[ii + 1] - 1;
	doublerealT ap = 0.0;
	// Распараллеливание почемуто тормозит очень сильно.
	//#pragma omp parallel for reduction(+:rsum)
	for (integer ii1 = is1; ii1 <= is2; ii1++)
	{
	//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
	integer ipos = Amat.j[ii1];

	if (ipos != Amat.i[ii1]) {
	if (this_is_C_node[ipos]) {
	if (Amat.abs_aij[ii1] > theta * threshold_quick_all[Amat.i[ii1]]) {
	rsum += Amat.abs_aij[ii1] * xP[ipos];
	rsum_ap += Amat.abs_aij[ii1];
	}
	}
	}
	//else {
	//ap = Amat.aij[ii1];
	//}
	}
	xP2[istr] += rsum;
	//xP2[istr] *= ap;// делим на диагональный элемент, т.к. в
	//диагонали хранится обратная величина для диагонального коэффициента.
	xP2[istr] /= rsum_ap;

	// Запускается только после первого раза сразу запускать нельзя
	// т.к. начальное значение не является адекватным.
	xP2[istr] = omega * xP2[istr] + (1.0 - omega) * rold; // this is SOR
	}
	}

	// find_max
	doublerealT xPmax = -1.0;
	doublerealT sumxP = 0.0;
	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
	{
	//integer	ii = P_copy[i_1].j;
	integer ii = Amat.j[s_1];

	{
	xP[ii] = xP2[ii];
	if (xP[ii] > xPmax) {
	xPmax = xP[ii];
	}
	sumxP += xP[ii];
	}
	}
	}
	// Просеивание.
	doublerealT sumxPnew = 0.0;
	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
	{
	//integer	ii = P_copy[i_1].j;
	integer ii = Amat.j[s_1];


	{
	if ((this_is_F_node[ii]) || (xP[ii] <= 0.2 * xPmax)) {
	xP[ii] = 0.0;
	}
	else {
	sumxPnew += xP[ii];
	}
	}
	}
	}
	if (sumxPnew > 1.0e-30) {

	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
	{
	//integer	ii = P_copy[i_1].j;
	integer ii = Amat.j[s_1];

	{
	xP[ii] *= sumxP / sumxPnew;
	if ((this_is_C_node[ii]) && (xP[ii] > 1.0e-30)) {
	P[k_1].aij = xP[ii];
	P[k_1].i = P_copy[i_1].i;
	P[k_1].j = ii;
	k_1++;
	if (k_1 >= nsizePR * n) {
	printf("memory error!!!\n");
	printf("not enough memory for the interpolation operator.\n");
	//system("PAUSE");
	//exit(1);
	deallocate_prolongation(nsizePR, n, P);
	}
	}
	}
	}
	}

	}
	else {
	printf("ERROR !!! zero string in matrix P\n");
	printf("string number=%lld\n", P_copy[i_1].i);
	system("PAUSE");
	//this_is_C_node[P_copy[i_1].i] = true;
	//this_is_F_node[P_copy[i_1].i] = false;
	//P[k_1].aij = 1.0;
	//P[k_1].i = P_copy[i_1].i;
	//P[k_1].j = P_copy[i_1].i;
	//k_1++;
	}

	bvisible[P_copy[i_1].i] = true;

	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); j_1++)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; s_1++)
	{
	//integer	ii = P_copy[i_1].j;
	integer ii = Amat.j[s_1];
	{
	xP[ii] = 0.0;
	xP2[ii] = 0.0;
	}
	}
	}
	}
	}
	}
	delete[] P_copy;
	delete[] bvisible;
	delete[] xP;
	delete[] xP2;

	icount1 = k_1;

	}
	*/

	//system("pause");

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif


} // my_interpolation_procedure_number3A_parallelx8





#endif
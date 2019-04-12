
#pragma once
#ifndef MY_AMG_INTERPOLATION_CPP 
#define MY_AMG_INTERPOLATION_CPP 1


// 3 ноября 2016.
// 17.02.2019
// Вызывается единожды лишь в случае нехватки выделеной памяти под оператор 
// интерполляции.
// Это происходит очень редко - на узком классе задач.
// Здесь происходит перевыделение оперативной памяти и увеличение значения nsizePR.
void deallocate_prolongation(integer &nsizePR, // Память под P в количествах n.
	integer n, // Количество неизвестных в СЛАУ на нулевом уровне.
	Ak1* &P // prolongation
)
{

	// Синтакис вызова:
	//deallocate_prolongation(nsizePR,n,R,P);

	// Для большинства реальных задач нам достаточно объёма памяти в nsizePR==12.
	// Однако встречаются задачи для котрыхнужно существенно большее количество памяти.
	// В этом случае я предлагаю увеличить её количество до nsizePR==35.
	if (nsizePR < 13) {

		Ak1* P_copy = NULL; // prolongation

		// Увеличиваем количество памяти до nsizePR==35.

		//Ak1* &P_copy = new Ak1[nsizePR*n + 1];

		//P = new Ak1[(integer)(35 * n) + 1]; // 3*nnz 2.4 // 35
		P_copy = (Ak1*)malloc(((nsizePR * n) + 1) * sizeof(Ak1));
		if (P_copy == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for P_copy matrix in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}

		if ((P_copy != NULL)) {
			for (integer i = 0; i <= (nsizePR * n); i++) {
				P_copy[i] = P[i];
			}
			if (P != NULL) {
				free(P);
			}
			P = NULL;
			integer nsizePR_old = nsizePR;
			nsizePR = 35; // Этого должно хватить с запасом.


			//P = new Ak1[(integer)(35 * n) + 1]; // 3*nnz 2.4 // 35
			P = (Ak1*)malloc(((nsizePR * n) + 1) * sizeof(Ak1));
			if (P == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for P matrix in my_agregat_amg.cpp algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}

			if ((P != NULL)) {
				for (integer i = 0; i <= (nsizePR_old * n); i++) {
					P[i] = P_copy[i];
				}
			}

			// Обязательно освобождаем оперативную память.
			if (P_copy != NULL) {
				free(P_copy);
				P_copy = NULL;
			}
		}
	}
	else {
#if doubleintprecision == 1
		printf("FATALL error!!! nsizePR=%lld\n", nsizePR);
#else
		printf("FATALL error!!! nsizePR=%d\n", nsizePR);
#endif

		printf("not enough memory for the interpolation operator.");
		printf("absolute fatal. see deallocate_prolongation in my_agregat_amg.cpp\n");
		system("PAUSE");
		exit(1);
	}

} // deallocate_prolongation


// Интерполяционная процедура №3.
// Протестированный вариант освобожденный от лишнего кода.
// За счёт освобождения от лишнего кода улучшена ясность кода.
// 12.05.2018
void my_interpolation_procedure_number3A(integer &the_number_of_neighbors_that_are_not_С_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak1* &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublereal RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublereal &theta, integer &n, Ak1* &R, integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublereal &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublereal &magic82, doublereal* &threshold_quick_all,
	doublereal* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	doublereal theta_strong_F = 0.23; // оптимальный выбор.
	theta_strong_F = theta83; // 3 июля 2016


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

	// the begining

	bool byes_add = false;
	// Быстрое добавление недостающих С узлов.
	if (1) {
		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//doublerealT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;		
			

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {


				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer icsosF = 0;

				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;

				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer iend_marker_position = row_startA[Amat[ii1].i + 1] - 1;

				if (bpositive_connections) {


#if doubleintprecision == 1
					//printf("i8=%lld n=%lld\n", i8, n_a[ilevel - 1]);
#else
					//printf("i8=%d n=%d\n", i8, n_a[ilevel - 1]);
#endif

				//getchar();


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
					doublereal maxelem_threshold = -1.0;
				
					//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
					
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
								//if (this_is_C_node[Amat[is0].j] == true) {
								if (fabs(Amat[is0].aij) > maxelem_threshold) {
									maxelem_threshold = fabs(Amat[is0].aij);
								}
								//}
							}
						}
					}
					else {
						maxelem_threshold = threshold_quick_all[Amat[ii1].i];
					}
					// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

					
					doublereal maxelem_threshold_theta = maxelem_threshold * theta;
					doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
						if (Amat[is0].j != Amat[ii1].i) {
							if (this_is_C_node[Amat[is0].j] == true) {
								//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
								if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
									sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat[is0].j] == true) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
										SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// Подсчитываем количество соседей которые не являются С узлами.
								the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
							}
						}
					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																			   // Следующий случай назван "безнадёжным".
																			   // Безнадёжный случай оставляется без изменений.
																			   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																			   // интерполляцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
					}
				}
				else {
					// only negative connections

					
						// Нам нужен разреженный оператор Галёркина.
						// 5 декабря 2015 года мы попробуем увеличить разреженность
						// оператора интерполляции а значит и оператора Галёркина.
						doublereal maxelem_threshold = -1.0;
						
						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
		
						if (!btreshold_on_new_vetv) {
							for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								if (Amat[is0].j != Amat[ii1].i) {
									// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
									//if (this_is_C_node[Amat[is0].j] == true) {
									if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold)) {
										maxelem_threshold = fabs(Amat[is0].aij);
									}
									//}
								}
							}
						}
						else {
							maxelem_threshold = threshold_quick_only_negative[Amat[ii1].i];
						}
						// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

						
						doublereal maxelem_threshold_theta = maxelem_threshold * theta;
						doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								if (this_is_C_node[Amat[is0].j] == true) {
									//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
										sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
										icsos++;
									}
								}
								else {
									if (this_is_F_node[Amat[is0].j] == true) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
											SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
											icsosF++;
										}
									}
									// Подсчитываем количество соседей которые не являются С узлами.
									the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
								}
							}
						}
						if (icsos == 1) {
							number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																				   // Следующий случай назван "безнадёжным".
																				   // Безнадёжный случай оставляется без изменений.
																				   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																				   // интерполляцию на этих уровнях можно улучшить.
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


				// Мы делаем из F узла C узел только в том случае если у нас нету С связей.
				//integer iend_marker_position = row_startA[Amat[ii1].i + 1] - 1;
				for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
					if (this_is_F_node[Amat[is0].j] == true) {
						if (Amat[is0].j != Amat[ii1].i) {


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
							    if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // Были добавления узлов.
												 //exit(1);
												 // здесь нужна непрямая интерполляция.
							    }
							}


						}
					}
				}//end

			}

	}
	




	/*
	9.05.2018
	1. Определили thresholdы.
	2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
	3. icsos - количество strong C. icsosF - количество strong F.
	4. Если SumP ==0 то узел становится С узлом. iadditionalCstatistic++;
	0;
	47; 0;
	186; 0;
	267; 0;
	296; 0;
	276; 0;
	221; 0;
	257.6; 0;
	335.2; 0;
	11L 3min 10s
	*/

	if (!byes_add) {

		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//--->doublereal magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 июля 2016 для конвективных задач
		// линейный прогноз не дает улучшений
		// ни по быстрдействию ни по числу V циклов.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublereal magic = magic82; // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold) {
								maxelem_threshold = fabs(Amat[is0].aij);
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

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j] == true) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j] == true) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				 {

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat[is0].j] == true) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
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
											// здесь нужна непрямая интерполляция.
										}
										else {
											// Это как раз единственная Strong C связь. 
											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = fabs(Amat[is0].aij) / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n, R, P);
												}
											}

										}
									}

								}
								else
									if (this_is_F_node[Amat[is0].j] == true) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
											// Рассматриваем Strong F связи.

											// Смысл :
											// 



											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat[is0].j;
											if (fabs(SumPall) < 1.0e-33) {
												printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
												//getchar();
												system("PAUSE");
												exit(1);
											}
											doublereal multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublereal maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat[ii1_loc].i;
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc) {
													// Важно что threshold должен быть определен только по С узлам.
													// Иначе потеря сходимости.
													if (this_is_C_node[Amat[is0_loc].j] == true) {
														if (Amat[is0_loc].j != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
														}
													}
												}
											}

											doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublereal sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat[is0_loc].j] == true) {

														if (Amat[is0_loc].j != istr_etalon) {

															sumP_loc += fabs(Amat[is0_loc].aij); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else {

														//if (Amat[is0_loc].j != istr_etalon) {
														// Подсчитываем количество соседей которые не являются С узлами.
														//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
														//}
													}

												}
											}

											doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// В общем это двойная интерполляционная сумма 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


												// Внедиагональный элемент из множества С узлов.

												// Данная вставка должна существенно сохранять 
												// разреженность оператора Галёркина на глубоких 
												// сеточных уровнях.
												// Модификация 5 декабря 2015.
												if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {
													
													if (this_is_C_node[Amat[is0_loc].j] == true) {
														if (Amat[is0_loc].j != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat[is0_loc].j];
															//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
															if (fabs(sumP_loc) < 1.0e-33) {
																printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																//getchar();
																system("PAUSE");
																exit(1);
															}
															P[icount1].aij = multiplyer_nu * fabs(Amat[is0_loc].aij) / sumP_loc;
															icount1++;
															if (icount1 >= nsizePR * n) {
																printf("memory error!!!\n");
																printf("not enough memory for the interpolation operator.\n");
																//system("PAUSE");
																//exit(1);
																deallocate_prolongation(nsizePR, n, R, P);
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
			////
		}
		else {
			// only negative connections.

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold)) {
								maxelem_threshold = fabs(Amat[is0].aij);
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

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j] == true) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j] == true) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				 {

					if (1) {
					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat[is0].j] == true) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
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
											// здесь нужна непрямая интерполляция.
										}
										else {
											// Это как раз единственная Strong C связь. 
											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < 1.0e-23) {
													printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = fabs(Amat[is0].aij) / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n, R, P);
												}
											}

										}
									}

								}
								else
									if (this_is_F_node[Amat[is0].j] == true) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
											// Рассматриваем Strong F связи.

											// Смысл :
											// 



											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat[is0].j;
											doublereal multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublereal maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat[ii1_loc].i;
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat[is0_loc].j] == true) {
														if (Amat[is0_loc].j != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
														}
													}
												}
											}

											doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublereal sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
													
													if (this_is_C_node[Amat[is0_loc].j] == true) {

														if (Amat[is0_loc].j != istr_etalon) {

															sumP_loc += fabs(Amat[is0_loc].aij); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else {

														//if (Amat[is0_loc].j != istr_etalon) {
														// Подсчитываем количество соседей которые не являются С узлами.
														//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
														//}
													}

												}
											}

											doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// В общем это двойная интерполляционная сумма 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


												// Внедиагональный элемент из множества С узлов.

												// Данная вставка должна существенно сохранять 
												// разреженность оператора Галёркина на глубоких 
												// сеточных уровнях.
												// Модификация 5 декабря 2015.
												if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
													
													if (this_is_C_node[Amat[is0_loc].j] == true) {
														if (Amat[is0_loc].j != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat[is0_loc].j];
															//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
															if (fabs(sumP_loc) < 1.0e-23) {
																printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																//getchar();
																system("PAUSE");
																exit(1);
															}
															P[icount1].aij = multiplyer_nu * fabs(Amat[is0_loc].aij) / sumP_loc;
															icount1++;
															if (icount1 >= nsizePR * n) {
																printf("memory error!!!\n");
																printf("not enough memory for the interpolation operator.\n");
																//system("PAUSE");
																//exit(1);
																deallocate_prolongation(nsizePR, n, R, P);
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
			} // end only negative connections

		}

	}



	//system("pause");


} // my_interpolation_procedure_number3A

// Интерполяционная процедура №3.
// Протестированный вариант освобожденный от лишнего кода.
// За счёт освобождения от лишнего кода улучшена ясность кода.
// 12.05.2018
void my_interpolation_procedure_number3A(integer &the_number_of_neighbors_that_are_not_С_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak2 &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublereal RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublereal &theta, integer &n,  integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublereal &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublereal &magic82, doublereal* &threshold_quick_all,
	doublereal* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	doublereal theta_strong_F = 0.23; // оптимальный выбор.
	theta_strong_F = theta83; // 3 июля 2016


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

	// the begining

	bool byes_add = false;
	// Быстрое добавление недостающих С узлов.
	if (1) {
		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//doublerealT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


		// Построение пролонгации для узлов которые составляют F nodes.
		// Каждый F-nodes окружён C-nodes.
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {


			// Найти соседей данного F-node которые C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
			// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
			doublereal sumP = 0.0;
			doublereal SumPall = 0.0;
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
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
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


				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]] == true) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
					// Следующий случай назван "безнадёжным".
					// Безнадёжный случай оставляется без изменений.
					// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
					// интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}
			}
			else {
				// only negative connections


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
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


				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]] == true) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
					// Следующий случай назван "безнадёжным".
					// Безнадёжный случай оставляется без изменений.
					// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
					// интерполляцию на этих уровнях можно улучшить.
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


			// Мы делаем из F узла C узел только в том случае если у нас нету С связей.
			//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
			for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				if (this_is_F_node[Amat.j[is0]] == true) {
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
								if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // Были добавления узлов.
								//exit(1);
								// здесь нужна непрямая интерполляция.
							}
						}


					}
				}
			}//end

		}

	}





	/*
	9.05.2018
	1. Определили thresholdы.
	2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
	3. icsos - количество strong C. icsosF - количество strong F.
	4. Если SumP ==0 то узел становится С узлом. iadditionalCstatistic++;
	0;
	47; 0;
	186; 0;
	267; 0;
	296; 0;
	276; 0;
	221; 0;
	257.6; 0;
	335.2; 0;
	11L 3min 10s
	*/

	if (!byes_add) {

		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//--->doublereal magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 июля 2016 для конвективных задач
		// линейный прогноз не дает улучшений
		// ни по быстрдействию ни по числу V циклов.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublereal magic = magic82; // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
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

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]] == true) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
					// Следующий случай назван "безнадёжным".
					// Безнадёжный случай оставляется без изменений.
					// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
					// интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat.j[is0]] == true) {
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
											// здесь нужна непрямая интерполляция.
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
								if (this_is_F_node[Amat.j[is0]] == true) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
										// Рассматриваем Strong F связи.

										// Смысл :
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
										doublereal multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
										// сканируем всех соседей iFpointeger 
										// чтобы найти среди них С узлы.

										// Идентифицируем локальный порог.
										doublereal maxelem_threshold_loc = -1.0;
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
												if (this_is_C_node[Amat.j[is0_loc]] == true) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

										// Найти соседей данного F-node которые C-node.
										integer icsos_loc = 0;

										// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
										// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
										doublereal sumP_loc = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// Будем отсекать сразу наибольшее число вариантов.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]] == true) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// Подсчитываем количество соседей которые не являются С узлами.
													//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
													//}
												}

											}
										}

										doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										// В общем это двойная интерполляционная сумма 
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]] == true) {
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
														P[icount1].aij = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_loc;
														icount1++;
														if (icount1 >= nsizePR * n) {
															printf("memory error!!!\n");
															printf("not enough memory for the interpolation operator.\n");
															//system("PAUSE");
															//exit(1);
															deallocate_prolongation(nsizePR, n,  P);
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
			////
		}
		else {
			// only negative connections.

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
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

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]] == true) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
					// Следующий случай назван "безнадёжным".
					// Безнадёжный случай оставляется без изменений.
					// К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
					// интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat.j[is0]] == true) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
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
											// здесь нужна непрямая интерполляция.
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
											if ((Amat.aij[is0] <0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
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
								if (this_is_F_node[Amat.j[is0]] == true) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
										// Рассматриваем Strong F связи.

										// Смысл :
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// Это лишнее условие, ранее уже было проверено что мы имеем дело
										// с сильной F связью.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];
										doublereal multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
										// сканируем всех соседей iFpointeger 
										// чтобы найти среди них С узлы.

										// Идентифицируем локальный порог.
										doublereal maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
											iend_for = is0_loc;
											if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												if (this_is_C_node[Amat.j[is0_loc]] == true) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

										// Найти соседей данного F-node которые C-node.
										integer icsos_loc = 0;

										// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
										// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
										doublereal sumP_loc = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// Будем отсекать сразу наибольшее число вариантов.
											if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]] == true) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// Подсчитываем количество соседей которые не являются С узлами.
													//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
													//}
												}

											}
										}

										doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										// В общем это двойная интерполляционная сумма 
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]] == true) {
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
															deallocate_prolongation(nsizePR, n,  P);
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
			} // end only negative connections

		}

	}



	//system("pause");


} // my_interpolation_procedure_number3A


  // Интерполяционная процедура №3.
  // Протестированный вариант освобожденный от лишнего кода.
  // За счёт освобождения от лишнего кода улучшена ясность кода.
  // 12.05.2018
void my_interpolation_procedure_number3B(integer &the_number_of_neighbors_that_are_not_С_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak1* &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublereal RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublereal &theta, integer &n, Ak1* &R, integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublereal &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublereal &magic82, doublereal* &threshold_quick_all,
	doublereal* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	doublereal theta_strong_F = 0.23; // оптимальный выбор.
	theta_strong_F = theta83; // 3 июля 2016

	// хеш таблица
	bool* this_is_Strong_C_node = new bool[n_a[ilevel-1]+1];
	for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
		// инициализация хеш таблицы.
		this_is_Strong_C_node[i_1] = false;
	}


	// четвёртая версия интерполяции.
	// интерполляция Руге и Стубена.
	integer ioneStrongC_and_0_StrongF = 0;

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

	// the begining

	bool byes_add = false;
	// Быстрое добавление недостающих С узлов.
	if (1) {
		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//doublerealT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


		// Построение пролонгации для узлов которые составляют F nodes.
		// Каждый F-nodes окружён C-nodes.
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {


			bool badd_amg1r5 = false;

			// Найти соседей данного F-node которые C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
			// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
			doublereal sumP = 0.0;
			doublereal SumPall = 0.0;
			integer icount_StronglyF = 0;

			//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
			integer ii1 = row_startA[i8];
			integer iend_marker_position = row_startA[Amat[ii1].i + 1] - 1;

			if (bpositive_connections) {


#if doubleintprecision == 1
				//printf("i8=%lld n=%lld\n", i8, n_a[ilevel - 1]);
#else
				//printf("i8=%d n=%d\n", i8, n_a[ilevel - 1]);
#endif

				//getchar();


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat[is0].j != Amat[ii1].i) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold) {
								maxelem_threshold = fabs(Amat[is0].aij);
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[Amat[ii1].i];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// Инициализация хеш таблицы.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								// заносим всех сильных С соседей.
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}

				
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j] == true) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
									
									integer ii2 = row_startA[Amat[is0].j];
									integer iend_marker_position2 = row_startA[Amat[ii2].i + 1] - 1;
									bool bfound = false;
									for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
										if (Amat[is2].j != Amat[ii2].i) {
											if (this_is_Strong_C_node[Amat[is2].j]) {
												bfound = true;
											}
										}
									}

									
										if (bfound) {
											// Заносим только тех сильных F соседей которые имеют в соседях(любых не обязательно сильных) первоначальных сильных С соселей.

											SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
												icount_StronglyF++;
												icsosF++;
										}
										else {
											badd_amg1r5 = true;
										}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}

				// Сброс - очистка хеш таблицы.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat[is0].j] = false;
							}
						}
					}
				}

				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}
			}
			else {
				// only negative connections


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat[is0].j != Amat[ii1].i) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold)) {
								maxelem_threshold = fabs(Amat[is0].aij);
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_only_negative[Amat[ii1].i];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.


				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				
				// Инициализация хеш таблицы.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}

				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j] == true) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
									
									
									integer ii2 = row_startA[Amat[is0].j];
									integer iend_marker_position2 = row_startA[Amat[ii2].i + 1] - 1;
									bool bfound = false;
									for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
										if (Amat[is2].j != Amat[ii2].i) {
											if (this_is_Strong_C_node[Amat[is2].j]) {												
												bfound = true;
											}
										}
									}
									
									if (bfound) {

										// Заносим только тех сильных F соседей которые имеют в соседях(любых не обязательно сильных) первоначальных сильных С соселей.
										SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
											icsosF++;
									}
									else {
										badd_amg1r5 = true;
									}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}

				// Сброс - очистка хеш таблицы.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = false;
							}
						}
					}
				}

				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
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


			// Мы делаем из F узла C узел только в том случае если у нас нету С связей.
			//integer iend_marker_position = row_startA[Amat[ii1].i + 1] - 1;
			for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				if (this_is_F_node[Amat[is0].j] == true) {
					if (Amat[is0].j != Amat[ii1].i) {


						// 20 jan 2016.
						// Можно не добавлять узел в случае если есть не менее трёх сильных F соседей.


						if ((fabs(sumP) < RealZERO)||(badd_amg1r5)) {
							// Это случай когда вообще нет сильных С соседей.

							// if (icsosF==0) полная потеря сходимости.
							//if (icsosF<2)// наблюдается замедление скорости сходимости и времени расчёта.
							{
								//printf("error interpolation zero diagonal sumP.\n");
								//printf("Fnode all sosed is F");
								//system("pause");
								//printf("i8 is Dirichlet node\n");
								if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // Были добавления узлов.
												 //exit(1);
												 // здесь нужна непрямая интерполляция.
							}
						}


					}
				}
			}//end





		}

	}





	/*
	9.05.2018
	1. Определили thresholdы.
	2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
	3. icsos - количество strong C. icsosF - количество strong F.
	4. Если SumP ==0 то узел становится С узлом. iadditionalCstatistic++;
	0;
	47; 0;
	186; 0;
	267; 0;
	296; 0;
	276; 0;
	221; 0;
	257.6; 0;
	335.2; 0;
	11L 3min 10s
	*/

	if (!byes_add) {

		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//---->doublereal magic = 0.4; // 0.4 optimum
								//magic = 0.3; // 3 июля 2016 для конвективных задач
								// линейный прогноз не дает улучшений
								// ни по быстрдействию ни по числу V циклов.
								//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublereal magic = magic82; // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold) {
								maxelem_threshold = fabs(Amat[is0].aij);
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

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}



				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j] == true) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j] == true) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {

									integer ii2 = row_startA[Amat[is0].j];
									integer iend_marker_position2 = row_startA[Amat[ii2].i + 1] - 1;
									bool bfound = false;
									for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
										if (Amat[is2].j != Amat[ii2].i) {
											if (this_is_Strong_C_node[Amat[is2].j]) {
												bfound = true;
											}
										}
									}
									if (bfound) {


										SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
											icount_StronglyF++;
										icsosF++;
									}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat[is0].j] == true) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
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
											// здесь нужна непрямая интерполляция.
										}
										else {
											// Это как раз единственная Strong C связь. 
											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = fabs(Amat[is0].aij) / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n, R, P);
												}
											}

										}
									}

								}
								else
									if (this_is_F_node[Amat[is0].j] == true) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
											// Рассматриваем Strong F связи.

											// Смысл :
											// 

											integer ii2 = row_startA[Amat[is0].j];
											integer iend_marker_position2 = row_startA[Amat[ii2].i + 1] - 1;
											bool bfound = false;
											for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
												if (Amat[is2].j != Amat[ii2].i) {
													if (this_is_Strong_C_node[Amat[is2].j]) {
														bfound = true;
													}
												}
											}
											if (bfound) {

												//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
												// Это лишнее условие, ранее уже было проверено что мы имеем дело
												// с сильной F связью.
												//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat[is0].j;
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												doublereal multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
												// сканируем всех соседей iFpointeger 
												// чтобы найти среди них С узлы.

												// Идентифицируем локальный порог.
												doublereal maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												integer ii1_loc = row_startA[iFpoint];
												integer istr_etalon = Amat[ii1_loc].i;
												integer iend_for = -1;
												integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
													iend_for = is0_loc;
													//if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc) {
														// Важно что threshold должен быть определен только по С узлам.
														// Иначе потеря сходимости.
														if (this_is_C_node[Amat[is0_loc].j] == true) {
															if (this_is_Strong_C_node[Amat[is0_loc].j]) {
																if (Amat[is0_loc].j != istr_etalon) {
																	maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
																}
															}
														}
													//}
												}

												doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

												// Найти соседей данного F-node которые C-node.
												integer icsos_loc = 0;

												// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
												// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
												doublereal sumP_loc = 0.0;
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// Будем отсекать сразу наибольшее число вариантов.
													//if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {

													if (this_is_Strong_C_node[Amat[is0_loc].j]) {
														if (this_is_C_node[Amat[is0_loc].j] == true) {

															if (Amat[is0_loc].j != istr_etalon) {

																sumP_loc += fabs(Amat[is0_loc].aij); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																icsos_loc++;
															}

														}
														else {

															//if (Amat[is0_loc].j != istr_etalon) {
															// Подсчитываем количество соседей которые не являются С узлами.
															//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
															//}
														}
													}

													//}
												}

												doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												// В общем это двойная интерполляционная сумма 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


													// Внедиагональный элемент из множества С узлов.

													// Данная вставка должна существенно сохранять 
													// разреженность оператора Галёркина на глубоких 
													// сеточных уровнях.
													// Модификация 5 декабря 2015.
													//if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {
													if (this_is_Strong_C_node[Amat[is0_loc].j]) {

														if (this_is_C_node[Amat[is0_loc].j] == true) {
															if (Amat[is0_loc].j != istr_etalon) {

																P[icount1].j = i8;
																P[icount1].i = C_numerate[Amat[is0_loc].j];
																//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
																if (fabs(sumP_loc) < 1.0e-33) {
																	printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																P[icount1].aij = multiplyer_nu * fabs(Amat[is0_loc].aij) / sumP_loc;
																icount1++;
																if (icount1 >= nsizePR * n) {
																	printf("memory error!!!\n");
																	printf("not enough memory for the interpolation operator.\n");
																	//system("PAUSE");
																	//exit(1);
																	deallocate_prolongation(nsizePR, n, R, P);
																}

															}
														}
													}
													//}
												}
											}

										}
									}
							}

						}
					}

				}

				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat[is0].j] = false;
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
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold)) {
								maxelem_threshold = fabs(Amat[is0].aij);
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

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j] == true) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								sumP += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j] == true) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
									
									integer ii2 = row_startA[Amat[is0].j];
									integer iend_marker_position2 = row_startA[Amat[ii2].i + 1] - 1;
									bool bfound = false;
									for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
										if (Amat[is2].j != Amat[ii2].i) {
											if (this_is_Strong_C_node[Amat[is2].j]) {
												bfound = true;
											}
										}
									}
									if (bfound) {


										SumPall += fabs(Amat[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat[is0].j] == true) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
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
											// здесь нужна непрямая интерполляция.
										}
										else {
											// Это как раз единственная Strong C связь. 
											// Внедиагональный элемент из множества С узлов.

											// Данная вставка должна существенно сохранять 
											// разреженность оператора Галёркина на глубоких 
											// сеточных уровнях.
											// Модификация 5 декабря 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < 1.0e-23) {
													printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = fabs(Amat[is0].aij) / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n, R, P);
												}
											}

										}
									}

								}
								else
									if (this_is_F_node[Amat[is0].j] == true) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
											// Рассматриваем Strong F связи.

											// Смысл :
											// 



											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat[is0].j;
											doublereal multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublereal maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat[ii1_loc].i;
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												//if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat[is0_loc].j] == true) {
														if (this_is_Strong_C_node[Amat[is0_loc].j]) {
															if (Amat[is0_loc].j != istr_etalon) {
																maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
															}
														}
													}
												//}
											}

											doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublereal sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												//if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat[is0_loc].j]) {


													if (this_is_C_node[Amat[is0_loc].j] == true) {

														if (Amat[is0_loc].j != istr_etalon) {

															sumP_loc += fabs(Amat[is0_loc].aij); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else {

														//if (Amat[is0_loc].j != istr_etalon) {
														// Подсчитываем количество соседей которые не являются С узлами.
														//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
														//}
													}
												}
												//}
											}

											doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// В общем это двойная интерполляционная сумма 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


												// Внедиагональный элемент из множества С узлов.

												// Данная вставка должна существенно сохранять 
												// разреженность оператора Галёркина на глубоких 
												// сеточных уровнях.
												// Модификация 5 декабря 2015.
												//if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat[is0_loc].j]) {


													if (this_is_C_node[Amat[is0_loc].j] == true) {
														if (Amat[is0_loc].j != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat[is0_loc].j];
															//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
															if (fabs(sumP_loc) < 1.0e-23) {
																printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																//getchar();
																system("PAUSE");
																exit(1);
															}
															P[icount1].aij = multiplyer_nu * fabs(Amat[is0_loc].aij) / sumP_loc;
															icount1++;
															if (icount1 >= nsizePR * n) {
																printf("memory error!!!\n");
																printf("not enough memory for the interpolation operator.\n");
																//system("PAUSE");
																//exit(1);
																deallocate_prolongation(nsizePR, n, R, P);
															}

														}
													}
												}
											//}
											}


										}
									}
							}

						}
					}

				}

				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j] == true) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = false;
							}
						}
					}
				}

			} // end only negative connections

		}

	}

#if doubleintprecision == 1
	//printf("one Strong C and 0 Strong F=%lld additional all %lld", ioneStrongC_and_0_StrongF, iadditionalCstatistic);
#else
	//printf("one Strong C and 0 Strong F=%d additional all %d", ioneStrongC_and_0_StrongF, iadditionalCstatistic);
#endif

	//system("pause");

	delete[] this_is_Strong_C_node;
	this_is_Strong_C_node = NULL;

} // my_interpolation_procedure_number3B

  // Интерполяционная процедура №3.
  // Протестированный вариант освобожденный от лишнего кода.
  // За счёт освобождения от лишнего кода улучшена ясность кода.
  // 12.05.2018
// 17-18 марта 2019 года подвергнута оптимизации. Быстродействие всего приложения возросло на 1.6%.
void my_interpolation_procedure_number3B(integer &the_number_of_neighbors_that_are_not_С_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak2 &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublereal RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublereal &theta, integer &n,  integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublereal &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublereal &magic82, doublereal* &threshold_quick_all,
	doublereal* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	doublereal theta_strong_F = 0.23; // оптимальный выбор.
	theta_strong_F = theta83; // 3 июля 2016

	// хеш таблица
	bool* this_is_Strong_C_node = new bool[n_a[ilevel - 1] + 1];
	for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
		// инициализация хеш таблицы.
		this_is_Strong_C_node[i_1] = false;
	}


	// четвёртая версия интерполяции.
	// интерполляция Руге и Стубена amg1r5.
	integer ioneStrongC_and_0_StrongF = 0;

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

	// the begining

	integer imax_count_sosed = -1;
	for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode] == true) {
		integer ii1 = row_startA[inode];
		//integer iend_1 = row_startA[Amat.i[ii1] + 1] - 1;
		integer iend_1 = row_startA[inode + 1] - 1;
		if (iend_1 - ii1 > imax_count_sosed) {
			imax_count_sosed = iend_1 - ii1;
		}
	}
	integer *iStrongC_nodes_ind = new integer[imax_count_sosed + 1];

	bool byes_add = false;
	// Быстрое добавление недостающих С узлов.
	if (1) {
		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//doublerealT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;
		integer icount_bad_string = 0;

		

		// Построение пролонгации для узлов которые составляют F nodes.
		// Каждый F-nodes окружён C-nodes.
		for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode] == true) {


			bool badd_amg1r5 = false;
			bool bfound_vneDiag_F_node = false;

			// Найти соседей данного F-node которые C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
			// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
			doublereal sumP = 0.0;
			doublereal SumPall = 0.0;
			integer icount_StronglyF = 0;

			//integer ii1 = BinarySearchAi(Amat, inode, 1 + iadd, nnz_a[ilevel - 1] + iadd);
			integer ii1 = row_startA[inode];
			//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
			//integer iend_1 = row_startA[Amat.i[ii1] + 1] - 1;
			integer iend_marker_position = row_startA[inode + 1] - 1;
			integer iend_1 = iend_marker_position;//row_startA[inode + 1] - 1;

			

			integer isize_iStrongC_nodes_ind = 0;

			if (bpositive_connections) {


#if doubleintprecision == 1
				//printf("inode=%lld n=%lld\n", inode, n_a[ilevel - 1]);
#else
				//printf("inode=%d n=%d\n", inode, n_a[ilevel - 1]);
#endif

				//getchar();


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != inode) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
							//doublereal mcand = fabs(Amat.aij[is0]);
							doublereal mcand = Amat.abs_aij[is0];
							if (mcand > maxelem_threshold) {
								maxelem_threshold = mcand;
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[Amat.i[ii1]];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				
				
				// Раньше было только вычисление this_is_Strong_C_node
				
				for (integer is0 = ii1; (is0 <= iend_1); is0++) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
					// Диагональный узел по умолчанию F типа
						if (this_is_C_node[Amat.j[is0]] == true) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
							if (Amat.abs_aij[is0] > maxelem_threshold_theta) {

								// заносим всех сильных С соседей.
								// Инициализация хеш таблицы.
								this_is_Strong_C_node[Amat.j[is0]] = true;
								sumP += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;
							}
						}
					//}
				}
				
				//if (1&&(!((fabs(sumP) < RealZERO) || (badd_amg1r5)))) 
				{


					for (integer is0 = ii1; (is0 <= iend_1); is0++) {
						
							//if (this_is_C_node[Amat.j[is0]] == false)
							{
								if (this_is_F_node[Amat.j[is0]] == true)
								{

									if (Amat.j[is0] != inode) 
									{

										the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 

										bfound_vneDiag_F_node = true;
										// F node
										// Подсчитываем количество соседей которые не являются С узлами.
									
										if (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F) {

											integer ii2 = row_startA[Amat.j[is0]];
											integer iend_marker_position2 = row_startA[Amat.i[ii2] + 1] - 1;
											bool bfound = false;
											for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
												if (Amat.j[is2] != Amat.i[ii2]) {
													if (this_is_Strong_C_node[Amat.j[is2]]) {
														bfound = true;
													}
												}
											}

											if (bfound) {
												// Заносим только тех сильных F соседей которые имеют в соседях(любых не обязательно сильных) первоначальных сильных С соселей.

												//SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
												SumPall += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
												icount_StronglyF++;
												icsosF++;
											}
											else {
												badd_amg1r5 = true;
											}
										}
									}
								}
							}
					}




					// Сброс - очистка хеш таблицы.
						/*//253 такта
					for (integer is0 = ii1; (is0 <= iend_1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (this_is_C_node[Amat.j[is0]] == true) {
								///if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
									this_is_Strong_C_node[Amat.j[is0]] = false;
								}
							}
						}
					}
					*/
					// Сброс - очистка хеш таблицы.
					for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
					}
					

					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																			   // Следующий случай назван "безнадёжным".
																			   // Безнадёжный случай оставляется без изменений.
																			   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																			   // интерполляцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
					}
				}
			}
			else {
				// only negative connections


				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != inode) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
							if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold)) {
								//maxelem_threshold = fabs(Amat.aij[is0]);
								maxelem_threshold = Amat.abs_aij[is0];
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_only_negative[Amat.i[ii1]];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.


				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// Инициализация хеш таблицы.
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					//if (Amat.j[is0] != inode) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat.j[is0]] = true;
								//sumP += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								sumP += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;
							}
						}
					//}
				}

				//if (1 && (!((fabs(sumP) < RealZERO) || (badd_amg1r5)))) 
				{

					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						
						if (this_is_F_node[Amat.j[is0]] == true) {

							if (Amat.j[is0] != inode) {

								bfound_vneDiag_F_node = true;

								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F)) {


									integer ii2 = row_startA[Amat.j[is0]];
									integer iend_marker_position2 = row_startA[Amat.i[ii2] + 1] - 1;
									bool bfound = false;
									for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
										if (Amat.j[is2] != Amat.i[ii2]) {
											if (this_is_Strong_C_node[Amat.j[is2]]) {
												bfound = true;
											}
										}
									}

									if (bfound) {

										// Заносим только тех сильных F соседей которые имеют в соседях(любых не обязательно сильных) первоначальных сильных С соселей.
										//SumPall += fabs(Amat.aij[is0]); // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										SumPall += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
									else {
										badd_amg1r5 = true;
									}
								}
								// Подсчитываем количество соседей которые не являются С узлами.
								the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
							}


						}
					}

					// Сброс - очистка хеш таблицы.
					/*
					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (this_is_C_node[Amat.j[is0]] == true) {
								if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
									this_is_Strong_C_node[Amat.j[is0]] = false;
								}
							}
						}
					}
					*/
					// Сброс - очистка хеш таблицы.
					for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
					}
					

					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																			   // Следующий случай назван "безнадёжным".
																			   // Безнадёжный случай оставляется без изменений.
																			   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																			   // интерполляцию на этих уровнях можно улучшить.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
					}
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


			// Мы делаем из F узла C узел только в том случае если у нас нету С связей.
			
			// 20 jan 2016.
			// Можно не добавлять узел в случае если есть не менее трёх сильных F соседей.
			/// Код 
			//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				//if (this_is_F_node[Amat.j[is0]] == true) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
						// Есть внедиагональный F узел.
			
					//}
				//}
			//}//end
			// Эквивалентен if (bfound_vneDiag_F_node&&((

			//236 тактов
			if (bfound_vneDiag_F_node&&((fabs(sumP) < RealZERO) || (badd_amg1r5))) {
				// Это случай когда вообще нет сильных С соседей.

				// if (icsosF==0) полная потеря сходимости.
				//if (icsosF<2)// наблюдается замедление скорости сходимости и времени расчёта.
				{
					//printf("error interpolation zero diagonal sumP.\n");
					//printf("Fnode all sosed is F");
					//system("pause");
					//printf("inode is Dirichlet node\n");
					if (this_is_C_node[inode] == false) iadditionalCstatistic++;
						this_is_F_node[inode] = false; // Этот узел Дирихле станет С нодом.
						this_is_C_node[inode] = true;
						bweSholdbeContinue = true;
						byes_add = true; // Были добавления узлов.
						 //exit(1);
						 // здесь нужна непрямая интерполляция.
						icount_bad_string++;
						
				}
			}
			
			

			



		}
		//if (icount_bad_string>0) {
			//printf("Algorithm s vozvratom %lld %lld\n", icount_bad_string, n_a[ilevel - 1]);
			//system("PAUSE");
		//}

		

	}





	/*
	9.05.2018
	1. Определили thresholdы.
	2. Собрали суммы. SumP - для strong C, SumPall - для strong F.
	3. icsos - количество strong C. icsosF - количество strong F.
	4. Если SumP ==0 то узел становится С узлом. iadditionalCstatistic++;
	0;
	47; 0;
	186; 0;
	267; 0;
	296; 0;
	276; 0;
	221; 0;
	257.6; 0;
	335.2; 0;
	11L 3min 10s
	*/

	if (!byes_add) {

		// с параметром 0.4 четвертая интерполяция обгоняет более сильную шестую.
		//--->doublereal magic = 0.4; // 0.4 optimum
								//magic = 0.3; // 3 июля 2016 для конвективных задач
								// линейный прогноз не дает улучшений
								// ни по быстрдействию ни по числу V циклов.
								//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublereal magic = magic82;// 0.4 is recomended.



		the_number_of_neighbors_that_are_not_С_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		

		if (bpositive_connections) {

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode] == true) {

				integer isize_iStrongC_nodes_ind = 0;

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, inode, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[inode];
				integer istr_etalon1 = inode;// Amat.i[ii1];
				//integer iend_for1 = -1;
				integer iend_1 = row_startA[istr_etalon1 + 1] - 1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_1); is0++) {
						//iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold) {
								//maxelem_threshold = fabs(Amat.aij[is0]);
							//}
							if (Amat.abs_aij[is0] > maxelem_threshold) {
								maxelem_threshold = Amat.abs_aij[is0];
							}
							//}
						}
					}
				}
				else {
					//for (integer is0 = ii1; (is0 <= iend_1); is0++) {
						//iend_for1 = is0;
					//}
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				for (integer is0 = ii1; (is0 <= iend_1); is0++) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
					// мы знаем априорно что диагональный узел у нас F node
						if (this_is_C_node[Amat.j[is0]] == true) {
							if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat.j[is0]] = true;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;
								sumP += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
					//}
				}

				/*
				if (fabs(sumP) < RealZERO) {
					//printf("error interpolation zero diagonal sumP.\n");
					printf("Fnode all sosed is F");
					system("pause");
					//	printf("inode is Dirichlet node\n");
					if (this_is_C_node[inode] == false) iadditionalCstatistic++;
					this_is_F_node[inode] = false; // Этот узел Дирихле станет С нодом.
					this_is_C_node[inode] = true;
					bweSholdbeContinue = true;
					//exit(1);
					// здесь нужна непрямая интерполляция.
				}
				*/

				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_1; is0++) {
					if (this_is_F_node[Amat.j[is0]] == true) {				

						//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
						if (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F) {
							// мы знаем априорно что диагональный узел у нас F node
							if (Amat.j[is0] != inode) {

								integer ii2 = row_startA[Amat.j[is0]];									
								integer inode8 = Amat.i[ii2];
								integer iend_marker_position2 = row_startA[inode8 + 1] - 1;
								bool bfound = false;
								for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
									if (Amat.j[is2] != inode8) {
										if (this_is_Strong_C_node[Amat.j[is2]]) {
												bfound = true;
											}
										}
									}
									if (bfound) {


										SumPall += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}

									// Подсчитываем количество соседей которые не являются С узлами.
									the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
								}								
							}						
						
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным  С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						
						//for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
							//	integer is0 = iStrongC_nodes_ind[is00];						
						//}

						for (integer is0 = ii1; (is0 <= iend_1); is0++) {
							
								
								// Диагональный узел F типа поэтому проверка на внедиагональность не производится
								if (this_is_C_node[Amat.j[is0]] == true) {
									
									// Нас интересуют только Strong связи.
									if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
										
										// Это как раз единственная Strong C связь. 
										// Внедиагональный элемент из множества С узлов.

										// Данная вставка должна существенно сохранять 
										// разреженность оператора Галёркина на глубоких 
										// сеточных уровнях.
										// Модификация 5 декабря 2015.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
										//if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
											P[icount1].j = inode;
											P[icount1].i = C_numerate[Amat.j[is0]];
											//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
											if (fabs(SumPall) < 1.0e-33) {
												printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
												//getchar();
												system("PAUSE");
												exit(1);
											}
											P[icount1].aij = Amat.abs_aij[is0] / SumPall;
											icount1++;
											if (icount1 >= nsizePR * n) {
												printf("memory error!!!\n");
												printf("not enough memory for the interpolation operator.\n");
												//system("PAUSE");
												//exit(1);
												deallocate_prolongation(nsizePR, n, P);
											}
										//}
									}
									
								}
								else
									if (this_is_F_node[Amat.j[is0]] == true) {
										

										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F) {

											if (Amat.j[is0] != inode) {
											// Рассматриваем Strong F связи.

											// Смысл :
											// 

											integer ii2 = row_startA[Amat.j[is0]];
											integer inode88 = Amat.i[ii2];
											integer iend_marker_position2 = row_startA[inode88 + 1] - 1;
											bool bfound = false;
											for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
												if (Amat.j[is2] != inode88) {
													if (this_is_Strong_C_node[Amat.j[is2]]) {
														bfound = true;
													}
												}
											}
											if (bfound) {

												//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
												// Это лишнее условие, ранее уже было проверено что мы имеем дело
												// с сильной F связью.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

												// Диагональный именно F узел. 
												//Если нам нужен С узел то можно не проверять на внедиагональность.
												integer iFpoint = Amat.j[is0];
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												doublereal multiplyer_nu = Amat.abs_aij[is0] / SumPall;
												// сканируем всех соседей iFpointeger 
												// чтобы найти среди них С узлы.

												// Идентифицируем локальный порог.
												doublereal maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												//integer ii1_loc = ii2;// row_startA[iFpoint];
												//integer istr_etalon = iFpoint;//  Amat.i[ii1_loc];
												//integer iend_for = iend_marker_position2;// -1;

												// Найти соседей данного F-node которые C-node.
												integer icsos_loc = 0;

												// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
												// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
												doublereal sumP_loc = 0.0;

												//integer iend_marker_position_loc = iend_marker_position2;// row_startA[iFpoint + 1] - 1;
												//for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												for (integer is0_loc = ii2; (is0_loc <= iend_marker_position2); is0_loc++) {
													//iend_for = is0_loc;
													//if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
													// Важно что threshold должен быть определен только по С узлам.
													// Иначе потеря сходимости.
													if (this_is_C_node[Amat.j[is0_loc]] == true) {
														if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
															// На диагонали F узел. Теперь можно не проверять на внедиагональность. 
															//if (Amat.j[is0_loc] != iFpoint) {
															if (Amat.abs_aij[is0_loc] > maxelem_threshold_loc) {
																maxelem_threshold_loc = Amat.abs_aij[is0_loc];
															}
																//sumP_loc += Amat.abs_aij[is0_loc]; // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																//icsos_loc++;
															//}
														}
													}
													//}
												}

												doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

												
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii2; is0_loc <= iend_marker_position2; is0_loc++) {

													// Будем отсекать сразу наибольшее число вариантов.
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

													if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
														if (this_is_C_node[Amat.j[is0_loc]] == true) {

															// На диагонали F узел. Теперь можно не проверять на внедиагональность.
															//if (Amat.j[is0_loc] != iFpoint) {

																sumP_loc += Amat.abs_aij[is0_loc]; // сумма модулей внедиагональных элементов которые принадлежат С узлам.
																icsos_loc++;
															//}

														}
														else {

															//if (Amat.j[is0_loc] != iFpoint) {
															// Подсчитываем количество соседей которые не являются С узлами.
															//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
															//}
														}
													}

													}
												}
												

												//doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												// В общем это двойная интерполляционная сумма 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii2; is0_loc <= iend_marker_position2; is0_loc++) {

													// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


													// Внедиагональный элемент из множества С узлов.

													// Данная вставка должна существенно сохранять 
													// разреженность оператора Галёркина на глубоких 
													// сеточных уровнях.
													// Модификация 5 декабря 2015.
													if (Amat.abs_aij[is0_loc] > maxelem_threshold_loc_magic) 
													{// нельзя это убирать - сходимость портится
													if (this_is_Strong_C_node[Amat.j[is0_loc]]) 
													{ // Обязательная проверка

														//if (this_is_C_node[Amat.j[is0_loc]] == true)
														{// можно убирать - сходимость НЕ портится
															//if (Amat.j[is0_loc] != iFpoint) {

																P[icount1].j = inode;
																P[icount1].i = C_numerate[Amat.j[is0_loc]];
																//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																if (fabs(sumP_loc) < 1.0e-33) {
																	printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																P[icount1].aij = multiplyer_nu * Amat.abs_aij[is0_loc] / sumP_loc;
																icount1++;
																if (icount1 >= nsizePR * n) {
																	printf("memory error!!!\n");
																	printf("not enough memory for the interpolation operator.\n");
																	//system("PAUSE");
																	//exit(1);
																	deallocate_prolongation(nsizePR, n,  P);
																}

															//}
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

				/*
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat.j[is0]] = false;
							}
						}
					}
				}
				*/
				// Сброс - очистка хеш таблицы.
				for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
					this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
				}


			}
			////
		}
		else {
			// only negative connections.

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode] == true) {

				// Нам нужен разреженный оператор Галёркина.
				// 5 декабря 2015 года мы попробуем увеличить разреженность
				// оператора интерполляции а значит и оператора Галёркина.
				doublereal maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, inode, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[inode];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = row_startA[istr_etalon1 + 1] - 1;// -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							// Если закоментировано то смотрится максимальный внедиагональный элемент в строке.
							//if (this_is_C_node[Amat.j[is0]] == true) {
							if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold)) {
								maxelem_threshold = Amat.abs_aij[is0];
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_only_negative[istr_etalon1];
				}
				// Здесь maxelem_threshold это модуль максимального внедиагонального элемента в строке среди С соседей.

				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublereal maxelem_threshold_theta = maxelem_threshold * theta;
				doublereal maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
					// Диагональный узел по умолчанию F типа
						if (this_is_C_node[Amat.j[is0]] == true) {
							if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat.j[is0]] = true;
							}
						}
					//}
				}


				// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
				// Для каждого такого члена суммы увличиваем счётчик iscos. По идее iscos должно быть 2 и более.
				doublereal sumP = 0.0;
				doublereal SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								sumP += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly С узлам.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]] == true) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F)) {

									integer ii2 = row_startA[Amat.j[is0]];
									integer iend_marker_position2 = row_startA[Amat.i[ii2] + 1] - 1;
									bool bfound = false;
									for (integer is2 = ii2; (is2 <= iend_marker_position2); is2++) {
										if (Amat.j[is2] != Amat.i[ii2]) {
											if (this_is_Strong_C_node[Amat.j[is2]]) {
												bfound = true;
											}
										}
									}
									if (bfound) {


										SumPall += Amat.abs_aij[is0]; // сумма модулей внедиагональных элементов которые принадлежат Strongly F узлам.
										icount_StronglyF++;
										icsosF++;
									}
								}
							}
							// Подсчитываем количество соседей которые не являются С узлами.
							the_number_of_neighbors_that_are_not_С_nodes++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
																		   // Следующий случай назван "безнадёжным".
																		   // Безнадёжный случай оставляется без изменений.
																		   // К счастью на глубоких уровнях сетки безнадёжный случай составляет лишь небольшой процент и 
																		   // интерполляцию на этих уровнях можно улучшить.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // Количество F узлов с одним единственным сильным  C соседом и к тому-же не имеющих сильных F соседов.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// Случай одного Strong C соседа в добавок к которому есть покрайней мере один Strong F сосед.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// Нас интересуют только Strong связи.


								if (this_is_C_node[Amat.j[is0]] == true) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("inode is Dirichlet node\n");
											if (this_is_C_node[inode] == false) iadditionalCstatistic++;
											this_is_F_node[inode] = false; // Этот узел Дирихле станет С нодом.
											this_is_C_node[inode] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// здесь нужна непрямая интерполляция.
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
											if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
												P[icount1].j = inode;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < 1.0e-23) {
													printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = Amat.abs_aij[is0] / SumPall;
												icount1++;
												if (icount1 >= nsizePR * n) {
													printf("memory error!!!\n");
													printf("not enough memory for the interpolation operator.\n");
													//system("PAUSE");
													//exit(1);
													deallocate_prolongation(nsizePR, n,  P);
												}
											}

										}
									}

								}
								else
									if (this_is_F_node[Amat.j[is0]] == true) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F)) {
											// Рассматриваем Strong F связи.

											// Смысл :
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// Это лишнее условие, ранее уже было проверено что мы имеем дело
											// с сильной F связью.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											doublereal multiplyer_nu = Amat.abs_aij[is0] / SumPall;
											// сканируем всех соседей iFpointeger 
											// чтобы найти среди них С узлы.

											// Идентифицируем локальный порог.
											doublereal maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												//if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												if (this_is_C_node[Amat.j[is0_loc]] == true) {
													if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = Amat.abs_aij[is0_loc];
														}
													}
												}
												//}
											}

											doublereal maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// Здесь maxelem_threshold_loc это модуль максимального внедиагонального элемента в строке среди С соседей локально.

											// Найти соседей данного F-node которые C-node.
											integer icsos_loc = 0;

											// Суммируем модули внедиагональных элементов из С узлов которые больше порога.
											// Для каждого такого члена суммы увличиваем счётчик iscos_loc. По идее iscos_loc должно быть 2 и более.
											doublereal sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// Будем отсекать сразу наибольшее число вариантов.
												//if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat.j[is0_loc]]) {


													if (this_is_C_node[Amat.j[is0_loc]] == true) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += Amat.abs_aij[is0_loc]; // сумма модулей внедиагональных элементов которые принадлежат С узлам.
															icsos_loc++;
														}

													}
													else {

														//if (Amat.j[is0_loc] != istr_etalon) {
														// Подсчитываем количество соседей которые не являются С узлами.
														//the_number_of_neighbors_that_are_not_С_nodes_loc++; // подсчитываем проблемы интерполяции 
														//}
													}
												}
												//}
											}

											doublereal maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// В общем это двойная интерполляционная сумма 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// При многоуровневом сравнении в первую очередь надо отсекать наибольшее число вариантов.


												// Внедиагональный элемент из множества С узлов.

												// Данная вставка должна существенно сохранять 
												// разреженность оператора Галёркина на глубоких 
												// сеточных уровнях.
												// Модификация 5 декабря 2015.
												//if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat.j[is0_loc]]) {


													if (this_is_C_node[Amat.j[is0_loc]] == true) {
														if (Amat.j[is0_loc] != istr_etalon) {

															P[icount1].j = inode;
															P[icount1].i = C_numerate[Amat.j[is0_loc]];
															//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
															if (fabs(sumP_loc) < 1.0e-23) {
																printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																//getchar();
																system("PAUSE");
																exit(1);
															}
															P[icount1].aij = multiplyer_nu * Amat.abs_aij[is0_loc] / sumP_loc;
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
												//}
											}


										}
									}
							}

						}
					}

				}

				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]] == true) {
							if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat.j[is0]] = false;
							}
						}
					}
				}

			} // end only negative connections

		}

	}

#if doubleintprecision == 1
	//printf("one Strong C and 0 Strong F=%lld additional all %lld", ioneStrongC_and_0_StrongF, iadditionalCstatistic);
#else
	//printf("one Strong C and 0 Strong F=%d additional all %d", ioneStrongC_and_0_StrongF, iadditionalCstatistic);
#endif

	//system("pause");

	delete[] iStrongC_nodes_ind;
	iStrongC_nodes_ind = NULL;

	delete[] this_is_Strong_C_node;
	this_is_Strong_C_node = NULL;

} // my_interpolation_procedure_number3B

#endif /*MY_AMG_INTERPOLATION_CPP*/
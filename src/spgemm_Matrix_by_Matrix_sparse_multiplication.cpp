// Алгоритмы умножения разреженной матрицы на разреженную матрицу.
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _SPGEMM_MATRIX_BY_MATRIX_SPARSE_MULTIPLICATION_CPP_
#define _SPGEMM_MATRIX_BY_MATRIX_SPARSE_MULTIPLICATION_CPP_ 1


// Обёртка для сишного realloc().
// В данной функции осуществляется проверка на nullptr.
// Копирование результата умножения двух разреженных матриц в исходную матрицу А 
// Может вызвать переполнение хранилища отводимого под матрицу А.
// В случае обнаружения переполнения делается попытка увеличить размер хранилища
// отводимого под хранение матрицы А.
// 07.08.2020
template <typename myARRT>
void my_realloc_memory1(myARRT*& x, integer n) {
	if (x != nullptr) {
		myARRT* x_temp = nullptr;
		x_temp = (myARRT*)realloc(x, ((n) * sizeof(myARRT)));
		if (x_temp == nullptr) {
			std::cout << "ERROR!!! Amat index overflow in AP parallel 8 multiplyer.\n" << std::endl;
			std::cout << "application crash for bx array in procedure my_realloc_memory1." << std::endl;
			std::cout << "Recomended: Increase the size of the memory allocated for storage of matrix A..." << std::endl;
			system("pause");
			exit(1);
		}
		else {
			x = x_temp;
		}
	}
}

// Разреженное матричное умножение двух разреженных матриц.
// Алгоритм Ф. Густавсона IBM 1978.
template <typename doublerealT>
void my_parallel8_sparse_matrix_by_matrix_multiplication_RA(Ak2& Amat,
	Ak1*& P, integer& istartAnew, integer*& istartAnew_m,
	integer*& row_ind_SR,
	integer*& row_ind_ER,
	integer*& row_startA,
	integer numberofcoarcenodes,
	const int iKnumber_thread,
	bool**& hash_table_m,
	integer**& index_visit_m,
	doublerealT**& vector_sum_m,
	integer*& index_size_m,
	integer& nsizeA,
	integer n, integer nnz, integer ilevel,
	bool bprint_mesage_diagnostic,
	integer*& n_a, doublereal AccumulqtorA_m_SIZE8,
	Ak1**& AccumulqtorA_m
	//, integer istartAnew_mem
) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	

	integer inum_core = number_cores();
	omp_set_num_threads(inum_core);
	

	const integer iSIZE =  n + 2;

	// Данные используемые для частичного формирователя суммы.

	for (int i_9 = 0; i_9 < iKnumber_thread; i_9++) {

#pragma omp parallel for schedule (static)
		for (integer i_91 = 0; i_91 < iSIZE; i_91++) {
			hash_table_m[i_9][i_91] = false;// inicialization
		}
		index_size_m[i_9] = 0;
		istartAnew_m[i_9] = 0;
	}

	const integer LIMIT_ANEW_m = (integer)(0.125 * AccumulqtorA_m_SIZE8 * nnz + 1);

	omp_set_num_threads(8); // только 8 потоков.

	// Сканируем первый операнд построчно.
	// глобальные переменные не перечисляются.
#pragma omp parallel for schedule (static) 
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

#ifdef _OPENMP 
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif

		// на основе хеш-таблицы. 

		// Сканируем текущую i-ую строку поэлементно
		for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
			integer col_ind = P[ii].j;
			// Сканируем col_ind строку второго операнда

			// Общую переменную объявим на уровень выше.
			doublerealT left_operand = P[ii].aij;

			integer i_1start = row_startA[col_ind];
			integer i_1end = row_startA[col_ind + 1] - 1;
			for (integer i_1 = i_1start; i_1 <= i_1end; i_1++) {

				doublerealT right_operand = Amat.aij[i_1];
				integer iaddind = Amat.j[i_1];

				if (hash_table_m[tid][iaddind]) {

					vector_sum_m[tid][iaddind] += left_operand * right_operand;
				}
				else {
					// Первое добавление.

					index_size_m[tid]++;

					index_visit_m[tid][index_size_m[tid]] = iaddind;

					hash_table_m[tid][iaddind] = true;

					vector_sum_m[tid][iaddind] = left_operand * right_operand;

				}
			}
		}
		
		

		integer ism = index_size_m[tid];
		if (nsizeA > istartAnew + ism) {
			if (istartAnew_m[tid] + ism < LIMIT_ANEW_m) {
				for (integer i_6 = 1; i_6 <= ism; i_6++) {
					integer jstr = index_visit_m[tid][i_6];

					hash_table_m[tid][jstr] = false; // инициализируем хеш-таблицу для следующих проходов.

					doublerealT vs1 = vector_sum_m[tid][jstr];

					// 7 ноября 2016 игнорируем чистые нули:
					if (fabs(vs1) > 1.0e-37) {
						// Мы не записываем в матрицу чистый ноль.
						Ak1 Atemp;
						Atemp.aij = vs1;
						Atemp.i = istr;
						Atemp.j = jstr;


						AccumulqtorA_m[tid][istartAnew_m[tid]++] = Atemp;
					}
				}
			}
			else {
				printf("AccumulqtorA_m index overflow\n");
				system("pause");
			}
		}
		else {
			// Слишком мало памяти выделено под матрицу А и в неё не умещаются все данные.
			// Нужно увеличить объём выделяемой памяти для А и перезапустить приложение.
			printf("Amat lack of memory\n");
			printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			printf("please, press any key to exit.\n");
			system("pause");
			exit(1);
		}



		index_size_m[tid] = 0; // Сброс индикатора, строка обработана.			

	}


	omp_set_num_threads(inum_core);
	

	//if (bprint_mesage_diagnostic) {
		//printf("oK. Counting Sort start.\n");
	//}
	for (int i_9 = 0; i_9 < iKnumber_thread; i_9++)
	{
		if (istartAnew_m[i_9] - 1 < (integer)(0.125 * AccumulqtorA_m_SIZE8 * nnz + 1)) {
			if (nsizeA > istartAnew + istartAnew_m[i_9] + 1) {
#pragma omp parallel for
				for (integer i_92 = 0; i_92 < istartAnew_m[i_9]; i_92++) {
					integer iz = istartAnew + i_92;
					Amat.aij[iz] = AccumulqtorA_m[i_9][i_92].aij;
					Amat.i[iz] = AccumulqtorA_m[i_9][i_92].i;
					Amat.j[iz] = AccumulqtorA_m[i_9][i_92].j;
				}
				istartAnew += istartAnew_m[i_9];
			}
		    else {
			   // 07.08.2020
			   // Нехватка памяти под хранение матрицы А.
			   my_realloc_memory1<real_mix_precision>(Amat.aij, istartAnew + 20 * istartAnew_m[i_9] + 1);
			   my_realloc_memory1<integer_mix_precision>(Amat.i, istartAnew + 20 * istartAnew_m[i_9] + 1);
			   my_realloc_memory1<integer_mix_precision>(Amat.j, istartAnew + 20 * istartAnew_m[i_9] + 1);
			   nsizeA = istartAnew + 20 * istartAnew_m[i_9];
			   //printf("ERROR!!! Amat index overflow in RA parallel 8 multiplyer.\n");
			   //system("PAUSE");
			   //exit(1);
#pragma omp parallel for
			   for (integer i_92 = 0; i_92 < istartAnew_m[i_9]; i_92++) {
				   integer iz = istartAnew + i_92;
				   Amat.aij[iz] = AccumulqtorA_m[i_9][i_92].aij;
				   Amat.i[iz] = AccumulqtorA_m[i_9][i_92].i;
				   Amat.j[iz] = AccumulqtorA_m[i_9][i_92].j;
			   }
			   istartAnew += istartAnew_m[i_9];

		   }
			
		}
		else {
			printf("AccumulqtorA_m index overflow\n");
			system("pause");
		}
	}

	// Сортировка не нужна т.к. используется schedule (static)
	//Counting_Sort(Amat, istartAnew_mem, istartAnew - 1, false, n_a[ilevel - 1]);
	//if (bprint_mesage_diagnostic) {
		//printf("Counting Sort End. \n");
	//}

	//getchar();

	omp_set_num_threads(i_my_num_core_parallelesation);

#endif

}

// Разреженное матричное умножение двух разреженных матриц.
// Алгоритм Ф. Густавсона IBM 1978.
template <typename doublerealT>
void my_parallel8_sparse_matrix_by_matrix_multiplication_AP(Ak2& Amat,
	Ak1 const *const P, integer& istartAnew, integer*& istartAnew_m,
	integer const *const row_ind_AS,
	integer const *const row_ind_AE,
	integer const *const row_ind_PS,
	integer const *const row_ind_PE,
	integer numberofcoarcenodes,
	const int iKnumber_thread,
	bool**& hash_table_m,
	integer**& index_visit_m,
	doublerealT**& vector_sum_m,
	integer*& index_size_m,
	integer& nsizeA,
	integer n, integer nnz, integer ilevel,
	bool bprint_mesage_diagnostic,
	integer*& n_a, doublereal AccumulqtorA_m_SIZE8,
	Ak1**& AccumulqtorA_m, integer& istartAnew2,
	integer& ibsp_length,
	integer i_bsp_LIMIT,
	BAD_STRING_PATCHING*& bsp,
	integer iVar,
	bool bcontinue_global)
{

#ifdef _OPENMP

	int i_my_num_core_parallelesation = omp_get_max_threads();

	integer inum_core = number_cores();
	omp_set_num_threads(inum_core);
	
	// Данные используемые для частичного формирователя суммы.

	for (int i_9 = 0; i_9 < iKnumber_thread; i_9++) {

		integer ISIZE =  n + 2;

#pragma omp parallel for
		for (integer i_91 = 0; i_91 < ISIZE; i_91++) {
			hash_table_m[i_9][i_91] = false;// inicialization
		}
		index_size_m[i_9] = 0;
		istartAnew_m[i_9] = 0;
	}


	// Мы будем сканировать левый операнд построчно, а
	// после окончания обработки одной строки левого операнда
	// получать готовую строку результата.

	const integer LIMIT_ANEW_m = (integer)(0.125 * AccumulqtorA_m_SIZE8 * nnz + 1);


	
	omp_set_num_threads(8); // только 8 потоков.

	// Сканируем первый операнд построчно.
	// глобальные переменные не перечисляются.
#pragma omp parallel for schedule (static) 
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

#ifdef _OPENMP 
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif

		// На основе хеш-таблицы.
		// сканируем все элементы строки левого операнда.
		for (integer ii1 = row_ind_AS[istr]; ii1 <= row_ind_AE[istr]; ii1++) {

			integer col_ind = Amat.j[ii1];
			doublerealT left_operand = Amat.aij[ii1];

			// Сканируем col_ind строку правого операнда накапливая сумму.
			integer i_1start = row_ind_PS[col_ind];
			integer i_1end = row_ind_PE[col_ind];

			
			
			for (integer ii2 = i_1start; ii2 <= i_1end; ii2++) {


				doublerealT right_operand = P[ii2].aij;

				integer iaddind = P[ii2].i;

				if (hash_table_m[tid][iaddind]) {

					vector_sum_m[tid][iaddind] += left_operand * right_operand;

				}
				else {

					// Первое добавление.
					index_size_m[tid]++;
					index_visit_m[tid][index_size_m[tid]] = iaddind;

					// Мгновенная вставка в hash table за O(1).
					hash_table_m[tid][iaddind] = true;

					//ifoundind = index_size;
					//vector_sum[index_visit[ifoundind]] = left_operand*right_operand;
					vector_sum_m[tid][iaddind] = left_operand * right_operand;


				}
				// требуется реализовать следующую логику:
				// 1. поиск элемента по ключу 
				// 2. если элемент не найден то добавление нового узла со значением ключа сохраняя балансировку.
				// Если элемент найден то нужно просто изменить foundnow на true. 
				// Т.е. достаточно просто поиска и вставки.
				// 3. В конце дерево необходимо ликвидировать.
				// Тип данных целочисленный ключ.


				//vector_sum[P[ii2].i] += rleft*rright;
			}
			
		}

		if (my_amg_manager.bdiagonal_dominant > 0) {
			doublereal additional_to_diagonal = 0.0;

			for (integer i_6 = 1; i_6 <= index_size_m[tid]; i_6++) {
				integer jstr = index_visit_m[tid][i_6];
				if ((istr != jstr) && (vector_sum_m[tid][jstr] >= 0.0)) {
					additional_to_diagonal += vector_sum_m[tid][jstr];
					vector_sum_m[tid][jstr] = 0.0;
					hash_table_m[tid][jstr] = false;
				}
			}
			vector_sum_m[tid][istr] += additional_to_diagonal;
		}

		//doublerealT maxth = -1.0;


		// huck: 16.04.2017
		integer ism = index_size_m[tid];
		for (integer i_61 = 1; i_61 <= ism; i_61++) {

			integer jstr61 = index_visit_m[tid][i_61];

			hash_table_m[tid][jstr61] = false; // initialization hash.

			doublerealT vs161 = vector_sum_m[tid][jstr61];

			//if (istr != jstr61) {
				// 14 января 2016 года.
				// Правильно определить барьер только по внедиагональным элементам.
				//if (fabs(vector_sum_m[tid][jstr61]) > maxth) maxth = fabs(vector_sum_m[tid][jstr61]);
			//}

	//#if doubleintprecision == 1
			//printf("i=%lld j=%lld aij=%e\n", istr, jstr61, vs161);
	//#else
			//printf("i=%d j=%d aij=%e\n", istr, jstr61, vs161);
	//#endif


			if ((istr == jstr61) && (vs161 < 1.0e-20)) {
				// отрицательный элемент на диагонали.

				printf("Negative diagonal coefficient found. No panic. Upwind patching. 16.04.2017. \n");

				printf("bad string: \n");
				for (integer i_63 = 1; i_63 <= index_size_m[tid]; i_63++) {

					integer jstr63 = index_visit_m[tid][i_63];
					doublerealT vs163 = vector_sum_m[tid][jstr63];
#if doubleintprecision == 1
					printf("i=%lld j=%lld aij=%e\n", istr, jstr63, vs163);
#else
					printf("i=%d j=%d aij=%e\n", istr, jstr63, vs163);
#endif

				}


				// Адаптированные три правила бак-труба:
				// Amat. Диагонали присваиваем сумма модулей только отрицательных внедиагональных коэффициентов +
				// вычитаем из этого отрицательную диагональ. Потом умножаем на два.
				// B. Удвоение отрицательных внедиагональных коэффициентов.
				// C. Полное зануление положительных внедиагональных коэффициентов (игнорирование).
				printf("patching string 16.04.2017: \n");
				/*
				for (integer i_62 = 1; i_62 <= index_size_m[tid]; i_62++) {
					integer jstr62 = index_visit_m[tid][i_62];
					if (istr != jstr62) {
						if (vector_sum_m[tid][jstr62] > 0.0) {
							index_visit_m[tid][i_62] = -1; // не существует такого элемента (игнорирование).
							vector_sum_m[tid][jstr61] += vector_sum_m[tid][jstr62];
							vector_sum_m[tid][jstr62] = 0.0;
						}
					}
				}

				if (vector_sum_m[tid][jstr61] < 0.0) {
					vector_sum_m[tid][jstr61] = 0.0;
					for (integer i_62 = 1; i_62 <= index_size_m[tid]; i_62++) {
						integer jstr62 = index_visit_m[tid][i_62];
						if (jstr62 > -1) {
							if (istr != jstr62) {
								if (vector_sum_m[tid][jstr62] < 0.0) {
									vector_sum_m[tid][jstr61] += fabs(vector_sum_m[tid][jstr62]);
								}
							}
						}
					}
				}
				*/
				// Denis Demidov recomendation
				// Всю строку домножаем на минус один, так чтобы диагональ стала положительна.
				for (integer i_62 = 1; i_62 <= index_size_m[tid]; i_62++) {
					integer jstr62 = index_visit_m[tid][i_62];
					vector_sum_m[tid][jstr62] *= -1.0;
				}
				// Запоминаем строку с отрицательной диагональю.
				if (ibsp_length < i_bsp_LIMIT) {
					bsp[ibsp_length].ilevel = ilevel;
					bsp[ibsp_length].istring_number = istr;
					ibsp_length++;
				}
				else {
					printf("ERROR!!! ibsp_length>=i_bsp_LIMIT\n");
					system("pause");
					exit(1);
				}

				// Выход из цикла for по переменной i_61.
				break;
			}


		}


		bool bCheck_ok = false; // проверяет наличие диагонали в строке матрицы.
								// Нам нужен чрезвычайно быстрый код поэтому мы принебрегаем проверками.

		

		if (nsizeA > istartAnew2 + ism) {
			if (istartAnew_m[tid]+ism < LIMIT_ANEW_m) {
				for (integer i_6 = 1; i_6 <= ism; i_6++) {

					integer jstr = index_visit_m[tid][i_6];
					//if (jstr > -1)
					{
						doublerealT vs1 = vector_sum_m[tid][jstr];
						//if (fabs(vs1) < 1.0e-37) {
#if doubleintprecision == 1
					//printf("zero vs1=%e, i==%lld j==%lld\n",vs1,istr,jstr);
#else
					//printf("zero vs1=%e, i==%d j==%d\n",vs1,istr,jstr);
#endif

					//}
					// 7 ноября 2016 игнорируем чистые нули:
						if (fabs(vs1) > 1.0e-37) {
							// Мы игнорируем чистые нули. 
							// Но вообще говоря непонятно почему они появляются.


								// алгебраический мультигрид Галёркина.
								// 22_10_2016.
							Ak1 Atemp;
							Atemp.aij = vs1;
							Atemp.i = istr;
							Atemp.j = jstr;

							if (istr == jstr) {

								bCheck_ok = true;

								if (vs1 < 1.0e-20) {
									// Ошибка проявляется в отсутствии диагонального элемента в результирующей матрице первого
									// произведения Галеркина. Надо смотреть ситуацию выше по коду.
									// 22737
									// сканируем все элементы строки левого операнда.
									//for (integer ii1_8 = row_ind_AS[istr]; ii1_8 <= row_ind_AE[istr]; ii1_8++) {
									//if (Amat.i[ii1_8] == 22737) {
#if doubleintprecision == 1
								//printf("i=%lld j=%lld aij=%e\n", Amat.i[ii1_8], Amat.j[ii1_8], Amat.aij[ii1_8]);
#else
								//printf("i=%d j=%d aij=%e\n", Amat.i[ii1_8], Amat.j[ii1_8], Amat.aij[ii1_8]);
#endif

								//}
								//integer col_ind = Amat.j[ii1_8];
								//}
#if doubleintprecision == 1
									printf("bad string %lld\n", istr);
#else
									printf("bad string %d\n", istr);
#endif

									printf("error: diagonal element is negative...\n");
									switch (iVar) {
									case PAM: printf("PAM equation\n"); break;
									case VELOCITY_X_COMPONENT: printf("VELOCITY_X_COMPONENT equation\n"); break;
									case VELOCITY_Y_COMPONENT: printf("VELOCITY_Y_COMPONENT equation\n"); break;
									case VELOCITY_Z_COMPONENT: printf("VELOCITY_Z_COMPONENT equation\n"); break;
									case NUSHA: printf("NU equation\n");  break;
									case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY equation\n");  break;
									case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA equation\n");  break;
									case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS equation\n");  break;
									case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS equation\n");  break;
									case TEMP: printf("TEMP equation\n"); break;
									case TOTALDEFORMATIONVAR: printf("STRESS system equation\n"); break;
									}
#if doubleintprecision == 1
									//for (integer ii1_8 = row_ind_AS[istr]; ii1_8 <= row_ind_AE[istr]; ii1_8++) {
									//printf("i=%lld j=%lld aij=%e\n", Amat.i[ii1_8], Amat.j[ii1_8], Amat.aij[ii1_8]);
									//}
#else
									//for (integer ii1_8 = row_ind_AS[istr]; ii1_8 <= row_ind_AE[istr]; ii1_8++) {
									//printf("i=%d j=%d aij=%e\n", Amat.i[ii1_8], Amat.j[ii1_8], Amat.aij[ii1_8]);
									//}
#endif

									for (integer i_61 = 1; i_61 <= index_size_m[tid]; i_61++) {

										integer jstr61 = index_visit_m[tid][i_61];
										doublerealT vs161 = vector_sum_m[tid][jstr61];
#if doubleintprecision == 1
										printf("i=%lld j=%lld aij=%e\n", istr, jstr61, vs161);
#else
										printf("i=%d j=%d aij=%e\n", istr, jstr61, vs161);
#endif

									}
									//getchar();
									system("pause");
									// прекращаем строить иерархию уровней.
									bcontinue_global = false;
									//goto BAD_STRING_MARKER;
									printf("fatall error bad string: goto BAD_STRING_MARKER;\n");
									system("pause");
									exit(1);

									doublerealT sum_dia = 0.0;
									for (integer i_8 = 1; i_8 <= index_size_m[tid]; i_8++) {
										if (i_8 != i_6) {
											integer jstr_8 = index_visit_m[tid][i_8];
											doublerealT vs1_8 = vector_sum_m[tid][jstr_8];
											sum_dia += fabs(vs1_8);
										}
									}
									// принудительное сильнейшее усиление диагонали.
									Atemp.aij = sum_dia;
									// ошибка признана не являющейся фатальной.
									// 22 декабря 2016
									//system("pause");
								}
							}



							AccumulqtorA_m[tid][istartAnew_m[tid]++] = Atemp;

						}
					}
				}
			}
			else {
				printf("AccumulqtorA_m index overflow\n");
				system("pause");
			}
		}
		else {
			// Слишком мало памяти выделено под матрицу А и в неё не умещаются все данные.
			// Нужно увеличить объём выделяемой памяти для А и перезапустить приложение.
			printf("Amat lack of memory\n");
			printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			printf("please, press any key to exit.\n");
			system("pause");
			exit(1);
		}

		if (!bCheck_ok) {
#if doubleintprecision == 1
			printf("bad string %lld\n", istr);
#else
			printf("bad string %d\n", istr);
#endif

			// прекращаем строить иерархию уровней.
			bcontinue_global = false;
			//goto BAD_STRING_MARKER;
			printf("fatall error bad string: goto BAD_STRING_MARKER 2;\n");
			system("pause");
			exit(1);

			for (integer i_6 = 1; i_6 <= index_size_m[tid]; i_6++) {

				integer jstr = index_visit_m[tid][i_6];
				doublerealT vs1 = vector_sum_m[tid][jstr];
#if doubleintprecision == 1
				printf("%lld %lld %e\n", istr, jstr, vs1);
#else
				printf("%d %d %e\n", istr, jstr, vs1);
#endif

			}
			system("pause");
		}

		index_size_m[tid] = 0;

	}

	omp_set_num_threads(inum_core);
	

	//integer istartAnew_mem2 = istartAnew2; // for Counting_Sort
	//if (bprint_mesage_diagnostic) {
		//printf("oK2. Counting Sort start.\n");
	//}
	for (int i_9 = 0; i_9 < iKnumber_thread; i_9++)
	{
		if (istartAnew_m[i_9] - 1 < (integer)(0.125 * AccumulqtorA_m_SIZE8 * nnz + 1)) {
			if (nsizeA > istartAnew2 + istartAnew_m[i_9] + 1) {
#pragma omp parallel for
				for (integer i_92 = 0; i_92 < istartAnew_m[i_9]; i_92++) {
					integer iz = i_92 + istartAnew2;
					Amat.aij[iz] = AccumulqtorA_m[i_9][i_92].aij;
					Amat.i[iz] = AccumulqtorA_m[i_9][i_92].i;
					Amat.j[iz] = AccumulqtorA_m[i_9][i_92].j;					
				}
				istartAnew2 += istartAnew_m[i_9];
			}
			else {
				// 07.08.2020
				// Нехватка памяти под хранение матрицы А.
				my_realloc_memory1<real_mix_precision>(Amat.aij, istartAnew2 + 20 * istartAnew_m[i_9] + 1);
				my_realloc_memory1<integer_mix_precision>(Amat.i, istartAnew2 + 20 * istartAnew_m[i_9] + 1);
				my_realloc_memory1<integer_mix_precision>(Amat.j, istartAnew2 + 20 * istartAnew_m[i_9] + 1);
				nsizeA = istartAnew2 + 20 * istartAnew_m[i_9];
				//printf("ERROR!!! Amat index overflow in AP parallel 8 multiplyer.\n");
				//system("PAUSE");
				//exit(1);

#pragma omp parallel for
				for (integer i_92 = 0; i_92 < istartAnew_m[i_9]; i_92++) {
					integer iz = i_92 + istartAnew2;
					Amat.aij[iz] = AccumulqtorA_m[i_9][i_92].aij;
					Amat.i[iz] = AccumulqtorA_m[i_9][i_92].i;
					Amat.j[iz] = AccumulqtorA_m[i_9][i_92].j;
				}
				istartAnew2 += istartAnew_m[i_9];
			}
		}
		else {
			printf("AccumulqtorA_m index overflow\n");
			system("pause");
		}
	}

	// Сортировка не нужна т.к. используется schedule (static)
	//Counting_Sort(Amat, istartAnew_mem2, istartAnew2 - 1, false, n_a[ilevel - 1]);
	//if (bprint_mesage_diagnostic) {
		//printf("Counting Sort End. \n");
	//}

	omp_set_num_threads(i_my_num_core_parallelesation);

#endif

}

// Разреженное матричное умножение двух разреженных матриц.
// Алгоритм Ф. Густавсона IBM 1978.
template <typename doublerealT>
void my_sparse_matrix_by_matrix_multiplication_RA_1(Ak2& Amat,
	Ak1*& P, integer& istartAnew,
	integer*& row_ind_SR,
	integer*& row_ind_ER,
	integer*& row_startA,
	integer numberofcoarcenodes,
	bool*& hash_table,
	integer*& index_visit,
	doublerealT*& vector_sum,
	integer& nsizeA) {

	integer index_size = 0;

	// Сканируем первый операнд построчно.
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

		// Начинаем обрабатывать новую строку.
		// Сброс формирователя суммы в ноль.

		node_AVL_Gus* root_Gus = 0;


		// на основе хеш-таблицы. 

		// Сканируем текущую i-ую строку поэлементно
		for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
			integer col_ind = P[ii].j;
			// Сканируем col_ind строку второго операнда

			// Общую переменную объявим на уровень выше.
			doublerealT left_operand = P[ii].aij;

			integer i_1start = row_startA[col_ind];
			integer i_1end = row_startA[col_ind + 1] - 1;
			for (integer i_1 = i_1start; i_1 <= i_1end; i_1++) {

				doublerealT right_operand = Amat.aij[i_1];
				integer iaddind = Amat.j[i_1];
				bool foundnow = false;

				// поиск .
				foundnow = hash_table[iaddind];

				if (foundnow) {
					vector_sum[iaddind] += left_operand * right_operand;
				}
				else {
					// Первое добавление.
					index_size++;
					index_visit[index_size] = iaddind;
					// Вставка 							
					hash_table[iaddind] = true;
					vector_sum[iaddind] = left_operand * right_operand;
				}


			}
		}



		doublerealT maxth = -1.0;
		// 22 октября 2016 Мы искоренили барьер из части P*Amat.
		for (integer i_6 = 1; i_6 <= index_size; i_6++) {
			integer jstr = index_visit[i_6];
			hash_table[jstr] = false; // инициализируем хеш-таблицу для следующих проходов.
		}

		// Нам нужен чрезвычайно быстрый код поэтому мы принебрегаем
		// проверками.
		bool bCheck_ok = false; // проверяет наличие диагонали в строке матрицы.
		integer istartAnew_8 = istartAnew; // запоминаем для вылечивания строки.

		if (nsizeA > istartAnew + index_size) {
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];

				doublerealT vs1 = vector_sum[jstr];

				if ((istr == jstr) && (vs1 > 1.0e-20)) {
					bCheck_ok = true;
				}

				// 22 октября 2016. Полностью искоренён барьер из части P*Amat произведения.


				// 7 ноября 2016 игнорируем чистые нули:
				if (fabs(vs1) > 1.0e-37) {
					// Мы не записываем в матрицу чистый ноль.
					Ak1 Atemp;
					Atemp.aij = vs1;
					Atemp.i = istr;
					Atemp.j = jstr;
					Amat.aij[istartAnew] = Atemp.aij;
					Amat.i[istartAnew] = Atemp.i;
					Amat.j[istartAnew] = Atemp.j;
					istartAnew++;
				}

			}
		}
		else {

			// 07.08.2020
			// Нехватка памяти под хранение матрицы А.
			my_realloc_memory1<real_mix_precision>(Amat.aij, istartAnew + 20 * index_size + 1);
			my_realloc_memory1<integer_mix_precision>(Amat.i, istartAnew + 20 * index_size + 1);
			my_realloc_memory1<integer_mix_precision>(Amat.j, istartAnew + 20 * index_size + 1);
			nsizeA = istartAnew + 20 * index_size;
			//printf("ERROR!!! Amat index overflow in AP parallel 8 multiplyer.\n");
			//system("PAUSE");
			//exit(1);

			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];

				doublerealT vs1 = vector_sum[jstr];

				if ((istr == jstr) && (vs1 > 1.0e-20)) {
					bCheck_ok = true;
				}

				// 22 октября 2016. Полностью искоренён барьер из части P*Amat произведения.


				// 7 ноября 2016 игнорируем чистые нули:
				if (fabs(vs1) > 1.0e-37) {
					// Мы не записываем в матрицу чистый ноль.
					Ak1 Atemp;
					Atemp.aij = vs1;
					Atemp.i = istr;
					Atemp.j = jstr;
					Amat.aij[istartAnew] = Atemp.aij;
					Amat.i[istartAnew] = Atemp.i;
					Amat.j[istartAnew] = Atemp.j;
					istartAnew++;
				}

			}


			// Слишком мало памяти выделено под матрицу А и в неё не умещаются все данные.
			// Нужно увеличить объём выделяемой памяти для А и перезапустить приложение.
			//printf("Amat lack of memory\n");
			//printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			//printf("please, press any key to exit.\n");
			//system("pause");
			//exit(1);
		}

		index_size = 0; // Сброс индикатора, строка обработана.
		clear_AVL_Gus(root_Gus);
		root_Gus = 0;

	}

} // my_sparse_multiplication


// Разреженное матричное умножение двух разреженных матриц.
// Алгоритм Ф. Густавсона IBM 1978.
template <typename doublerealT>
void my_sparse_matrix_by_matrix_multiplication_RA(Ak2& Amat,
	Ak1 const *const P, integer& istartAnew,
	integer const *const row_ind_SR,
	integer const *const row_ind_ER,
	integer const *const row_startA,
	integer numberofcoarcenodes,
	bool*& hash_table,
	integer*& index_visit,
	doublerealT*& vector_sum,
	integer& nsizeA) {

	integer index_size = 0;

	// Сканируем первый операнд построчно.
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

		// Начинаем обрабатывать новую строку.
		// Сброс формирователя суммы в ноль.

		// на основе хеш-таблицы. 

		// Сканируем текущую i-ую строку поэлементно
		for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
			integer col_ind = P[ii].j;
			// Сканируем col_ind строку второго операнда

			// Общую переменную объявим на уровень выше.
			doublerealT left_operand = P[ii].aij;

			integer i_1start = row_startA[col_ind];
			integer i_1end = row_startA[col_ind + 1] - 1;
			for (integer i_1 = i_1start; i_1 <= i_1end; i_1++) {

				doublerealT right_operand = Amat.aij[i_1];
				integer iaddind = Amat.j[i_1];
				bool foundnow = false;

				// поиск .
				foundnow = hash_table[iaddind];

				if (foundnow) {
					vector_sum[iaddind] += left_operand * right_operand;
				}
				else {
					// Первое добавление.
					index_size++;
					index_visit[index_size] = iaddind;
					// Вставка 							
					hash_table[iaddind] = true;
					vector_sum[iaddind] = left_operand * right_operand;
				}
			}
		}

		// Нам нужен чрезвычайно быстрый код поэтому мы принебрегаем
		// проверками.		

		if (nsizeA > istartAnew + index_size) {
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];
				hash_table[jstr] = false; // инициализируем хеш-таблицу для следующих проходов.


				doublerealT vs1 = vector_sum[jstr];

				// 22 октября 2016. Полностью искоренён барьер из части R*Amat произведения.


				// 7 ноября 2016 игнорируем чистые нули:
				if (fabs(vs1) > 1.0e-37) {
					// Мы не записываем в матрицу чистый ноль.
					Amat.aij[istartAnew] = vs1;
					Amat.i[istartAnew] = istr;
					Amat.j[istartAnew] = jstr;
					istartAnew++;
				}

			}
		}
		else {

			my_realloc_memory1<real_mix_precision>(Amat.aij, istartAnew + 20 * index_size + 1);
			my_realloc_memory1<integer_mix_precision>(Amat.i, istartAnew + 20 * index_size + 1);
			my_realloc_memory1<integer_mix_precision>(Amat.j, istartAnew + 20 * index_size + 1);
			nsizeA = istartAnew + 20 * index_size;


			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];
				hash_table[jstr] = false; // инициализируем хеш-таблицу для следующих проходов.


				doublerealT vs1 = vector_sum[jstr];

				// 22 октября 2016. Полностью искоренён барьер из части R*Amat произведения.


				// 7 ноября 2016 игнорируем чистые нули:
				if (fabs(vs1) > 1.0e-37) {
					// Мы не записываем в матрицу чистый ноль.
					Amat.aij[istartAnew] = vs1;
					Amat.i[istartAnew] = istr;
					Amat.j[istartAnew] = jstr;
					istartAnew++;
				}

			}


			// Слишком мало памяти выделено под матрицу А и в неё не умещаются все данные.
			// Нужно увеличить объём выделяемой памяти для А и перезапустить приложение.
			//printf("Amat lack of memory\n");
			//printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			//printf("please, press any key to exit.\n");
			//system("pause");
			//exit(1);
		}

		index_size = 0; // Сброс индикатора, строка обработана.		

	}

} // my_sparse_multiplication


// Разреженное матричное умножение двух разреженных матриц.
// Алгоритм Ф. Густавсона IBM 1978.
template <typename doublerealT>
void my_sparse_matrix_by_matrix_multiplication_AP(Ak2& Amat,
	Ak1*& P, integer& istartAnew2,
	integer*& row_ind_AS,
	integer*& row_ind_AE,
	integer*& row_ind_PS,
	integer*& row_ind_PE,
	integer numberofcoarcenodes,
	bool*& hash_table,
	integer*& index_visit,
	doublerealT*& vector_sum,
	integer& nsizeA,
	integer& ibsp_length,
	integer i_bsp_LIMIT,
	integer ilevel,
	BAD_STRING_PATCHING*& bsp)
{
	integer index_size = 0;

	// Левый операнд сканируется построчно, а
	// после окончания обработки одной строки левого операнда
	// получается готовая строка результата.

	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

		node_AVL_Gus* root_Gus = 0;

		// На основе хеш-таблицы.
		// сканируем все элементы строки левого операнда.
		integer row_ind = Amat.i[row_ind_AS[istr]];// номер строки.
		for (integer ii1 = row_ind_AS[istr]; ii1 <= row_ind_AE[istr]; ii1++) {

			integer col_ind = Amat.j[ii1];
			doublerealT left_operand = Amat.aij[ii1];

			// Сканируем col_ind строку правого операнда накапливая сумму.
			for (integer ii2 = row_ind_PS[col_ind]; ii2 <= row_ind_PE[col_ind]; ii2++) {

				doublerealT right_operand = P[ii2].aij;

				integer iaddind = P[ii2].i;
				bool foundnow = false;

				// мгновенный поиск за O(1).
				foundnow = hash_table[iaddind];

				if (foundnow) {

					// Диагональный элемент
					vector_sum[iaddind] += left_operand * right_operand;


				}
				else {
					// Первое добавление.


					index_size++;
					index_visit[index_size] = iaddind;
					// Вставка
					// Мгновенная вставка в хеш-таблицу за O(1).
					hash_table[iaddind] = true;

					vector_sum[iaddind] = left_operand * right_operand;




				}
				// требуется реализовать следующую логику:
				// 1. поиск элемента по ключу 
				// 2. если элемент не найден то добавление нового узла со значением ключа сохраняя балансировку.
				// Если элемент найден то нужно просто изменить foundnow на true. 
				// Т.е. достаточно просто поиска и вставки.
				// 3. В конце дерево необходимо ликвидировать.
				// Тип данных целочисленный ключ.

			}
		}

		if (my_amg_manager.bdiagonal_dominant > 0) {
			doublereal additional_to_diagonal = 0.0;
			integer jstr_diag = -1;
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];
				if ((row_ind != jstr) && (vector_sum[jstr] >= 0.0)) {
					additional_to_diagonal += vector_sum[jstr];
					vector_sum[jstr] = 0.0;
					hash_table[jstr] = false;
				}
				if (jstr == row_ind) jstr_diag = jstr;
			}
			vector_sum[jstr_diag] += additional_to_diagonal;
		}

		doublerealT maxth = -1.0;
		for (integer i_6 = 1; i_6 <= index_size; i_6++) {
			integer jstr = index_visit[i_6];
			hash_table[jstr] = false; // initialization hash.
			if (istr != jstr) {
				// 14 января 2016 года.
				// Правильно определить барьер только по внедиагональным элементам.
				if (fabs(vector_sum[jstr]) > maxth) maxth = fabs(vector_sum[jstr]);
			}
		}



		// huck: 16.04.2017

		for (integer i_61 = 1; i_61 <= index_size; i_61++) {

			integer jstr61 = index_visit[i_61];
			doublerealT vs161 = vector_sum[jstr61];



			if ((istr == jstr61) && (vs161 < 1.0e-20)) {
				// отрицательный элемент на диагонали.

				printf("Negative diagonal coefficient found. No panic. Upwind patching. 16.04.2017. \n");

				printf("bad string: \n");
				for (integer i_63 = 1; i_63 <= index_size; i_63++) {

					integer jstr63 = index_visit[i_63];
					doublerealT vs163 = vector_sum[jstr63];
#if doubleintprecision == 1
					printf("i=%lld j=%lld aij=%e\n", istr, jstr63, vs163);
#else
					printf("i=%d j=%d aij=%e\n", istr, jstr63, vs163);
#endif

				}


				// Адаптированные три правила бак-труба:
				// Amat. Диагонали присваиваем сумма модулей только отрицательных внедиагональных коэффициентов +
				// вычитаем из этого отрицательную диагональ. Потом умножаем на два.
				// B. Удвоение отрицательных внедиагональных коэффициентов.
				// C. Полное обнуление положительных внедиагональных коэффициентов (игнорирование).
				printf("patching string 16.04.2017: \n");
				/*
				for (integer i_62 = 1; i_62 <= index_size; i_62++) {
					integer jstr62 = index_visit[i_62];
					if (istr != jstr62) {
						if (vector_sum[jstr62] > 0.0) {
							index_visit[i_62] = -1; // не существует такого элемента (игнорирование).
							vector_sum[jstr61] += vector_sum[jstr62];
							vector_sum[jstr62] = 0.0;
						}
					}
				}

				if (vector_sum[jstr61] < 0.0) {
					vector_sum[jstr61] = 0.0;
					for (integer i_62 = 1; i_62 <= index_size; i_62++) {
						integer jstr62 = index_visit[i_62];
						if (jstr62 > -1) {
							if (istr != jstr62) {
								if (vector_sum[jstr62] < 0.0) {
									vector_sum[jstr61] += fabs(vector_sum[jstr62]);
								}
							}
						}
					}
				}
				*/

				// Denis Demidov recomendation
				// Всю строку домножаем на минус один, так чтобы диагональ стала положительна.
				for (integer i_62 = 1; i_62 <= index_size; i_62++) {
					integer jstr62 = index_visit[i_62];
					vector_sum[jstr62] *= -1.0;
				}
				if (ibsp_length < i_bsp_LIMIT) {
					// Запоминаем строку с отрицательной диагональю.
					bsp[ibsp_length].ilevel = ilevel;
					bsp[ibsp_length].istring_number = istr;
					ibsp_length++;
				}
				else {
					printf("ERROR!!! ibsp_length>=i_bsp_LIMIT\n");
					system("pause");
					exit(1);
				}
				// Выход из цикла for по переменной i_61.
				break;
			}


		}


		if (nsizeA > istartAnew2 + index_size) {
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {

				integer jstr = index_visit[i_6];
				doublerealT vs1 = vector_sum[jstr];

				// 7 ноября 2016 игнорируем чистые нули:
				if ((jstr > -1) && (fabs(vs1) > 1.0e-37)) {
					// Мы игнорируем чистые нули. 
					// Но вообще говоря непонятно почему они появляются.


						// алгебраический мультигрид Галёркина.
						// 22_10_2016.
					Ak1 Atemp;
					Atemp.aij = vs1;
					Atemp.i = istr;
					Atemp.j = jstr;

					Amat.aij[istartAnew2] = Atemp.aij;
					Amat.i[istartAnew2] = Atemp.i;
					Amat.j[istartAnew2] = Atemp.j;
					istartAnew2++;

				}
			}
		}
		else {

			my_realloc_memory1<real_mix_precision>(Amat.aij, istartAnew2 + 20 * index_size + 1);
			my_realloc_memory1<integer_mix_precision>(Amat.i, istartAnew2 + 20 * index_size + 1);
			my_realloc_memory1<integer_mix_precision>(Amat.j, istartAnew2 + 20 * index_size + 1);
			nsizeA = istartAnew2 + 20 * index_size;


			for (integer i_6 = 1; i_6 <= index_size; i_6++) {

				integer jstr = index_visit[i_6];
				doublerealT vs1 = vector_sum[jstr];

				// 7 ноября 2016 игнорируем чистые нули:
				if ((jstr > -1) && (fabs(vs1) > 1.0e-37)) {
					// Мы игнорируем чистые нули. 
					// Но вообще говоря непонятно почему они появляются.


					// алгебраический мультигрид Галёркина.
					// 22_10_2016.
					Ak1 Atemp;
					Atemp.aij = vs1;
					Atemp.i = istr;
					Atemp.j = jstr;

					Amat.aij[istartAnew2] = Atemp.aij;
					Amat.i[istartAnew2] = Atemp.i;
					Amat.j[istartAnew2] = Atemp.j;
					istartAnew2++;

				}
			}

			// Слишком мало памяти выделено под матрицу А и в неё не умещаются все данные.
			// Нужно увеличить объём выделяемой памяти для А и перезапустить приложение.
			//printf("Amat lack of memory\n");
			//printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			//printf("please, press any key to exit.\n");
			//system("pause");
			//exit(1);
		}

		index_size = 0;
		// Освобождение памяти из под АВЛ дерева.
		clear_AVL_Gus(root_Gus);
		root_Gus = 0;
	}

}

#endif

#pragma once
#ifndef BASIC_FUNCTIONS_MY_AGREGAT_AMG_RELAXATION_RUNGE_KUTT3or5_CPP
#define BASIC_FUNCTIONS_MY_AGREGAT_AMG_RELAXATION_RUNGE_KUTT3or5_CPP 1

// smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
// Достоинство в том что оптимальные параметры известны из литературы и не надо
// ничего вручную подбирать.
// Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
// что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer iorder)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		//doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy


		if (x_jacoby_buffer == nullptr) {
			printf("ERROR: x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}



		for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {

				integer istr = ii - iadd;
				doublerealT rold = x_jacoby_buffer[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// игнорирование positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

					//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
					//}
					//else {
					// не рабтает.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;


				x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
				//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

				//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
				//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
			}

		}

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5

  // smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
  // Достоинство в том что оптимальные параметры известны из литературы и не надо
  // ничего вручную подбирать.
  // Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
  // что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start,
	integer*& row_ptr_end, integer iadd, integer iorder, bool*& F_false_C_true, integer idirect)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		//doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy





		for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

			if (idirect == 1) {

				// Восходящая ветвь: сначала F потом C.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

					}
				}

			}
			else {
				// Нисходящая ветвь: сначала C потом F.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

					}
				}


			}

		}

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5


// smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
// Достоинство в том что оптимальные параметры известны из литературы и не надо
// ничего вручную подбирать.
// Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
// что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer iorder, doublerealT*& diag_minus_one, integer ibsp_length)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		//doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy


		if (x_jacoby_buffer == nullptr) {
			printf("ERROR: x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}

		if (ibsp_length == 0) {

			for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					//x[istr] = diag_minus_one[istr] * b[istr];
					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

						//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;


					x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
					//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

					//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
				}

			}
		}
		else {


			for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = diag_minus_one[istr] * b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

						//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;


					x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
					//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

					//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
				}

			}
		}
	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5

  // smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
  // Достоинство в том что оптимальные параметры известны из литературы и не надо
  // ничего вручную подбирать.
  // Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
  // что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5/*_oldV*/(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start,
	integer*& row_ptr_end, integer iadd, integer iorder, bool*& F_false_C_true, integer idirect, doublerealT*& diag_minus_one, integer ibsp_length)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		//doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy


		if (ibsp_length == 0) {
			for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

				if (idirect == 1) {

					// Восходящая ветвь: сначала F потом C.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					// marker 1
					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}


					// marker2
					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}

				}
				else {
					// Нисходящая ветвь: сначала C потом F.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}


				}

			}

		}
		else {


			for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

				if (idirect == 1) {

					// Восходящая ветвь: сначала F потом C.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}

				}
				else {
					// Нисходящая ветвь: сначала C потом F.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}


				}

			}
		}
	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5

// smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
  // Достоинство в том что оптимальные параметры известны из литературы и не надо
  // ничего вручную подбирать.
  // Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
  // что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
// Разворачивание цикла помогает использовать конвейер процессора более эффективно. 19.06.2021.
template <typename doublerealT>
void Runge_Kutt_3or5_optimize(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start,
	integer*& row_ptr_end, integer iadd, integer iorder, bool*& F_false_C_true, integer idirect, doublerealT*& diag_minus_one, integer ibsp_length)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy


		/*integer* id_in = new integer[endpos - startpos + 3];
		integer* id_out = new integer[endpos - startpos + 3];
		integer i_id_in = 0;
		integer i_id_out = 0;
		for (integer ii = startpos; ii <= endpos; ++ii) {
			if (F_false_C_true[ii] == false) {
				id_in[i_id_in] = ii;
				++i_id_in;
			}
			else {
				id_out[i_id_out] = ii;
				++i_id_out;
			}
		}*/

		if (ibsp_length == 0) {
			for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

				if (idirect == 1) {

					// Восходящая ветвь: сначала F потом C.


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}




					// marker 1
					//#pragma loop(hint_parallel(8))
					// //for (integer ii_left = 0; ii_left < i_id_in; ++ii_left)
					// 
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii)
					{
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false))
						{
							//integer ii = id_in[ii_left];							
							integer istr = ii - iadd;

							doublerealT rold = x_jacoby_buffer[istr];


							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];


							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];


							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)

							if (is2 - is1 + 1 >= 8) {
								integer ii1;
								doublerealT rsum1 = 0.0;
								doublerealT rsum2 = 0.0;
								doublerealT rsum3 = 0.0;
								for (ii1 = is1; ii1 + 3 <= is2; ii1 += 4)
								{
									//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									integer ipos1 = Amat.j[ii1 + 1];
									integer ipos2 = Amat.j[ii1 + 2];
									integer ipos3 = Amat.j[ii1 + 3];

									// 13.07.2016
									// игнорирование positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
									rsum1 += -Amat.aij[ii1 + 1] * x_jacoby_buffer[ipos1];
									rsum2 += -Amat.aij[ii1 + 2] * x_jacoby_buffer[ipos2];
									rsum3 += -Amat.aij[ii1 + 3] * x_jacoby_buffer[ipos3];

									//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
									//}
									//else {
									// не рабтает.
									//	ap_now += Amat.aij[ii1];
									//}
								}

								while (ii1 <= is2) {
									//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									// 13.07.2016
									// игнорирование positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

									//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
									//}
									//else {
									// не рабтает.
									//	ap_now += Amat.aij[ii1];
									//}
									++ii1;
								}
								rsum += rsum1 + rsum2 + rsum3;
							}
							else {
								for (integer ii1 = is1; ii1 <= is2; ++ii1)
								{
									//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									// 13.07.2016
									// игнорирование positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

									//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
									//}
									//else {
									// не рабтает.
									//	ap_now += Amat.aij[ii1];
									//}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}




					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}


					// marker2
					//#pragma loop(hint_parallel(8))
					// for (integer ii_left = 0; ii_left < i_id_out; ++ii_left)
					//   
#pragma omp parallel for					
					for (integer ii = startpos; ii <= endpos; ++ii)
					{
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii]))
						{

							// integer ii = id_out[ii_left];
							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)


							if (is2 - is1 + 1 >= 8) {
								integer ii1;
								doublerealT rsum1 = 0.0;
								doublerealT rsum2 = 0.0;
								doublerealT rsum3 = 0.0;
								for (ii1 = is1; ii1 + 3 <= is2; ii1 += 4)
								{
									//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									integer ipos1 = Amat.j[ii1 + 1];
									integer ipos2 = Amat.j[ii1 + 2];
									integer ipos3 = Amat.j[ii1 + 3];
									// 13.07.2016
									// игнорирование positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
									rsum1 += -Amat.aij[ii1 + 1] * x_jacoby_buffer[ipos1];
									rsum2 += -Amat.aij[ii1 + 2] * x_jacoby_buffer[ipos2];
									rsum3 += -Amat.aij[ii1 + 3] * x_jacoby_buffer[ipos3];

									//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
									//}
									//else {
									// не рабтает.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								while (ii1 <= is2) {
									//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									// 13.07.2016
									// игнорирование positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

									//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
									//}
									//else {
									// не рабтает.
									//	ap_now += Amat.aij[ii1];
									//}

									++ii1;
								}
								rsum += rsum1 + rsum2 + rsum3;
							}
							else {
								for (integer ii1 = is1; ii1 <= is2; ++ii1)
								{
									//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									// 13.07.2016
									// игнорирование positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

									//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
									//}
									//else {
									// не рабтает.
									//	ap_now += Amat.aij[ii1];
									//}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}

				}
				else {
					// Нисходящая ветвь: сначала C потом F.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							//x[istr] = diag_minus_one[istr] * b[istr];
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}


				}

			}

		}
		else {


			for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; ++inumber_step_Runge_Kutt) {

				if (idirect == 1) {

					// Восходящая ветвь: сначала F потом C.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}

				}
				else {
					// Нисходящая ветвь: сначала C потом F.

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						}
					}


					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						integer istr = ii - iadd;
						x_jacoby_buffer[istr] = x[istr];
					}

					//#pragma loop(hint_parallel(8))
#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) {

							integer istr = ii - iadd;
							doublerealT rold = x_jacoby_buffer[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// Распараллеливание почемуто тормозит очень сильно.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

								//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
								//}
								//else {
								// не рабтает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;


							x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
							//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

							//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
							//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

						}
					}


				}

			}
		}

		//delete[] id_in;
		//delete[] id_out;

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5



#endif 

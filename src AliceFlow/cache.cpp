// перенумерация неизвестных в соответствии с алгоритмом Катхила Маки 
// в надежде лучше попадать в кэш процессора.
// Очередь обход графа в ширину.
// 24.02.2022

#pragma once
#ifndef CACHE_CPP
#define CACHE_CPP 1

// Функционирует только для структурированной расчётной сетки.
// Гипотеза: данная нумерация улучшает локальность данных и способствует лучшему попаданию в кэш.
// Реализует алгоритм Катхила Маки.
// 26.02.2022
void renumerate_setup(equation3D*& sl, equation3D_bon*& slb, int maxelm, int maxbound,
	int*& new_number, int*& new_number_internal, int*& new_number_bound, // Прямая нумерация.
	int*& rev_number, bool bactive_ren) // Реверсированная нумерация.
{

	if (bactive_ren) {

		// Изменение порядка следования (нумерации) узлов для лучшего попадания в кэш устройства или процессора.

		// Внутри функции происходит выделение оперативной памяти под массивы:
		// new_number, new_number_internal, new_number_bound
		// rev_number, rev_number_internal, rev_number_bound

		const int size0 = maxelm + maxbound;

		//std::cout << "apriory size0=" << size0 << " maxelm=" << maxelm << " maxbound=" << maxbound << " " << std::endl;

		//for (int i = 0; i < maxelm; ++i) {
			//std::cout << sl[i].iP << " " << i << std::endl;
		//}
		//for (int i = 0; i < maxbound; ++i) {
			//std::cout << slb[i].iW << " " << i+maxelm << std::endl;
		//}
		//getchar();

		// Очередь номеров узлов.
		int* queue = new int[7 * static_cast<integer>(maxelm) + 2 * static_cast<integer>(maxbound) + 1];

		int pop = 0, push = 0;
		// индекс pop вытолкнуть из очереди.
		// индекс push вставить в очередь.



		new_number = new int[size0]; // Новая нумерация в x и векторе правой части.


		bool* in_queue = new bool[size0];
		bool* in_out_queue = new bool[size0]; // покинули очередь.
		int* in_queue_id = new int[size0]; // индекс в очережи.

#pragma omp parallel for
		for (int i = 0; i < size0; ++i)
		{
			new_number[i] = -1;
			in_queue[i] = false;
			in_out_queue[i] = false;
			in_queue_id[i] = -1; // не в очереди.
		}


		int iP = iP_perefirier_start;
		int ic = 0;

		int iPi = 0; // internal
		int iPb = 0; // boundary

		new_number_internal = new int[maxelm];
		new_number_bound = new int[maxbound];

		new_number[ic] = iP;
		++ic;
		new_number_internal[iPi] = iP;
		++iPi;
		in_queue[iP] = true;
		in_out_queue[iP] = true; // Покинул очередь и вошел очередным номером в new_number

		queue[push] = sl[iP].iW;
		in_queue_id[sl[iP].iW] = push;
		++push;
		in_queue[sl[iP].iW] = true;

		queue[push] = sl[iP].iE;
		in_queue_id[sl[iP].iE] = push;
		++push;
		in_queue[sl[iP].iE] = true;

		queue[push] = sl[iP].iS;
		in_queue_id[sl[iP].iS] = push;
		++push;
		in_queue[sl[iP].iS] = true;

		queue[push] = sl[iP].iN;
		in_queue_id[sl[iP].iN] = push;
		++push;
		in_queue[sl[iP].iN] = true;

		queue[push] = sl[iP].iB;
		in_queue_id[sl[iP].iB] = push;
		++push;
		in_queue[sl[iP].iB] = true;

		queue[push] = sl[iP].iT;
		in_queue_id[sl[iP].iT] = push;
		++push;
		in_queue[sl[iP].iT] = true;



		while (pop < push) {

			iP = queue[pop];
			++pop;

			if (ic < size0) {
				new_number[ic] = iP;
				++ic;
			}
			else {
				std::cout << "buffer new_number overflow.\n";
				system("PAUSE");
				exit(1);
			}

			if (iP < maxelm) {
				if (iPi < maxelm) {
					new_number_internal[iPi] = iP;
					++iPi;
				}
				else {
					std::cout << "buffer new_number_internal overflow.\n";
					system("PAUSE");
					exit(1);
				}
			}
			else {
				if (iPb < maxbound) {
					new_number_bound[iPb] = iP - maxelm;
					++iPb;
				}
				else {
					std::cout << "buffer new_number_bound overflow.\n";
					system("PAUSE");
					exit(1);
				}
			}

			if (iP < size0) {
			    in_out_queue[iP] = true; // Покинул очередь и вошел очередным номером в new_number
		    }
			else {
				std::cout << "buffer in_out_queue overflow.\n";
				system("PAUSE");
				exit(1);
			}

			// Добавляем в очередь соседей узла iP только в том случае если это 
			// внутренний узел.

			if (iP < maxelm) {

				// Сортировка сначала идут узлы с минимальной степенью, т.е. граничные.

				if (sl[iP].iW >= maxelm) {
					if (sl[iP].iW < size0) {
						if (in_out_queue[sl[iP].iW] == false)
						{// Еще не покидал очередь и не входил в список номеров new_number. 
							if (in_queue[sl[iP].iW] == false) {
								queue[push] = sl[iP].iW;
								in_queue_id[sl[iP].iW] = push;
								++push;
								in_queue[sl[iP].iW] = true;
							}

						}
					}
					else {
						std::cout << "buffer in_queue[iW] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iE >= maxelm) {
					if (sl[iP].iE < size0) {
						if (in_out_queue[sl[iP].iE] == false)
						{// Еще не покидал очередь и не входил в список номеров new_number. 
							if (in_queue[sl[iP].iE] == false) {
								queue[push] = sl[iP].iE;
								in_queue_id[sl[iP].iE] = push;
								++push;
								in_queue[sl[iP].iE] = true;
							}

						}
					}
					else {
						std::cout << "buffer in_queue[iE] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iS >= maxelm) {
					if (sl[iP].iS < size0) {
						if (in_out_queue[sl[iP].iS] == false)
						{// Еще не покидал очередь и не входил в список номеров new_number. 
							if (in_queue[sl[iP].iS] == false) {
								queue[push] = sl[iP].iS;
								in_queue_id[sl[iP].iS] = push;
								++push;
								in_queue[sl[iP].iS] = true;
							}

						}
					}
					else {
						std::cout << "buffer in_queue[iS] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iN >= maxelm) {
					if (sl[iP].iN < size0) {
						if (in_out_queue[sl[iP].iN] == false)
						{// Еще не покидал очередь и не входил в список номеров new_number. 
							if (in_queue[sl[iP].iN] == false) {
								queue[push] = sl[iP].iN;
								in_queue_id[sl[iP].iN] = push;
								++push;
								in_queue[sl[iP].iN] = true;
							}

						}
					}
					else {
						std::cout << "buffer in_queue[iN] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iB >= maxelm) {
					if (sl[iP].iB < size0) {
						if (in_out_queue[sl[iP].iB] == false)
						{// Еще не покидал очередь и не входил в список номеров new_number. 
							if (in_queue[sl[iP].iB] == false) {
								queue[push] = sl[iP].iB;
								in_queue_id[sl[iP].iB] = push;
								++push;
								in_queue[sl[iP].iB] = true;
							}

						}
					}
					else {
						std::cout << "buffer in_queue[iB] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iT >= maxelm) {
					if (sl[iP].iT < size0) {
					if (in_out_queue[sl[iP].iT] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iT] == false) {
							queue[push] = sl[iP].iT;
							in_queue_id[sl[iP].iT] = push;
							++push;
							in_queue[sl[iP].iT] = true;
						}

					}
					}
					else {
						std::cout << "buffer in_queue[iT] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}



				if (sl[iP].iW < maxelm) {
					if (sl[iP].iW < size0) {
					if (in_out_queue[sl[iP].iW] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iW] == false) {
							queue[push] = sl[iP].iW;
							in_queue_id[sl[iP].iW] = push;
							++push;
							in_queue[sl[iP].iW] = true;
						}

					}
					}
					else {
						std::cout << "tail buffer in_queue[iW] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iE < maxelm) {
					if (sl[iP].iE < size0) {
					if (in_out_queue[sl[iP].iE] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iE] == false) {
							queue[push] = sl[iP].iE;
							in_queue_id[sl[iP].iE] = push;
							++push;
							in_queue[sl[iP].iE] = true;
						}

					}
					}
					else {
						std::cout << "tail buffer in_queue[iE] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iS < maxelm) {
					if (sl[iP].iS < size0) {
					if (in_out_queue[sl[iP].iS] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iS] == false) {
							queue[push] = sl[iP].iS;
							in_queue_id[sl[iP].iS] = push;
							++push;
							in_queue[sl[iP].iS] = true;
						}

					}
					}
					else {
						std::cout << "tail buffer in_queue[iS] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iN < maxelm) {
					if (sl[iP].iN < size0) {
					if (in_out_queue[sl[iP].iN] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iN] == false) {
							queue[push] = sl[iP].iN;
							in_queue_id[sl[iP].iN] = push;
							++push;
							in_queue[sl[iP].iN] = true;
						}

					}
					}
					else {
						std::cout << "tail buffer in_queue[iN] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iB < maxelm) {
					if (sl[iP].iB < size0) {
					if (in_out_queue[sl[iP].iB] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iB] == false) {
							queue[push] = sl[iP].iB;
							in_queue_id[sl[iP].iB] = push;
							++push;
							in_queue[sl[iP].iB] = true;
						}

					}
					}
					else {
						std::cout << "tail buffer in_queue[iB] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}

				if (sl[iP].iT < maxelm) {
					if (sl[iP].iT < size0) {
					if (in_out_queue[sl[iP].iT] == false)
					{// Еще не покидал очередь и не входил в список номеров new_number. 
						if (in_queue[sl[iP].iT] == false) {
							queue[push] = sl[iP].iT;
							in_queue_id[sl[iP].iT] = push;
							++push;
							in_queue[sl[iP].iT] = true;
						}

					}
					}
					else {
						std::cout << "tail buffer in_queue[iT] overflow.\n";
						system("PAUSE");
						exit(1);
					}
				}
				


			}

		}



		//for (int i = 0; i < size0; ++i)
		//{
			//std::cout << new_number[i] << " \n";
			//if (i % 100 == 0) getchar();
		//}
		rev_number = new int[size0]; // Реверсированная (обратная нумерация).

		int* new_number_shadow = new int[size0];
		for (int i = 0; i < size0; ++i)
		{
			// 32 стал 1.
			new_number_shadow[new_number[i]] = i;
			//rev_number[i] = new_number[i];
		}

		for (int i = 0; i < size0; ++i)
		{
			new_number[i] = new_number_shadow[i];
		}

		delete[] new_number_shadow;

		delete[] queue;
		delete[] in_queue;
		delete[] in_out_queue;
		delete[] in_queue_id;



		//rev_number_internal = new int[size0]; // Реверсированная (обратная нумерация).
		//rev_number_bound = new int[size0]; // Реверсированная (обратная нумерация).

		// i_1 - счётчик внутренних узлов.
		// i_2 - счётчик граничных узлов.
		//int i_1 = 0, i_2 = 0;

		for (int i = 0; i < size0; ++i) {
			//new_number 32 ->1


			rev_number[new_number[i]] = i;
			/*if (i < maxelm) {

				// Внутренний
				rev_number_internal[new_number[i]] = i_1;
				++i_1;
			}
			else {

				// Граничный
				rev_number_bound[new_number[i]] = i_2;
				++i_2;
			}*/
		}

	}

} // renumerate

// Функционирует только для структурированной расчётной сетки.
// Гипотеза: данная нумерация улучшает локальность данных и способствует лучшему попаданию в кэш.
// Первая версия.
// 26.02.2022
void renumerate_setup1(equation3D* &sl, equation3D_bon* &slb, int maxelm, int maxbound,
	int* &new_number, int* &new_number_internal, int* &new_number_bound, // Прямая нумерация.
	int* &rev_number, bool bactive_ren) // Реверсированная нумерация.
{
	
	if (bactive_ren) {

		// Изменение порядка следования (нумерации) узлов для лучшего попадания в кэш устройства или процессора.

		// Внутри функции происходит выделение оперативной памяти под массивы:
		// new_number, new_number_internal, new_number_bound
		// rev_number, rev_number_internal, rev_number_bound

		const int size0 = maxelm + maxbound;

		//std::cout << "apriory size0=" << size0 << " maxelm=" << maxelm << " maxbound=" << maxbound << " " << std::endl;

		//for (int i = 0; i < maxelm; ++i) {
			//std::cout << sl[i].iP << " " << i << std::endl;
		//}
		//for (int i = 0; i < maxbound; ++i) {
			//std::cout << slb[i].iW << " " << i+maxelm << std::endl;
		//}
		//getchar();

		// Очередь номеров узлов.
		int* queue = new int[7 * static_cast<integer>(maxelm) + 2 * static_cast<integer>(maxbound) + 1];

		int pop = 0, push = 0;
		// индекс pop вытолкнуть из очереди.
		// индекс push вставить в очередь.

		

		new_number = new int[size0]; // Новая нумерация в x и векторе правой части.


		bool* in_queue = new bool[size0];
		bool* in_out_queue = new bool[size0]; // покинули очередь.
		int* in_queue_id = new int[size0]; // индекс в очережи.

#pragma omp parallel for
		for (int i = 0; i < size0; ++i)
		{
			new_number[i] = -1;
			in_queue[i] = false;
			in_out_queue[i] = false;
			in_queue_id[i] = -1; // не в очереди.
		}


		int iP = iP_perefirier_start;
		int ic = 0;

		int iPi = 0; // internal
		int iPb = 0; // boundary

		new_number_internal = new int[maxelm];
		new_number_bound = new int[maxbound];

		new_number[ic] = iP;
		++ic;
		new_number_internal[iPi] = iP;
		++iPi;
		in_queue[iP] = true;
		in_out_queue[iP] = true; // Покинул очередь и вошел очередным номером в new_number

		queue[push] = sl[iP].iW;
		in_queue_id[sl[iP].iW] = push;
		++push;
		in_queue[sl[iP].iW] = true;

		queue[push] = sl[iP].iE;
		in_queue_id[sl[iP].iE] = push;
		++push;
		in_queue[sl[iP].iE] = true;

		queue[push] = sl[iP].iS;
		in_queue_id[sl[iP].iS] = push;
		++push;
		in_queue[sl[iP].iS] = true;

		queue[push] = sl[iP].iN;
		in_queue_id[sl[iP].iN] = push;
		++push;
		in_queue[sl[iP].iN] = true;

		queue[push] = sl[iP].iB;
		in_queue_id[sl[iP].iB] = push;
		++push;
		in_queue[sl[iP].iB] = true;

		queue[push] = sl[iP].iT;
		in_queue_id[sl[iP].iT] = push;
		++push;
		in_queue[sl[iP].iT] = true;



		while (pop < push) {

			iP = queue[pop];
			++pop;

			new_number[ic] = iP;
			++ic;

			if (iP < maxelm) {
				new_number_internal[iPi] = iP;
				++iPi;
			}
			else {
				new_number_bound[iPb] = iP - maxelm;
				++iPb;
			}

			in_out_queue[iP] = true; // Покинул очередь и вошел очередным номером в new_number


			// Добавляем в очередь соседей узла iP только в том случае если это 
			// внутренний узел.

			if (iP < maxelm) {

				if (in_out_queue[sl[iP].iW] == false)
				{// Еще не покидал очередь и не входил в список номеров new_number. 
					if (in_queue[sl[iP].iW] == false) {
						queue[push] = sl[iP].iW;
						in_queue_id[sl[iP].iW] = push;
						++push;
						in_queue[sl[iP].iW] = true;
					}
					else {
						//swap pop <-> in_queue_id[sl[iP].iW]

						// Меняем повторно и даже троекратно встретившийся узел с первым в очереди на выход узел.
						// Повторно встретившийся при формировании очереди узел войдёт первым кандидатом в new_number.


						int tmp = queue[pop];
						int iP1 = queue[pop];
						int iP2 = queue[in_queue_id[sl[iP].iW]];
						queue[pop] = queue[in_queue_id[sl[iP].iW]];
						queue[in_queue_id[sl[iP].iW]] = tmp;

						in_queue_id[iP1] = in_queue_id[sl[iP].iW];
						in_queue_id[iP2] = pop;

					}
				}

				if (in_out_queue[sl[iP].iE] == false)
				{// Еще не покидал очередь и не входил в список номеров new_number. 
					if (in_queue[sl[iP].iE] == false) {
						queue[push] = sl[iP].iE;
						in_queue_id[sl[iP].iE] = push;
						++push;
						in_queue[sl[iP].iE] = true;
					}
					else {

						//swap pop <-> in_queue_id[sl[iP].iE]

						// Меняем повторно и даже троекратно встретившийся узел с первым в очереди на выход узел.
						// Повторно встретившийся при формировании очереди узел войдёт первым кандидатом в new_number.


						int tmp = queue[pop];
						int iP1 = queue[pop];
						int iP2 = queue[in_queue_id[sl[iP].iE]];
						queue[pop] = queue[in_queue_id[sl[iP].iE]];
						queue[in_queue_id[sl[iP].iE]] = tmp;

						in_queue_id[iP1] = in_queue_id[sl[iP].iE];
						in_queue_id[iP2] = pop;
					}
				}

				if (in_out_queue[sl[iP].iS] == false)
				{// Еще не покидал очередь и не входил в список номеров new_number. 
					if (in_queue[sl[iP].iS] == false) {
						queue[push] = sl[iP].iS;
						in_queue_id[sl[iP].iS] = push;
						++push;
						in_queue[sl[iP].iS] = true;
					}
					else {

						//swap pop <-> in_queue_id[sl[iP].iS]

						// Меняем повторно и даже троекратно встретившийся узел с первым в очереди на выход узел.
						// Повторно встретившийся при формировании очереди узел войдёт первым кандидатом в new_number.


						int tmp = queue[pop];
						int iP1 = queue[pop];
						int iP2 = queue[in_queue_id[sl[iP].iS]];
						queue[pop] = queue[in_queue_id[sl[iP].iS]];
						queue[in_queue_id[sl[iP].iS]] = tmp;

						in_queue_id[iP1] = in_queue_id[sl[iP].iS];
						in_queue_id[iP2] = pop;
					}
				}

				if (in_out_queue[sl[iP].iN] == false)
				{// Еще не покидал очередь и не входил в список номеров new_number. 
					if (in_queue[sl[iP].iN] == false) {
						queue[push] = sl[iP].iN;
						in_queue_id[sl[iP].iN] = push;
						++push;
						in_queue[sl[iP].iN] = true;
					}
					else {

						//swap pop <-> in_queue_id[sl[iP].iN]

						// Меняем повторно и даже троекратно встретившийся узел с первым в очереди на выход узел.
						// Повторно встретившийся при формировании очереди узел войдёт первым кандидатом в new_number.


						int tmp = queue[pop];
						int iP1 = queue[pop];
						int iP2 = queue[in_queue_id[sl[iP].iN]];
						queue[pop] = queue[in_queue_id[sl[iP].iN]];
						queue[in_queue_id[sl[iP].iN]] = tmp;

						in_queue_id[iP1] = in_queue_id[sl[iP].iN];
						in_queue_id[iP2] = pop;
					}
				}

				if (in_out_queue[sl[iP].iB] == false)
				{// Еще не покидал очередь и не входил в список номеров new_number. 
					if (in_queue[sl[iP].iB] == false) {
						queue[push] = sl[iP].iB;
						in_queue_id[sl[iP].iB] = push;
						++push;
						in_queue[sl[iP].iB] = true;
					}
					else {

						//swap pop <-> in_queue_id[sl[iP].iB]

						// Меняем повторно и даже троекратно встретившийся узел с первым в очереди на выход узел.
						// Повторно встретившийся при формировании очереди узел войдёт первым кандидатом в new_number.


						int tmp = queue[pop];
						int iP1 = queue[pop];
						int iP2 = queue[in_queue_id[sl[iP].iB]];
						queue[pop] = queue[in_queue_id[sl[iP].iB]];
						queue[in_queue_id[sl[iP].iB]] = tmp;

						in_queue_id[iP1] = in_queue_id[sl[iP].iB];
						in_queue_id[iP2] = pop;
					}
				}

				if (in_out_queue[sl[iP].iT] == false)
				{// Еще не покидал очередь и не входил в список номеров new_number. 
					if (in_queue[sl[iP].iT] == false) {
						queue[push] = sl[iP].iT;
						in_queue_id[sl[iP].iT] = push;
						++push;
						in_queue[sl[iP].iT] = true;
					}
					else {

						//swap pop <-> in_queue_id[sl[iP].iT]

						// Меняем повторно и даже троекратно встретившийся узел с первым в очереди на выход узел.
						// Повторно встретившийся при формировании очереди узел войдёт первым кандидатом в new_number.


						int tmp = queue[pop];
						int iP1 = queue[pop];
						int iP2 = queue[in_queue_id[sl[iP].iT]];
						queue[pop] = queue[in_queue_id[sl[iP].iT]];
						queue[in_queue_id[sl[iP].iT]] = tmp;

						in_queue_id[iP1] = in_queue_id[sl[iP].iT];
						in_queue_id[iP2] = pop;
					}
				}
			}

		}



		
		rev_number = new int[size0]; // Реверсированная (обратная нумерация).

		int* new_number_shadow = new int[size0];
		for (int i = 0; i < size0; ++i)
		{
			// 32 стал 1.
			new_number_shadow[new_number[i]] = i;
			//rev_number[i] = new_number[i];
		}

		for (int i = 0; i < size0; ++i)
		{
			new_number[i] = new_number_shadow[i];
		}

		delete[] new_number_shadow;

		

		delete[] queue;
		delete[] in_queue;
		delete[] in_out_queue;
		delete[] in_queue_id;

		
		
		//rev_number_internal = new int[size0]; // Реверсированная (обратная нумерация).
		//rev_number_bound = new int[size0]; // Реверсированная (обратная нумерация).

		// i_1 - счётчик внутренних узлов.
		// i_2 - счётчик граничных узлов.
		//int i_1 = 0, i_2 = 0;

		for (int i = 0; i < size0; ++i) {
			//new_number 32 ->1


			rev_number[new_number[i]] = i;
			/*if (i < maxelm) {

				// Внутренний
				rev_number_internal[new_number[i]] = i_1;
				++i_1;
			}
			else {

				// Граничный
				rev_number_bound[new_number[i]] = i_2;
				++i_2;
			}*/
		}

	}

} // renumerate_setup1

// Делает нумерацию в соответствии с лучшим попаданием в кеш.
  // Функционирует только для структурированной расчётной сетки.
void renumerate_direct(equation3D* &sl, equation3D_bon* &slb, int maxelm, int maxbound,
	doublereal* &x, doublereal* &b,
	int* &new_number, int* &new_number_internal, int* &new_number_bound, // Прямая нумерация.
	int* &rev_number,
	bool bactive_ren) // Обратная нумерация
{

	if (bactive_ren) {

		const int size0 = maxelm + maxbound;

		equation3D* sl_copy = new equation3D[maxelm];
		doublereal* x_copy = new doublereal[size0];
		doublereal* b_copy = new doublereal[size0];

#pragma omp parallel for
		for (int i = 0; i < maxelm; ++i) {

			x_copy[new_number[i]] = x[i];
			b_copy[new_number[i]] = b[i];

			sl_copy[i] = sl[new_number_internal[i]];
			

			sl_copy[i].iP = new_number[sl_copy[i].iP];
			sl_copy[i].iW = new_number[sl_copy[i].iW];
			sl_copy[i].iE = new_number[sl_copy[i].iE];
			sl_copy[i].iS = new_number[sl_copy[i].iS];
			sl_copy[i].iN = new_number[sl_copy[i].iN];
			sl_copy[i].iB = new_number[sl_copy[i].iB];
			sl_copy[i].iT = new_number[sl_copy[i].iT];

		}


		equation3D_bon* slb_copy = new equation3D_bon[maxbound];

#pragma omp parallel for
		for (int j = 0; j < maxbound; ++j) {

			int i = j + maxelm;

			x_copy[new_number[i]] = x[i];
			b_copy[new_number[i]] = b[i];

			slb_copy[j] = slb[new_number_bound[j]];
			

			if (slb_copy[j].iI > -1) slb_copy[j].iI = new_number[slb_copy[j].iI];
			slb_copy[j].iW = new_number[slb_copy[j].iW];

		}

#pragma omp parallel for
		for (int i = 0; i < maxelm; ++i) {
			x[i] = x_copy[i];
			b[i] = b_copy[i];

			sl[i] = sl_copy[i];
		}

#pragma omp parallel for
		for (int j = 0; j < maxbound; ++j) {

			int i = j + maxelm;

			x[i] = x_copy[i];
			b[i] = b_copy[i];

			

			slb[j] = slb_copy[j];

		}

		delete[] sl_copy;
		delete[] x_copy;
		delete[] b_copy;
		delete[] slb_copy;


		//std::cout << "apostoriory size0=" << size0 << " maxelm=" << maxelm << " maxbound=" << maxbound << " " << std::endl;

		//for (int i = 0; i < maxelm; ++i) {
			//std::cout << sl[i].iP << " " << i << std::endl;
		//}
		//for (int i = 0; i < maxbound; ++i) {
			//std::cout << slb[i].iW << " " << i + maxelm << std::endl;
		//}
		//getchar();

	}

}

// Восстанавливает прямую нумерацию в векторе x.
// 24.02.2022
void renumerate_reverse(doublereal* &x, int* &rev_number, int maxelm, int maxbound, bool bactive_ren) // Обратная нумерация
{
	if (bactive_ren) {

		const int size0 = maxelm + maxbound;

		// Создаём временный вектор y в котором будем хранить нужную нумерацию.
		doublereal* y = new doublereal[size0];

#pragma omp parallel for
		for (int i = 0; i < size0; ++i) 
		{			
			y[rev_number[i]] = x[i];
		}

		// copy
#pragma omp parallel for
		for (int i = 0; i < size0; ++i) 
		{
			x[i] = y[i];
		}

		delete[] y;

	}

}

#endif // CACHE_CPP
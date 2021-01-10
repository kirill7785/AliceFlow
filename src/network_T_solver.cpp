// Файл network_T_solver.cpp
// Графовый метод решения уравнений кондукционной теплопередачи.
// begin 24.06.2020 - начало разработки стационарной версии решателя.
// 27.06.2020 - стационарная версия решателя на АЛИС сетке (3262 строки).
// end **.**.**

#pragma once
#ifndef MY_NETWORK_T_SOLVER_CPP
#define MY_NETWORK_T_SOLVER_CPP 1

void Seidel_network(integer n, integer maxelm,
	doublereal*& rthdsd, doublereal*& potent,
	doublereal*& val, integer*& col_ind, integer*& row_ptr,
	bool b_nonlinear_network, integer*& id, WALL*& w)
{
	// Seidel start
	for (integer i = 0; i < n; i++) {
		if (i < maxelm) {
			doublereal sum = rthdsd[i];
			for (integer j = row_ptr[i] + 1; j < row_ptr[i + 1]; j++) {
				sum += -(val[j] * potent[col_ind[j]]);
			}
			// делим на диагональ.
			if (b_nonlinear_network) {
				// Верхняя релаксация.
				//printf("nonlinear");
				potent[col_ind[row_ptr[i]]] = potent[col_ind[row_ptr[i]]] + 1.85 * ((sum / val[row_ptr[i]]) - potent[col_ind[row_ptr[i]]]);
			}
			else {
				potent[col_ind[row_ptr[i]]] = sum / val[row_ptr[i]];
			}
		}
		else {
			if ((w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
				(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
				(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)) {
				doublereal sum = rthdsd[i];
				for (integer j = row_ptr[i] + 1; j < row_ptr[i + 1]; j++) {
					sum += -(val[j] * potent[col_ind[j]]);
				}
				potent[col_ind[row_ptr[i]]] = potent[col_ind[row_ptr[i]]] + 1.85 * ((sum / val[row_ptr[i]]) - potent[col_ind[row_ptr[i]]]);

			}
		}
	}
	//Seidel end
}


doublereal residual_network(integer n, integer maxelm,
	doublereal*& rthdsd, doublereal*& potent,
	doublereal*& val, integer*& col_ind, integer*& row_ptr,
	bool b_nonlinear_network, integer*& id, WALL*& w)
{
	// Residual start
	doublereal s = 0.0;
	for (integer i = 0; i < n; i++) {
		if (i < maxelm) {
			doublereal sum = rthdsd[i];
			for (integer j = row_ptr[i] + 1; j < row_ptr[i + 1]; j++) {
				sum += -(val[j] * potent[col_ind[j]]);
			}
			// делим на диагональ.
			s += (val[row_ptr[i]] * potent[col_ind[row_ptr[i]]] - sum)*(val[row_ptr[i]] * potent[col_ind[row_ptr[i]]] - sum);			
		}
		else {
			if ((w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
				(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY)||
				(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)) {
				doublereal sum = rthdsd[i];
				for (integer j = row_ptr[i] + 1; j < row_ptr[i + 1]; j++) {
					sum += -(val[j] * potent[col_ind[j]]);
				}
				s+=(val[row_ptr[i]]*potent[col_ind[row_ptr[i]]]  - sum)*(val[row_ptr[i]] * potent[col_ind[row_ptr[i]]] - sum);

			}
		}
	}
	return sqrt(s/n);
	//Residual end
}


void visible_CRS_Matrix(integer n, integer nnz, doublereal*& val,
	integer*& col_ind, integer*& row_ptr,
	doublereal*& rthdsd, doublereal*& potent,
	BLOCK*& b, integer& lb, integer*& id,
	WALL*& w, integer& lw, integer& maxelm,
	bool &b_first_start_matrix_print)
{
	if (b_first_start_matrix_print) {

		
		FILE* fp_inicialization_data = NULL;

#ifdef MINGW_COMPILLER
		int err_inicialization_data = 0;
		fp_inicialization_data = fopen64("graph.txt", "w");
		if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
		errno_t err_inicialization_data = 0;
		err_inicialization_data = fopen_s(&fp_inicialization_data, "graph.txt", "w");
#endif

		if (err_inicialization_data == 0) {

			for (integer i_1 = 0; i_1 < maxelm; i_1++) {
				fprintf(fp_inicialization_data, "%s ", b[id[i_1]].name);
			}
			for (integer i_1 = maxelm; i_1 < n; i_1++) {
				fprintf(fp_inicialization_data, "%s ", w[id[i_1]].name);
			}
			fprintf(fp_inicialization_data, "\n");


			integer* is = new integer[n];
			for (integer i_1 = 0; i_1 < n; i_1++) {
				if (i_1 < maxelm) {
					fprintf(fp_inicialization_data, "%s ", b[id[i_1]].name);
				}
				else {
					fprintf(fp_inicialization_data, "%s ", w[id[i_1]].name);
				}

				// Обработка одной строки.
				for (integer j_1 = 0; j_1 < n; j_1++) {
					is[j_1] = 0;
				}
				for (integer j = row_ptr[i_1] + 1; j < row_ptr[i_1 + 1]; j++) {
					is[col_ind[j]] = 1;
				}

				for (integer j_1 = 0; j_1 < n; j_1++) {

					if (j_1 < n - 1) {
						/*if (i_1 == j_1) {
							if (i_1 < maxelm) {
								std::cout << b[id[i_1]].name << ",";
							}
							else {
								std::cout << w[id[i_1]].name << ",";
							}
						}
						else*/ {
							fprintf(fp_inicialization_data, "%lld ", is[j_1]);
						}
					}
					else {
						/*if (i_1 == j_1) {
							if (i_1 < maxelm) {
								std::cout << b[id[i_1]].name;
							}
							else {
								std::cout << w[id[i_1]].name;
							}
						}
						else*/ {
							fprintf(fp_inicialization_data, "%lld ", is[j_1]);
						}
					}
				}
				//std::cout << std::endl;
				fprintf(fp_inicialization_data, "\n");
			}
			delete[] is;


			fclose(fp_inicialization_data);

		}
		std::cout << "Matrix is print..." << std::endl;
		system("PAUSE");
		b_first_start_matrix_print = false;

		
	}
}

void visible_CRS_Matrix2(integer n, integer nnz, doublereal*& val,
	integer*& col_ind, integer*& row_ptr,
	doublereal*& rthdsd, doublereal*& potent,
	BLOCK*& b, integer& lb, integer*& id,
	WALL*& w, integer& lw, integer& maxelm,
	bool &b_first_start_matrix_print)
{
	if (b_first_start_matrix_print) {
		integer* is = new integer[n];
		for (integer i_1 = 0; i_1 < n; i_1++) {
			// Обработка одной строки.
			for (integer j_1 = 0; j_1 < n; j_1++) {
				is[j_1] = 0;
			}
			for (integer j = row_ptr[i_1] + 1; j < row_ptr[i_1 + 1]; j++) {
				is[col_ind[j]] = 1;
			}
			std::cout << std::endl;
			for (integer j_1 = 0; j_1 < n; j_1++) {

				if (j_1 < n - 1) {
					/*if (i_1 == j_1) {
						if (i_1 < maxelm) {
							std::cout << b[id[i_1]].name << ",";
						}
						else {
							std::cout << w[id[i_1]].name << ",";
						}
					}
					else*/ {
						std::cout << is[j_1] << ",";
					}
				}
				else {
					/*if (i_1 == j_1) {
						if (i_1 < maxelm) {
							std::cout << b[id[i_1]].name;
						}
						else {
							std::cout << w[id[i_1]].name;
						}
					}
					else*/ {
						std::cout << is[j_1];
					}
				}
			}
			//std::cout << std::endl;
		}
		delete[] is;
		std::cout << std::endl;
		for (integer i_1 = 0; i_1 < maxelm; i_1++) {
			std::cout << b[id[i_1]].name << std::endl;
		}
		for (integer i_1 = maxelm; i_1 < n; i_1++) {
			std::cout << w[id[i_1]].name << std::endl;
		}
		std::cout << "Matrix is print..." << std::endl;
		system("PAUSE");
		b_first_start_matrix_print = false;
	}
}

// Вывод на консоль задачи линейной алгебры для проверки.  
void check_CRS_matrix(integer n, integer nnz, doublereal*& val,
	integer*& col_ind, integer*& row_ptr,
	doublereal*& rthdsd, doublereal*& potent,
	BLOCK*& b, integer& lb, integer* &id,
	WALL* &w, integer &lw, integer &maxelm)
{
	//std::cout << n << std::endl;
	for (integer i_1 = 0; i_1 < n; i_1++) {
		/*
		std::cout << i_1 << " " << val[i_1] << " " << col_ind[i_1] << " ";
		if (i_1 <= n) {
			std::cout << row_ptr[i_1] << " ";
		}
		if (i_1 < n) {
			std::cout << rthdsd[i_1] << " " << potent[i_1] << std::endl;
		}
		else {
			std::cout << std::endl;
		}
		*/
		if (val[row_ptr[i_1]] <= 0.0) {
			if (i_1 < maxelm) {
				if ((id[i_1] > 0) && (id[i_1] < lb)) {
					std::cout << b[id[i_1]].name << " ";
				}
				else {
					printf("error check_CRS_matrix id[i_1] is unknown:  ");
				}
			}
			else {
				std::cout << w[id[i_1]].name << " ";
			}
			printf("i_1==%lld %e\n", i_1, val[row_ptr[i_1]]);
			system("PAUSE");
		}

	}
}// check_CRS_matrix


// Вывод на консоль задачи линейной алгебры для проверки.  
void print_CRS_matrix(integer n, integer nnz, doublereal*& val,
	integer*& col_ind, integer*& row_ptr,
	doublereal*& rthdsd, doublereal*& potent)
{
	std::cout << n << std::endl;
	for (integer i_1 = 0; i_1 < nnz; i_1++) {
		std::cout << i_1 << " " << val[i_1] << " " << col_ind[i_1] << " ";
		if (i_1 <= n) {
			std::cout << row_ptr[i_1] << " ";
		}
		if (i_1 < n) {
			std::cout << rthdsd[i_1] << " " << potent[i_1] << std::endl;
		}
		else {
			std::cout << std::endl;
		}
		system("PAUSE");
	}
}// print_CRS_matrix

// 24.06.2020 - начало разработки.
// Графовый метод решения стационарного уравнения теплопередачи.
// 27.06.2020 Поддерживается АЛИС.
// 11.07.2020 начало тестирования.
void calculate_Network_T(TEMPER& t,
	BLOCK*& b, integer& lb,
	WALL*& w, integer& lw,
	integer& ls, TPROP* &matlist)
{
	bool b_first_start_matrix_print = true;

	bool* block_is_active = new bool[lb];// только активные блоки (true). Если один блок полностью перекрыт другим то он (false). cabinet всегда false.
	for (integer i_4 = 0; i_4 < lb; i_4++) {
		block_is_active[i_4] = false;
	}
	for (integer i_4 = 0; i_4 < t.maxelm; i_4++) {
		block_is_active[t.whot_is_block[i_4]] = true;
	}
	block_is_active[0] = false;

	integer maxelm = 0;
	for (integer i = 1; i < lb; i++) {
		if ((b[i].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[i])) maxelm++;
	}
	integer* id = new integer[maxelm + lw];
	integer* id_reverse = new integer[lb + lw];
	integer* wall2block_link = new integer[lw];
	integer ic = 0;
	for (integer i = 1; i < lb; i++) {
		if ((b[i].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[i])) {
			id[ic] = i;// идентификатор блока.
			id_reverse[i] = ic;
			ic++;
		}
	}
	for (integer i = 0; i < lw; i++) {
		id[ic] = i; // идентификатор стенки.
		id_reverse[lb + i] = ic;
		ic++;
		wall2block_link[i] = -1; // инициализация.
	}

	bool* hash = new bool[lb];
	bool* hash_wall = new bool[lw];
	integer** ilink = new integer * [maxelm + lw];
	integer** ilink_reverse = new integer * [maxelm + lw];
	doublereal** dS = new doublereal * [maxelm + lw];
	integer* inumber_neighbour = new integer[maxelm + lw];
	integer* inumber_neighbour_only_body = new integer[maxelm + lw];
	for (integer i = 0; i < maxelm + lw; i++) {
		inumber_neighbour[i] = 0;// нет соседей.
		inumber_neighbour_only_body[i] = 0;// нет соседей среди соседних блоков.
		ilink[i] = nullptr;
		ilink_reverse[i] = nullptr;
		dS[i] = nullptr;
	}

	{ // делаем у переменной ic локальную область видимости.
		ic = 0;
		for (integer i = 1; i < lb; i++) {
		if ((b[i].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[i])) {
			for (integer j_1 = 0; j_1 < lb; j_1++) {
				hash[j_1] = false;
			}
			for (integer j_1 = 0; j_1 < lw; j_1++) {
				hash_wall[j_1] = false;
			}
			hash[0] = true;
			hash[i] = true;
			integer ic1 = 0;
			// Находим количество соседей блока i, в массиве id он имеет номер ic.
			for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {
				if (t.whot_is_block[j_1] == i) {
					if ((t.neighbors_for_the_internal_node[E_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ic1++;
							}
						}
					}


				}
			}
			integer ic_block = ic1;
			//printf("Number of count block 2 block neighbours %lld.\n",ic1);
			//system("PAUSE");
			inumber_neighbour_only_body[ic] = ic_block;
			// Стенки с условием Дирихле, Ньютона-Рихмана или Стефана-Больцмана.
			for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {
				if (t.whot_is_block[j_1] == i) {
					//стенки 
					if (t.neighbors_for_the_internal_node[E_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][0][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][0][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][0][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if (t.neighbors_for_the_internal_node[T_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if (t.neighbors_for_the_internal_node[B_SIDE][0][j_1] >= t.maxelm) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ic1++;
							}
						}
					}

				}
			}

			//printf("Number of count block 2 block and wall neighbours %lld.\n", ic1);
			//system("PAUSE");
			ilink[ic] = new integer[ic1];
			for (integer j_1 = 0; j_1 < ic1; j_1++) {
				ilink[ic][j_1] = -1; // инициализация несуществующим индексом.
			}
			// обратное преобразование требует проиндексировать все блоки и все стенки.
			ilink_reverse[ic] = new integer[lb + lw + 1];
			for (integer j_1 = 0; j_1 < lb + lw + 1; j_1++) {
				ilink_reverse[ic][j_1] = -1;// инициализация несуществующим индексом.
			}

			inumber_neighbour[ic] = ic1; // количество блоков соседей.
			for (integer j_1 = 0; j_1 < lb; j_1++) {
				hash[j_1] = false;
			}

			hash[0] = true;
			hash[i] = true;
			ic1 = 0;
			// Записываем всех соседей блока i, в массиве id он имеет номер ic.
			for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {
				if (t.whot_is_block[j_1] == i) {

					if ((t.neighbors_for_the_internal_node[E_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][0][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][1][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][2][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][3][j_1];
						if (hash[t.whot_is_block[i_1]] == false) {
							if ((b[t.whot_is_block[i_1]].itype != PHYSICS_TYPE_IN_BODY::HOLLOW) && (block_is_active[t.whot_is_block[i_1]])) {
								hash[t.whot_is_block[i_1]] = true;
								ilink[ic][ic1] = t.whot_is_block[i_1];// Соседи блока id[ic].
								ilink_reverse[ic][t.whot_is_block[i_1]] = ic1;
								ic1++;
							}
						}
					}
				}
			}

			//printf("Additional internal neighbourhuuds on block %lld.\n", ic1);
			//system("PAUSE");

			for (integer j_1 = 0; j_1 < lw; j_1++) {
				hash_wall[j_1] = false;
			}

			ic_block = ic1;
			// Стенки с условием Дирихле, Ньютона-Рихмана или Стефана-Больцмана.
			for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {
				if (t.whot_is_block[j_1] == i) {
					//стенки 
					if (t.neighbors_for_the_internal_node[E_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if (t.neighbors_for_the_internal_node[W_SIDE][0][j_1] >= t.maxelm) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if (t.neighbors_for_the_internal_node[N_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if (t.neighbors_for_the_internal_node[S_SIDE][0][j_1] >= t.maxelm) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if (t.neighbors_for_the_internal_node[T_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							// стенка с условием Дирихле.
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле которая еще не встречалась.
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][0][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] >= t.maxelm)) {
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))) {
							if (hash_wall[t.border_neighbor[inumber].MCB - ls] == false) {
								// Это стенка с условием Дирихле
								hash_wall[t.border_neighbor[inumber].MCB - ls] = true;
								ilink[ic][ic1] = t.border_neighbor[inumber].MCB - ls;
								ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls] = ic1;
								wall2block_link[t.border_neighbor[inumber].MCB - ls] = i;// идентификатор блока с которым связана стенка.
								ic1++;
							}
						}
					}

				}
			}

			//printf("Additional internal and boundary neighbourhuuds on block %lld.\n", ic1);
			//system("PAUSE");

			dS[ic] = new doublereal[ic1];
			for (integer j_1 = 0; j_1 < ic1; j_1++) {
				dS[ic][j_1] = 0.0; // инициализация
			}

			// Вычисление общей площади контакта у блоков.
			for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {
				if (t.whot_is_block[j_1] == i) {
					if ((t.neighbors_for_the_internal_node[E_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][0][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "E_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "E_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][1][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "E_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "E_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][2][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "E_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "E_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][3][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "E_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "E_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][0][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "W_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "W_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr) && 
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][1][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "W_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "W_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][2][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "W_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "W_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][3][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "W_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "W_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hy * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][0][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "N_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "N_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][1][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "N_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "N_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][2][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "N_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "N_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][3][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "N_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "N_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][0][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "S_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "S_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][1][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "S_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "S_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][2][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "S_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "S_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][3][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "S_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "S_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hz;
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][0][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "T_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "T_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][1][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "T_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "T_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][2][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "T_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "T_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][3][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "T_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "T_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][0][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][0][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][0][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "B_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "B_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][1][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "B_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "B_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][2][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "B_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "B_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] > -1) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] < t.maxelm))
					{
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][3][j_1];
						if ((t.whot_is_block[i_1] != i) && (t.whot_is_block[i_1] != 0) && (block_is_active[t.whot_is_block[i_1]])) {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							if (ilink_reverse[ic][t.whot_is_block[i_1]] == -1)
							{
								std::cout << "B_SIDE problem ilink_reverse block2block" << std::endl;
								system("PAUSE");
							}
							if (ilink_reverse[ic][t.whot_is_block[i_1]] > ic1) {
								std::cout << "B_SIDE problem ilink_reverse block2block ic1 limit" << std::endl;
								system("PAUSE");
							}
							dS[ic][ilink_reverse[ic][t.whot_is_block[i_1]]] += hx * hy;
						}
					}

					// Стенки.

					if (t.neighbors_for_the_internal_node[E_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[E_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[E_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if (t.neighbors_for_the_internal_node[W_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[W_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[W_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hy * hz;
						}

					}
					if (t.neighbors_for_the_internal_node[N_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[N_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[N_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if (t.neighbors_for_the_internal_node[S_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							//std::cout << "ic_block = " << ic_block << "; t.border_neighbor[inumber].MCB=" << t.border_neighbor[inumber].MCB << "; ls=" << ls << std::endl;
							//std::cout << "inumber_neighbour[ic]=" << inumber_neighbour[ic]<<std::endl;
							//system("PAUSE");
							//dS[ic][ilink_reverse[ic][ic_block + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr) &&
						(t.neighbors_for_the_internal_node[S_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							//std::cout << "ic_block = " << ic_block << "; t.border_neighbor[inumber].MCB=" << t.border_neighbor[inumber].MCB << "; ls=" << ls << std::endl;
							//std::cout << "inumber_neighbour[ic]=" << inumber_neighbour[ic]<<std::endl;
							//system("PAUSE");
							//dS[ic][ilink_reverse[ic][ic_block + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							//std::cout << "ic_block = " << ic_block << "; t.border_neighbor[inumber].MCB=" << t.border_neighbor[inumber].MCB << "; ls=" << ls << std::endl;
							//std::cout << "inumber_neighbour[ic]=" << inumber_neighbour[ic]<<std::endl;
							//system("PAUSE");
							//dS[ic][ilink_reverse[ic][ic_block + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[S_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[S_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							//std::cout << "ic_block = " << ic_block << "; t.border_neighbor[inumber].MCB=" << t.border_neighbor[inumber].MCB << "; ls=" << ls << std::endl;
							//std::cout << "inumber_neighbour[ic]=" << inumber_neighbour[ic]<<std::endl;
							//system("PAUSE");
							//dS[ic][ilink_reverse[ic][ic_block + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hz;
						}

					}
					if (t.neighbors_for_the_internal_node[T_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr)&&
						(t.neighbors_for_the_internal_node[T_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[T_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if (t.neighbors_for_the_internal_node[B_SIDE][0][j_1] >= t.maxelm) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][0][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][1][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][1][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr)&&
						(t.neighbors_for_the_internal_node[B_SIDE][2][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][2][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}
					if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr) &&
						(t.neighbors_for_the_internal_node[B_SIDE][3][j_1] >= t.maxelm)) {
						// граничный узел
						integer i_1 = t.neighbors_for_the_internal_node[B_SIDE][3][j_1];
						integer inumber = i_1 - t.maxelm;
						// идентификатор граничного узла.
						if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
							(t.border_neighbor[inumber].MCB >= ls) &&
							((w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
								(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)))
						{
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);
							dS[ic][ilink_reverse[ic][lb + t.border_neighbor[inumber].MCB - ls]] += hx * hy;
						}

					}

				}
			}

			//printf("Square calculate Ok.");
			//system("PAUSE");

			ic++;
		}
	}
	}
	doublereal* potent = new doublereal[maxelm + lw]; // вектор решения.
	doublereal* rthdsd = new doublereal[maxelm + lw]; // правая часть.
	
	//unsteady TODO
	
	for (integer i = 0; i < maxelm + lw; i++) {
		potent[i] = t.operatingtemperature;
		if (i < maxelm) {
			// Внутренний блок.
			rthdsd[i] = (b[id[i]].g.xE - b[id[i]].g.xS) *
				(b[id[i]].g.yE - b[id[i]].g.yS) *
				(b[id[i]].g.zE - b[id[i]].g.zS) *
				get_power(b[id[i]].n_Sc, b[id[i]].temp_Sc, b[id[i]].arr_Sc, potent[i]);
		}
		else {
			// стенка.
			if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
				rthdsd[i] = w[id[i]].Tamb; // Только условия Дирихле.
			}
			else {
				// Нелинейное граничное условие.
				rthdsd[i] = 0.0;
			}
		}
	}
	// declarate matrix.

	// Расчёт числа ненулевых значений в матрице СЛАУ.
	integer nnz = maxelm + lw;
	for (integer i = 0; i < maxelm + lw; i++) {
		nnz += inumber_neighbour[i]; // количество связей с соседями.
		if (i >= maxelm) {
			if ((w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
				(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY)||
				(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY))
			{
				nnz++; // двухточечное нелинейное граничное условие.
			}
		}
	}
	doublereal* val = new doublereal[nnz];
	integer* col_ind = new integer[nnz];
	integer* row_ptr = new integer[maxelm + lw + 1];

	// Установка нелинейного флага.
	bool b_nonlinear_network = false;
	bool b_Newton_Richman = false;
	bool b_Stefan_Bolcman = false;
	for (integer i = maxelm; i < maxelm + lw; i++) {
		if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) {
			b_Newton_Richman = true; // Нелинейность Ньютона-Рихмана.
			b_nonlinear_network = true; // Задача нелинейна. Нужно применять нижнюю релаксацию.
		}
	}

	for (integer i = maxelm; i < maxelm + lw; i++) {
		if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
			b_Stefan_Bolcman = true; // Нелинейность Стефана - Больцмана, нужна более сильная релаксация.
			b_Newton_Richman = false;
			b_nonlinear_network = true; // Задача нелинейна. Нужно применять нижнюю релаксацию.
		}
	}


	printf("Assemble matrix Ok.");
	//system("PAUSE");

	// Хранит тепловые мощности в вт для каждой ячейки от радиационных потоков.
	//doublereal *rthdsd_radiation_global = new doublereal[t.maxelm];
	doublereal *rthdsd_radiation_loc = new doublereal[maxelm + lw];
	doublereal *rthdsd_radiation_loc_relax = new doublereal[maxelm + lw];

	for (integer i = 0; i < maxelm + lw; i++) {
		rthdsd_radiation_loc_relax[i] = 0.0;
	}



	doublereal tmax_old = -1.0e30;

	doublereal* potent_old = new doublereal[maxelm + lw];
	for (integer i25 = 0; i25 < maxelm + lw; i25++) potent_old[i25] = potent[i25];

	integer iter = -1;

	doublereal r1 = 1.0e-30;
	doublereal r2 = 1.0e30;

	// solve SLAU
	// Seidel method.
	for ( ;  ; ) {


		// Нулевая тепловая мощность от радиационных потоков.
		//#pragma omp parallel for
		//				for (integer i = 0; i < t.maxelm; i++) {
		//				rthdsd_radiation_global[i] = 0.0;
		//		}

		for (integer i = 0; i < maxelm + lw; i++) {
			rthdsd_radiation_loc[i] = 0.0;
		}

		iter++;
		if (iter > 300000) break;

		//r1 = 1.0e-30;

		doublereal tmax = -1.0e30;
		doublereal tmin = 1.0e30;
		for (integer i = 0; i < maxelm + lw; i++) {
			if (potent[i] > tmax) tmax = potent[i];
			if (potent[i] < tmin) tmin = potent[i];
		}

		if (fabs(tmax - tmin) < 1.0e-30) {
			printf("%4lld %e\n", iter, 10000.0);
		}
		else {
			//if (fabs(tmax_old - tmax) < 0.00002 * fabs(tmax - tmin)) {
				//break; // досрочный выход из цикла for.
			//}
			if (iter % 100 == 0) {
				printf("%4lld %e residual=%e min =%e max=%e\n", iter, fabs(tmax_old - tmax) / fabs(tmax - tmin), r2, tmin, tmax);
			}
		}
		if ((fabs(r2) < 1.0e-2) && (fabs(tmax_old - tmax) < 0.0005)) {
			std::cout << "break bPhysics_stop, dres<1e-2 && (fabs(maxnew - maxold) < 0.0005)" << std::endl;
			getchar();
			break;
		}
		if ((fabs(r2) < 1.4e-2) && (fabs(tmax_old - tmax) < 0.0005)) {
			std::cout << "break bPhysics_stop, dres<1.3e-2 && (fabs(maxnew - maxold) < 0.0005)" << std::endl;
			getchar();
			break;
		}


		bool bvacuum_Prism123 = false;
		for (integer i237 = 1; i237 < lb; i237++) {
			if (b[i237].radiation.binternalRadiation) {
				bvacuum_Prism123 = true;
			}
		}
		if (bvacuum_Prism123) {

			//if (iter > 249) break; // 02.11.2020

			//if ((fabs(r2) < 1.4e-2) && (fabs(tmax_old - tmax) < 0.01))
			//if ((fabs(tmax_old - tmax) < 0.0005))
			//{
			// Сразу досрочный выход из цикла решения но до истинного решения еще далеко.
			//std::cout << "break bPhysics_stop, NO dres<1e-2 && (fabs(maxnew - maxold) < 0.0005)" << std::endl;
			//getchar();
			//break;
			//}
		}



		r1 = r2;

		tmax_old = tmax;

		// Update matrix:
		// initializate matrix.
		for (integer i = 0; i < nnz; i++) {
			val[i] = 0.0;
			col_ind[i] = -1;
		}
		for (integer i = 0; i < maxelm + lw + 1; i++) {
			row_ptr[i] = nnz;
		}
		row_ptr[0] = 0;

		if (bvacuum_Prism123) {
			// Вычисление лучистых потоков на границе вакуумных промежутков.
			// Вычисление осреднённых температур в К на границах вакуумных промежутков:
			for (integer i = 0; i < lb; i++) {
				update_avg_temperatures(t.potent, b[i]);
			}

			// Вычисление плотностей радиационных тепловых потоков:
			for (integer i = 0; i < lb; i++) {
				calculation_density_radiation_heat_flux(b[i]);
			}
		}

		if (0) {
			/*
			radiosity_patch_for_vacuum_Prism_Object_(rthdsd_radiation_global, b, lb, t.maxelm, t.whot_is_block);

			for (integer i_4 = 0; i_4 < t.maxelm; i_4++) {
			rthdsd_radiation_loc[id_reverse[t.whot_is_block[i_4]]] += rthdsd_radiation_global[i_4];
			}

			for (integer i = 0; i < maxelm + lw; i++) {
			doublereal alphaR = 1.0;
			rthdsd_radiation_loc[i] = alphaR*rthdsd_radiation_loc[i] + (1.0 - alphaR)*rthdsd_radiation_loc_relax[i];
			}

			for (integer i = 0; i < maxelm + lw; i++) {
			rthdsd_radiation_loc_relax[i] = rthdsd_radiation_loc[i];
			}
			*/
		}
		else {
			for (integer i = 0; i < maxelm; i++) {
				if (b[id[i]].radiation.binternalRadiation) {
					if ((b[id[i]].radiation.nodelistW != nullptr) &&
						(b[id[i]].radiation.nodelistE != nullptr) &&
						(b[id[i]].radiation.nodelistS != nullptr) &&
						(b[id[i]].radiation.nodelistN != nullptr) &&
						(b[id[i]].radiation.nodelistB != nullptr) &&
						(b[id[i]].radiation.nodelistT != nullptr))
					{
						for (integer j_1 = 0; j_1 < inumber_neighbour[i]; j_1++) {
							if (j_1 < inumber_neighbour_only_body[i]) {
								// блок id[i] к блоку ilink[i][j_1].
								if ((b[ilink[i][j_1]].g.itypegeom == PRISM) && (!b[ilink[i][j_1]].radiation.binternalRadiation)) {
									if (fabs(b[id[i]].g.xS - b[ilink[i][j_1]].g.xE) < 1.0e-20) {
										// W ilink[i][j_1]] ---E bid[i] 
										rthdsd_radiation_loc[id_reverse[ilink[i][j_1]]] += -(b[id[i]].radiation.JW - b[id[i]].radiation.JE)*dS[i][j_1];
									}

									if (fabs(b[id[i]].g.xE - b[ilink[i][j_1]].g.xS) < 1.0e-20) {
										// E ilink[i][j_1]] ---W bid[i] 
										rthdsd_radiation_loc[id_reverse[ilink[i][j_1]]] += (b[id[i]].radiation.JW - b[id[i]].radiation.JE)*dS[i][j_1];
									}

									if (fabs(b[id[i]].g.yS - b[ilink[i][j_1]].g.yE) < 1.0e-20) {
										// S ilink[i][j_1]] ---N bid[i] 
										rthdsd_radiation_loc[id_reverse[ilink[i][j_1]]] += -(b[id[i]].radiation.JS - b[id[i]].radiation.JN)*dS[i][j_1];
									}

									if (fabs(b[id[i]].g.yE - b[ilink[i][j_1]].g.yS) < 1.0e-20) {
										// N ilink[i][j_1]] ---S bid[i] 
										rthdsd_radiation_loc[id_reverse[ilink[i][j_1]]] += (b[id[i]].radiation.JS - b[id[i]].radiation.JN)*dS[i][j_1];
									}

									if (fabs(b[id[i]].g.zS - b[ilink[i][j_1]].g.zE) < 1.0e-20) {
										// B ilink[i][j_1]] ---T bid[i] 
										rthdsd_radiation_loc[id_reverse[ilink[i][j_1]]] += -(b[id[i]].radiation.JB - b[id[i]].radiation.JT)*dS[i][j_1];
									}

									if (fabs(b[id[i]].g.zE - b[ilink[i][j_1]].g.zS) < 1.0e-20) {
										// T ilink[i][j_1]] ---B bid[i] 
										rthdsd_radiation_loc[id_reverse[ilink[i][j_1]]] += (b[id[i]].radiation.JB - b[id[i]].radiation.JT)*dS[i][j_1];
									}
								}

							}
						}
					}
				}
			}

			for (integer i = 0; i < maxelm + lw; i++) {
				doublereal alphaR = 1.0;
				rthdsd_radiation_loc[i] = alphaR*rthdsd_radiation_loc[i] + (1.0 - alphaR)*rthdsd_radiation_loc_relax[i];
			}

			for (integer i = 0; i < maxelm + lw; i++) {
				rthdsd_radiation_loc_relax[i] = rthdsd_radiation_loc[i];
			}
		}

		for (integer i = 0; i < maxelm + lw; i++) {
			if (i < maxelm) {
				// Внутренний блок.
				if (fabs(rthdsd_radiation_loc[i]) > 1.0e-300) {
					//std::cout << b[id[i]].name << " " << rthdsd_radiation_loc[i] << " W\n";

				}
			}
		}
		//system("PAUSE");


		for (integer i = 0; i < maxelm + lw; i++) {
			if (i < maxelm) {
				// Внутренний блок.
				rthdsd[i] = (b[id[i]].g.xE - b[id[i]].g.xS) *
					(b[id[i]].g.yE - b[id[i]].g.yS) *
					(b[id[i]].g.zE - b[id[i]].g.zS) *
					get_power(b[id[i]].n_Sc, b[id[i]].temp_Sc, b[id[i]].arr_Sc, potent[i])+
					rthdsd_radiation_loc[i]; // Тепловая мощность в Вт из за излучения, принадлежащая блоку id[i].
			}
			else {
				// стенка.
				if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
					rthdsd[i] = w[id[i]].Tamb; // Только условия Дирихле.
				}
				else {
					// Нелинейное граничное условие.
					rthdsd[i] = 0.0;
				}
			}
		}

		// assemble heat transfer martrix
		integer im = 1;
		integer idiag = 0;
		for (integer i = 0; i < maxelm + lw; i++) {
			if (i < maxelm) {
				// Внутренний блок.
				doublereal sum = 0.0;
				//printf("i==%lld inumber_neighbour[%lld]=%lld inumber_neighbour_only_body[%lld]=%lld\n", i, i, inumber_neighbour[i], i, inumber_neighbour_only_body[i]);
				//system("PAUSE");
				for (integer j = 0; j < inumber_neighbour[i]; j++) {

					bool ortho_k1 = false;
					bool ortho_k2 = false;
					doublereal ortho_m1 = 1.0;
					doublereal ortho_m2 = 1.0;

					// блок id[i] граничит с блоком ilink[i][j].
					doublereal distance = 0.0;
					if (j < inumber_neighbour_only_body[i]) {
						// блок id[i] к блоку ilink[i][j].
						TOCHKA bp0, bp1;
						if (b[id[i]].g.itypegeom == PRISM) {
							bp0.x = 0.5 * (b[id[i]].g.xS + b[id[i]].g.xE);
							bp0.y = 0.5 * (b[id[i]].g.yS + b[id[i]].g.yE);
							bp0.z = 0.5 * (b[id[i]].g.zS + b[id[i]].g.zE);
						}
						else if (b[id[i]].g.itypegeom == CYLINDER) {
							switch (b[id[i]].g.iPlane) {
							case XY_PLANE:
								bp0.x = b[id[i]].g.xC;
								bp0.y = b[id[i]].g.yC;
								bp0.z = b[id[i]].g.zC+0.5* b[id[i]].g.Hcyl;
								break;
							case XZ_PLANE:
								bp0.x = b[id[i]].g.xC;
								bp0.z = b[id[i]].g.zC;
								bp0.y = b[id[i]].g.yC + 0.5 * b[id[i]].g.Hcyl;
								break;
							case YZ_PLANE:
								bp0.y = b[id[i]].g.yC;
								bp0.z = b[id[i]].g.zC;
								bp0.x = b[id[i]].g.xC + 0.5 * b[id[i]].g.Hcyl;
								break;
							}
						}
						else if (b[id[i]].g.itypegeom == POLYGON) {
							// Вычисляем геометрический центр полигона.
							integer iscan = 0;
							switch (b[id[i]].g.iPlane_obj2) {
							case XY_PLANE:
								bp0.x = 0.0; bp0.y = 0.0;
								for (iscan = 0; iscan < b[id[i]].g.nsizei; iscan++) {
									bp0.x += b[id[i]].g.xi[iscan];
									bp0.y += b[id[i]].g.yi[iscan];
								}
								bp0.x /= 1.0 * b[id[i]].g.nsizei;
								bp0.y /= 1.0 * b[id[i]].g.nsizei;
								bp0.z = b[id[i]].g.zi[0] + 0.5 * b[id[i]].g.hi[0];
								break;
							case XZ_PLANE:
								bp0.x = 0.0; bp0.z = 0.0;
								for (iscan = 0; iscan < b[id[i]].g.nsizei; iscan++) {
									bp0.x += b[id[i]].g.xi[iscan];
									bp0.z += b[id[i]].g.zi[iscan];
								}
								bp0.x /= 1.0 * b[id[i]].g.nsizei;
								bp0.z /= 1.0 * b[id[i]].g.nsizei;
								bp0.y = b[id[i]].g.yi[0] + 0.5 * b[id[i]].g.hi[0];
								break;
							case YZ_PLANE:
								bp0.y = 0.0; bp0.z = 0.0;
								for (iscan = 0; iscan < b[id[i]].g.nsizei; iscan++) {
									bp0.y += b[id[i]].g.yi[iscan];
									bp0.z += b[id[i]].g.zi[iscan];
								}
								bp0.y /= 1.0 * b[id[i]].g.nsizei;
								bp0.z /= 1.0 * b[id[i]].g.nsizei;
								bp0.x = b[id[i]].g.xi[0] + 0.5 * b[id[i]].g.hi[0];
								break;
							}
						}
						else {
							bp0.x = 0.5 * (b[id[i]].g.xS + b[id[i]].g.xE);
							bp0.y = 0.5 * (b[id[i]].g.yS + b[id[i]].g.yE);
							bp0.z = 0.5 * (b[id[i]].g.zS + b[id[i]].g.zE);
						}
						if (b[ilink[i][j]].g.itypegeom == PRISM) {
							bp1.x = 0.5 * (b[ilink[i][j]].g.xS + b[ilink[i][j]].g.xE);
							bp1.y = 0.5 * (b[ilink[i][j]].g.yS + b[ilink[i][j]].g.yE);
							bp1.z = 0.5 * (b[ilink[i][j]].g.zS + b[ilink[i][j]].g.zE);
						}
						else if (b[ilink[i][j]].g.itypegeom == CYLINDER) {
							switch (b[ilink[i][j]].g.iPlane) {
							case XY_PLANE:
								bp1.x = b[ilink[i][j]].g.xC;
								bp1.y = b[ilink[i][j]].g.yC;
								bp1.z = b[ilink[i][j]].g.zC + 0.5 * b[ilink[i][j]].g.Hcyl;
								break;
							case XZ_PLANE:
								bp1.x = b[ilink[i][j]].g.xC;
								bp1.z = b[ilink[i][j]].g.zC;
								bp1.y = b[ilink[i][j]].g.yC + 0.5 * b[ilink[i][j]].g.Hcyl;
								break;
							case YZ_PLANE:
								bp1.y = b[ilink[i][j]].g.yC;
								bp1.z = b[ilink[i][j]].g.zC;
								bp1.x = b[ilink[i][j]].g.xC + 0.5 * b[ilink[i][j]].g.Hcyl;
								break;
							}
						}
						else if (b[ilink[i][j]].g.itypegeom == POLYGON) {
							// Вычисляем геометрический центр полигона.
							integer iscan = 0;
							switch (b[ilink[i][j]].g.iPlane_obj2) {
							case XY_PLANE:
								bp1.x = 0.0; bp1.y = 0.0;
								for (iscan = 0; iscan < b[ilink[i][j]].g.nsizei; iscan++) {
									bp1.x += b[ilink[i][j]].g.xi[iscan];
									bp1.y += b[ilink[i][j]].g.yi[iscan];
								}
								bp1.x /= 1.0 * b[ilink[i][j]].g.nsizei;
								bp1.y /= 1.0 * b[ilink[i][j]].g.nsizei;
								bp1.z = b[ilink[i][j]].g.zi[0] + 0.5 * b[ilink[i][j]].g.hi[0];
								break;
							case XZ_PLANE:
								bp1.x = 0.0; bp1.z = 0.0;
								for (iscan = 0; iscan < b[ilink[i][j]].g.nsizei; iscan++) {
									bp1.x += b[ilink[i][j]].g.xi[iscan];
									bp1.z += b[ilink[i][j]].g.zi[iscan];
								}
								bp1.x /= 1.0 * b[ilink[i][j]].g.nsizei;
								bp1.z /= 1.0 * b[ilink[i][j]].g.nsizei;
								bp1.y = b[ilink[i][j]].g.yi[0] + 0.5 * b[ilink[i][j]].g.hi[0];
								break;
							case YZ_PLANE:
								bp1.y = 0.0; bp1.z = 0.0;
								for (iscan = 0; iscan < b[ilink[i][j]].g.nsizei; iscan++) {
									bp1.y += b[ilink[i][j]].g.yi[iscan];
									bp1.z += b[ilink[i][j]].g.zi[iscan];
								}
								bp1.y /= 1.0 * b[ilink[i][j]].g.nsizei;
								bp1.z /= 1.0 * b[ilink[i][j]].g.nsizei;
								bp1.x = b[ilink[i][j]].g.xi[0] + 0.5 * b[ilink[i][j]].g.hi[0];
								break;
							}
						}
						else {
							bp1.x = 0.5 * (b[ilink[i][j]].g.xS + b[ilink[i][j]].g.xE);
							bp1.y = 0.5 * (b[ilink[i][j]].g.yS + b[ilink[i][j]].g.yE);
							bp1.z = 0.5 * (b[ilink[i][j]].g.zS + b[ilink[i][j]].g.zE);
						}

						distance = sqrt((bp0.x - bp1.x)*(bp0.x - bp1.x) +
							(bp0.y - bp1.y)*(bp0.y - bp1.y) +
							(bp0.z - bp1.z)*(bp0.z - bp1.z));
						if (distance < 1.0e-12) {
							// Защита от деления на ноль.
							// Берется длина минимального ребра из двух блоков.
							// Центры блоков случайным образом совпали
							distance = fmin(fmin(fabs((b[id[i]].g.xS - b[id[i]].g.xE)),
								fmin(fabs((b[id[i]].g.yS - b[id[i]].g.yE)),
									fabs((b[id[i]].g.zS - b[id[i]].g.zE)))),
								fmin(fabs(b[ilink[i][j]].g.xS - b[ilink[i][j]].g.xE),
									fmin((fabs(b[ilink[i][j]].g.yS - b[ilink[i][j]].g.yE)),
										(fabs(b[ilink[i][j]].g.zS - b[ilink[i][j]].g.zE)))));
						}
					}
					else {
						// блок id[i] к стенке MCB-ls==ilink[i][j].
						//printf("ilink[i][j] = %lld, i==%lld j==%lld\n", ilink[i][j],i,j);
						//system("PAUSE");

						TOCHKA bp0;
						if (b[id[i]].g.itypegeom == PRISM) {
							bp0.x = 0.5 * (b[id[i]].g.xS + b[id[i]].g.xE);
							bp0.y = 0.5 * (b[id[i]].g.yS + b[id[i]].g.yE);
							bp0.z = 0.5 * (b[id[i]].g.zS + b[id[i]].g.zE);
						}
						else if (b[id[i]].g.itypegeom == CYLINDER) {
							switch (b[id[i]].g.iPlane) {
							case XY_PLANE:
								bp0.x = b[id[i]].g.xC;
								bp0.y = b[id[i]].g.yC;
								bp0.z = b[id[i]].g.zC + 0.5 * b[id[i]].g.Hcyl;
								break;
							case XZ_PLANE:
								bp0.x = b[id[i]].g.xC;
								bp0.z = b[id[i]].g.zC;
								bp0.y = b[id[i]].g.yC + 0.5 * b[id[i]].g.Hcyl;
								break;
							case YZ_PLANE:
								bp0.y = b[id[i]].g.yC;
								bp0.z = b[id[i]].g.zC;
								bp0.x = b[id[i]].g.xC + 0.5 * b[id[i]].g.Hcyl;
								break;
							}
						}
						else if (b[id[i]].g.itypegeom == POLYGON) {
							// Вычисляем геометрический центр полигона.
							integer iscan = 0;
							switch (b[id[i]].g.iPlane_obj2) {
							case XY_PLANE:
								bp0.x = 0.0; bp0.y = 0.0;
								for (iscan = 0; iscan < b[id[i]].g.nsizei; iscan++) {
									bp0.x += b[id[i]].g.xi[iscan];
									bp0.y += b[id[i]].g.yi[iscan];
								}
								bp0.x /= 1.0 * b[id[i]].g.nsizei;
								bp0.y /= 1.0 * b[id[i]].g.nsizei;
								bp0.z = b[id[i]].g.zi[0] + 0.5 * b[id[i]].g.hi[0];
								break;
							case XZ_PLANE:
								bp0.x = 0.0; bp0.z = 0.0;
								for (iscan = 0; iscan < b[id[i]].g.nsizei; iscan++) {
									bp0.x += b[id[i]].g.xi[iscan];
									bp0.z += b[id[i]].g.zi[iscan];
								}
								bp0.x /= 1.0 * b[id[i]].g.nsizei;
								bp0.z /= 1.0 * b[id[i]].g.nsizei;
								bp0.y = b[id[i]].g.yi[0] + 0.5 * b[id[i]].g.hi[0];
								break;
							case YZ_PLANE:
								bp0.y = 0.0; bp0.z = 0.0;
								for (iscan = 0; iscan < b[id[i]].g.nsizei; iscan++) {
									bp0.y += b[id[i]].g.yi[iscan];
									bp0.z += b[id[i]].g.zi[iscan];
								}
								bp0.y /= 1.0 * b[id[i]].g.nsizei;
								bp0.z /= 1.0 * b[id[i]].g.nsizei;
								bp0.x = b[id[i]].g.xi[0] + 0.5 * b[id[i]].g.hi[0];
								break;
							}
						}
						else {
							bp0.x = 0.5 * (b[id[i]].g.xS + b[id[i]].g.xE);
							bp0.y = 0.5 * (b[id[i]].g.yS + b[id[i]].g.yE);
							bp0.z = 0.5 * (b[id[i]].g.zS + b[id[i]].g.zE);
						}

						distance = sqrt((bp0.x - 0.5 * (w[ilink[i][j]].g.xS + w[ilink[i][j]].g.xE)) *
							(bp0.x - 0.5 * (w[ilink[i][j]].g.xS + w[ilink[i][j]].g.xE)) +
							(bp0.y - 0.5 * (w[ilink[i][j]].g.yS + w[ilink[i][j]].g.yE)) *
							(bp0.y - 0.5 * (w[ilink[i][j]].g.yS + w[ilink[i][j]].g.yE)) +
							(bp0.z - 0.5 * (w[ilink[i][j]].g.zS + w[ilink[i][j]].g.zE)) *
							(bp0.z - 0.5 * (w[ilink[i][j]].g.zS + w[ilink[i][j]].g.zE)));
					}
					doublereal rho, cp, lam;
					rho = 1.1614; cp = 1005; lam = 0.025; // инициализация default  dry air 300K 1atm properties
					if (matlist[b[id[i]].imatid].blibmat == 1) {
						// библиотечный, находящийся внутри программы AliceFlow материал.
						if (b[id[i]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
							my_solid_properties(potent[i], rho, cp, lam, matlist[b[id[i]].imatid].ilibident);
							// проверка на допустимость температур.
							diagnostic_critical_temperature(potent[i], f, t, b, lb);
						} // SOLID
						if (b[id[i]].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							doublereal mu, beta_t; // значения не используются но требуются.
							doublereal pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).

							my_fluid_properties(potent[i], pressure, rho, cp, lam, mu, beta_t, matlist[b[id[i]].imatid].ilibident);
						} // FLUID
					}
					else if (matlist[b[id[i]].imatid].blibmat == 0) {
						// материал определённый пользователем:
						// постоянные свойства.
						if (!((fabs(matlist[b[id[i]].imatid].orthotropy_multiplyer_x - 1.0) < 0.0001) &&
							(fabs(matlist[b[id[i]].imatid].orthotropy_multiplyer_y - 1.0) < 0.0001) &&
							(fabs(matlist[b[id[i]].imatid].orthotropy_multiplyer_z - 1.0) < 0.0001))) {
							ortho_k1 = true;
							ortho_m1 = fmax(matlist[b[id[i]].imatid].orthotropy_multiplyer_x, fmax(
								matlist[b[id[i]].imatid].orthotropy_multiplyer_y,
								matlist[b[id[i]].imatid].orthotropy_multiplyer_z));
						}
						rho = matlist[b[id[i]].imatid].rho;
						//cp=matlist[b[ib].imatid].cp;
						//lam=matlist[b[ib].imatid].lam;
						cp = get_cp(matlist[b[id[i]].imatid].n_cp, matlist[b[id[i]].imatid].temp_cp, matlist[b[id[i]].imatid].arr_cp, potent[i]);
						lam = get_lam(matlist[b[id[i]].imatid].n_lam, matlist[b[id[i]].imatid].temp_lam, matlist[b[id[i]].imatid].arr_lam, potent[i]);

					}
					doublereal lambda_G = lam;
					if (j < inumber_neighbour_only_body[i]) {
						// Блок граничит с блоком.
						if (matlist[b[ilink[i][j]].imatid].blibmat == 1) {
							// библиотечный, находящийся внутри программы AliceFlow материал.
							if (b[ilink[i][j]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								my_solid_properties(potent[i], rho, cp, lam, matlist[b[ilink[i][j]].imatid].ilibident);
								// проверка на допустимость температур.
								diagnostic_critical_temperature(potent[i], f, t, b, lb);
							} // SOLID
							if (b[ilink[i][j]].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
								doublereal mu, beta_t; // значения не используются но требуются.
								doublereal pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).

								my_fluid_properties(potent[i], pressure, rho, cp, lam, mu, beta_t, matlist[b[ilink[i][j]].imatid].ilibident);
							} // FLUID
						}
						else if (matlist[b[ilink[i][j]].imatid].blibmat == 0) {
							// материал определённый пользователем:
							// постоянные свойства.
							if (!((fabs(matlist[b[ilink[i][j]].imatid].orthotropy_multiplyer_x - 1.0) < 0.0001) &&
								(fabs(matlist[b[ilink[i][j]].imatid].orthotropy_multiplyer_y - 1.0) < 0.0001) &&
								(fabs(matlist[b[ilink[i][j]].imatid].orthotropy_multiplyer_z - 1.0) < 0.0001))) {
								ortho_k2 = true;
								ortho_m2 = fmax(matlist[b[ilink[i][j]].imatid].orthotropy_multiplyer_x, fmax(
									matlist[b[ilink[i][j]].imatid].orthotropy_multiplyer_y,
									matlist[b[ilink[i][j]].imatid].orthotropy_multiplyer_z));
							}
							rho = matlist[b[ilink[i][j]].imatid].rho;
							//cp=matlist[b[ilink[i][j]].imatid].cp;
							//lam=matlist[b[ilink[i][j]].imatid].lam;
							cp = get_cp(matlist[b[ilink[i][j]].imatid].n_cp, matlist[b[ilink[i][j]].imatid].temp_cp, matlist[b[ilink[i][j]].imatid].arr_cp, potent[i]);
							lam = get_lam(matlist[b[ilink[i][j]].imatid].n_lam, matlist[b[ilink[i][j]].imatid].temp_lam, matlist[b[ilink[i][j]].imatid].arr_lam, potent[i]);

						}
						lambda_G = 2.0 * lambda_G * lam / (lambda_G + lam);
						if (ortho_k1 && ortho_k2) {// плата к плате.
							// учет ортотропности коэффициента теплопроводности.
							lambda_G *= 2.0 * ortho_m1 * ortho_m2 / (ortho_m1 + ortho_m2);
						}
					}
					val[im] = -((lambda_G * dS[i][j]) / distance);
					if (j < inumber_neighbour_only_body[i]) {
						// блок id[i] к блоку ilink[i][j].
						col_ind[im] = id_reverse[ilink[i][j]]; // номер столбца.
					}
					else {
						col_ind[im] = id_reverse[lb + ilink[i][j]];
					}
					sum += ((lambda_G * dS[i][j]) / distance);
					im++;
				}
				val[idiag] = sum; // диагональное преобладание.
				col_ind[idiag] = i;
				idiag = im;
				row_ptr[i + 1] = im;
				im++;
			}
			else {
				// стенка
				//((w[id[i]].ifamily == DIRICHLET_FAMILY) ||
					//(w[id[i]].ifamily == NEWTON_RICHMAN_FAMILY) ||
					//(w[id[i]].ifamily == STEFAN_BOLCMAN_FAMILY)))

				if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
					val[idiag] = 1.0; // Только условия Дирихле.
					col_ind[idiag] = i;
					row_ptr[i + 1] = idiag + 1;
					idiag++;
				}
				else if ((w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
					(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY)||
					(w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)) {
					// Находим номер блока ib с которым контактирует стенка.
					integer ib = wall2block_link[id[i]];
					integer ic = id_reverse[ib];

					TOCHKA bp0;
					if (b[ib].g.itypegeom == PRISM) {
						bp0.x = 0.5 * (b[ib].g.xS + b[ib].g.xE);
						bp0.y = 0.5 * (b[ib].g.yS + b[ib].g.yE);
						bp0.z = 0.5 * (b[ib].g.zS + b[ib].g.zE);
					}
					else if (b[ib].g.itypegeom == CYLINDER) {
						switch (b[ib].g.iPlane) {
						case XY_PLANE:
							bp0.x = b[ib].g.xC;
							bp0.y = b[ib].g.yC;
							bp0.z = b[ib].g.zC + 0.5 * b[ib].g.Hcyl;
							break;
						case XZ_PLANE:
							bp0.x = b[ib].g.xC;
							bp0.z = b[ib].g.zC;
							bp0.y = b[ib].g.yC + 0.5 * b[ib].g.Hcyl;
							break;
						case YZ_PLANE:
							bp0.y = b[ib].g.yC;
							bp0.z = b[ib].g.zC;
							bp0.x = b[ib].g.xC + 0.5 * b[ib].g.Hcyl;
							break;
						}
					}
					else if (b[ib].g.itypegeom == POLYGON) {
						// Вычисляем геометрический центр полигона.
						integer iscan = 0;
						switch (b[ib].g.iPlane_obj2) {
						case XY_PLANE:
							bp0.x = 0.0; bp0.y = 0.0;
							for (iscan = 0; iscan < b[ib].g.nsizei; iscan++) {
								bp0.x += b[ib].g.xi[iscan];
								bp0.y += b[ib].g.yi[iscan];
							}
							bp0.x /= 1.0 * b[ib].g.nsizei;
							bp0.y /= 1.0 * b[ib].g.nsizei;
							bp0.z = b[ib].g.zi[0] + 0.5 * b[ib].g.hi[0];
							break;
						case XZ_PLANE:
							bp0.x = 0.0; bp0.z = 0.0;
							for (iscan = 0; iscan < b[ib].g.nsizei; iscan++) {
								bp0.x += b[ib].g.xi[iscan];
								bp0.z += b[ib].g.zi[iscan];
							}
							bp0.x /= 1.0 * b[ib].g.nsizei;
							bp0.z /= 1.0 * b[ib].g.nsizei;
							bp0.y = b[ib].g.yi[0] + 0.5 * b[ib].g.hi[0];
							break;
						case YZ_PLANE:
							bp0.y = 0.0; bp0.z = 0.0;
							for (iscan = 0; iscan < b[ib].g.nsizei; iscan++) {
								bp0.y += b[ib].g.yi[iscan];
								bp0.z += b[ib].g.zi[iscan];
							}
							bp0.y /= 1.0 * b[ib].g.nsizei;
							bp0.z /= 1.0 * b[ib].g.nsizei;
							bp0.x = b[ib].g.xi[0] + 0.5 * b[ib].g.hi[0];
							break;
						}
					}
					else {
						bp0.x = 0.5 * (b[ib].g.xS + b[ib].g.xE);
						bp0.y = 0.5 * (b[ib].g.yS + b[ib].g.yE);
						bp0.z = 0.5 * (b[ib].g.zS + b[ib].g.zE);
					}

					doublereal distance = sqrt((bp0.x - 0.5 * (w[id[i]].g.xS + w[id[i]].g.xE)) *
						(bp0.x - 0.5 * (w[id[i]].g.xS + w[id[i]].g.xE)) +
						(bp0.y - 0.5 * (w[id[i]].g.yS + w[id[i]].g.yE)) *
						(bp0.y - 0.5 * (w[id[i]].g.yS + w[id[i]].g.yE)) +
						(bp0.z - 0.5 * (w[id[i]].g.zS + w[id[i]].g.zE)) *
						(bp0.z - 0.5 * (w[id[i]].g.zS + w[id[i]].g.zE)));
					doublereal rho, cp, lam;
					rho = 1.1614; cp = 1005; lam = 0.025; // инициализация default  dry air 300K 1atm properties
					if (matlist[b[ib].imatid].blibmat == 1) {
						// библиотечный, находящийся внутри программы AliceFlow материал.
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
							my_solid_properties(potent[ic], rho, cp, lam, matlist[b[ib].imatid].ilibident);
							// проверка на допустимость температур.
							diagnostic_critical_temperature(potent[ic], f, t, b, lb);
						} // SOLID
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							doublereal mu, beta_t; // значения не используются но требуются.
							doublereal pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).

							my_fluid_properties(potent[ic], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
						} // FLUID
					}
					else if (matlist[b[ib].imatid].blibmat == 0) {
						// материал определённый пользователем:
						// постоянные свойства.
						rho = matlist[b[ib].imatid].rho;
						//cp=matlist[b[ib].imatid].cp;
						//lam=matlist[b[ib].imatid].lam;
						cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, potent[ic]);
						lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, potent[ic]);

					}
					doublereal lambda_G = lam;
					// Гипотеза!!! Только одна стенка и она записана в конце списка.
					//((lambda_G*dS[ic][inumber_neighbour[ic]-1]) / distance);
					val[idiag] = ((lambda_G) / distance);
					col_ind[idiag] = i;//стенка
					idiag++;
					val[idiag] = -((lambda_G) / distance);
					col_ind[idiag] = ic;//блок
					idiag++;
					row_ptr[i + 1] = idiag;
					im = idiag + 1;
					if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) {

						//printf("nonlinear:  potent[ic]=%e  potent[i]=%e Tamb=%e\n", potent[ic], potent[i], w[id[i]].Tamb);
						//printf("lambda_G=%e distance=%e w[id[i]].film_coefficient=%e\n", lambda_G, distance, w[id[i]].film_coefficient);
						//system("PAUSE");
						if (potent[i] > w[id[i]].Tamb) {
							rthdsd[i] = -w[id[i]].film_coefficient * (potent[i] - w[id[i]].Tamb); // Условия Ньютона-Рихмана.
						}
						else {
							rthdsd[i] = 0.0;
							//rthdsd[i] = -w[id[i]].film_coefficient * (2.0); // Условия Ньютона-Рихмана.
						}
					}
					if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
						// Условие Стефана - Больцмана.
						rthdsd[i] = -w[id[i]].emissivity *w[id[i]].ViewFactor * STEFAN_BOLCMAN_CONST * 
							((273.15 + potent[i]) * (273.15 + potent[i]) * (273.15 + potent[i]) * (273.15 + potent[i]) -
							(273.15 + w[id[i]].Tamb) * (273.15 + w[id[i]].Tamb) * (273.15 + w[id[i]].Tamb) * (273.15 + w[id[i]].Tamb));
					}
					if (w[id[i]].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY) {
						// Однородное условие Неймана
						rthdsd[i] = 0.0;
					}
				}
			}
		}

		integer n = maxelm + lw;

		// One iteration
		r2 = residual_network(n, maxelm, rthdsd, potent, val, col_ind, row_ptr, b_nonlinear_network, id, w);
		//Seidel_network(n, maxelm, rthdsd, potent, val, col_ind, row_ptr, b_nonlinear_network, id, w);

		//printf("nnz=%lld row_ptr[n]=%lld\n",nnz, row_ptr[n]);

		/*
		// Нормировка
		for (integer i_1 = 0; i_1 < n; i_1++) {
			for (integer j_1 = row_ptr[i_1] + 1; j_1 < row_ptr[i_1 + 1]; j_1++) {
				val[j_1] /= val[row_ptr[i_1]];
			}
			rthdsd[i_1] /= val[row_ptr[i_1]];
			val[row_ptr[i_1]] = 1.0;
		}
		*/

		//  b_first_start_matrix_print - защита от повторного срабатывания.
		visible_CRS_Matrix(n, nnz, val, col_ind, row_ptr, rthdsd, potent, b, lb, id, w, lw, maxelm, b_first_start_matrix_print);
		check_CRS_matrix(n, nnz, val, col_ind, row_ptr, rthdsd, potent,b,lb,id,w,lw,maxelm);
		//print_CRS_matrix(n, nnz, val, col_ind, row_ptr, rthdsd, potent);
		if (b_nonlinear_network) {
			// граничное условие Ньютона - Рихмана 


			doublereal alpha = 0.98; // нижняя релаксация.

			if (b_Newton_Richman) {
				// граничное условие Ньютона - Рихмана
				// Ok релаксация достаточна.
			}

			if (b_Stefan_Bolcman) {
				// Нелинейность Стефана - Больцмана, нужна более сильная релаксация.
				alpha = 0.2; // было 0.1
			}

			// Если задача нелинейна.
			if ((b_Newton_Richman)||(b_Stefan_Bolcman)) {
				// нижняя релаксация введённая в матрицу СЛАУ.
				doublereal alphaA = 0.9;// было 0.4
				for (integer i_1 = 0; i_1 < maxelm; i_1++) {					
					// Это нужно чтобы сошелся солвер решения СЛАУ.
					rthdsd[i_1] += (1 - alphaA) * (val[row_ptr[i_1]] / alphaA) * potent_old[i_1];
					val[row_ptr[i_1]] = val[row_ptr[i_1]] / alphaA;
				}
				for (integer i_1 = maxelm; i_1 < maxelm + lw; i_1++) {
					if ((w[id[i_1]].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) ||
						(w[id[i_1]].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
						(w[id[i_1]].ifamily == WALL_BOUNDARY_CONDITION::NEIMAN_FAMILY)) {
						// К нелинейному граничному условию применим также нижнюю релаксацию.
						// Это нужно чтобы сошелся солвер решения СЛАУ.
						rthdsd[i_1] += (1 - alphaA) * (val[row_ptr[i_1]] / alphaA) * potent_old[i_1];
						val[row_ptr[i_1]] = val[row_ptr[i_1]] / alphaA;
					}
				}
			}		


			//Bi_CGStabCRS(n, val, col_ind, row_ptr, rthdsd, potent, 2000);
			bool worked_successfully;
			amg_loc_memory_networkT(val, col_ind, row_ptr, n, nnz,
				rthdsd, potent, 1.0, true, 0, worked_successfully,
				b, lb, ls, maxelm, false);

			// Решатель Дениса Демидова.
			//amgcl_networkT_solver(val, col_ind, row_ptr, n, rthdsd, potent, false);
				
			for (integer i25 = 0; i25 < n; i25++) {
				potent[i25] = potent_old[i25] + alpha * (potent[i25] - potent_old[i25]);
				//if (potent[i25] < t.operatingtemperature) potent[i25] = t.operatingtemperature;
			}

			for (integer i25 = 0; i25 < n; i25++) potent_old[i25] = potent[i25];
			//system("PAUSE");
		}
		else {
			//Bi_CGStabCRS(n, val, col_ind, row_ptr, rthdsd, potent, 2000);
			bool worked_successfully;
			amg_loc_memory_networkT(val, col_ind, row_ptr, n, nnz,
				rthdsd, potent, 1.0, true, 0, worked_successfully,
				b, lb, ls, maxelm,false);
		}
		// system("PAUSE");

	}

	// print report
	doublereal tmax = -1.0e30;
	doublereal tmin = 1.0e30;
	for (integer i = 0; i < maxelm + lw; i++) {
		if (potent[i] > tmax) tmax = potent[i];
		if (potent[i] < tmin) tmin = potent[i];
		if (i < maxelm) {
			std::cout << "Temperature block " << b[id[i]].name << "is equal =" << potent[i] << " gradus Celsius." << std::endl;
		}
		else {
			std::cout << "Temperature wall " << w[id[i]].name << "is equal =" << potent[i] << " gradus Celsius." << std::endl;
		}
	}
	std::cout << "maximum temperature is equal =" << tmax << " gradus Celsius." << std::endl;
	std::cout << "minimum temperature is equal =" << tmin << " gradus Celsius." << std::endl;

	// Сохранение температуры на сетке.
	for (integer i = 0; i < t.maxelm; i++) {
		if (block_is_active[t.whot_is_block[i]]) {
			t.potent[i] = potent[id_reverse[t.whot_is_block[i]]];
		}
	}
	for (integer i = 0; i < t.maxbound; i++) {
		// Копируем температуру из ближайшего внутреннего узла.
		if ((t.border_neighbor[i].MCB >= ls) && (t.border_neighbor[i].MCB < ls + lw)) {
			// Твердая стенка. 
			t.potent[t.maxelm + i] = potent[id_reverse[lb+ t.border_neighbor[i].MCB-ls]];
		}
		else {
			// Температуры во внутреннем и граничных узлах равны.
			t.potent[t.maxelm + i] = t.potent[t.border_neighbor[i].iI];
		}
	}

	bool* bw = new bool[lw];
	for (integer i = 0; i < lw; i++) {
		bw[i] = false;
	}


	for (integer j52=0; j52<lw; j52++) {
		doublereal pdiss = 0.0;
		for (integer i = 0; i < t.maxbound; i++) {
			// Копируем температуру из ближайшего внутреннего узла.
			if ((t.border_neighbor[i].MCB >= ls) && (t.border_neighbor[i].MCB < ls + lw)) {
				if (t.border_neighbor[i].MCB - ls == j52) {
					// Твердая стенка. 
					bw[t.border_neighbor[i].MCB - ls] = true;
					doublereal dS52 = 0.0;
					
					integer ib = t.whot_is_block[t.border_neighbor[i].iI];

					doublereal dx = 0.0, dy = 0.0, dz = 0.0;
					volume3D(t.border_neighbor[i].iI, t.nvtx, t.pa, dx, dy, dz);

					switch (w[t.border_neighbor[i].MCB - ls].iPlane) {
					case XY_PLANE:
						dS52 = dx*dy;
						break;
					case XZ_PLANE:
						dS52 = dx*dz;
						break;
					case YZ_PLANE:
						dS52 = dy*dz;
						break;
					}

					TOCHKA bp0;
					if (b[ib].g.itypegeom == PRISM) {
						bp0.x = 0.5 * (b[ib].g.xS + b[ib].g.xE);
						bp0.y = 0.5 * (b[ib].g.yS + b[ib].g.yE);
						bp0.z = 0.5 * (b[ib].g.zS + b[ib].g.zE);
					}
					else if (b[ib].g.itypegeom == CYLINDER) {
						switch (b[ib].g.iPlane) {
						case XY_PLANE:
							bp0.x = b[ib].g.xC;
							bp0.y = b[ib].g.yC;
							bp0.z = b[ib].g.zC + 0.5 * b[ib].g.Hcyl;
							break;
						case XZ_PLANE:
							bp0.x = b[ib].g.xC;
							bp0.z = b[ib].g.zC;
							bp0.y = b[ib].g.yC + 0.5 * b[ib].g.Hcyl;
							break;
						case YZ_PLANE:
							bp0.y = b[ib].g.yC;
							bp0.z = b[ib].g.zC;
							bp0.x = b[ib].g.xC + 0.5 * b[ib].g.Hcyl;
							break;
						}
					}
					else if (b[ib].g.itypegeom == POLYGON) {
						// Вычисляем геометрический центр полигона.
						integer iscan = 0;
						switch (b[ib].g.iPlane_obj2) {
						case XY_PLANE:
							bp0.x = 0.0; bp0.y = 0.0;
							for (iscan = 0; iscan < b[ib].g.nsizei; iscan++) {
								bp0.x += b[ib].g.xi[iscan];
								bp0.y += b[ib].g.yi[iscan];
							}
							bp0.x /= 1.0 * b[ib].g.nsizei;
							bp0.y /= 1.0 * b[ib].g.nsizei;
							bp0.z = b[ib].g.zi[0] + 0.5 * b[ib].g.hi[0];
							break;
						case XZ_PLANE:
							bp0.x = 0.0; bp0.z = 0.0;
							for (iscan = 0; iscan < b[ib].g.nsizei; iscan++) {
								bp0.x += b[ib].g.xi[iscan];
								bp0.z += b[ib].g.zi[iscan];
							}
							bp0.x /= 1.0 * b[ib].g.nsizei;
							bp0.z /= 1.0 * b[ib].g.nsizei;
							bp0.y = b[ib].g.yi[0] + 0.5 * b[ib].g.hi[0];
							break;
						case YZ_PLANE:
							bp0.y = 0.0; bp0.z = 0.0;
							for (iscan = 0; iscan < b[ib].g.nsizei; iscan++) {
								bp0.y += b[ib].g.yi[iscan];
								bp0.z += b[ib].g.zi[iscan];
							}
							bp0.y /= 1.0 * b[ib].g.nsizei;
							bp0.z /= 1.0 * b[ib].g.nsizei;
							bp0.x = b[ib].g.xi[0] + 0.5 * b[ib].g.hi[0];
							break;
						}
					}
					else {
						bp0.x = 0.5 * (b[ib].g.xS + b[ib].g.xE);
						bp0.y = 0.5 * (b[ib].g.yS + b[ib].g.yE);
						bp0.z = 0.5 * (b[ib].g.zS + b[ib].g.zE);
					}

					doublereal distance = sqrt((bp0.x - 0.5 * (w[t.border_neighbor[i].MCB - ls].g.xS + w[t.border_neighbor[i].MCB - ls].g.xE)) *
						(bp0.x - 0.5 * (w[t.border_neighbor[i].MCB - ls].g.xS + w[t.border_neighbor[i].MCB - ls].g.xE)) +
						(bp0.y - 0.5 * (w[t.border_neighbor[i].MCB - ls].g.yS + w[t.border_neighbor[i].MCB - ls].g.yE)) *
						(bp0.y - 0.5 * (w[t.border_neighbor[i].MCB - ls].g.yS + w[t.border_neighbor[i].MCB - ls].g.yE)) +
						(bp0.z - 0.5 * (w[t.border_neighbor[i].MCB - ls].g.zS + w[t.border_neighbor[i].MCB - ls].g.zE)) *
						(bp0.z - 0.5 * (w[t.border_neighbor[i].MCB - ls].g.zS + w[t.border_neighbor[i].MCB - ls].g.zE)));
					doublereal rho, cp, lam;
					rho = 1.1614; cp = 1005; lam = 0.025; // инициализация default  dry air 300K 1atm properties
					if (matlist[b[ib].imatid].blibmat == 1) {
						// библиотечный, находящийся внутри программы AliceFlow материал.
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
							my_solid_properties(potent[ic], rho, cp, lam, matlist[b[ib].imatid].ilibident);
							// проверка на допустимость температур.
							diagnostic_critical_temperature(potent[ic], f, t, b, lb);
						} // SOLID
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
							doublereal mu, beta_t; // значения не используются но требуются.
							doublereal pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).

							my_fluid_properties(potent[ic], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
						} // FLUID
					}
					else if (matlist[b[ib].imatid].blibmat == 0) {
						// материал определённый пользователем:
						// постоянные свойства.
						rho = matlist[b[ib].imatid].rho;
						//cp=matlist[b[ib].imatid].cp;
						//lam=matlist[b[ib].imatid].lam;
						cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, potent[ic]);
						lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, potent[ic]);

					}
					doublereal lambda_G = lam;

					pdiss += dS52*lambda_G*fabs(potent[id_reverse[ib]] - potent[id_reverse[lb + t.border_neighbor[i].MCB - ls]]) / distance;

					
				}
			}


		}

		printf("%s power is = %e\n", w[j52].name, pdiss);
	}

	delete[] bw;

	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;

	delete[] potent; // вектор решения.
	delete[] potent_old;
	delete[] rthdsd; // правая часть.

	delete[] id;
	delete[] id_reverse;
	delete[] hash;
	delete[] hash_wall;
	delete[] wall2block_link;
	delete[] block_is_active;
	delete[] inumber_neighbour;
	delete[] inumber_neighbour_only_body;
	for (integer i = 0; i < maxelm + lw; i++) {
		delete[] dS[i];
		delete[] ilink[i];
		delete[] ilink_reverse[i];
	}
	delete[] dS;
	delete[] ilink;
	delete[] ilink_reverse;

	//system("PAUSE");

} // calculate_Network_T

#endif
// ��������� ��������� ����������� ������� �� ����������� �������.
// ������������ � �������������� ������������� ������.
#pragma once
#ifndef _SPGEMM_MATRIX_BY_MATRIX_SPARSE_MULTIPLICATION_CPP_
#define _SPGEMM_MATRIX_BY_MATRIX_SPARSE_MULTIPLICATION_CPP_ 1


// ������ ��� ������� realloc().
// � ������ ������� �������������� �������� �� nullptr.
// ����������� ���������� ��������� ���� ����������� ������ � �������� ������� � 
// ����� ������� ������������ ��������� ���������� ��� ������� �.
// � ������ ����������� ������������ �������� ������� ��������� ������ ���������
// ���������� ��� �������� ������� �.
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

// ����������� ��������� ��������� ���� ����������� ������.
// �������� �. ���������� IBM 1978.
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

	// ������ ������������ ��� ���������� ������������� �����.

	for (int i_9 = 0; i_9 < iKnumber_thread; i_9++) {

#pragma omp parallel for schedule (static)
		for (integer i_91 = 0; i_91 < iSIZE; i_91++) {
			hash_table_m[i_9][i_91] = false;// inicialization
		}
		index_size_m[i_9] = 0;
		istartAnew_m[i_9] = 0;
	}

	const integer LIMIT_ANEW_m = (integer)(0.125 * AccumulqtorA_m_SIZE8 * nnz + 1);

	omp_set_num_threads(8); // ������ 8 �������.

	// ��������� ������ ������� ���������.
	// ���������� ���������� �� �������������.
#pragma omp parallel for schedule (static) 
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

#ifdef _OPENMP 
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif

		// �� ������ ���-�������. 

		// ��������� ������� i-�� ������ �����������
		for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
			integer col_ind = P[ii].j;
			// ��������� col_ind ������ ������� ��������

			// ����� ���������� ������� �� ������� ����.
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
					// ������ ����������.

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

					hash_table_m[tid][jstr] = false; // �������������� ���-������� ��� ��������� ��������.

					doublerealT vs1 = vector_sum_m[tid][jstr];

					// 7 ������ 2016 ���������� ������ ����:
					if (fabs(vs1) > 1.0e-37) {
						// �� �� ���������� � ������� ������ ����.
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
			// ������� ���� ������ �������� ��� ������� � � � �� �� ��������� ��� ������.
			// ����� ��������� ����� ���������� ������ ��� � � ������������� ����������.
			printf("Amat lack of memory\n");
			printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			printf("please, press any key to exit.\n");
			system("pause");
			exit(1);
		}



		index_size_m[tid] = 0; // ����� ����������, ������ ����������.			

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
			   // �������� ������ ��� �������� ������� �.
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

	// ���������� �� ����� �.�. ������������ schedule (static)
	//Counting_Sort(Amat, istartAnew_mem, istartAnew - 1, false, n_a[ilevel - 1]);
	//if (bprint_mesage_diagnostic) {
		//printf("Counting Sort End. \n");
	//}

	//getchar();

	omp_set_num_threads(i_my_num_core_parallelesation);

#endif

}

// ����������� ��������� ��������� ���� ����������� ������.
// �������� �. ���������� IBM 1978.
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
	
	// ������ ������������ ��� ���������� ������������� �����.

	for (int i_9 = 0; i_9 < iKnumber_thread; i_9++) {

		integer ISIZE =  n + 2;

#pragma omp parallel for
		for (integer i_91 = 0; i_91 < ISIZE; i_91++) {
			hash_table_m[i_9][i_91] = false;// inicialization
		}
		index_size_m[i_9] = 0;
		istartAnew_m[i_9] = 0;
	}


	// �� ����� ����������� ����� ������� ���������, �
	// ����� ��������� ��������� ����� ������ ������ ��������
	// �������� ������� ������ ����������.

	const integer LIMIT_ANEW_m = (integer)(0.125 * AccumulqtorA_m_SIZE8 * nnz + 1);


	
	omp_set_num_threads(8); // ������ 8 �������.

	// ��������� ������ ������� ���������.
	// ���������� ���������� �� �������������.
#pragma omp parallel for schedule (static) 
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

#ifdef _OPENMP 
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif

		// �� ������ ���-�������.
		// ��������� ��� �������� ������ ������ ��������.
		for (integer ii1 = row_ind_AS[istr]; ii1 <= row_ind_AE[istr]; ii1++) {

			integer col_ind = Amat.j[ii1];
			doublerealT left_operand = Amat.aij[ii1];

			// ��������� col_ind ������ ������� �������� ���������� �����.
			integer i_1start = row_ind_PS[col_ind];
			integer i_1end = row_ind_PE[col_ind];

			
			
			for (integer ii2 = i_1start; ii2 <= i_1end; ii2++) {


				doublerealT right_operand = P[ii2].aij;

				integer iaddind = P[ii2].i;

				if (hash_table_m[tid][iaddind]) {

					vector_sum_m[tid][iaddind] += left_operand * right_operand;

				}
				else {

					// ������ ����������.
					index_size_m[tid]++;
					index_visit_m[tid][index_size_m[tid]] = iaddind;

					// ���������� ������� � hash table �� O(1).
					hash_table_m[tid][iaddind] = true;

					//ifoundind = index_size;
					//vector_sum[index_visit[ifoundind]] = left_operand*right_operand;
					vector_sum_m[tid][iaddind] = left_operand * right_operand;


				}
				// ��������� ����������� ��������� ������:
				// 1. ����� �������� �� ����� 
				// 2. ���� ������� �� ������ �� ���������� ������ ���� �� ��������� ����� �������� ������������.
				// ���� ������� ������ �� ����� ������ �������� foundnow �� true. 
				// �.�. ���������� ������ ������ � �������.
				// 3. � ����� ������ ���������� �������������.
				// ��� ������ ������������� ����.


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
				// 14 ������ 2016 ����.
				// ��������� ���������� ������ ������ �� ��������������� ���������.
				//if (fabs(vector_sum_m[tid][jstr61]) > maxth) maxth = fabs(vector_sum_m[tid][jstr61]);
			//}

	//#if doubleintprecision == 1
			//printf("i=%lld j=%lld aij=%e\n", istr, jstr61, vs161);
	//#else
			//printf("i=%d j=%d aij=%e\n", istr, jstr61, vs161);
	//#endif


			if ((istr == jstr61) && (vs161 < 1.0e-20)) {
				// ������������� ������� �� ���������.

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


				// �������������� ��� ������� ���-�����:
				// Amat. ��������� ����������� ����� ������� ������ ������������� ��������������� ������������� +
				// �������� �� ����� ������������� ���������. ����� �������� �� ���.
				// B. �������� ������������� ��������������� �������������.
				// C. ������ ��������� ������������� ��������������� ������������� (�������������).
				printf("patching string 16.04.2017: \n");
				/*
				for (integer i_62 = 1; i_62 <= index_size_m[tid]; i_62++) {
					integer jstr62 = index_visit_m[tid][i_62];
					if (istr != jstr62) {
						if (vector_sum_m[tid][jstr62] > 0.0) {
							index_visit_m[tid][i_62] = -1; // �� ���������� ������ �������� (�������������).
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
				// ��� ������ ��������� �� ����� ����, ��� ����� ��������� ����� ������������.
				for (integer i_62 = 1; i_62 <= index_size_m[tid]; i_62++) {
					integer jstr62 = index_visit_m[tid][i_62];
					vector_sum_m[tid][jstr62] *= -1.0;
				}
				// ���������� ������ � ������������� ����������.
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

				// ����� �� ����� for �� ���������� i_61.
				break;
			}


		}


		bool bCheck_ok = false; // ��������� ������� ��������� � ������ �������.
								// ��� ����� ����������� ������� ��� ������� �� ������������ ����������.

		

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
					// 7 ������ 2016 ���������� ������ ����:
						if (fabs(vs1) > 1.0e-37) {
							// �� ���������� ������ ����. 
							// �� ������ ������ ��������� ������ ��� ����������.


								// �������������� ���������� ��������.
								// 22_10_2016.
							Ak1 Atemp;
							Atemp.aij = vs1;
							Atemp.i = istr;
							Atemp.j = jstr;

							if (istr == jstr) {

								bCheck_ok = true;

								if (vs1 < 1.0e-20) {
									// ������ ����������� � ���������� ������������� �������� � �������������� ������� �������
									// ������������ ���������. ���� �������� �������� ���� �� ����.
									// 22737
									// ��������� ��� �������� ������ ������ ��������.
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
									// ���������� ������� �������� �������.
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
									// �������������� ���������� �������� ���������.
									Atemp.aij = sum_dia;
									// ������ �������� �� ���������� ���������.
									// 22 ������� 2016
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
			// ������� ���� ������ �������� ��� ������� � � � �� �� ��������� ��� ������.
			// ����� ��������� ����� ���������� ������ ��� � � ������������� ����������.
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

			// ���������� ������� �������� �������.
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
				// �������� ������ ��� �������� ������� �.
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

	// ���������� �� ����� �.�. ������������ schedule (static)
	//Counting_Sort(Amat, istartAnew_mem2, istartAnew2 - 1, false, n_a[ilevel - 1]);
	//if (bprint_mesage_diagnostic) {
		//printf("Counting Sort End. \n");
	//}

	omp_set_num_threads(i_my_num_core_parallelesation);

#endif

}

// ����������� ��������� ��������� ���� ����������� ������.
// �������� �. ���������� IBM 1978.
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

	// ��������� ������ ������� ���������.
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

		// �������� ������������ ����� ������.
		// ����� ������������� ����� � ����.

		node_AVL_Gus* root_Gus = 0;


		// �� ������ ���-�������. 

		// ��������� ������� i-�� ������ �����������
		for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
			integer col_ind = P[ii].j;
			// ��������� col_ind ������ ������� ��������

			// ����� ���������� ������� �� ������� ����.
			doublerealT left_operand = P[ii].aij;

			integer i_1start = row_startA[col_ind];
			integer i_1end = row_startA[col_ind + 1] - 1;
			for (integer i_1 = i_1start; i_1 <= i_1end; i_1++) {

				doublerealT right_operand = Amat.aij[i_1];
				integer iaddind = Amat.j[i_1];
				bool foundnow = false;

				// ����� .
				foundnow = hash_table[iaddind];

				if (foundnow) {
					vector_sum[iaddind] += left_operand * right_operand;
				}
				else {
					// ������ ����������.
					index_size++;
					index_visit[index_size] = iaddind;
					// ������� 							
					hash_table[iaddind] = true;
					vector_sum[iaddind] = left_operand * right_operand;
				}


			}
		}



		doublerealT maxth = -1.0;
		// 22 ������� 2016 �� ���������� ������ �� ����� P*Amat.
		for (integer i_6 = 1; i_6 <= index_size; i_6++) {
			integer jstr = index_visit[i_6];
			hash_table[jstr] = false; // �������������� ���-������� ��� ��������� ��������.
		}

		// ��� ����� ����������� ������� ��� ������� �� ������������
		// ����������.
		bool bCheck_ok = false; // ��������� ������� ��������� � ������ �������.
		integer istartAnew_8 = istartAnew; // ���������� ��� ����������� ������.

		if (nsizeA > istartAnew + index_size) {
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];

				doublerealT vs1 = vector_sum[jstr];

				if ((istr == jstr) && (vs1 > 1.0e-20)) {
					bCheck_ok = true;
				}

				// 22 ������� 2016. ��������� �������� ������ �� ����� P*Amat ������������.


				// 7 ������ 2016 ���������� ������ ����:
				if (fabs(vs1) > 1.0e-37) {
					// �� �� ���������� � ������� ������ ����.
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
			// �������� ������ ��� �������� ������� �.
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

				// 22 ������� 2016. ��������� �������� ������ �� ����� P*Amat ������������.


				// 7 ������ 2016 ���������� ������ ����:
				if (fabs(vs1) > 1.0e-37) {
					// �� �� ���������� � ������� ������ ����.
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


			// ������� ���� ������ �������� ��� ������� � � � �� �� ��������� ��� ������.
			// ����� ��������� ����� ���������� ������ ��� � � ������������� ����������.
			//printf("Amat lack of memory\n");
			//printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			//printf("please, press any key to exit.\n");
			//system("pause");
			//exit(1);
		}

		index_size = 0; // ����� ����������, ������ ����������.
		clear_AVL_Gus(root_Gus);
		root_Gus = 0;

	}

} // my_sparse_multiplication


// ����������� ��������� ��������� ���� ����������� ������.
// �������� �. ���������� IBM 1978.
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

	// ��������� ������ ������� ���������.
	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

		// �������� ������������ ����� ������.
		// ����� ������������� ����� � ����.

		// �� ������ ���-�������. 

		// ��������� ������� i-�� ������ �����������
		for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
			integer col_ind = P[ii].j;
			// ��������� col_ind ������ ������� ��������

			// ����� ���������� ������� �� ������� ����.
			doublerealT left_operand = P[ii].aij;

			integer i_1start = row_startA[col_ind];
			integer i_1end = row_startA[col_ind + 1] - 1;
			for (integer i_1 = i_1start; i_1 <= i_1end; i_1++) {

				doublerealT right_operand = Amat.aij[i_1];
				integer iaddind = Amat.j[i_1];
				bool foundnow = false;

				// ����� .
				foundnow = hash_table[iaddind];

				if (foundnow) {
					vector_sum[iaddind] += left_operand * right_operand;
				}
				else {
					// ������ ����������.
					index_size++;
					index_visit[index_size] = iaddind;
					// ������� 							
					hash_table[iaddind] = true;
					vector_sum[iaddind] = left_operand * right_operand;
				}
			}
		}

		// ��� ����� ����������� ������� ��� ������� �� ������������
		// ����������.		

		if (nsizeA > istartAnew + index_size) {
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];
				hash_table[jstr] = false; // �������������� ���-������� ��� ��������� ��������.


				doublerealT vs1 = vector_sum[jstr];

				// 22 ������� 2016. ��������� �������� ������ �� ����� R*Amat ������������.


				// 7 ������ 2016 ���������� ������ ����:
				if (fabs(vs1) > 1.0e-37) {
					// �� �� ���������� � ������� ������ ����.
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
				hash_table[jstr] = false; // �������������� ���-������� ��� ��������� ��������.


				doublerealT vs1 = vector_sum[jstr];

				// 22 ������� 2016. ��������� �������� ������ �� ����� R*Amat ������������.


				// 7 ������ 2016 ���������� ������ ����:
				if (fabs(vs1) > 1.0e-37) {
					// �� �� ���������� � ������� ������ ����.
					Amat.aij[istartAnew] = vs1;
					Amat.i[istartAnew] = istr;
					Amat.j[istartAnew] = jstr;
					istartAnew++;
				}

			}


			// ������� ���� ������ �������� ��� ������� � � � �� �� ��������� ��� ������.
			// ����� ��������� ����� ���������� ������ ��� � � ������������� ����������.
			//printf("Amat lack of memory\n");
			//printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			//printf("please, press any key to exit.\n");
			//system("pause");
			//exit(1);
		}

		index_size = 0; // ����� ����������, ������ ����������.		

	}

} // my_sparse_multiplication


// ����������� ��������� ��������� ���� ����������� ������.
// �������� �. ���������� IBM 1978.
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

	// ����� ������� ����������� ���������, �
	// ����� ��������� ��������� ����� ������ ������ ��������
	// ���������� ������� ������ ����������.

	for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

		node_AVL_Gus* root_Gus = 0;

		// �� ������ ���-�������.
		// ��������� ��� �������� ������ ������ ��������.
		integer row_ind = Amat.i[row_ind_AS[istr]];// ����� ������.
		for (integer ii1 = row_ind_AS[istr]; ii1 <= row_ind_AE[istr]; ii1++) {

			integer col_ind = Amat.j[ii1];
			doublerealT left_operand = Amat.aij[ii1];

			// ��������� col_ind ������ ������� �������� ���������� �����.
			for (integer ii2 = row_ind_PS[col_ind]; ii2 <= row_ind_PE[col_ind]; ii2++) {

				doublerealT right_operand = P[ii2].aij;

				integer iaddind = P[ii2].i;
				bool foundnow = false;

				// ���������� ����� �� O(1).
				foundnow = hash_table[iaddind];

				if (foundnow) {

					// ������������ �������
					vector_sum[iaddind] += left_operand * right_operand;


				}
				else {
					// ������ ����������.


					index_size++;
					index_visit[index_size] = iaddind;
					// �������
					// ���������� ������� � ���-������� �� O(1).
					hash_table[iaddind] = true;

					vector_sum[iaddind] = left_operand * right_operand;




				}
				// ��������� ����������� ��������� ������:
				// 1. ����� �������� �� ����� 
				// 2. ���� ������� �� ������ �� ���������� ������ ���� �� ��������� ����� �������� ������������.
				// ���� ������� ������ �� ����� ������ �������� foundnow �� true. 
				// �.�. ���������� ������ ������ � �������.
				// 3. � ����� ������ ���������� �������������.
				// ��� ������ ������������� ����.

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
				// 14 ������ 2016 ����.
				// ��������� ���������� ������ ������ �� ��������������� ���������.
				if (fabs(vector_sum[jstr]) > maxth) maxth = fabs(vector_sum[jstr]);
			}
		}



		// huck: 16.04.2017

		for (integer i_61 = 1; i_61 <= index_size; i_61++) {

			integer jstr61 = index_visit[i_61];
			doublerealT vs161 = vector_sum[jstr61];



			if ((istr == jstr61) && (vs161 < 1.0e-20)) {
				// ������������� ������� �� ���������.

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


				// �������������� ��� ������� ���-�����:
				// Amat. ��������� ����������� ����� ������� ������ ������������� ��������������� ������������� +
				// �������� �� ����� ������������� ���������. ����� �������� �� ���.
				// B. �������� ������������� ��������������� �������������.
				// C. ������ ��������� ������������� ��������������� ������������� (�������������).
				printf("patching string 16.04.2017: \n");
				/*
				for (integer i_62 = 1; i_62 <= index_size; i_62++) {
					integer jstr62 = index_visit[i_62];
					if (istr != jstr62) {
						if (vector_sum[jstr62] > 0.0) {
							index_visit[i_62] = -1; // �� ���������� ������ �������� (�������������).
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
				// ��� ������ ��������� �� ����� ����, ��� ����� ��������� ����� ������������.
				for (integer i_62 = 1; i_62 <= index_size; i_62++) {
					integer jstr62 = index_visit[i_62];
					vector_sum[jstr62] *= -1.0;
				}
				if (ibsp_length < i_bsp_LIMIT) {
					// ���������� ������ � ������������� ����������.
					bsp[ibsp_length].ilevel = ilevel;
					bsp[ibsp_length].istring_number = istr;
					ibsp_length++;
				}
				else {
					printf("ERROR!!! ibsp_length>=i_bsp_LIMIT\n");
					system("pause");
					exit(1);
				}
				// ����� �� ����� for �� ���������� i_61.
				break;
			}


		}


		if (nsizeA > istartAnew2 + index_size) {
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {

				integer jstr = index_visit[i_6];
				doublerealT vs1 = vector_sum[jstr];

				// 7 ������ 2016 ���������� ������ ����:
				if ((jstr > -1) && (fabs(vs1) > 1.0e-37)) {
					// �� ���������� ������ ����. 
					// �� ������ ������ ��������� ������ ��� ����������.


						// �������������� ���������� ��������.
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

				// 7 ������ 2016 ���������� ������ ����:
				if ((jstr > -1) && (fabs(vs1) > 1.0e-37)) {
					// �� ���������� ������ ����. 
					// �� ������ ������ ��������� ������ ��� ����������.


					// �������������� ���������� ��������.
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

			// ������� ���� ������ �������� ��� ������� � � � �� �� ��������� ��� ������.
			// ����� ��������� ����� ���������� ������ ��� � � ������������� ����������.
			//printf("Amat lack of memory\n");
			//printf("yuo mast increase the size of the matrix Amat and restart solver\n");
			//printf("please, press any key to exit.\n");
			//system("pause");
			//exit(1);
		}

		index_size = 0;
		// ������������ ������ �� ��� ��� ������.
		clear_AVL_Gus(root_Gus);
		root_Gus = 0;
	}

}

#endif
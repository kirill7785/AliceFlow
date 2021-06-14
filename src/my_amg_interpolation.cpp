
#pragma once
#ifndef MY_AMG_INTERPOLATION_CPP 
#define MY_AMG_INTERPOLATION_CPP 1


// 3 ������ 2016.
// 17.02.2019
// ���������� �������� ���� � ������ �������� ���������� ������ ��� �������� 
// ������������.
// ��� ���������� ����� ����� - �� ����� ������ �����.
// ����� ���������� ������������� ����������� ������ � ���������� �������� nsizePR.
void deallocate_prolongation(integer &nsizePR, // ������ ��� P � ����������� n.
	integer n, // ���������� ����������� � ���� �� ������� ������.
	Ak1* &P // prolongation
)
{

	// ��������� ������:
	//deallocate_prolongation(nsizePR,n,R,P);

	// ��� ����������� �������� ����� ��� ���������� ������ ������ � nsizePR==12.
	// ������ ����������� ������ ��� ������� ����� ����������� ������� ���������� ������.
	// � ���� ������ � ��������� ��������� � ���������� �� nsizePR==35.
	if (nsizePR < 13) {

		Ak1* P_copy = NULL; // prolongation

		// ����������� ���������� ������ �� nsizePR==35.

		//Ak1* &P_copy = new Ak1[nsizePR*n + 1];

		//P = new Ak1[(integer)(35 * n) + 1]; // 3*nnz 2.4 // 35
		P_copy = (Ak1*)malloc(((nsizePR * n) + 1) * sizeof(Ak1));
		if (P_copy == NULL) {
			// ������������ ������ �� ������ ������������.
			printf("Problem: not enough memory on your equipment for P_copy matrix in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}

		if ((P_copy != NULL)) {
			
			if (P != NULL) {
				for (integer i = 0; i <= (nsizePR * n); i++) {
					P_copy[i] = P[i];
				}
				free(P);
			}
			else {
				printf("error dealocate_prolongation function P is null. my_amg_interpolation.cpp\n");
				system("pause");
				exit(1);
			}
			P = NULL;
			integer nsizePR_old = nsizePR;
			nsizePR = 35; // ����� ������ ������� � �������.


			//P = new Ak1[(integer)(35 * n) + 1]; // 3*nnz 2.4 // 35
			P = (Ak1*)malloc(((nsizePR * n) + 1) * sizeof(Ak1));
			if (P == NULL) {
				// ������������ ������ �� ������ ������������.
				printf("Problem: not enough memory on your equipment for P matrix in my_agregat_amg.cpp algorithm...\n");
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

			// ����������� ����������� ����������� ������.
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


// ������������ ������������. 18.09.2020.
#include "my_long_range_amg_interpolation.cpp"

// ���������������� ��������� �3.
// ���������������� ������� ������������� �� ������� ����.
// �� ���� ������������ �� ������� ���� �������� ������� ����.
// 12.05.2018
template <typename doublerealT>
void my_interpolation_procedure_number3A(integer &the_number_of_neighbors_that_are_not_C_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak1* &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublerealT RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublerealT &theta, integer &n, Ak1* &R, integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublerealT &magic82, doublerealT* &threshold_quick_all,
	doublerealT* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = theta83; // 3 ���� 2016


	// �������� ������ ������������.

	// ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
	// ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
	// ������� F �����.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	// the begining

	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if (1) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;		
			

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {


				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
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


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;
				
					//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
					
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
								//if (this_is_C_node[Amat[is0].j]  ) {
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
					// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

					
					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
						if (Amat[is0].j != Amat[ii1].i) {
							if (this_is_C_node[Amat[is0].j]  ) {
								//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
								if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
									sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat[is0].j]  ) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
										SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}
						}
					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}
				}
				else {
					// only negative connections

					
						// ��� ����� ����������� �������� ��������.
						// 5 ������� 2015 ���� �� ��������� ��������� �������������
						// ��������� ������������ � ������ � ��������� ��������.
						doublerealT maxelem_threshold = -1.0;
						
						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
		
						if (!btreshold_on_new_vetv) {
							for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								if (Amat[is0].j != Amat[ii1].i) {
									// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
									//if (this_is_C_node[Amat[is0].j]  ) {
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
						// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

						
						doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
						doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								if (this_is_C_node[Amat[is0].j]  ) {
									//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
										sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
										icsos++;
									}
								}
								else {
									if (this_is_F_node[Amat[is0].j]  ) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
											SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
											icount_StronglyF++;
											icsosF++;
										}
									}
									// ������������ ���������� ������� ������� �� �������� � ������.
									the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
								}
							}
						}
						if (icsos == 1) {
							number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																				   // ��������� ������ ������ "����������".
																				   // ���������� ������ ����������� ��� ���������.
																				   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																				   // ������������ �� ���� ������� ����� ��������.
							if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
						}

				}

				/*
				9.05.2018
				1. ���������� threshold�.
				2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
				3. icsos - ���������� strong C. icsosF - ���������� strong F.
				*/

				// 1 ������ 2016 ���� ����� ��� ������������.
				// ������� � ������ ������ ������ ������ ����� ���� � �����.


				// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
				//integer iend_marker_position = row_startA[Amat[ii1].i + 1] - 1;
				for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
					if (this_is_F_node[Amat[is0].j]  ) {
						if (Amat[is0].j != Amat[ii1].i) {


							// 20 jan 2016.
							// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.


							if (fabs(sumP) < RealZERO) {
								// ��� ������ ����� ������ ��� ������� � �������.

								// if (icsosF==0) ������ ������ ����������.
								//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
								{
								//printf("error interpolation zero diagonal sumP.\n");
								//printf("Fnode all sosed is F");
								//system("pause");
								//printf("i8 is Dirichlet node\n");
							    if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // ���� ���������� �����.
												 //exit(1);
												 // ����� ����� �������� ������������.
							    }
							}


						}
					}
				}//end

			}

	}
	




	/*
	9.05.2018
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
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

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
		// �������� ������� �� ���� ���������
		// �� �� �������������� �� �� ����� V ������.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82; // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat[is0].j]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j]  ) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j]  ) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				 {

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat[is0].j]  ) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;											
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat[is0].j]  ) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat[is0].j;
											if (fabs(SumPall) < RealZERO) {
												printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
												//getchar();
												system("PAUSE");
												exit(1);
											}
											doublerealT multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat[ii1_loc].i;
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc) {
													// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
													// ����� ������ ����������.
													if (this_is_C_node[Amat[is0_loc].j]  ) {
														if (Amat[is0_loc].j != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
														}
													}
												}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat[is0_loc].j]  ) {

														if (Amat[is0_loc].j != istr_etalon) {

															sumP_loc += fabs(Amat[is0_loc].aij); // ����� ������� ��������������� ��������� ������� ����������� � �����.
															icsos_loc++;
														}

													}
													else {

														//if (Amat[is0_loc].j != istr_etalon) {
														// ������������ ���������� ������� ������� �� �������� � ������.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
														//}
													}

												}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {
													
													if (this_is_C_node[Amat[is0_loc].j]  ) {
														if (Amat[is0_loc].j != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat[is0_loc].j];
															//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
															if (fabs(sumP_loc) < RealZERO) {
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

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat[is0].j]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j]  ) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j]  ) {
								//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
								if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				 {

					if (1) {
					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat[is0].j]  ) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ �-�����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat[is0].j]  ) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat[is0].j;
											doublerealT multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat[ii1_loc].i;
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat[is0_loc].j]  ) {
														if (Amat[is0_loc].j != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
														}
													}
												}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
													
													if (this_is_C_node[Amat[is0_loc].j]  ) {

														if (Amat[is0_loc].j != istr_etalon) {

															sumP_loc += fabs(Amat[is0_loc].aij); // ����� ������� ��������������� ��������� ������� ����������� � �����.
															icsos_loc++;
														}

													}
													else {

														//if (Amat[is0_loc].j != istr_etalon) {
														// ������������ ���������� ������� ������� �� �������� � ������.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
														//}
													}

												}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
													
													if (this_is_C_node[Amat[is0_loc].j]  ) {
														if (Amat[is0_loc].j != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat[is0_loc].j];
															//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
															if (fabs(sumP_loc) < RealZERO) {
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

// ���������������� ��������� �3.
// ���������������� ������� ������������� �� ������� ����.
// �� ���� ������������ �� ������� ���� �������� ������� ����.
// 12.05.2018
template <typename doublerealT>
void my_interpolation_procedure_number3A(integer &the_number_of_neighbors_that_are_not_C_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak2 &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublerealT RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublerealT &theta, integer &n,  integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublerealT &magic82, doublerealT* &threshold_quick_all,
	doublerealT* &threshold_quick_only_negative, bool &bsuffix_work)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = theta83; // 3 ���� 2016


	// �������� ������ ������������.

	// ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
	// ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
	// ������� F �����.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	// the begining

	// ����� � ���� ����������� ����� ��������� ����� ���������� �����������
	// ���������� �������. 359 V ������ ������ 27. �.�. ���������� 
	// ������������� ����� ��� �� �������.
	

	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if ((bsuffix_work)&&(1)) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


		// ���������� ����������� ��� ����� ������� ���������� F nodes.
		// ������ F-nodes ������ C-nodes.
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {


			// ����� ������� ������� F-node ������� C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
			// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
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


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
					// ��������� ������ ������ "����������".
					// ���������� ������ ����������� ��� ���������.
					// � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
					// ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}
			}
			else {
				// only negative connections


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
					// ��������� ������ ������ "����������".
					// ���������� ������ ����������� ��� ���������.
					// � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
					// ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}

			}

			/*
			9.05.2018
			1. ���������� threshold�.
			2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
			3. icsos - ���������� strong C. icsosF - ���������� strong F.
			*/

			// 1 ������ 2016 ���� ����� ��� ������������.
			// ������� � ������ ������ ������ ������ ����� ���� � �����.


			// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
			//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
			//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				//if (this_is_F_node[Amat.j[is0]]  ) {
					//if (Amat.j[is0] != Amat.i[ii1]) {


						// 20 jan 2016.
						// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.
			if ((iend_marker_position > ii1)&&
				(fabs(sumP) < RealZERO)) {
				// ���� ���� ���� ���� �����.

				
				// ��� ������ ����� ������ ��� ������� � �������.

				// if (icsosF==0) ������ ������ ����������.
				//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
				{
					//printf("error interpolation zero diagonal sumP.\n");
					//printf("Fnode all sosed is F");
					//system("pause");
					//printf("i8 is Dirichlet node\n");
					if (this_is_C_node[i8] == false) iadditionalCstatistic++;
					this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
					this_is_C_node[i8] = true;
					bweSholdbeContinue = true;
					byes_add = true; // ���� ���������� �����.
					//exit(1);
					// ����� ����� �������� ������������.
				}
				
			}


					//}
				//}
			//}//end

		}

	}





	/*
	9.05.2018
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
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

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
		// �������� ������� �� ���� ���������
		// �� �� �������������� �� �� ����� V ������.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82; // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F-nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
					// ��������� ������ ������ "����������".
					// ���������� ������ ����������� ��� ���������.
					// � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
					// ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
								if (this_is_F_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
										// ������������� Strong F �����.

										// �����:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
										// � ������� F ������.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];
										if (fabs(SumPall) < RealZERO) {
											printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
											//getchar();
											system("PAUSE");
											exit(1);
										}
										doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
										// ��������� ���� ������� iFpointeger 
										// ����� ����� ����� ��� � ����.

										// �������������� ��������� �����.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
											iend_for = is0_loc;
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
												// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
												// ����� ������ ����������.
												if (this_is_C_node[Amat.j[is0_loc]]  ) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

										// ����� ������� ������� F-node ������� C-node.
										integer icsos_loc = 0;

										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
										doublerealT sumP_loc = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// ����� �������� ����� ���������� ����� ���������.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]]  ) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// ������������ ���������� ������� ������� �� �������� � ������.
													//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
													//}
												}

											}
										}

										doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										// � ����� ��� ������� ���������������� ����� 
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]]  ) {
													if (Amat.j[is0_loc] != istr_etalon) {

														P[icount1].j = i8;
														P[icount1].i = C_numerate[Amat.j[is0_loc]];
														//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
														if (fabs(sumP_loc) < RealZERO) {
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

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
					// ��������� ������ ������ "����������".
					// ���������� ������ ����������� ��� ���������.
					// � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
					// ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if ((Amat.aij[is0] <0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
								if (this_is_F_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
										// ������������� Strong F �����.

										// �����:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
										// � ������� F ������.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];
										doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
										// ��������� ���� ������� iFpointeger 
										// ����� ����� ����� ��� � ����.

										// �������������� ��������� �����.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
											iend_for = is0_loc;
											if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												if (this_is_C_node[Amat.j[is0_loc]]  ) {
													if (Amat.j[is0_loc] != istr_etalon) {
														maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
													}
												}
											}
										}

										doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
										// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

										// ����� ������� ������� F-node ������� C-node.
										integer icsos_loc = 0;

										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
										doublerealT sumP_loc = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// ����� �������� ����� ���������� ����� ���������.
											if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]]  ) {

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_loc += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
														icsos_loc++;
													}

												}
												else {

													//if (Amat.j[is0_loc] != istr_etalon) {
													// ������������ ���������� ������� ������� �� �������� � ������.
													//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
													//}
												}

											}
										}

										doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										// � ����� ��� ������� ���������������� ����� 
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

											// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]]  ) {
													if (Amat.j[is0_loc] != istr_etalon) {

														P[icount1].j = i8;
														P[icount1].i = C_numerate[Amat.j[is0_loc]];
														//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
														if (fabs(sumP_loc) < RealZERO) {
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




// ��� PMIS ������������.
template <typename doublerealT>
doublerealT activation_function_for_thrteshold(doublerealT threshold, integer ilevel, bool bmagic)
{

	// ���� ������������� ����� ������� ��������� �� ��� ���� ����� ���������� ��������� ���������
	// ����� ����������� ���������� ����� ����������� �� �������� <50
	// ���������� ������� �������� threshold � �� ����� ��� ����������,
	// ��������� �� �������� ���� � ������� ������ �� ������ ������ �� �������� 0.25.

	// ���������� theta �������� � ����� �������� ���������� ��������� ������������,
	// ������� ����� ���������� P ���������� ������ ������� ��������� �������� 
	// ���������������� ������.

	// ��� �������� �������� ���������� 64 ������ ������ � ���� (12��� �����������).

	if (bmagic) {
		//0.4; 0.38.
		if (ilevel == 1) {
			return (doublerealT)(0.55);// 0.55;//0.42; 
		}
		else if (ilevel == 2) {
			return (doublerealT)(0.55);// 0.55;
		}
		else if (ilevel == 3) {
			return (doublerealT)(0.21); // 0.18; 0.21 optimum
		}
		else {
			if (10 < ilevel) {
				return (doublerealT)(pow(0.38, pow(10.0, 0.55)));//0.52 optimum; 0.54; 0.55; 0.58; 0.61
			}
			else {
				return (doublerealT)(pow(0.38, pow(1.0*ilevel, 0.55)));//0.52 optimum; 0.54; 0.55; 0.58; 0.61
			}
		}
	}
	else {

		// ������������ iwamura ������ ������������� ������.
		if (ilevel <= 2) {
			return (doublerealT)(0.98); // iwamura
		}
		else if (ilevel == 3) {
			return (doublerealT)(0.75); // iwamura
		}
		else {
			//return threshold;
			//return pow(threshold, min(5, ilevel));
			// ����������� ��� ����������� ��������� ����� � �������� �����
			// ����� ���� �������� �� 10 ������.
			//return pow(0.24, pow(min(10, ilevel), 0.51));
			if (10 < ilevel) {
				return (doublerealT)(pow(0.24, pow(10.0, 0.64)));//0.51; 0.64 optimum
			}
			else {
				return (doublerealT)(pow(0.24, pow(1.0*ilevel, 0.64)));//0.51; 0.64 optimum
			}
			//return pow(threshold, pow(ilevel, 0.5));
			//return pow(threshold, pow(min(10, ilevel), d_my_optimetric1_6_12_2019));//optimetric
			//0.51 optimum.
		}
	}

} // activation_function_for_thrteshold

  // ���������������� ��������� �3 ��� PMIS.
  // ���������������� ������� ������������� �� ������� ����.
  // �� ���� ������������ �� ������� ���� �������� ������� ����.
  // 12.05.2018
template <typename doublerealT>
void my_interpolation_procedure_number3A_PMIS0(integer &the_number_of_neighbors_that_are_not_C_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak2 &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublerealT RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublerealT &theta_gl, integer &n, integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublerealT &magic82, doublerealT* &threshold_quick_all,
	doublerealT* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = activation_function_for_thrteshold(theta83, ilevel,false); // 3 ���� 2016
	doublerealT theta = activation_function_for_thrteshold(theta_gl, ilevel,false);


										 // �������� ������ ������������.

										 // ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
										 // ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
										 // ������� F �����.

										 //6interpolation 0.4 6.77 11 26 28.355
										 //6interpolation 0.45 6.6 10 27 28.151
										 //6interpolation 0.5 6.42 12 32 28.735
										 //4interpolation 0.4 3.7  52 24.736 // best
										 //4interpolation 0.3 3.78 13 59 27.525
										 //4interpolation 0.5 3.61 12 55 25.533
										 //4interpolation 0.45 3.65 10 63 30.24

										 // the begining

	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if (1) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


		// ���������� ����������� ��� ����� ������� ���������� F nodes.
		// ������ F-nodes ������ C-nodes.
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {


			// ����� ������� ������� F-node ������� C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
			// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
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


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}
			}
			else {
				// only negative connections


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}

			}

			/*
			9.05.2018
			1. ���������� threshold�.
			2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
			3. icsos - ���������� strong C. icsosF - ���������� strong F.
			*/


			// 1 ������ 2016 ���� ����� ��� ������������.
			// ������� � ������ ������ ������ ������ ����� ���� � �����.


			// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
			//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
			for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				if (this_is_F_node[Amat.j[is0]]  ) {
					if (Amat.j[is0] != Amat.i[ii1]) {


						// 20 jan 2016.
						// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.

						//if ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO)) 
						if (fabs(sumP) < RealZERO) 
						{
							// ��� ������ ����� ������ ��� ������� � �������.

							// if (icsosF==0) ������ ������ ����������.
							//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
							{
								//printf("error interpolation zero diagonal sumP.\n");
								//printf("Fnode all sosed is F");
								//system("pause");
								//printf("i8 is Dirichlet node\n");
								if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // ���� ���������� �����.

								// ���� ����������� ����������� ������������� ������ ��������:
								//printf("error!!! ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO))\n");
								//system("pause");
								//exit(1);
												 // ����� ����� �������� ������������.
							}
						}


					}
				}
			}//end

		}

	}


	if (byes_add) {
		// ���� ����������� ����������� ������������� ������ ��������:
		printf("number (fabs(sumP) < RealZERO) = %lld; n_a=%lld.\n", iadditionalCstatistic, n_a[ilevel - 1]);
		//system("pause");
		//exit(1);
	}


	/*
	9.05.2018
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
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

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
		// �������� ������� �� ���� ���������
		// �� �� �������������� �� �� ����� V ������.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = activation_function_for_thrteshold(magic82,ilevel,true); // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
										/*
										// �� ����� ����� ������ � ����. � ����������� ���� ����� ����� ������� F ���� �� � ����.
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else */
										{
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat.j[is0]]  ) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											if (fabs(SumPall) < RealZERO) {
												printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
												//getchar();
												system("PAUSE");
												exit(1);
											}
											doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.
											
											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if (Amat.j[is0_loc] != istr_etalon) {
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
														// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
														// ����� ������ ����������.
														if (this_is_C_node[Amat.j[is0_loc]]  )
														{
															if (Amat.j[is0_loc] != istr_etalon) {
																maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
															}
														}
													}
												}
											}

											//doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
															icsos_loc++;
														}

													}
													else {

														//if (Amat.j[is0_loc] != istr_etalon) {
														// ������������ ���������� ������� ������� �� �������� � ������.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
														//}
													}

												}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (Amat.j[is0_loc] != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat.j[is0_loc]];
															//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
															if (fabs(sumP_loc) < RealZERO) {
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
			////
		}
		else {
			// only negative connections.

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
								sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
									SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
										/*
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else*/
										{
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if ((Amat.aij[is0] <0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat.j[is0]]  ) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if ((Amat.aij[is0]<0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											doublerealT multiplyer_nu = fabs(Amat.aij[is0]) / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = fabs(Amat.aij[is0_loc]);
														}
													}
												}
											}

											//doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
															icsos_loc++;
														}

													}
													else {

														//if (Amat.j[is0_loc] != istr_etalon) {
														// ������������ ���������� ������� ������� �� �������� � ������.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
														//}
													}

												}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (Amat.j[is0_loc] != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat.j[is0_loc]];
															//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
															if (fabs(sumP_loc) < RealZERO) {
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


} // my_interpolation_procedure_number3A_PMIS


 // ���������������� ��������� �3 ��� PMIS.
  // ���������������� ������� ������������� �� ������� ����.
  // �� ���� ������������ �� ������� ���� �������� ������� ����.
  // 12.05.2018
template <typename doublerealT>
void my_interpolation_procedure_number3A_PMIS(integer& the_number_of_neighbors_that_are_not_C_nodes,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative)
{
	// ���������� ���������� ����� � ����� 01.10.2020
	// true - ����������� ����� ��������, �� ��������� ����������� ����������� ������.
	bool weakening_the_addition_of_new_C_nodes = true;// true;

	bool bprint_mesage_diagnostic = true;
	if (my_amg_manager.iprint_log == 0) {
		bprint_mesage_diagnostic = false;
	}

	if (!b_on_adaptive_local_refinement_mesh) {
		// �� ����������������� ����� ���������� ������ ������ � ����.
		// ������� ��������� ������ �����������.
		// �� ���� ����� ������ ����������� �������� �� 10% ����������� ��������� ��� ������ 
		// �� ������� �������. ������� ����������� �������� � ���������� ����� �������� �� ����������.
		weakening_the_addition_of_new_C_nodes = false;// ����������� ����� ��������� ���� �����.
		//if (ilevel < 4) weakening_the_addition_of_new_C_nodes = false;// ������ 4, ��������� ������ ����.
	}

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = activation_function_for_thrteshold(theta83, ilevel, false); // 3 ���� 2016
	doublerealT theta = activation_function_for_thrteshold(theta_gl, ilevel, false);

	// ���-�������
	bool* this_is_Strong_C_node = new bool[n_a[ilevel - 1] + 1];
	for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
		// ������������� ���-�������.
		this_is_Strong_C_node[i_1] = false;
	}

	integer ioneStrongC_and_0_StrongF = 0;

	// �������� ������ ������������.

	// ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
	// ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
	// ������� F �����.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	// the begining
    // ��������� ������������ ����� ������� � ������ �� ���� �������. 
	integer imax_count_sosed = -1;
	for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode]  ) {
		integer ii1 = row_startA[inode];
		//integer iend_1 = row_startA[Amat.i[ii1] + 1] - 1;
		integer iend_1 = row_startA[inode + 1] - 1;
		if (iend_1 - ii1 > imax_count_sosed) {
			imax_count_sosed = iend_1 - ii1;
		}
	}
	integer* iStrongC_nodes_ind = new integer[imax_count_sosed + 1];



	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if (1) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;
		integer icount_bad_string = 0;


		if (bpositive_connections) {
			// ���������� ����������� ��� ����� ������� ���������� F nodes.
		// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				bool badd_amg1r5 = false;
				bool bfound_vneDiag_F_node = false;

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;

				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
				integer iend_1 = iend_marker_position;

				integer isize_iStrongC_nodes_ind = 0;





#if doubleintprecision == 1
				//printf("i8=%lld n=%lld\n", i8, n_a[ilevel - 1]);
#else
				//printf("i8=%d n=%d\n", i8, n_a[ilevel - 1]);
#endif

				//getchar();


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
							if (aij_abs > maxelem_threshold) {
								maxelem_threshold = aij_abs;
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[Amat.i[ii1]];
				}
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						//doublerealT aij_abs = fabs(Amat.aij[is0]);
						doublerealT aij_abs = Amat.abs_aij[is0];
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (aij_abs > maxelem_threshold*barjer) {
							//if (aij_abs > maxelem_threshold*theta) {
							if (aij_abs > maxelem_threshold_theta) {
								sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.

								// ������� ���� ������� � �������.
								// ������������� ���-�������.
								this_is_Strong_C_node[Amat.j[is0]] = true;

								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;

								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (aij_abs > maxelem_threshold*theta_strong_F) {
								if (aij_abs > maxelem_threshold_theta_strong_F) {
									SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}

				// ����� � ������������ amg1r5 ���� ��������� ���������� ����. �� �� ����
				// ������������ ������ ��� ������� � ������������ ����������� 
				// ��� FF1 ������������.
				// ���� � ���� ������ ��������� FF ���� �� ���� ������ � ������ �� ����� ����������� ���������� ������ ����.
				// � PMIS � ��� ����� ���� ��� ������ ��������� FF ���� �� ������� ������ � ������.
				// � ���� ������ ������ ������������ ����������� ��
				// distance 2.

				// ����� - ������� ���-�������.
				for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
					this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
				}

				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				/*
				9.05.2018
				1. ���������� threshold�.
				2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
				3. icsos - ���������� strong C. icsosF - ���������� strong F.
				*/

				// 1 ������ 2016 ���� ����� ��� ������������.
				// ������� � ������ ������ ������ ������ ����� ���� � �����.


				
				bool baddC = true;
				for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
					if (this_is_C_node[Amat.j[is0]]) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (weakening_the_addition_of_new_C_nodes) {
								// ���� ��� �� ������ �������� � ������
								// �� ��� ������� ���������� �� �����������
								// ����� �������� ����� F ���� ���� ������ � �����.
								// ����� � ������ ������� � F ���� ��� ������� � �������
								// ������ ������ ������� � ������ (������ �������� �� ������)
								// �� �� ��������������� �� �����
								// ������� � ������ � �� ��������� ����� � �����.
								// ��� �������� �������� ���������� � ��� ����.

								baddC = false;
							}
						}
					}
				}

				if (baddC) {

					if (sumP < RealZERO)
					{
						// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
						//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
						//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							//if (this_is_F_node[Amat.j[is0]]  ) {
								//if (Amat.j[is0] != Amat.i[ii1]) {

						if (iend_marker_position > ii1) {
							// ���� ���� ���� ���� ����� F ����.


							// 20 jan 2016.
							// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.

							//if ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO)) 

								// ��� ������ ����� ������ ��� ������� � �������.

								// if (icsosF==0) ������ ������ ����������.
								//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
							{
								//printf("error interpolation zero diagonal sumP.\n");
								//printf("Fnode all sosed is F");
								//system("pause");
								//printf("i8 is Dirichlet node\n");
								if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // ���� ���������� �����.
								//continue;

								// ���� ����������� ����������� ������������� ������ ��������:
								//printf("error!!! ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO))\n");
								//system("pause");
								//exit(1);
								 // ����� ����� �������� ������������.
							}


							//}
						//}
					//}//end
						}
					}
				}

			}
	    }
		else {
			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				bool badd_amg1r5 = false;
				bool bfound_vneDiag_F_node = false;

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;

				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
				integer iend_1 = iend_marker_position;

				integer isize_iStrongC_nodes_ind = 0;


			
				// only negative connections


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
							if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold)) {
								maxelem_threshold = aij_abs;
							}
							//}
						}
					}
				}
				else {
					maxelem_threshold = threshold_quick_only_negative[Amat.i[ii1]];
				}
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					if (Amat.j[is0] != Amat.i[ii1]) {
						//doublerealT aij_abs = fabs(Amat.aij[is0]);
						doublerealT aij_abs = Amat.abs_aij[is0];
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (aij_abs > maxelem_threshold*barjer) {
							//if (aij_abs > maxelem_threshold*theta) {
							if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
								sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
								this_is_Strong_C_node[Amat.j[is0]] = true;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (aij_abs > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta_strong_F)) {
									SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}

				// ����� - ������� ���-�������.
				for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
					this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
				}

				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}

			

				/*
				9.05.2018
				1. ���������� threshold�.
				2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
				3. icsos - ���������� strong C. icsosF - ���������� strong F.
				*/


				// 1 ������ 2016 ���� ����� ��� ������������.
				// ������� � ������ ������ ������ ������ ����� ���� � �����.

				bool baddC = true;
				for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
					if (this_is_C_node[Amat.j[is0]]) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (weakening_the_addition_of_new_C_nodes) {
								// ���� ��� �� ������ �������� � ������
								// �� ��� ������� ���������� �� �����������
								// ����� �������� ����� F ���� ���� ������ � �����.
								// ����� � ������ ������� � F ���� ��� ������� � �������
								// ������ ������ ������� � ������ (������ �������� �� ������)
								// �� �� ��������������� �� �����
								// ������� � ������ � �� ��������� ����� � �����.
								// ��� �������� �������� ���������� � ��� ����.

								baddC = false;
							}
						}
					}
				}

				if (baddC) {
					// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
					if (sumP < RealZERO)
					{
						if (iend_marker_position > ii1) {
							// ���� ���� ���� ���� ����� F ����.

							//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
							//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								//if (this_is_F_node[Amat.j[is0]]  ) {
									//if (Amat.j[is0] != Amat.i[ii1]) {


										// 20 jan 2016.
										// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.

										//if ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO)) 

											// ��� ������ ����� ������ ��� ������� � �������.

											// if (icsosF==0) ������ ������ ����������.
											//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
							{
								//printf("error interpolation zero diagonal sumP.\n");
								//printf("Fnode all sosed is F");
								//system("pause");
								//printf("i8 is Dirichlet node\n");
								if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // ���� ���������� �����.
								//continue;

								// ���� ����������� ����������� ������������� ������ ��������:
								//printf("error!!! ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO))\n");
								//system("pause");
								//exit(1);
											 // ����� ����� �������� ������������.
							}



							//}
						//}
					//}//end
						}
					}
				}
			}
		}		

	}


	if (byes_add) {
		if (bprint_mesage_diagnostic) {
			// ���� ����������� ����������� ������������� ������ ��������:
			printf("number (fabs(sumP) < RealZERO) = %lld; n_a=%lld.\n", iadditionalCstatistic, n_a[ilevel - 1]);
		}
		//system("pause");
		//exit(1);
	}


	/*
	9.05.2018
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
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

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
		// �������� ������� �� ���� ���������
		// �� �� �������������� �� �� ����� V ������.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = activation_function_for_thrteshold(magic82, ilevel, true); // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				integer isize_iStrongC_nodes_ind = 0;

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
							if (aij_abs > maxelem_threshold) {
								maxelem_threshold = aij_abs;
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						//doublerealT aij_abs = fabs(Amat.aij[is0]);
						doublerealT aij_abs = Amat.abs_aij[is0];
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (aij_abs > maxelem_threshold*barjer) {
							//if (aij_abs > maxelem_threshold*theta) {
							if (aij_abs > maxelem_threshold_theta) {
								sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;

								this_is_Strong_C_node[Amat.j[is0]] = true;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;

							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (aij_abs > maxelem_threshold*theta_strong_F) {
								if (aij_abs > maxelem_threshold_theta_strong_F) {
									SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
					bool ivariant=0; // 0 - ������� �������(�����), 1 - �������������� ������.
					bool baddC = true;
					doublerealT Cmax = -1.0; // �������� ������� � ����� �� ������ � ������.
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
					if (!baddC)
					{
						maxelem_threshold_theta = 0.98*Cmax;
						if (ivariant == 1) {
							maxelem_threshold_theta = -1.0;
						}
						for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
							if (this_is_C_node[Amat.j[is0]]) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									doublerealT aij_abs = Amat.abs_aij[is0];
									if ((ivariant==1)||(aij_abs > 0.98*Cmax))
									{
										// ���� (ivariant==1) �� ������� ��������� � �� ��������� ��� � ����.

										sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly C �����.
										icsos++;
										this_is_Strong_C_node[Amat.j[is0]] = true;
										iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
										isize_iStrongC_nodes_ind++;
									}
								}
							}
						}
					}

				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.
								//doublerealT aij_abs = fabs(Amat.aij[is0]);
								doublerealT aij_abs = Amat.abs_aij[is0];

								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (aij_abs > maxelem_threshold*theta) {
									if (aij_abs > maxelem_threshold_theta) {
										/*
										// �� ����� ����� ������ � ����. � ����������� ���� ����� ����� ������� F ���� �� � ����.
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else */
										{
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (aij_abs > maxelem_threshold*barjer) {
											//if (aij_abs > maxelem_threshold*theta) {
											//if (aij_abs > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = aij_abs / sumP;
												if (SumPall < RealZERO) {
													printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = aij_abs / SumPall;
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

								}
								else
									if (this_is_F_node[Amat.j[is0]]  ) {
										//if (aij_abs > maxelem_threshold*theta_strong_F) {
										if (aij_abs > maxelem_threshold_theta_strong_F) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (aij_abs > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (aij_abs > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											if (SumPall < RealZERO) {
												printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
												//getchar();
												system("PAUSE");
												exit(1);
											}
											doublerealT multiplyer_nu = aij_abs / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												iend_for = is0_loc;
												if (Amat.j[is0_loc] != istr_etalon) {
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
													if (aij_loc_abs > maxelem_threshold_loc) {
														// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
														// ����� ������ ����������.
														if (this_is_C_node[Amat.j[is0_loc]]  )
														{
															if (Amat.j[is0_loc] != istr_etalon) {
																maxelem_threshold_loc = aij_loc_abs;
															}
														}
													}
												}
											}

											//doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											bool bfound_amg1r5 = false;
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												if (Amat.abs_aij[is0_loc] > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {

														if (Amat.j[is0_loc] != istr_etalon) {
															if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
																bfound_amg1r5 = true;
															}
														}
													}
												}
											}

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;

											if (bfound_amg1r5) {
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
													// ����� �������� ����� ���������� ����� ���������.
													if (aij_loc_abs > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {

															if (Amat.j[is0_loc] != istr_etalon) {
																if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
																	sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																	icsos_loc++;
																}
															}

														}
														else {

															//if (Amat.j[is0_loc] != istr_etalon) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}
														}

													}
												}
											}
											else {

												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
													// ����� �������� ����� ���������� ����� ���������.
													if (aij_loc_abs > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {

															if (Amat.j[is0_loc] != istr_etalon) {

																sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																icsos_loc++;
																break; // ���������� ������ ��������� ����������.
															}

														}
														else {

															//if (Amat.j[is0_loc] != istr_etalon) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}
														}

													}
												}
											}  

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											if (bfound_amg1r5) {
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

													// ��������������� ������� �� ��������� � �����.

													// ������ ������� ������ ����������� ��������� 
													// ������������� ��������� �������� �� �������� 
													// �������� �������.
													// ����������� 5 ������� 2015.
													if (aij_loc_abs > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {
															if (Amat.j[is0_loc] != istr_etalon) {

																if (this_is_Strong_C_node[Amat.j[is0_loc]]) {

																	P[icount1].j = i8;
																	P[icount1].i = C_numerate[Amat.j[is0_loc]];
																	//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																	if (sumP_loc < RealZERO) {
																		printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//getchar();
																		system("PAUSE");
																		exit(1);
																	}
																	P[icount1].aij = multiplyer_nu * aij_loc_abs / sumP_loc;
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
											else {
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

													// ��������������� ������� �� ��������� � �����.

													// ������ ������� ������ ����������� ��������� 
													// ������������� ��������� �������� �� �������� 
													// �������� �������.
													// ����������� 5 ������� 2015.
													if (aij_loc_abs > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {
															if (Amat.j[is0_loc] != istr_etalon) {

																P[icount1].j = i8;
																P[icount1].i = C_numerate[Amat.j[is0_loc]];
																//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																if (sumP_loc < RealZERO) {
																	printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																P[icount1].aij = multiplyer_nu * aij_loc_abs / sumP_loc;
																icount1++;
																if (icount1 >= nsizePR * n) {
																	printf("memory error!!!\n");
																	printf("not enough memory for the interpolation operator.\n");
																	//system("PAUSE");
																	//exit(1);
																	deallocate_prolongation(nsizePR, n, P);
																}
																break;// �� ����� ������ ��������� �����������.

															}
														}
													}
												}
												//end
											}

										}
									}
							
									
                            }

						}
					}

				}

				// ����� - ������� ���-�������.
				for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
					this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
				}

			}
			////
		}
		else {
			// only negative connections.

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				integer isize_iStrongC_nodes_ind = 0;

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
							if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold)) {
								maxelem_threshold = aij_abs;
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						//doublerealT aij_abs = fabs(Amat.aij[is0]);
						doublerealT aij_abs = Amat.abs_aij[is0];
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (aij_abs > maxelem_threshold*barjer) {
							//if (aij_abs > maxelem_threshold*theta) {
							if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
								sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;

								this_is_Strong_C_node[Amat.j[is0]] = true;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
								//if (aij_abs > maxelem_threshold*theta_strong_F) {
								if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta_strong_F)) {
									SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
									icount_StronglyF++;
									icsosF++;
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}

				
				if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
					bool ivariant = 0; // 0 - ������� �������(�����), 1 - �������������� ������.
					bool baddC = true;
					doublerealT Cmax = -1.0; // �������� ������� � ����� �� ������ � ������.
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
					if (!baddC)
					{
						maxelem_threshold_theta = 0.98*Cmax;
						if (ivariant == 1) {
							maxelem_threshold_theta = -1.0;
						}
						for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
							if (this_is_C_node[Amat.j[is0]]) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									doublerealT aij_abs = Amat.abs_aij[is0];
									if ((ivariant == 1) || (aij_abs > 0.98*Cmax))
									{
										// ���� (ivariant==1) �� ������� ��������� � �� ��������� ��� � ����.

										sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly C �����.
										icsos++;
										this_is_Strong_C_node[Amat.j[is0]] = true;
										iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
										isize_iStrongC_nodes_ind++;
									}
								}
							}
						}
					}

				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {

							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];

							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.

								if (this_is_C_node[Amat.j[is0]]  ) {
									
									//if (aij_abs > maxelem_threshold*theta) {
									if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
										/*
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else*/
										{
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (aij_abs > maxelem_threshold*barjer) {
											//if (aij_abs > maxelem_threshold*theta) {
											if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = aij_abs / sumP;
												if (SumPall < RealZERO) {
													printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												P[icount1].aij = aij_abs / SumPall;
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
									if (this_is_F_node[Amat.j[is0]]  ) {
										//if (aij_abs > maxelem_threshold*theta_strong_F) {
										if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta_strong_F)) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (aij_abs > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (aij_abs > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											doublerealT multiplyer_nu = aij_abs / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
												doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
												iend_for = is0_loc;
												if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = aij_loc_abs;
														}
													}
												}
											}

											//doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											bool bfound_amg1r5 = false;
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
												//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
												doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
												// ����� �������� ����� ���������� ����� ���������.
												if (aij_loc_abs > maxelem_threshold_loc_magic) {

													if (this_is_C_node[Amat.j[is0_loc]]  ) {

														if (Amat.j[is0_loc] != istr_etalon) {
															if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
																bfound_amg1r5 = true;
															}
														}
													}
												}
											}


											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											if (bfound_amg1r5) {
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
													// ����� �������� ����� ���������� ����� ���������.
													if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {

															if (Amat.j[is0_loc] != istr_etalon) {
																if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
																	sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																	icsos_loc++;
																}
															}
														}
														else {

															//if (Amat.j[is0_loc] != istr_etalon) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}
														}

													}
												}
											}
											else {
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

													// ����� �������� ����� ���������� ����� ���������.
													if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {

															if (Amat.j[is0_loc] != istr_etalon) {

																sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																icsos_loc++;
																break; // ���������� ������ ��������� ����������.
															}

														}
														else {

															//if (Amat.j[is0_loc] != istr_etalon) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}
														}

													}
												}
											}
											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											
											if (bfound_amg1r5) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.

													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

													// ��������������� ������� �� ��������� � �����.

													// ������ ������� ������ ����������� ��������� 
													// ������������� ��������� �������� �� �������� 
													// �������� �������.
													// ����������� 5 ������� 2015.
													if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {
															if (Amat.j[is0_loc] != istr_etalon) {
																if (this_is_Strong_C_node[Amat.j[is0_loc]]) {

																	P[icount1].j = i8;
																	P[icount1].i = C_numerate[Amat.j[is0_loc]];
																	//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																	if (sumP_loc < RealZERO) {
																		printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//getchar();
																		system("PAUSE");
																		exit(1);
																	}
																	P[icount1].aij = multiplyer_nu * aij_loc_abs / sumP_loc;
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
											else {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.

													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

													// ��������������� ������� �� ��������� � �����.

													// ������ ������� ������ ����������� ��������� 
													// ������������� ��������� �������� �� �������� 
													// �������� �������.
													// ����������� 5 ������� 2015.
													if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]  ) {
															if (Amat.j[is0_loc] != istr_etalon) {

																P[icount1].j = i8;
																P[icount1].i = C_numerate[Amat.j[is0_loc]];
																//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																if (sumP_loc < RealZERO) {
																	printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																	//getchar();
																	system("PAUSE");
																	exit(1);
																}
																P[icount1].aij = multiplyer_nu * aij_loc_abs / sumP_loc;
																icount1++;
																if (icount1 >= nsizePR * n) {
																	printf("memory error!!!\n");
																	printf("not enough memory for the interpolation operator.\n");
																	//system("PAUSE");
																	//exit(1);
																	deallocate_prolongation(nsizePR, n, P);
																}
																break; // ���������� ������ ��������� ����������.

															}
														}
													}
												}
											}


											//end
										}
									}
							}

						}
					}

				}

				// ����� - ������� ���-�������.
				for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
					this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
				}

			} // end only negative connections

		}

	}

	delete[] iStrongC_nodes_ind;
	delete[] this_is_Strong_C_node;

	//system("pause");


} // my_interpolation_procedure_number3A_PMIS


  // ���������������� ��������� �3 ��� PMIS.
  // ���������������� ������� ������������� �� ������� ����.
  // �� ���� ������������ �� ������� ���� �������� ������� ����.
  // 12.05.2018
template <typename doublerealT>
void my_interpolation_procedure_number3A_PMIS_parallel8(integer& the_number_of_neighbors_that_are_not_C_nodes_gl,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor_gl,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic_gl,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF_gl,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative)
{
	// ���������� ���������� ����� � ����� 01.10.2020
	// true - ����������� ����� ��������, �� ��������� ����������� ����������� ������.
	bool weakening_the_addition_of_new_C_nodes = true;// true;

	bool bprint_mesage_diagnostic = true;
	if (my_amg_manager.iprint_log == 0) {
		bprint_mesage_diagnostic = false;
	}

	if (!b_on_adaptive_local_refinement_mesh) {
		// �� ����������������� ����� ���������� ������ ������ � ����.
		// ������� ��������� ������ �����������.
		// �� ���� ����� ������ ����������� �������� �� 10% ����������� ��������� ��� ������ 
		// �� ������� �������. ������� ����������� �������� � ���������� ����� �������� �� ����������.
		weakening_the_addition_of_new_C_nodes = false;// ����������� ����� ��������� ���� �����.
													  //if (ilevel < 4) weakening_the_addition_of_new_C_nodes = false;// ������ 4, ��������� ������ ����.
	}

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(8); // ���������� 8 �������, 10 ������� ��� �������� �� �������.
#endif

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = activation_function_for_thrteshold(theta83, ilevel, false); // 3 ���� 2016
	doublerealT theta = activation_function_for_thrteshold(theta_gl, ilevel, false);

	// ���-�������
	//bool* this_is_Strong_C_node = new bool[n_a[ilevel - 1] + 1];

//#pragma omp parallel for
	//for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
		// ������������� ���-�������.
		//this_is_Strong_C_node[i_1] = false;
	//}

	integer ioneStrongC_and_0_StrongF = 0;

	// �������� ������ ������������.

	// ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
	// ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
	// ������� F �����.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	// the begining
	// ��������� ������������ ����� ������� � ������ �� ���� �������. 
	integer imax_count_sosed = -1;

//#pragma omp parallel for firstprivate(imax_count_sosed) lastprivate(imax_count_sosed)
	// �� � ���� ������ �� ���������������� � ������ ����� �.�. ���������� ������.
	// �������� �������� �������� ���������� imax_count_sosed � ���������� ������ ����������
	// ������������ � ������ ������� �� ������� ������� iStrongC_nodes_ind ��� ��������
	// � ������ ����������� ����� � �����. 06.10.2020.
#pragma omp parallel
	{
		integer imax_count_sosed_loc = -1;

#pragma omp for
		for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) {
			if (this_is_F_node[inode]) {
				integer ii1 = row_startA[inode];
				//integer iend_1 = row_startA[Amat.i[ii1] + 1] - 1;
				integer iend_1 = row_startA[inode + 1] - 1;

				if (iend_1 - ii1 > imax_count_sosed_loc) {
					imax_count_sosed_loc = iend_1 - ii1;
				}
			}
		}

#pragma omp critical 
		{
			if (imax_count_sosed_loc > imax_count_sosed) {
				imax_count_sosed = imax_count_sosed_loc;
			}
		}
	}
	//integer* iStrongC_nodes_ind = new integer[imax_count_sosed + 1];

	integer iadditionalCstatistic = 0;
	integer number_of_F_nodes_with_one_single_strong_C_neighborF = 0;

	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if (1) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		integer the_number_of_neighbors_that_are_not_C_nodes = 0;
		integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;
		//integer icount_bad_string = 0;


		if (bpositive_connections) {
			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.

#pragma omp parallel for reduction(+: number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes, iadditionalCstatistic, number_of_F_nodes_with_one_single_strong_C_neighborF) 
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) {
				if (this_is_F_node[i8]) {

					bool badd_amg1r5 = false;
					bool bfound_vneDiag_F_node = false;

					// ����� ������� ������� F-node ������� C-node.
					integer icsos = 0;
					integer icsosF = 0;

					// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
					// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;

					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
					integer iend_1 = iend_marker_position;

					//integer isize_iStrongC_nodes_ind = 0;





#if doubleintprecision == 1
					//printf("i8=%lld n=%lld\n", i8, n_a[ilevel - 1]);
#else
					//printf("i8=%d n=%d\n", i8, n_a[ilevel - 1]);
#endif

				//getchar();


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;

					//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								//doublerealT aij_abs = fabs(Amat.aij[is0]);
								doublerealT aij_abs = Amat.abs_aij[is0];
								// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
								//if (this_is_C_node[Amat.j[is0]]  ) {
								if (aij_abs > maxelem_threshold) {
									maxelem_threshold = aij_abs;
								}
								//}
							}
						}
					}
					else {
						maxelem_threshold = threshold_quick_all[Amat.i[ii1]];
					}
					// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


					//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (aij_abs > maxelem_threshold*barjer) {
								//if (aij_abs > maxelem_threshold*theta) {
								if (aij_abs > maxelem_threshold_theta) {
									sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.

													 // ������� ���� ������� � �������.
													 // ������������� ���-�������.
									//this_is_Strong_C_node[Amat.j[is0]] = true;

									//iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
									//isize_iStrongC_nodes_ind++;

									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (aij_abs > maxelem_threshold*theta_strong_F) {
									if (aij_abs > maxelem_threshold_theta_strong_F) {
										SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}
						}
					}

					// ����� � ������������ amg1r5 ���� ��������� ���������� ����. �� �� ����
					// ������������ ������ ��� ������� � ������������ ����������� 
					// ��� FF1 ������������.
					// ���� � ���� ������ ��������� FF ���� �� ���� ������ � ������ �� ����� ����������� ���������� ������ ����.
					// � PMIS � ��� ����� ���� ��� ������ ��������� FF ���� �� ������� ������ � ������.
					// � ���� ������ ������ ������������ ����������� ��
					// distance 2.

					// ����� - ������� ���-�������.
					//for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						//this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
					//}

					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}



					/*
					9.05.2018
					1. ���������� threshold�.
					2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
					3. icsos - ���������� strong C. icsosF - ���������� strong F.
					*/

					// 1 ������ 2016 ���� ����� ��� ������������.
					// ������� � ������ ������ ������ ������ ����� ���� � �����.



					bool baddC = true;
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (this_is_C_node[Amat.j[is0]]) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (weakening_the_addition_of_new_C_nodes) {
									// ���� ��� �� ������ �������� � ������
									// �� ��� ������� ���������� �� �����������
									// ����� �������� ����� F ���� ���� ������ � �����.
									// ����� � ������ ������� � F ���� ��� ������� � �������
									// ������ ������ ������� � ������ (������ �������� �� ������)
									// �� �� ��������������� �� �����
									// ������� � ������ � �� ��������� ����� � �����.
									// ��� �������� �������� ���������� � ��� ����.

									baddC = false;
								}
							}
						}
					}

					if (baddC) {

						if (sumP < RealZERO)
						{
							// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
							//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
							//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							//if (this_is_F_node[Amat.j[is0]]  ) {
							//if (Amat.j[is0] != Amat.i[ii1]) {

							if (iend_marker_position > ii1) {
								// ���� ���� ���� ���� ����� F ����.


								// 20 jan 2016.
								// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.

								//if ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO)) 

								// ��� ������ ����� ������ ��� ������� � �������.

								// if (icsosF==0) ������ ������ ����������.
								//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
								{
									//printf("error interpolation zero diagonal sumP.\n");
									//printf("Fnode all sosed is F");
									//system("pause");
									//printf("i8 is Dirichlet node\n");
									if (this_is_C_node[i8] == false) iadditionalCstatistic++;
									this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
									this_is_C_node[i8] = true;
									
									
#pragma omp critical
									{// ����� ������ ������ ���, ��� ��� ������ Acomp
										if (!byes_add) {
											byes_add = true;
										}
									}
									

													 //continue;

													 // ���� ����������� ����������� ������������� ������ ��������:
													 //printf("error!!! ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO))\n");
													 //system("pause");
													 //exit(1);
													 // ����� ����� �������� ������������.
								}


								//}
								//}
								//}//end
							}
						}
					}

				}
			}
		}
		else {
			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
#pragma omp parallel for reduction(+:  number_of_F_nodes_with_one_single_strong_C_neighbor, the_number_of_neighbors_that_are_not_C_nodes, iadditionalCstatistic, number_of_F_nodes_with_one_single_strong_C_neighborF) 
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) {
				if (this_is_F_node[i8]) {

					bool badd_amg1r5 = false;
					bool bfound_vneDiag_F_node = false;

					// ����� ������� ������� F-node ������� C-node.
					integer icsos = 0;
					integer icsosF = 0;

					// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
					// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;

					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
					integer iend_1 = iend_marker_position;

					//integer isize_iStrongC_nodes_ind = 0;



					// only negative connections


					// ��� ����� ����������� �������� ��������.
					// 5 ������� 2015 ���� �� ��������� ��������� �������������
					// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;

					//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								//doublerealT aij_abs = fabs(Amat.aij[is0]);
								doublerealT aij_abs = Amat.abs_aij[is0];
								// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
								//if (this_is_C_node[Amat.j[is0]]  ) {
								if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold)) {
									maxelem_threshold = aij_abs;
								}
								//}
							}
						}
					}
					else {
						maxelem_threshold = threshold_quick_only_negative[Amat.i[ii1]];
					}
					// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


					//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (aij_abs > maxelem_threshold*barjer) {
								//if (aij_abs > maxelem_threshold*theta) {
								if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
									sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
									icsos++;
									//this_is_Strong_C_node[Amat.j[is0]] = true;
									//iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
									//isize_iStrongC_nodes_ind++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (aij_abs > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta_strong_F)) {
										SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}
						}
					}

					// ����� - ������� ���-�������.
					//for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						//this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
					//}

					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}



					/*
					9.05.2018
					1. ���������� threshold�.
					2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
					3. icsos - ���������� strong C. icsosF - ���������� strong F.
					*/


					// 1 ������ 2016 ���� ����� ��� ������������.
					// ������� � ������ ������ ������ ������ ����� ���� � �����.

					bool baddC = true;
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (this_is_C_node[Amat.j[is0]]) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (weakening_the_addition_of_new_C_nodes) {
									// ���� ��� �� ������ �������� � ������
									// �� ��� ������� ���������� �� �����������
									// ����� �������� ����� F ���� ���� ������ � �����.
									// ����� � ������ ������� � F ���� ��� ������� � �������
									// ������ ������ ������� � ������ (������ �������� �� ������)
									// �� �� ��������������� �� �����
									// ������� � ������ � �� ��������� ����� � �����.
									// ��� �������� �������� ���������� � ��� ����.

									baddC = false;
								}
							}
						}
					}

					if (baddC) {
						// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
						if (sumP < RealZERO)
						{
							if (iend_marker_position > ii1) {
								// ���� ���� ���� ���� ����� F ����.

								//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
								//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
								//if (this_is_F_node[Amat.j[is0]]  ) {
								//if (Amat.j[is0] != Amat.i[ii1]) {


								// 20 jan 2016.
								// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.

								//if ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO)) 

								// ��� ������ ����� ������ ��� ������� � �������.

								// if (icsosF==0) ������ ������ ����������.
								//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
								{
									//printf("error interpolation zero diagonal sumP.\n");
									//printf("Fnode all sosed is F");
									//system("pause");
									//printf("i8 is Dirichlet node\n");
									// ��� F node �� ��������� ��� ���� �� ���� ��� ������ i8.

									if (this_is_C_node[i8] == false) iadditionalCstatistic++;
									this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
									this_is_C_node[i8] = true;
									
#pragma omp critical
									{// ����� ������ ������ ���, ��� ��� ������ Acomp
										if (!byes_add) {
											byes_add = true;
										}
									}
									
													 //continue;

													 // ���� ����������� ����������� ������������� ������ ��������:
													 //printf("error!!! ((fabs(sumP) < RealZERO)&&(fabs(SumPall) < RealZERO))\n");
													 //system("pause");
													 //exit(1);
													 // ����� ����� �������� ������������.
								}



								//}
								//}
								//}//end
							}
						}
					}
				}
			}
		}

		number_of_F_nodes_with_one_single_strong_C_neighbor_gl=number_of_F_nodes_with_one_single_strong_C_neighbor;
		the_number_of_neighbors_that_are_not_C_nodes_gl = the_number_of_neighbors_that_are_not_C_nodes;
	}

	//if (iadditionalCstatistic > 0) byes_add = true;
	bweSholdbeContinue = bweSholdbeContinue || byes_add;

	iadditionalCstatistic_gl += iadditionalCstatistic;
	iadditionalCstatistic = 0;

	number_of_F_nodes_with_one_single_strong_C_neighborF_gl += number_of_F_nodes_with_one_single_strong_C_neighborF;
	number_of_F_nodes_with_one_single_strong_C_neighborF = 0;

	if (byes_add) {
		if (bprint_mesage_diagnostic) {
			// ���� ����������� ����������� ������������� ������ ��������:
			printf("number (fabs(sumP) < RealZERO) = %lld; n_a=%lld.\n", iadditionalCstatistic_gl, n_a[ilevel - 1]);
		}
		//system("pause");
		//exit(1);
	}


	/*
	9.05.2018
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic_gl++;
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

		integer LOC_P_LIM_SIZE = -1;
		double dsize = 2.2*0.125*nsizePR*n_a[0];
		LOC_P_LIM_SIZE = (integer)(dsize);

		Ak1** P_loc = new Ak1*[8];
		integer* isize_P_loc = new integer[8];
		for (int itid = 0; itid < 8; itid++) {
			isize_P_loc[itid] = 0;
			
			P_loc[itid] = new Ak1[LOC_P_LIM_SIZE];
		}


		// ���-�������
		bool** this_is_Strong_C_node = new bool*[8];
		for (int itid = 0; itid < 8; itid++) {
			this_is_Strong_C_node[itid] = new bool[n_a[ilevel - 1] + 1];
		}

		for (int itid = 0; itid < 8; itid++) {

#pragma omp parallel for
			for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
				// ������������� ���-�������.
				this_is_Strong_C_node[itid][i_1] = false;
			}
		}

		// ������ ������� ������� ����������������.
		integer** iStrongC_nodes_ind = new integer*[8];
		for (int itid = 0; itid < 8; itid++) {
			iStrongC_nodes_ind[itid] = new integer[imax_count_sosed + 1];
		}


		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
		// �������� ������� �� ���� ���������
		// �� �� �������������� �� �� ����� V ������.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = activation_function_for_thrteshold(magic82, ilevel, true); // 0.4 is recomended.



		integer the_number_of_neighbors_that_are_not_C_nodes = 0;
		integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
#pragma omp parallel for reduction(+ : number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes,  number_of_F_nodes_with_one_single_strong_C_neighborF)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) {
#ifdef _OPENMP 
				int itid = omp_get_thread_num();
#else
				int itid = 0;
#endif

				if (this_is_F_node[i8]) {

					

					integer isize_iStrongC_nodes_ind = 0;

					// ��� ����� ����������� �������� ��������.
					// 5 ������� 2015 ���� �� ��������� ��������� �������������
					// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;
					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer istr_etalon1 = Amat.i[ii1];
					integer iend_for1 = -1;
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
							iend_for1 = is0;
							if (Amat.j[is0] != istr_etalon1) {
								//doublerealT aij_abs = fabs(Amat.aij[is0]);
								doublerealT aij_abs = Amat.abs_aij[is0];
								// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
								//if (this_is_C_node[Amat.j[is0]]  ) {
								if (aij_abs > maxelem_threshold) {
									maxelem_threshold = aij_abs;
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
					// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

					// ����� ������� ������� F-node ������� C-node.
					integer icsos = 0;
					integer icsosF = 0;

					//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

					// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
					// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
					for (integer is0 = ii1; is0 <= iend_for1; is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (aij_abs > maxelem_threshold*barjer) {
								//if (aij_abs > maxelem_threshold*theta) {
								if (aij_abs > maxelem_threshold_theta) {
									sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
									icsos++;

									this_is_Strong_C_node[itid][Amat.j[is0]] = true;
									iStrongC_nodes_ind[itid][isize_iStrongC_nodes_ind] = is0;
									isize_iStrongC_nodes_ind++;

								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (aij_abs > maxelem_threshold*theta_strong_F) {
									if (aij_abs > maxelem_threshold_theta_strong_F) {
										SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}
						}
					}
					if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
						bool ivariant = 0; // 0 - ������� �������(�����), 1 - �������������� ������.
						bool baddC = true;
						doublerealT Cmax = -1.0; // �������� ������� � ����� �� ������ � ������.
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
						if (!baddC)
						{
							maxelem_threshold_theta = 0.98*Cmax;
							if (ivariant == 1) {
								maxelem_threshold_theta = -1.0;
							}
							for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
								if (this_is_C_node[Amat.j[is0]]) {
									if (Amat.j[is0] != Amat.i[ii1]) {
										doublerealT aij_abs = Amat.abs_aij[is0];
										if ((ivariant == 1) || (aij_abs > 0.98*Cmax))
										{
											// ���� (ivariant==1) �� ������� ��������� � �� ��������� ��� � ����.

											sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly C �����.
											icsos++;
											this_is_Strong_C_node[itid][Amat.j[is0]] = true;
											iStrongC_nodes_ind[itid][isize_iStrongC_nodes_ind] = is0;
											isize_iStrongC_nodes_ind++;
										}
									}
								}
							}
						}

					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}



					{

						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						if (1) {
							//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
							// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
							//
							SumPall += sumP;

							for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// ��� ���������� ������ Strong �����.
									//doublerealT aij_abs = fabs(Amat.aij[is0]);
									doublerealT aij_abs = Amat.abs_aij[is0];

									if (this_is_C_node[Amat.j[is0]]) {
										//if (aij_abs > maxelem_threshold*theta) {
										if (aij_abs > maxelem_threshold_theta) {
											/*
											// �� ����� ����� ������ � ����. � ����������� ���� ����� ����� ������� F ���� �� � ����.
											if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
											}
											else */
											{
												// ��� ��� ��� ������������ Strong C �����. 
												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												//if (aij_abs > maxelem_threshold*barjer) {
												//if (aij_abs > maxelem_threshold*theta) {
												//if (aij_abs > maxelem_threshold_theta) {
													Ak1 Ptmp;
													Ptmp.j = i8;
													Ptmp.i = C_numerate[Amat.j[is0]];
													//Ptmp.aij = aij_abs / sumP;
													if (SumPall < RealZERO) {
														printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
														//getchar();
														system("PAUSE");
														exit(1);
													}
													Ptmp.aij = aij_abs / SumPall;

													if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

														P_loc[itid][isize_P_loc[itid]] = Ptmp;
														isize_P_loc[itid]++;
													}
													else {
														std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
														system("PAUSE");
													}

													//icount1++;

												//}

											}
										}

									}
									else
										if (this_is_F_node[Amat.j[is0]]) {
											//if (aij_abs > maxelem_threshold*theta_strong_F) {
											if (aij_abs > maxelem_threshold_theta_strong_F) {
												// ������������� Strong F �����.

												// �����:
												// 



												//if (aij_abs > maxelem_threshold*barjer) {
												// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
												// � ������� F ������.
												//if (aij_abs > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat.j[is0];
												if (SumPall < RealZERO) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												doublerealT multiplyer_nu = aij_abs / SumPall;
												// ��������� ���� ������� iFpointeger 
												// ����� ����� ����� ��� � ����.

												// �������������� ��������� �����.
												doublerealT maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												integer ii1_loc = row_startA[iFpoint];
												integer istr_etalon = Amat.i[ii1_loc];
												integer iend_for = -1;
												integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
													iend_for = is0_loc;
													if (Amat.j[is0_loc] != istr_etalon) {
														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
														if (aij_loc_abs > maxelem_threshold_loc) {
															// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
															// ����� ������ ����������.
															if (this_is_C_node[Amat.j[is0_loc]])
															{
																if (Amat.j[is0_loc] != istr_etalon) {
																	maxelem_threshold_loc = aij_loc_abs;
																}
															}
														}
													}
												}

												//doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

												// ����� ������� ������� F-node ������� C-node.
												integer icsos_loc = 0;

												bool bfound_amg1r5 = false;
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ����� �������� ����� ���������� ����� ���������.
													if (Amat.abs_aij[is0_loc] > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]) {

															if (Amat.j[is0_loc] != istr_etalon) {
																if (this_is_Strong_C_node[itid][Amat.j[is0_loc]]) {
																	bfound_amg1r5 = true;
																}
															}
														}
													}
												}

												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
												// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
												doublerealT sumP_loc = 0.0;

												if (bfound_amg1r5) {
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
														// ����� �������� ����� ���������� ����� ���������.
														if (aij_loc_abs > maxelem_threshold_loc_magic) {

															if (this_is_C_node[Amat.j[is0_loc]]) {

																if (Amat.j[is0_loc] != istr_etalon) {
																	if (this_is_Strong_C_node[itid][Amat.j[is0_loc]]) {
																		sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																		icsos_loc++;
																	}
																}

															}
															else {

																//if (Amat.j[is0_loc] != istr_etalon) {
																// ������������ ���������� ������� ������� �� �������� � ������.
																//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
																//}
															}

														}
													}
												}
												else {

													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
														// ����� �������� ����� ���������� ����� ���������.
														if (aij_loc_abs > maxelem_threshold_loc_magic) {

															if (this_is_C_node[Amat.j[is0_loc]]) {

																if (Amat.j[is0_loc] != istr_etalon) {

																	sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																	icsos_loc++;
																	break; // ���������� ������ ��������� ����������.
																}

															}
															else {

																//if (Amat.j[is0_loc] != istr_etalon) {
																// ������������ ���������� ������� ������� �� �������� � ������.
																//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
																//}
															}

														}
													}
												}

												doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												// � ����� ��� ������� ���������������� ����� 
												if (bfound_amg1r5) {
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.
														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

														// ��������������� ������� �� ��������� � �����.

														// ������ ������� ������ ����������� ��������� 
														// ������������� ��������� �������� �� �������� 
														// �������� �������.
														// ����������� 5 ������� 2015.
														if (aij_loc_abs > maxelem_threshold_loc_magic) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {

																	if (this_is_Strong_C_node[itid][Amat.j[is0_loc]]) {

																		//if (sumP_loc > RealZERO) {
																			// WARNING ������������� ���� ������ ���� �� ����.

																			Ak1 Ptmp;
																			Ptmp.j = i8;
																			Ptmp.i = C_numerate[Amat.j[is0_loc]];
																			//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																			if (sumP_loc < RealZERO) {
																				printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																				//getchar();
																				system("PAUSE");
																				exit(1);
																			}
																			Ptmp.aij = multiplyer_nu * aij_loc_abs / sumP_loc;

																			if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																				P_loc[itid][isize_P_loc[itid]] = Ptmp;
																				isize_P_loc[itid]++;
																			}
																			else {
																				std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																				system("PAUSE");
																			}
																		//}
																	}
																}
															}
														}
													}
												}
												else {
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.
														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

														// ��������������� ������� �� ��������� � �����.

														// ������ ������� ������ ����������� ��������� 
														// ������������� ��������� �������� �� �������� 
														// �������� �������.
														// ����������� 5 ������� 2015.
														if (aij_loc_abs > maxelem_threshold_loc_magic) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {

																	Ak1 Ptmp;

																	Ptmp.j = i8;
																	Ptmp.i = C_numerate[Amat.j[is0_loc]];
																	//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																	if (sumP_loc < RealZERO) {
																		printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//getchar();
																		system("PAUSE");
																		exit(1);
																	}
																	Ptmp.aij = multiplyer_nu * aij_loc_abs / sumP_loc;

																	if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																		P_loc[itid][isize_P_loc[itid]] = Ptmp;
																		isize_P_loc[itid]++;
																	}
																	else {
																		std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																		system("PAUSE");
																	}

																	break;// �� ����� ������ ��������� �����������.

																}
															}
														}
													}
													//end
												}

											}
										}


								}

							}
						}

					}

					// ����� - ������� ���-�������.
					for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						this_is_Strong_C_node[itid][Amat.j[iStrongC_nodes_ind[itid][is00]]] = false;
					}

				}
			}
			////

			
			

		}
		else {
			// only negative connections.

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
#pragma omp parallel for reduction(+ : number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes,  number_of_F_nodes_with_one_single_strong_C_neighborF)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) {

#ifdef _OPENMP 
				int itid = omp_get_thread_num();
#else
				int itid = 0;
#endif

				if (this_is_F_node[i8]) {

					integer isize_iStrongC_nodes_ind = 0;

					// ��� ����� ����������� �������� ��������.
					// 5 ������� 2015 ���� �� ��������� ��������� �������������
					// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;
					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer istr_etalon1 = Amat.i[ii1];
					integer iend_for1 = -1;
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
							iend_for1 = is0;
							if (Amat.j[is0] != istr_etalon1) {
								//doublerealT aij_abs = fabs(Amat.aij[is0]);
								doublerealT aij_abs = Amat.abs_aij[is0];
								// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
								//if (this_is_C_node[Amat.j[is0]]  ) {
								if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold)) {
									maxelem_threshold = aij_abs;
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
					// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

					// ����� ������� ������� F-node ������� C-node.
					integer icsos = 0;
					integer icsosF = 0;

					//doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					//doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
					doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
					doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

					// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
					// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
					doublerealT sumP = 0.0;
					doublerealT SumPall = 0.0;
					integer icount_StronglyF = 0;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
					for (integer is0 = ii1; is0 <= iend_for1; is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							//doublerealT aij_abs = fabs(Amat.aij[is0]);
							doublerealT aij_abs = Amat.abs_aij[is0];
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (aij_abs > maxelem_threshold*barjer) {
								//if (aij_abs > maxelem_threshold*theta) {
								if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
									sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
									icsos++;

									this_is_Strong_C_node[itid][Amat.j[is0]] = true;
									iStrongC_nodes_ind[itid][isize_iStrongC_nodes_ind] = is0;
									isize_iStrongC_nodes_ind++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (aij_abs > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta_strong_F)) {
										SumPall += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}
						}
					}


					if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
						bool ivariant = 0; // 0 - ������� �������(�����), 1 - �������������� ������.
						bool baddC = true;
						doublerealT Cmax = -1.0; // �������� ������� � ����� �� ������ � ������.
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
						if (!baddC)
						{
							maxelem_threshold_theta = 0.98*Cmax;
							if (ivariant == 1) {
								maxelem_threshold_theta = -1.0;
							}
							for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
								if (this_is_C_node[Amat.j[is0]]) {
									if (Amat.j[is0] != Amat.i[ii1]) {
										doublerealT aij_abs = Amat.abs_aij[is0];
										if ((ivariant == 1) || (aij_abs > 0.98*Cmax))
										{
											// ���� (ivariant==1) �� ������� ��������� � �� ��������� ��� � ����.

											sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly C �����.
											icsos++;
											this_is_Strong_C_node[itid][Amat.j[is0]] = true;
											iStrongC_nodes_ind[itid][isize_iStrongC_nodes_ind] = is0;
											isize_iStrongC_nodes_ind++;
										}
									}
								}
							}
						}

					}
					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}



					{

						if (1) {
							//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
							// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
							//
							SumPall += sumP;

							for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {

								//doublerealT aij_abs = fabs(Amat.aij[is0]);
								doublerealT aij_abs = Amat.abs_aij[is0];

								if (Amat.j[is0] != Amat.i[ii1]) {
									// ��� ���������� ������ Strong �����.

									if (this_is_C_node[Amat.j[is0]]) {

										//if (aij_abs > maxelem_threshold*theta) {
										if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {
											/*
											if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
											}
											else*/
											{
												// ��� ��� ��� ������������ Strong C �����. 
												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												//if (aij_abs > maxelem_threshold*barjer) {
												//if (aij_abs > maxelem_threshold*theta) {
												if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta)) {

													Ak1 Ptmp;

													Ptmp.j = i8;
													Ptmp.i = C_numerate[Amat.j[is0]];
													//Ptmp.aij = aij_abs / sumP;
													if (SumPall < RealZERO) {
														printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
														//getchar();
														system("PAUSE");
														exit(1);
													}
													Ptmp.aij = aij_abs / SumPall;

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
										if (this_is_F_node[Amat.j[is0]]) {
											//if (aij_abs > maxelem_threshold*theta_strong_F) {
											if ((Amat.aij[is0] < 0.0) && (aij_abs > maxelem_threshold_theta_strong_F)) {
												// ������������� Strong F �����.

												// �����:
												// 



												//if (aij_abs > maxelem_threshold*barjer) {
												// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
												// � ������� F ������.
												//if (aij_abs > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat.j[is0];
												doublerealT multiplyer_nu = aij_abs / SumPall;
												// ��������� ���� ������� iFpointeger 
												// ����� ����� ����� ��� � ����.

												// �������������� ��������� �����.
												doublerealT maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												integer ii1_loc = row_startA[iFpoint];
												integer istr_etalon = Amat.i[ii1_loc];
												integer iend_for = -1;
												integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
													iend_for = is0_loc;
													if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc)) {
														if (this_is_C_node[Amat.j[is0_loc]]) {
															if (Amat.j[is0_loc] != istr_etalon) {
																maxelem_threshold_loc = aij_loc_abs;
															}
														}
													}
												}

												//doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

												// ����� ������� ������� F-node ������� C-node.
												integer icsos_loc = 0;

												bool bfound_amg1r5 = false;
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
													//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
													doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
													// ����� �������� ����� ���������� ����� ���������.
													if (aij_loc_abs > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]) {

															if (Amat.j[is0_loc] != istr_etalon) {
																if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
																	bfound_amg1r5 = true;
																}
															}
														}
													}
												}


												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
												// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
												doublerealT sumP_loc = 0.0;
												if (bfound_amg1r5) {
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {
														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];
														// ����� �������� ����� ���������� ����� ���������.
														if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

															if (this_is_C_node[Amat.j[is0_loc]]) {

																if (Amat.j[is0_loc] != istr_etalon) {
																	if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
																		sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																		icsos_loc++;
																	}
																}
															}
															else {

																//if (Amat.j[is0_loc] != istr_etalon) {
																// ������������ ���������� ������� ������� �� �������� � ������.
																//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
																//}
															}

														}
													}
												}
												else {
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

														// ����� �������� ����� ���������� ����� ���������.
														if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

															if (this_is_C_node[Amat.j[is0_loc]]) {

																if (Amat.j[is0_loc] != istr_etalon) {

																	sumP_loc += aij_loc_abs; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																	icsos_loc++;
																	break; // ���������� ������ ��������� ����������.
																}

															}
															else {

																//if (Amat.j[is0_loc] != istr_etalon) {
																// ������������ ���������� ������� ������� �� �������� � ������.
																//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
																//}
															}

														}
													}
												}
												doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												// � ����� ��� ������� ���������������� ����� 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {

												if (bfound_amg1r5) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.

														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

														// ��������������� ������� �� ��������� � �����.

														// ������ ������� ������ ����������� ��������� 
														// ������������� ��������� �������� �� �������� 
														// �������� �������.
														// ����������� 5 ������� 2015.
														if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {
																	if (this_is_Strong_C_node[Amat.j[is0_loc]]) {

																		Ak1 Ptmp;

																		Ptmp.j = i8;
																		Ptmp.i = C_numerate[Amat.j[is0_loc]];
																		//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																		if (sumP_loc < RealZERO) {
																			printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																			//getchar();
																			system("PAUSE");
																			exit(1);
																		}
																		Ptmp.aij = multiplyer_nu * aij_loc_abs / sumP_loc;
																		
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
												else {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

														// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.

														//doublerealT aij_loc_abs = fabs(Amat.aij[is0_loc]);
														doublerealT aij_loc_abs = Amat.abs_aij[is0_loc];

														// ��������������� ������� �� ��������� � �����.

														// ������ ������� ������ ����������� ��������� 
														// ������������� ��������� �������� �� �������� 
														// �������� �������.
														// ����������� 5 ������� 2015.
														if ((Amat.aij[is0_loc] < 0.0) && (aij_loc_abs > maxelem_threshold_loc_magic)) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {

																	Ak1 Ptmp;

																	Ptmp.j = i8;
																	Ptmp.i = C_numerate[Amat.j[is0_loc]];
																	//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																	if (sumP_loc < RealZERO) {
																		printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//getchar();
																		system("PAUSE");
																		exit(1);
																	}
																	Ptmp.aij = multiplyer_nu * aij_loc_abs / sumP_loc;


																	if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																		P_loc[itid][isize_P_loc[itid]] = Ptmp;
																		isize_P_loc[itid]++;
																	}
																	else {
																		std::cout << "Error Increase size P_loc memory limit\n" << std::endl;
																		system("PAUSE");
																	}
																	break; // ���������� ������ ��������� ����������.

																}
															}
														}
													}
												}


												//end
											}
										}
								}

							}
						}

					}

					// ����� - ������� ���-�������.
					for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						this_is_Strong_C_node[itid][Amat.j[iStrongC_nodes_ind[itid][is00]]] = false;
					}

				} // end only negative connections
			}
		}


		// �������� ����������� ����������.
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
			delete[] P_loc[itid];
			delete[] this_is_Strong_C_node[itid];
			delete[] iStrongC_nodes_ind[itid];
		}
		delete[] isize_P_loc;
		delete[] P_loc;
		delete[] this_is_Strong_C_node;
		delete[] iStrongC_nodes_ind;

		number_of_F_nodes_with_one_single_strong_C_neighbor_gl = number_of_F_nodes_with_one_single_strong_C_neighbor;
		the_number_of_neighbors_that_are_not_C_nodes_gl = the_number_of_neighbors_that_are_not_C_nodes;
	}

	number_of_F_nodes_with_one_single_strong_C_neighborF_gl += number_of_F_nodes_with_one_single_strong_C_neighborF;

	
	

	//system("pause");
#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

} // my_interpolation_procedure_number3A_PMIS_parallel8

  // ���������������� ��������� �3.
  // ���������������� ������� ������������� �� ������� ����.
  // �� ���� ������������ �� ������� ���� �������� ������� ����.
  // 12.05.2018
template <typename doublerealT>
void my_interpolation_procedure_number3B(integer &the_number_of_neighbors_that_are_not_C_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak1* &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublerealT RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublerealT &theta, integer &n, Ak1* &R, integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublerealT &magic82, doublerealT* &threshold_quick_all,
	doublerealT* &threshold_quick_only_negative)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = theta83; // 3 ���� 2016

	// ���-�������
	bool* this_is_Strong_C_node = new bool[n_a[ilevel-1]+1];
	for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
		// ������������� ���-�������.
		this_is_Strong_C_node[i_1] = false;
	}


	// �������� ������ ������������.
	// ������������ ���� � �������.
	integer ioneStrongC_and_0_StrongF = 0;

	// ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
	// ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
	// ������� F �����.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	// the begining

	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if (1) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


		// ���������� ����������� ��� ����� ������� ���������� F nodes.
		// ������ F-nodes ������ C-nodes.
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {


			bool badd_amg1r5 = false;

			// ����� ������� ������� F-node ������� C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
			// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
			doublerealT sumP = 0.0;
			doublerealT SumPall = 0.0;
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


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat[is0].j != Amat[ii1].i) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat[is0].j]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// ������������� ���-�������.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								// ������� ���� ������� � �������.
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}

				
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j]  ) {
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
											// ������� ������ ��� ������� F ������� ������� ����� � �������(����� �� ����������� �������) �������������� ������� � �������.

											SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
												icount_StronglyF++;
												icsosF++;
										}
										else {
											badd_amg1r5 = true;
										}
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}

				// ����� - ������� ���-�������.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat[is0].j] = false;
							}
						}
					}
				}

				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}
			}
			else {
				// only negative connections


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat[is0].j != Amat[ii1].i) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat[is0].j]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
				
				// ������������� ���-�������.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}

				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j]  ) {
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

										// ������� ������ ��� ������� F ������� ������� ����� � �������(����� �� ����������� �������) �������������� ������� � �������.
										SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
											icount_StronglyF++;
											icsosF++;
									}
									else {
										badd_amg1r5 = true;
									}
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}

				// ����� - ������� ���-�������.
				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							if ((Amat[is0].aij < 0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = false;
							}
						}
					}
				}

				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}

			}

			/*
			9.05.2018
			1. ���������� threshold�.
			2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
			3. icsos - ���������� strong C. icsosF - ���������� strong F.
			*/

			// 1 ������ 2016 ���� ����� ��� ������������.
			// ������� � ������ ������ ������ ������ ����� ���� � �����.


			// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
			//integer iend_marker_position = row_startA[Amat[ii1].i + 1] - 1;
			for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				if (this_is_F_node[Amat[is0].j]  ) {
					if (Amat[is0].j != Amat[ii1].i) {


						// 20 jan 2016.
						// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.


						if ((fabs(sumP) < RealZERO)||(badd_amg1r5)) {
							// ��� ������ ����� ������ ��� ������� � �������.

							// if (icsosF==0) ������ ������ ����������.
							//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
							{
								//printf("error interpolation zero diagonal sumP.\n");
								//printf("Fnode all sosed is F");
								//system("pause");
								//printf("i8 is Dirichlet node\n");
								if (this_is_C_node[i8] == false) iadditionalCstatistic++;
								this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								byes_add = true; // ���� ���������� �����.
												 //exit(1);
												 // ����� ����� �������� ������������.
							}
						}


					}
				}
			}//end





		}

	}





	/*
	9.05.2018
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
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

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//---->doublerealT magic = 0.4; // 0.4 optimum
								//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
								// �������� ������� �� ���� ���������
								// �� �� �������������� �� �� ����� V ������.
								//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82; // 0.4 is recomended.



		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat[is0].j]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}



				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j]  ) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
								sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j]  ) {
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


										SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
											icount_StronglyF++;
										icsosF++;
									}
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat[is0].j]  ) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if (fabs(Amat[is0].aij) > maxelem_threshold_theta) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat[is0].j]  ) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F) {
											// ������������� Strong F �����.

											// �����:
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
												// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
												// � ������� F ������.
												//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat[is0].j;
												if (fabs(SumPall) < RealZERO) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												doublerealT multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
												// ��������� ���� ������� iFpointeger 
												// ����� ����� ����� ��� � ����.

												// �������������� ��������� �����.
												doublerealT maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												integer ii1_loc = row_startA[iFpoint];
												integer istr_etalon = Amat[ii1_loc].i;
												integer iend_for = -1;
												integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
													iend_for = is0_loc;
													//if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc) {
														// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
														// ����� ������ ����������.
														if (this_is_C_node[Amat[is0_loc].j]  ) {
															if (this_is_Strong_C_node[Amat[is0_loc].j]) {
																if (Amat[is0_loc].j != istr_etalon) {
																	maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
																}
															}
														}
													//}
												}

												doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

												// ����� ������� ������� F-node ������� C-node.
												integer icsos_loc = 0;

												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
												// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
												doublerealT sumP_loc = 0.0;
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ����� �������� ����� ���������� ����� ���������.
													//if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {

													if (this_is_Strong_C_node[Amat[is0_loc].j]) {
														if (this_is_C_node[Amat[is0_loc].j]  ) {

															if (Amat[is0_loc].j != istr_etalon) {

																sumP_loc += fabs(Amat[is0_loc].aij); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																icsos_loc++;
															}

														}
														else {

															//if (Amat[is0_loc].j != istr_etalon) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}
														}
													}

													//}
												}

												doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												// � ����� ��� ������� ���������������� ����� 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

													// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


													// ��������������� ������� �� ��������� � �����.

													// ������ ������� ������ ����������� ��������� 
													// ������������� ��������� �������� �� �������� 
													// �������� �������.
													// ����������� 5 ������� 2015.
													//if (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic) {
													if (this_is_Strong_C_node[Amat[is0_loc].j]) {

														if (this_is_C_node[Amat[is0_loc].j]  ) {
															if (Amat[is0_loc].j != istr_etalon) {

																P[icount1].j = i8;
																P[icount1].i = C_numerate[Amat[is0_loc].j];
																//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
																if (fabs(sumP_loc) < RealZERO) {
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
						if (this_is_C_node[Amat[is0].j]  ) {
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

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat[ii1].i;
				integer iend_for1 = -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); is0++) {
						iend_for1 = is0;
						if (Amat[is0].j != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat[is0].j]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
					if (Amat[is0].j != Amat[ii1].i) {
						if (this_is_C_node[Amat[is0].j]  ) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat[is0].j] = true;
							}
						}
					}
				}


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat[is0].i == Amat[ii1].i); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat[is0].j != istr_etalon1) {
						if (this_is_C_node[Amat[is0].j]  ) {
							//	if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
							//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
							if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
								sumP += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat[is0].j]  ) {
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


										SumPall += fabs(Amat[is0].aij); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat[ii1].i + 1] - 1); is0++) {
							if (Amat[is0].j != Amat[ii1].i) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat[is0].j]  ) {
									//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
									if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("i8 is Dirichlet node\n");
											if (this_is_C_node[i8] == false) iadditionalCstatistic++;
											this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[i8] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta) {
											if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta)) {
												P[icount1].j = i8;
												P[icount1].i = C_numerate[Amat[is0].j];
												//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat[is0].j]  ) {
										//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {
										if ((Amat[is0].aij<0.0) && (fabs(Amat[is0].aij) > maxelem_threshold_theta_strong_F)) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (fabs(Amat[is0].aij) > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (fabs(Amat[is0].aij) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat[is0].j;
											doublerealT multiplyer_nu = fabs(Amat[is0].aij) / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat[ii1_loc].i;
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												//if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc)) {
													if (this_is_C_node[Amat[is0_loc].j]  ) {
														if (this_is_Strong_C_node[Amat[is0_loc].j]) {
															if (Amat[is0_loc].j != istr_etalon) {
																maxelem_threshold_loc = fabs(Amat[is0_loc].aij);
															}
														}
													}
												//}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												//if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat[is0_loc].j]) {


													if (this_is_C_node[Amat[is0_loc].j]  ) {

														if (Amat[is0_loc].j != istr_etalon) {

															sumP_loc += fabs(Amat[is0_loc].aij); // ����� ������� ��������������� ��������� ������� ����������� � �����.
															icsos_loc++;
														}

													}
													else {

														//if (Amat[is0_loc].j != istr_etalon) {
														// ������������ ���������� ������� ������� �� �������� � ������.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
														//}
													}
												}
												//}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat[is0_loc].i == Amat[ii1_loc].i); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												//if ((Amat[is0_loc].aij<0.0) && (fabs(Amat[is0_loc].aij) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat[is0_loc].j]) {


													if (this_is_C_node[Amat[is0_loc].j]  ) {
														if (Amat[is0_loc].j != istr_etalon) {

															P[icount1].j = i8;
															P[icount1].i = C_numerate[Amat[is0_loc].j];
															//P[icount1].aij = fabs(Amat[is0].aij) / sumP;
															if (fabs(sumP_loc) < RealZERO) {
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
						if (this_is_C_node[Amat[is0].j]  ) {
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

  // ���������������� ��������� �3.
  // ���������������� ������� ������������� �� ������� ����.
  // �� ���� ������������ �� ������� ���� �������� ������� ����.
  // 12.05.2018
// 17-18 ����� 2019 ���� ����������� �����������. �������������� ����� ���������� �������� �� 1.6%.
template <typename doublerealT>
void my_interpolation_procedure_number3B(integer &the_number_of_neighbors_that_are_not_C_nodes,
	integer &number_of_F_nodes_with_one_single_strong_C_neighbor,
	integer* &n_a, bool* &this_is_F_node, integer* &row_startA,
	integer* &nnz_a, bool &bpositive_connections, Ak2 &Amat,
	bool &bweSholdbeContinue, bool* &this_is_C_node, integer &iadditionalCstatistic,
	const doublerealT RealZERO, integer &icount1, Ak1* &P, integer &nsizePR, integer &ilevel,
	integer &iadd, doublerealT &theta, integer &n,  integer* &C_numerate,
	integer &number_of_F_nodes_with_one_single_strong_C_neighborF,
	doublerealT &theta83, bool &btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool &from_re_operation_protection0, doublerealT &magic82, doublerealT* &threshold_quick_all,
	doublerealT* &threshold_quick_only_negative, bool& bsuffix_work)
{

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = theta83; // 3 ���� 2016

	// ���-�������
	bool* this_is_Strong_C_node = new bool[n_a[ilevel - 1] + 1];
	for (integer i_1 = 0; i_1 <= n_a[ilevel - 1]; i_1++) {
		// ������������� ���-�������.
		this_is_Strong_C_node[i_1] = false;
	}


	// �������� ������ ������������.
	// ������������ ����� ���� � ������ ������� amg1r5.
	integer ioneStrongC_and_0_StrongF = 0;

	// ����� ���� F �� ������� Strong � ������� ��� ���������� � �����.
	// ���� F ������� ������ Strong  � ������ �������������� � ������� ������� � ������� 
	// ������� F �����.

	//6interpolation 0.4 6.77 11 26 28.355
	//6interpolation 0.45 6.6 10 27 28.151
	//6interpolation 0.5 6.42 12 32 28.735
	//4interpolation 0.4 3.7  52 24.736 // best
	//4interpolation 0.3 3.78 13 59 27.525
	//4interpolation 0.5 3.61 12 55 25.533
	//4interpolation 0.45 3.65 10 63 30.24

	// the begining

	integer imax_count_sosed = -1;
	for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode]  ) {
		integer ii1 = row_startA[inode];
		//integer iend_1 = row_startA[Amat.i[ii1] + 1] - 1;
		integer iend_1 = row_startA[inode + 1] - 1;
		if (iend_1 - ii1 > imax_count_sosed) {
			imax_count_sosed = iend_1 - ii1;
		}
	}
	integer *iStrongC_nodes_ind = new integer[imax_count_sosed + 1];

	bool byes_add = false;
	// ������� ���������� ����������� � �����.
	if ((bsuffix_work)&&(1)) {
		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//doublerealTT magic = 0.4; // 0.4 optimum


		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;
		integer icount_bad_string = 0;

		

		// ���������� ����������� ��� ����� ������� ���������� F nodes.
		// ������ F-nodes ������ C-nodes.
		for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode]  ) {


			bool badd_amg1r5 = false;
			bool bfound_vneDiag_F_node = false;

			// ����� ������� ������� F-node ������� C-node.
			integer icsos = 0;
			integer icsosF = 0;

			// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
			// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
			doublerealT sumP = 0.0;
			doublerealT SumPall = 0.0;
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


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != inode) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
							//doublerealT mcand = fabs(Amat.aij[is0]);
							doublerealT mcand = Amat.abs_aij[is0];
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				
				
				// ������ ���� ������ ���������� this_is_Strong_C_node
				
				for (integer is0 = ii1; (is0 <= iend_1); is0++) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
					// ������������ ���� �� ��������� F ����
						if (this_is_C_node[Amat.j[is0]]  ) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
							if (Amat.abs_aij[is0] > maxelem_threshold_theta) {

								// ������� ���� ������� � �������.
								// ������������� ���-�������.
								this_is_Strong_C_node[Amat.j[is0]] = true;
								sumP += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
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
								if (this_is_F_node[Amat.j[is0]]  )
								{

									if (Amat.j[is0] != inode) 
									{

										the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 

										bfound_vneDiag_F_node = true;
										// F-node
										// ������������ ���������� ������� ������� �� �������� � ������.
									
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
												// ������� ������ ��� ������� F ������� ������� ����� � �������(����� �� ����������� �������) �������������� ������� � �������.

												//SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
												SumPall += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
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




					// ����� - ������� ���-�������.
						/*//253 �����
					for (integer is0 = ii1; (is0 <= iend_1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (this_is_C_node[Amat.j[is0]]  ) {
								///if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
									this_is_Strong_C_node[Amat.j[is0]] = false;
								}
							}
						}
					}
					*/
					// ����� - ������� ���-�������.
					for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
					}
					

					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}
				}
			}
			else {
				// only negative connections


				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;

				//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {

				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
						if (Amat.j[is0] != inode) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				// ������������� ���-�������.
				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					//if (Amat.j[is0] != inode) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat.j[is0]] = true;
								//sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								sumP += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
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
						
						if (this_is_F_node[Amat.j[is0]]  ) {

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

										// ������� ������ ��� ������� F ������� ������� ����� � �������(����� �� ����������� �������) �������������� ������� � �������.
										//SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										SumPall += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
									else {
										badd_amg1r5 = true;
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}


						}
					}

					// ����� - ������� ���-�������.
					/*
					for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
						if (Amat.j[is0] != Amat.i[ii1]) {
							if (this_is_C_node[Amat.j[is0]]  ) {
								if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
									this_is_Strong_C_node[Amat.j[is0]] = false;
								}
							}
						}
					}
					*/
					// ����� - ������� ���-�������.
					for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
						this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
					}
					

					if (icsos == 1) {
						number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																			   // ��������� ������ ������ "����������".
																			   // ���������� ������ ����������� ��� ���������.
																			   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																			   // ������������ �� ���� ������� ����� ��������.
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
					}
				}

			}

			/*
			9.05.2018
			1. ���������� threshold�.
			2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
			3. icsos - ���������� strong C. icsosF - ���������� strong F.
			*/

			// 1 ������ 2016 ���� ����� ��� ������������.
			// ������� � ������ ������ ������ ������ ����� ���� � �����.


			// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
			
			// 20 jan 2016.
			// ����� �� ��������� ���� � ������ ���� ���� �� ����� ��� ������� F �������.
			/// ��� 
			//for (integer is0 = ii1; (is0 <= iend_marker_position); is0++) {
				//if (this_is_F_node[Amat.j[is0]]  ) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
						// ���� ��������������� F ����.
			
					//}
				//}
			//}//end
			// ������������ if (bfound_vneDiag_F_node&&((

			//236 ������
			if (bfound_vneDiag_F_node&&((fabs(sumP) < RealZERO) || (badd_amg1r5))) {
				// ��� ������ ����� ������ ��� ������� � �������.

				// if (icsosF==0) ������ ������ ����������.
				//if (icsosF<2)// ����������� ���������� �������� ���������� � ������� �������.
				{
					//printf("error interpolation zero diagonal sumP.\n");
					//printf("Fnode all sosed is F");
					//system("pause");
					//printf("inode is Dirichlet node\n");
					if (this_is_C_node[inode] == false) iadditionalCstatistic++;
						this_is_F_node[inode] = false; // ���� ���� ������� ������ � �����.
						this_is_C_node[inode] = true;
						bweSholdbeContinue = true;
						byes_add = true; // ���� ���������� �����.
						 //exit(1);
						 // ����� ����� �������� ������������.
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
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
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

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
								//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
								// �������� ������� �� ���� ���������
								// �� �� �������������� �� �� ����� V ������.
								//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82;// 0.4 is recomended.



		the_number_of_neighbors_that_are_not_C_nodes = 0;
		number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F-nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode]  ) {

				integer isize_iStrongC_nodes_ind = 0;

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, inode, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[inode];
				integer istr_etalon1 = inode;// Amat.i[ii1];
				//integer iend_for1 = -1;
				integer iend_1 = row_startA[istr_etalon1 + 1] - 1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_1); is0++) {
						//iend_for1 = is0;
						if (Amat.j[is0] != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;

				for (integer is0 = ii1; (is0 <= iend_1); is0++) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
					// �� ����� �������� ��� ������������ ���� � ��� F-����.
						if (this_is_C_node[Amat.j[is0]]  ) {
							if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat.j[is0]] = true;
								iStrongC_nodes_ind[isize_iStrongC_nodes_ind] = is0;
								isize_iStrongC_nodes_ind++;
								sumP += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
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
					this_is_F_node[inode] = false; // ���� ���� ������� ������ �-�����.
					this_is_C_node[inode] = true;
					bweSholdbeContinue = true;
					//exit(1);
					// ����� ����� �������� ������������.
				}
				*/

				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_1; is0++) {
					if (this_is_F_node[Amat.j[is0]]  ) {				

						//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
						if (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F) {
							// �� ����� �������� ��� ������������ ���� � ��� F node
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


										SumPall += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}

									// ������������ ���������� ������� ������� �� �������� � ������.
									the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
								}								
							}						
						
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ �������  � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						
						//for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
							//	integer is0 = iStrongC_nodes_ind[is00];						
						//}

						for (integer is0 = ii1; (is0 <= iend_1); is0++) {
							
								
								// ������������ ���� F ���� ������� �������� �� ����������������� �� ������������
								if (this_is_C_node[Amat.j[is0]]  ) {
									
									// ��� ���������� ������ Strong �����.
									if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
										
										// ��� ��� ��� ������������ Strong C �����. 
										// ��������������� ������� �� ��������� � �����.

										// ������ ������� ������ ����������� ��������� 
										// ������������� ��������� �������� �� �������� 
										// �������� �������.
										// ����������� 5 ������� 2015.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
										//if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
											P[icount1].j = inode;
											P[icount1].i = C_numerate[Amat.j[is0]];
											//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
											if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat.j[is0]]  ) {
										

										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F) {

											if (Amat.j[is0] != inode) {
											// ������������� Strong F �����.

											// �����:
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
												// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
												// � ������� F ������.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

												// ������������ ������ F ����. 
												//���� ��� ����� � ���� �� ����� �� ��������� �� �����������������.
												integer iFpoint = Amat.j[is0];
												if (fabs(SumPall) < RealZERO) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//getchar();
													system("PAUSE");
													exit(1);
												}
												doublerealT multiplyer_nu = Amat.abs_aij[is0] / SumPall;
												// ��������� ���� ������� iFpointeger 
												// ����� ����� ����� ��� � ����.

												// �������������� ��������� �����.
												doublerealT maxelem_threshold_loc = -1.0;
												//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
												//integer ii1_loc = ii2;// row_startA[iFpoint];
												//integer istr_etalon = iFpoint;//  Amat.i[ii1_loc];
												//integer iend_for = iend_marker_position2;// -1;

												// ����� ������� ������� F-node ������� C-node.
												integer icsos_loc = 0;

												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
												// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
												doublerealT sumP_loc = 0.0;

												//integer iend_marker_position_loc = iend_marker_position2;// row_startA[iFpoint + 1] - 1;
												//for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); is0_loc++) {
												for (integer is0_loc = ii2; (is0_loc <= iend_marker_position2); is0_loc++) {
													//iend_for = is0_loc;
													//if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
													// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
													// ����� ������ ����������.
													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
															// �� ��������� F ����. ������ ����� �� ��������� �� �����������������. 
															//if (Amat.j[is0_loc] != iFpoint) {
															if (Amat.abs_aij[is0_loc] > maxelem_threshold_loc) {
																maxelem_threshold_loc = Amat.abs_aij[is0_loc];
															}
																//sumP_loc += Amat.abs_aij[is0_loc]; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																//icsos_loc++;
															//}
														}
													}
													//}
												}

												doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
												// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

												
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii2; is0_loc <= iend_marker_position2; is0_loc++) {

													// ����� �������� ����� ���������� ����� ���������.
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

													if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
														if (this_is_C_node[Amat.j[is0_loc]]  ) {

															// �� ��������� F ����. ������ ����� �� ��������� �� �����������������.
															//if (Amat.j[is0_loc] != iFpoint) {

																sumP_loc += Amat.abs_aij[is0_loc]; // ����� ������� ��������������� ��������� ������� ����������� � �����.
																icsos_loc++;
															//}

														}
														else {

															//if (Amat.j[is0_loc] != iFpoint) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}
														}
													}

													}
												}
												

												//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												// � ����� ��� ������� ���������������� ����� 
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
												for (integer is0_loc = ii2; is0_loc <= iend_marker_position2; is0_loc++) {

													// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


													// ��������������� ������� �� ��������� � �����.

													// ������ ������� ������ ����������� ��������� 
													// ������������� ��������� �������� �� �������� 
													// �������� �������.
													// ����������� 5 ������� 2015.
													if (Amat.abs_aij[is0_loc] > maxelem_threshold_loc_magic) 
													{// ������ ��� ������� - ���������� ��������
													if (this_is_Strong_C_node[Amat.j[is0_loc]]) 
													{ // ������������ ��������

														//if (this_is_C_node[Amat.j[is0_loc]]  )
														{// ����� ������� - ���������� �� ��������
															//if (Amat.j[is0_loc] != iFpoint) {

																P[icount1].j = inode;
																P[icount1].i = C_numerate[Amat.j[is0_loc]];
																//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
																if (fabs(sumP_loc) < RealZERO) {
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
						if (this_is_C_node[Amat.j[is0]]  ) {
							if (Amat.abs_aij[is0] > maxelem_threshold_theta) {
								this_is_Strong_C_node[Amat.j[is0]] = false;
							}
						}
					}
				}
				*/
				// ����� - ������� ���-�������.
				for (integer is00 = 0; is00 < isize_iStrongC_nodes_ind; is00++) {
					this_is_Strong_C_node[Amat.j[iStrongC_nodes_ind[is00]]] = false;
				}


			}
			////
		}
		else {
			// only negative connections.

			// ���������� ����������� ��� ����� ������� ���������� F-nodes.
			// ������ F-nodes ������ C-nodes.
			for (integer inode = 1; inode <= n_a[ilevel - 1]; inode++) if (this_is_F_node[inode]  ) {

				// ��� ����� ����������� �������� ��������.
				// 5 ������� 2015 ���� �� ��������� ��������� �������������
				// ��������� ������������ � ������ � ��������� ��������.
				doublerealT maxelem_threshold = -1.0;
				//integer ii1 = BinarySearchAi(Amat, inode, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer ii1 = row_startA[inode];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = row_startA[istr_etalon1 + 1] - 1;// -1;
				if (!btreshold_on_new_vetv) {
					for (integer is0 = ii1; (is0 <= iend_for1); is0++) {
						if (Amat.j[is0] != istr_etalon1) {
							// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
							//if (this_is_C_node[Amat.j[is0]]  ) {
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
				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.

				// ����� ������� ������� F-node ������� C-node.
				integer icsos = 0;
				integer icsosF = 0;

				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
					//if (Amat.j[is0] != Amat.i[ii1]) {
					// ������������ ���� �� ��������� F ����
						if (this_is_C_node[Amat.j[is0]]  ) {
							if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								this_is_Strong_C_node[Amat.j[is0]] = true;
							}
						}
					//}
				}


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumP = 0.0;
				doublerealT SumPall = 0.0;
				integer icount_StronglyF = 0;
				//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); is0++) {
				for (integer is0 = ii1; is0 <= iend_for1; is0++) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]  ) {
							//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
							//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
							if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
								sumP += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
								icsos++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]  ) {
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


										SumPall += Amat.abs_aij[is0]; // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
										icount_StronglyF++;
										icsosF++;
									}
								}
							}
							// ������������ ���������� ������� ������� �� �������� � ������.
							the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
						}
					}
				}
				if (icsos == 1) {
					number_of_F_nodes_with_one_single_strong_C_neighbor++; // ���������� F ����� � ����� ������������ ������� � �������.
																		   // ��������� ������ ������ "����������".
																		   // ���������� ������ ����������� ��� ���������.
																		   // � ������� �� �������� ������� ����� ���������� ������ ���������� ���� ��������� ������� � 
																		   // ������������ �� ���� ������� ����� ��������.
					if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C ������� � � ����-�� �� ������� ������� F �������.
				}



				{

					if (1) {
						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
						//
						SumPall += sumP;

						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); is0++) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								// ��� ���������� ������ Strong �����.


								if (this_is_C_node[Amat.j[is0]]  ) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
										if (fabs(sumP) < RealZERO) {
											//printf("error interpolation zero diagonal sumP.\n");
											//printf("Fnode all sosed is F");
											//system("pause");
											//	printf("inode is Dirichlet node\n");
											if (this_is_C_node[inode] == false) iadditionalCstatistic++;
											this_is_F_node[inode] = false; // ���� ���� ������� ������ � �����.
											this_is_C_node[inode] = true;
											bweSholdbeContinue = true;
											//exit(1);
											// ����� ����� �������� ������������.
										}
										else {
											// ��� ��� ��� ������������ Strong C �����. 
											// ��������������� ������� �� ��������� � �����.

											// ������ ������� ������ ����������� ��������� 
											// ������������� ��������� �������� �� �������� 
											// �������� �������.
											// ����������� 5 ������� 2015.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
											if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta)) {
												P[icount1].j = inode;
												P[icount1].i = C_numerate[Amat.j[is0]];
												//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
												if (fabs(SumPall) < RealZERO) {
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
									if (this_is_F_node[Amat.j[is0]]  ) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
										if ((Amat.aij[is0]<0.0) && (Amat.abs_aij[is0] > maxelem_threshold_theta_strong_F)) {
											// ������������� Strong F �����.

											// �����:
											// 



											//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
											// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
											// � ������� F ������.
											//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

											integer iFpoint = Amat.j[is0];
											doublerealT multiplyer_nu = Amat.abs_aij[is0] / SumPall;
											// ��������� ���� ������� iFpointeger 
											// ����� ����� ����� ��� � ����.

											// �������������� ��������� �����.
											doublerealT maxelem_threshold_loc = -1.0;
											//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
											integer ii1_loc = row_startA[iFpoint];
											integer istr_etalon = Amat.i[ii1_loc];
											integer iend_for = -1;
											integer iend_marker_position = row_startA[istr_etalon + 1] - 1;
											for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position); is0_loc++) {
												iend_for = is0_loc;
												//if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												if (this_is_C_node[Amat.j[is0_loc]]  ) {
													if (this_is_Strong_C_node[Amat.j[is0_loc]]) {
														if (Amat.j[is0_loc] != istr_etalon) {
															maxelem_threshold_loc = Amat.abs_aij[is0_loc];
														}
													}
												}
												//}
											}

											doublerealT maxelem_threshold_loc_magic = maxelem_threshold_loc * magic;
											// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

											// ����� ������� ������� F-node ������� C-node.
											integer icsos_loc = 0;

											// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
											// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
											doublerealT sumP_loc = 0.0;
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ����� �������� ����� ���������� ����� ���������.
												//if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat.j[is0_loc]]) {


													if (this_is_C_node[Amat.j[is0_loc]]  ) {

														if (Amat.j[is0_loc] != istr_etalon) {

															sumP_loc += Amat.abs_aij[is0_loc]; // ����� ������� ��������������� ��������� ������� ����������� � �����.
															icsos_loc++;
														}

													}
													else {

														//if (Amat.j[is0_loc] != istr_etalon) {
														// ������������ ���������� ������� ������� �� �������� � ������.
														//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
														//}
													}
												}
												//}
											}

											doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

											// � ����� ��� ������� ���������������� ����� 
											//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); is0_loc++) {
											for (integer is0_loc = ii1_loc; is0_loc <= iend_for; is0_loc++) {

												// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


												// ��������������� ������� �� ��������� � �����.

												// ������ ������� ������ ����������� ��������� 
												// ������������� ��������� �������� �� �������� 
												// �������� �������.
												// ����������� 5 ������� 2015.
												//if ((Amat.aij[is0_loc]<0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {
												if (this_is_Strong_C_node[Amat.j[is0_loc]]) {


													if (this_is_C_node[Amat.j[is0_loc]]  ) {
														if (Amat.j[is0_loc] != istr_etalon) {

															P[icount1].j = inode;
															P[icount1].i = C_numerate[Amat.j[is0_loc]];
															//P[icount1].aij = fabs(Amat.aij[is0]) / sumP;
															if (fabs(sumP_loc) < RealZERO) {
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
						if (this_is_C_node[Amat.j[is0]]  ) {
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
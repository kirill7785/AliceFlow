// �������� �� ������ �� ��������� � ���� ? 
// Future TODO: 19-20.10.2021 FFFFFC ������. ����� F ����� ������.
// CFFFFFFC
// �������� ��������� � PMIS ����������� �� �������� ������� � (������ �������).
template <typename doublerealT>
void my_long_rang_interpolation_procedure_FFFFFFC_number3A_parallelx8(integer& the_number_of_neighbors_that_are_not_C_nodes_gl,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor_gl,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic_gl,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF_gl,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative, bool bsuffix_work, bool baggressive_interpol)
{

#ifdef _OPENMP
	int inth = number_cores();
	if (inth > 1) {
		omp_set_num_threads(inth);
	}
	//omp_set_num_threads(8);
#else 
	// ���� �����.
	int inth = 1;
#endif

	bool old_def_threshold = true;

	static int ivprjsk = 0;
	++ivprjsk;


	const int ilimCadd = 4;

	// ������
	int iadditionalCstatistic = 0;
	int iadds0 = 0;

	int iOPORNAQ_C0 = 0;
	int iOPORNAQ_C1 = 0;
	int iOPORNAQ_C2 = 0;
	int iOPORNAQ_C3 = 0;
	int iOPORNAQ_C4 = 0;
	int iOPORNAQ_C5 = 0;

	if (ivprjsk % 2 != 0) {


#pragma omp parallel for reduction(+: iadditionalCstatistic, iadds0, iOPORNAQ_C0, iOPORNAQ_C1, iOPORNAQ_C2, iOPORNAQ_C3, iOPORNAQ_C4, iOPORNAQ_C5)
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) if (this_is_F_node[i8]) {

			bool bproblem = false;

			int iOPORNAQ_C0_1 = 0;
			int iOPORNAQ_C1_1 = 0;
			int iOPORNAQ_C2_1 = 0;
			int iOPORNAQ_C3_1 = 0;
			int iOPORNAQ_C4_1 = 0;
			int iOPORNAQ_C5_1 = 0;

			int i_na_odnoi_vetke_C2 = 0;

			integer ii1 = row_startA[i8];
			integer istr_etalon1 = Amat.i[ii1];
			integer iend_for1 = -1;
			for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
				iend_for1 = is0;
			}

			if (iend_for1 != ii1) {
				// ������ ����.
				doublerealT theta_strong_F = theta83; // 3 ���� 2016
				doublerealT theta = theta_gl;

				const doublerealT magic = magic82; // 0.4 is recomended.
				const doublerealT magic_FF = magic82; // 18.09.2020
				//const doublerealT magic_FFF = magic82; // 24.09.2021

				int iOPORNAQ_C = 0;

				// ��� ����� ����������� �������� ��������.
							// 5 ������� 2015 ���� �� ��������� ��������� �������������
							// ��������� ������������ � ������ � ��������� ��������.



				doublerealT maxelem_threshold = -1.0;

				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
					iend_for1 = is0;
				}
				if (old_def_threshold) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
						if ((Amat.j[is0] != istr_etalon1) && (Amat.aij[is0] < 0.0) && (-Amat.aij[is0] > maxelem_threshold)) maxelem_threshold = -Amat.aij[is0];
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}

				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumPC = 0.0;
				doublerealT sumPF = 0.0;
				doublerealT SumPall = 0.0;

				// ����� 1.1

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumPC += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
								iOPORNAQ_C++;
								iOPORNAQ_C0++;
								iOPORNAQ_C0_1++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

									// ������������� Strong F �����.

									sumPF += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
								}
							}
						}
					}
				}

				if (iOPORNAQ_C >= ilimCadd) {
					SumPall = sumPC;

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {

								}
							}
						}
					}

				}
				else {

					SumPall = sumPC + sumPF;


					// ����� 1.2

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {

								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {

									if (i_na_odnoi_vetke_C2 != 2) {
										i_na_odnoi_vetke_C2 = 0;
									}

									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

										// ������������� Strong F �����.

										integer i8_1 = Amat.j[is0];


										doublerealT nu1 = fabs(Amat.aij[is0]) / SumPall;

										// ����� ���������� ����� ��� ���� i8_1.

										//********************************************************

										doublerealT maxelem_threshold_1 = -1.0;

										integer ii1_1 = row_startA[i8_1];
										integer istr_etalon1_1 = Amat.i[ii1_1];
										integer iend_for1_1 = -1;
										for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
											iend_for1_1 = is0_1;
										}
										if (old_def_threshold) {
											for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1_1) && (Amat.aij[is0_1] < 0.0) && (-Amat.aij[is0_1] > maxelem_threshold_1)) maxelem_threshold_1 = -Amat.aij[is0_1];
											}
										}
										else {
											maxelem_threshold_1 = threshold_quick_all[istr_etalon1_1];
										}

										// ����� maxelem_threshold_1 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



										doublerealT maxelem_threshold_theta_1 = maxelem_threshold_1 * theta;
										doublerealT maxelem_threshold_theta_strong_F_1 = maxelem_threshold_1 * magic;


										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
										doublerealT sumPC_1 = 0.0;
										doublerealT sumPF_1 = 0.0;
										doublerealT SumPall_1 = 0.0;

										// ����� 2.1

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														sumPC_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
														iOPORNAQ_C++;
														iOPORNAQ_C1++;
														iOPORNAQ_C1_1++;
														i_na_odnoi_vetke_C2++;
													}
												}
												else {
													if (this_is_F_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

															// ������������� Strong F �����.

															sumPF_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
														}
													}
												}
											}
										}


										if (iOPORNAQ_C >= ilimCadd) {
											SumPall_1 = sumPC_1;

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {

														}
													}
												}
											}

										}
										else {

											SumPall_1 = sumPC_1 + sumPF_1;


											// ����� 2.2

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {

														}
													}
													else {
														if (this_is_F_node[Amat.j[is0_1]]) {

															if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

																// ������������� Strong F �����.

																integer i8_2 = Amat.j[is0_1];


																doublerealT nu2 = fabs(Amat.aij[is0_1]) / SumPall_1;

																// ����� ���������� ����� ��� ���� i8_2.


																//********************************************************

																doublerealT maxelem_threshold_2 = -1.0;

																integer ii1_2 = row_startA[i8_2];
																integer istr_etalon1_2 = Amat.i[ii1_2];
																integer iend_for1_2 = -1;
																for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																	iend_for1_2 = is0_2;
																}
																if (old_def_threshold) {
																	for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1_2) && (Amat.aij[is0_2] < 0.0) && (-Amat.aij[is0_2] > maxelem_threshold_2)) maxelem_threshold_2 = -Amat.aij[is0_2];
																	}
																}
																else {
																	maxelem_threshold_2 = threshold_quick_all[istr_etalon1_2];
																}

																// ����� maxelem_threshold_2 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																doublerealT maxelem_threshold_theta_2 = maxelem_threshold_2 * theta;
																doublerealT maxelem_threshold_theta_strong_F_2 = maxelem_threshold_2 * magic_FF;


																// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																doublerealT sumPC_2 = 0.0;
																doublerealT sumPF_2 = 0.0;
																doublerealT SumPall_2 = 0.0;

																// ����� 3.1

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				sumPC_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
																				iOPORNAQ_C++;
																				iOPORNAQ_C2++;
																				iOPORNAQ_C2_1++;
																				i_na_odnoi_vetke_C2++;
																			}
																		}
																		else {
																			if (this_is_F_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																					// ������������� Strong F �����.

																					sumPF_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																				}
																			}
																		}
																	}
																}

																if (iOPORNAQ_C >= ilimCadd) {
																	SumPall_2 = sumPC_2;

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {

																				}
																			}
																		}
																	}

																}
																else {

																	SumPall_2 = sumPC_2 + sumPF_2;


																	// ����� 3.2

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {

																				}
																			}
																			else {
																				if (this_is_F_node[Amat.j[is0_2]]) {

																					if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																						// ������������� Strong F �����.

																						integer i8_3 = Amat.j[is0_2];


																						doublerealT nu3 = fabs(Amat.aij[is0_2]) / SumPall_2;


																						///******************************************************************
																						// ����� ���������� ����� ��� ���� i8_3.

										                                               //********************************************************

																						doublerealT maxelem_threshold_3 = -1.0;

																						integer ii1_3 = row_startA[i8_3];
																						integer istr_etalon1_3 = Amat.i[ii1_3];
																						integer iend_for1_3 = -1;
																						for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																							iend_for1_3 = is0_3;
																						}
																						if (old_def_threshold) {
																							for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1_3) && (Amat.aij[is0_3] < 0.0) && (-Amat.aij[is0_3] > maxelem_threshold_3)) maxelem_threshold_3 = -Amat.aij[is0_3];
																							}
																						}
																						else {
																							maxelem_threshold_3 = threshold_quick_all[istr_etalon1_3];
																						}

																						// ����� maxelem_threshold_3 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																						doublerealT maxelem_threshold_theta_3 = maxelem_threshold_3 * theta;
																						doublerealT maxelem_threshold_theta_strong_F_3 = maxelem_threshold_3 * magic;


																						// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																						// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																						doublerealT sumPC_3 = 0.0;
																						doublerealT sumPF_3 = 0.0;
																						doublerealT SumPall_3 = 0.0;

																						// ����� 4.1

																						for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																								if (this_is_C_node[Amat.j[is0_3]]) {

																									if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										sumPC_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																										iOPORNAQ_C++;
																										iOPORNAQ_C3++;
																										iOPORNAQ_C3_1++;
																										i_na_odnoi_vetke_C2++;
																									}
																								}
																								else {
																									if (this_is_F_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																											// ������������� Strong F �����.

																											sumPF_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																										}
																									}
																								}
																							}
																						}


																						if (iOPORNAQ_C >= ilimCadd) {
																							SumPall_3 = sumPC_3;

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {

																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {

																										}
																									}
																								}
																							}

																						}
																						else {

																							SumPall_3 = sumPC_3 + sumPF_3;


																							// ����� 4.2

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {

																										}
																									}
																									else {
																										if (this_is_F_node[Amat.j[is0_3]]) {

																											if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																												// ������������� Strong F �����.

																												integer i8_4 = Amat.j[is0_3];


																												doublerealT nu4 = fabs(Amat.aij[is0_3]) / SumPall_3;



																												// ����� ���������� ����� ��� ���� i8_4.

										                                                                        //********************************************************

																												doublerealT maxelem_threshold_4 = -1.0;

																												integer ii1_4 = row_startA[i8_4];
																												integer istr_etalon1_4 = Amat.i[ii1_4];
																												integer iend_for1_4 = -1;
																												for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																													iend_for1_4 = is0_4;
																												}
																												if (old_def_threshold) {
																													for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1_4) && (Amat.aij[is0_4] < 0.0) && (-Amat.aij[is0_4] > maxelem_threshold_4)) maxelem_threshold_4 = -Amat.aij[is0_4];
																													}
																												}
																												else {
																													maxelem_threshold_4 = threshold_quick_all[istr_etalon1_4];
																												}

																												// ����� maxelem_threshold_4 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																												doublerealT maxelem_threshold_theta_4 = maxelem_threshold_4 * theta;
																												doublerealT maxelem_threshold_theta_strong_F_4 = maxelem_threshold_4 * magic;


																												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																												// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																												doublerealT sumPC_4 = 0.0;
																												doublerealT sumPF_4 = 0.0;
																												doublerealT SumPall_4 = 0.0;

																												// ����� 5.1

																												for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																													if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_4)) {
																														if (this_is_C_node[Amat.j[is0_4]]) {

																															if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																sumPC_4 += fabs(Amat.aij[is0_4]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																																iOPORNAQ_C++;
																																iOPORNAQ_C4++;
																																iOPORNAQ_C4_1++;
																																i_na_odnoi_vetke_C2++;
																															}
																														}
																														else {
																															if (this_is_F_node[Amat.j[is0_4]]) {

																																if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_strong_F_4) {

																																	// ������������� Strong F �����.

																																	sumPF_4 += fabs(Amat.aij[is0_4]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																																}
																															}
																														}
																													}
																												}


																												if (iOPORNAQ_C >= ilimCadd) {
																													SumPall_4 = sumPC_4;

																													for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																															if (this_is_C_node[Amat.j[is0_4]]) {

																																if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {

																																}
																															}
																														}
																													}

																												}
																												else {

																													SumPall_4 = sumPC_4 + sumPF_4;


																													// ����� 5.2

																													for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																															if (this_is_C_node[Amat.j[is0_4]]) {

																																if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {

																																}
																															}
																															else {
																																if (this_is_F_node[Amat.j[is0_4]]) {

																																	if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_strong_F_4) {

																																		// ������������� Strong F �����.

																																		integer i8_5 = Amat.j[is0_4];


																																		doublerealT nu5 = fabs(Amat.aij[is0_4]) / SumPall_4;

																																		


																																		// ����� ���������� ����� ��� ���� i8_5.




																																		//********************************************************

																																		doublerealT maxelem_threshold_5 = -1.0;

																																		integer ii1_5 = row_startA[i8_5];
																																		integer istr_etalon1_5 = Amat.i[ii1_5];
																																		integer iend_for1_5 = -1;
																																		for (integer is0_5 = ii1_5; (is0_5 <= row_startA[istr_etalon1_5 + 1] - 1); ++is0_5) {
																																			iend_for1_5 = is0_5;
																																		}
																																		if (old_def_threshold) {
																																			for (integer is0_5 = ii1_5; (is0_5 <= row_startA[istr_etalon1_5 + 1] - 1); ++is0_5) {
																																				if ((Amat.j[is0_5] != istr_etalon1_5) && (Amat.aij[is0_5] < 0.0) && (-Amat.aij[is0_5] > maxelem_threshold_5)) maxelem_threshold_5 = -Amat.aij[is0_5];
																																			}
																																		}
																																	else {
																																		maxelem_threshold_5 = threshold_quick_all[istr_etalon1_5];
																																	}

																																	// ����� maxelem_threshold_5 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																																	doublerealT maxelem_threshold_theta_5 = maxelem_threshold_5 * theta;
																																	//doublerealT maxelem_threshold_theta_strong_F_5 = maxelem_threshold_5 * magic_FFF;


																																	// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																																	// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																																	doublerealT sumPC_5 = 0.0;
																																	//doublerealT sumPF_5 = 0.0;
																																	//doublerealT SumPall_5 = 0.0;

																																	// ����� 6.1

																																	for (integer is0_5 = ii1_5; is0_5 <= iend_for1_5; ++is0_5) {
																																		if ((Amat.j[is0_5] != istr_etalon1) && (Amat.j[is0_5] != istr_etalon1_1) && (Amat.j[is0_5] != istr_etalon1_2) && (Amat.j[is0_5] != istr_etalon1_3) && (Amat.j[is0_5] != istr_etalon1_4) && (Amat.j[is0_5] != istr_etalon1_5)) {
																																			if (this_is_C_node[Amat.j[is0_5]]) {

																																				if (fabs(Amat.aij[is0_5]) > maxelem_threshold_theta_5) {
																																					sumPC_5 += fabs(Amat.aij[is0_5]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																																					iOPORNAQ_C++;
																																					iOPORNAQ_C5++;
																																					iOPORNAQ_C5_1++;
																																					i_na_odnoi_vetke_C2++;
																																				}
																																			}

																																		}
																																	}


																																	//SumPall_5 = sumPC_5;


																																	// ����� 6.2

																																	for (integer is0_5 = ii1_5; is0_5 <= iend_for1_5; ++is0_5) {
																																		if ((Amat.j[is0_5] != istr_etalon1) && (Amat.j[is0_5] != istr_etalon1_1) && (Amat.j[is0_5] != istr_etalon1_2) && (Amat.j[is0_5] != istr_etalon1_3) && (Amat.j[is0_5] != istr_etalon1_4) && (Amat.j[is0_5] != istr_etalon1_5)) {
																																			if (this_is_C_node[Amat.j[is0_5]]) {
		
																																				if (fabs(Amat.aij[is0_5]) > maxelem_threshold_theta_5) {

																																			    }
																																		    }

																																	    }
																																    }

																																//********************************************************

																																// ��������� ������ ���������� i8_5.
																																//********************************************************

																																	}
																																}
																															}
																														}
																													}
																						                        }

																												//********************************************************

																												// ��������� ������ ���������� i8_4.
																												//********************************************************

																											}
																										}
																									}
																								}
																							}
																						}

																						//*******************************************18.3

																					// ��������� ������ ���������� i8_3.

																					}
																				}
																			}
																		}
																	}


																}
																//********************************************************

																// ��������� ������ ���������� i8_2.



															}
														}
													}
												}
											}


										}
										//********************************************************
										// ��������� ������ ���������� i8_1.



									}
								}
							}
						}
					}

				}

				// ��� ����� ��������� �������� ����� �������� ������� �����.
				// ����������� ����� �������� ����� � ������������ � ����� ����������.
				if ((iOPORNAQ_C == 0) ||
					((iOPORNAQ_C == 1) && (iOPORNAQ_C0_1 == 0) && (iOPORNAQ_C1_1 == 0)) ||
					((iOPORNAQ_C == 2) && (i_na_odnoi_vetke_C2==2) && (iOPORNAQ_C0_1 == 0) && (iOPORNAQ_C1_1 == 0))
					) {

					bproblem = true;

				}
				// ���� ��� � ����� ������, ���� ���� ������ ���� � �����, �� �� ���������� ������ �� ���������� 2 � ����.
				// ���� ���� ������ ��� ���� ������ �� ����� ����� � ������ �������.
			}

			if (ii1 == (row_startA[istr_etalon1 + 1] - 1)) bproblem = false; // ������� ������� �� �� ��������� � � ����. 01.10.2021
			if (ii1 + 1 == (row_startA[istr_etalon1 + 1] - 1)) bproblem = false; // ���������� ������� ������� 01.10.2021
			//if (ii1 + 2 <= (row_startA[istr_etalon1 + 1] - 1)) {
				// ���� ���� ��� ������
				//if (iOPORNAQ_C == 0) {
					//� ����� ������� ���� �������� � ����.
				//    bproblem = false;
				//}
		//	}


			if (bproblem)
			{

				//int ir = rand() % 1000;
				//if (ir < 1000) {

					// ������ 10%.

					this_is_C_node[istr_etalon1] = true;
					this_is_F_node[istr_etalon1] = false;

					iadditionalCstatistic++;
				//}
			}

		}

		iadditionalCstatistic_gl += iadditionalCstatistic;
		iadds0 = iadditionalCstatistic;
		iadditionalCstatistic = 0;

		if (my_amg_manager.iprint_log) {

			std::cout << "iadditional C=" << iadds0 << std::endl;

			std::cout << iOPORNAQ_C0 << " " << iOPORNAQ_C1 << " " << iOPORNAQ_C2 << " " << iOPORNAQ_C3 << " " << iOPORNAQ_C4 << " " << iOPORNAQ_C5 << " " << std::endl;
		}

		bweSholdbeContinue = true;

	}
	else {
		bweSholdbeContinue = false;
	}

	//if (iadditionalCstatistic == 0) 
	{

		doublerealT theta_strong_F = theta83; // 3 ���� 2016
		doublerealT theta = theta_gl;

		const doublerealT magic = magic82; // 0.4 is recomended.
		const doublerealT magic_FF = magic82; // 18.09.2020
		//const doublerealT magic_FFF = magic82; // 24.09.2021


		//int itid = 0;

		integer LOC_P_LIM_SIZE = -1;
		double dsize = 2.2 * (1.0 / (1.0 * inth)) * nsizePR * n_a[0];
		LOC_P_LIM_SIZE = static_cast<integer>(dsize);

		Ak1** P_loc = new Ak1 * [inth];
		integer* isize_P_loc = new integer[inth];
		for (int itid = 0; itid < inth; ++itid)
		{
			isize_P_loc[itid] = 0;

			P_loc[itid] = new Ak1[LOC_P_LIM_SIZE];
		}


		//int istat[5] = { 0,0,0,0,0 };
		//istat[4] = iadds0;

		iOPORNAQ_C0 = 0;
		iOPORNAQ_C1 = 0;
		iOPORNAQ_C2 = 0;
		iOPORNAQ_C3 = 0;
		iOPORNAQ_C4 = 0;
		iOPORNAQ_C5 = 0;

#pragma omp parallel for reduction(+:  iOPORNAQ_C0, iOPORNAQ_C1, iOPORNAQ_C2, iOPORNAQ_C3, iOPORNAQ_C4, iOPORNAQ_C5)
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8)
		{
#ifdef _OPENMP
			int itid = omp_get_thread_num();
#else
			int itid = 0;
#endif

			if (this_is_F_node[i8]) {


				int iOPORNAQ_C = 0;

				// ��� ����� ����������� �������� ��������.
							// 5 ������� 2015 ���� �� ��������� ��������� �������������
							// ��������� ������������ � ������ � ��������� ��������.



				doublerealT maxelem_threshold = -1.0;

				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
					iend_for1 = is0;
				}
				if (old_def_threshold) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
						if ((Amat.j[is0] != istr_etalon1) && (Amat.aij[is0] < 0.0) && (-Amat.aij[is0] > maxelem_threshold)) maxelem_threshold = -Amat.aij[is0];
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}

				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumPC = 0.0;
				doublerealT sumPF = 0.0;
				doublerealT SumPall = 0.0;

				// ����� 1.1

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumPC += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
								iOPORNAQ_C++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

									// ������������� Strong F �����.

									sumPF += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
								}
							}
						}
					}
				}

				if (iOPORNAQ_C >= ilimCadd) {
					SumPall = sumPC;

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									Ak1 Ptmp;

									Ptmp.j = static_cast<integer_mix_precision>(i8);
									Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
									//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
									if (fabs(SumPall) < 1.0e-33) {
										printf("error 1.0.1 ! division by zero. SumPall =%e\n", SumPall);
										//system("pause");
										system("PAUSE");
										exit(1);
									}
									Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
				else {

					SumPall = sumPC + sumPF;


					// ����� 1.2

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									Ak1 Ptmp;

									Ptmp.j = static_cast<integer_mix_precision>(i8);
									Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
									//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
									if (fabs(SumPall) < 1.0e-33) {
										printf("error 1.0.2 ! division by zero. SumPall =%e\n", SumPall);
										//system("pause");
										system("PAUSE");
										exit(1);
									}
									Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
							else {
								if (this_is_F_node[Amat.j[is0]]) {

									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

										// ������������� Strong F �����.

										integer i8_1 = Amat.j[is0];


										doublerealT nu1 = fabs(Amat.aij[is0]) / SumPall;

										// ����� ���������� ����� ��� ���� i8_1.

										//********************************************************

										doublerealT maxelem_threshold_1 = -1.0;

										integer ii1_1 = row_startA[i8_1];
										integer istr_etalon1_1 = Amat.i[ii1_1];
										integer iend_for1_1 = -1;
										for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
											iend_for1_1 = is0_1;
										}
										if (old_def_threshold) {
											for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1_1) && (Amat.aij[is0_1] < 0.0) && (-Amat.aij[is0_1] > maxelem_threshold_1)) maxelem_threshold_1 = -Amat.aij[is0_1];
											}
										}
										else {
											maxelem_threshold_1 = threshold_quick_all[istr_etalon1_1];
										}

										// ����� maxelem_threshold_1 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



										doublerealT maxelem_threshold_theta_1 = maxelem_threshold_1 * theta;
										doublerealT maxelem_threshold_theta_strong_F_1 = maxelem_threshold_1 * magic;


										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
										doublerealT sumPC_1 = 0.0;
										doublerealT sumPF_1 = 0.0;
										doublerealT SumPall_1 = 0.0;

										// ����� 2.1

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														sumPC_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
														iOPORNAQ_C++;
													}
												}
												else {
													if (this_is_F_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

															// ������������� Strong F �����.

															sumPF_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
														}
													}
												}
											}
										}


										if (iOPORNAQ_C >= ilimCadd) {
											SumPall_1 = sumPC_1;

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															Ak1 Ptmp;

															Ptmp.j = static_cast<integer_mix_precision>(i8);
															Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_1]]);
															//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
															if (fabs(SumPall_1) < 1.0e-33) {
																printf("error 1.0.3 ! division by zero. SumPall =%e\n", SumPall_1);
																//system("pause");
																system("PAUSE");
																exit(1);
															}
															Ptmp.aij = static_cast<real_mix_precision>(nu1 * (fabs(Amat.aij[is0_1]) / SumPall_1));

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
										else {

											SumPall_1 = sumPC_1 + sumPF_1;


											// ����� 2.2

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															Ak1 Ptmp;

															Ptmp.j = static_cast<integer_mix_precision>(i8);
															Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_1]]);

															if (fabs(SumPall_1) < 1.0e-33) {
																printf("error 2.0 ! division by zero. SumPall_1 =%e\n", SumPall_1);
																//system("pause");
																system("PAUSE");
																exit(1);
															}
															Ptmp.aij = static_cast<real_mix_precision>(nu1 * (fabs(Amat.aij[is0_1]) / SumPall_1));

															if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																P_loc[itid][isize_P_loc[itid]] = Ptmp;
																isize_P_loc[itid]++;
															}
															else {
																std::cout << "Error 2.0 Increase size P_loc memory limit\n" << std::endl;
																system("PAUSE");
															}
														}
													}
													else {
														if (this_is_F_node[Amat.j[is0_1]]) {

															if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

																// ������������� Strong F �����.

																integer i8_2 = Amat.j[is0_1];


																doublerealT nu2 = fabs(Amat.aij[is0_1]) / SumPall_1;

																// ����� ���������� ����� ��� ���� i8_2.


																//********************************************************

																doublerealT maxelem_threshold_2 = -1.0;

																integer ii1_2 = row_startA[i8_2];
																integer istr_etalon1_2 = Amat.i[ii1_2];
																integer iend_for1_2 = -1;
																for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																	iend_for1_2 = is0_2;
																}
																if (old_def_threshold) {
																	for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1_2) && (Amat.aij[is0_2] < 0.0) && (-Amat.aij[is0_2] > maxelem_threshold_2)) maxelem_threshold_2 = -Amat.aij[is0_2];
																	}
																}
																else {
																	maxelem_threshold_2 = threshold_quick_all[istr_etalon1_2];
																}

																// ����� maxelem_threshold_2 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																doublerealT maxelem_threshold_theta_2 = maxelem_threshold_2 * theta;
																doublerealT maxelem_threshold_theta_strong_F_2 = maxelem_threshold_2 * magic_FF;


																// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																doublerealT sumPC_2 = 0.0;
																doublerealT sumPF_2 = 0.0;
																doublerealT SumPall_2 = 0.0;

																// ����� 3.1

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				sumPC_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
																				iOPORNAQ_C++;
																			}
																		}
																		else {
																			if (this_is_F_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																					// ������������� Strong F �����.

																					sumPF_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																				}
																			}
																		}
																	}
																}

																if (iOPORNAQ_C >= ilimCadd) {
																	SumPall_2 = sumPC_2;

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					Ak1 Ptmp;

																					Ptmp.j = static_cast<integer_mix_precision>(i8);
																					Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_2]]);
																					//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																					if (fabs(SumPall_2) < 1.0e-33) {
																						printf("error 1.0.4 ! division by zero. SumPall =%e\n", SumPall_2);
																						//system("pause");
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * (fabs(Amat.aij[is0_2]) / SumPall_2));

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
																else {

																	SumPall_2 = sumPC_2 + sumPF_2;


																	// ����� 3.2

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					Ak1 Ptmp;

																					Ptmp.j = static_cast<integer_mix_precision>(i8);
																					Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_2]]);

																					if (fabs(SumPall_2) < 1.0e-33) {
																						printf("error 3.0 ! division by zero. SumPall_2 =%e\n", SumPall_2);
																						//system("pause");
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * (fabs(Amat.aij[is0_2]) / SumPall_2));

																					if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																						P_loc[itid][isize_P_loc[itid]] = Ptmp;
																						isize_P_loc[itid]++;
																					}
																					else {
																						std::cout << "Error 3.0 Increase size P_loc memory limit\n" << std::endl;
																						system("PAUSE");
																					}
																				}
																			}
																			else {
																				if (this_is_F_node[Amat.j[is0_2]]) {

																					if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																						// ������������� Strong F �����.

																						integer i8_3 = Amat.j[is0_2];


																						doublerealT nu3 = fabs(Amat.aij[is0_2]) / SumPall_2;


																						///******************************************************************
																						// ����� ���������� ����� ��� ���� i8_3.

										                                                //********************************************************

																						doublerealT maxelem_threshold_3 = -1.0;

																						integer ii1_3 = row_startA[i8_3];
																						integer istr_etalon1_3 = Amat.i[ii1_3];
																						integer iend_for1_3 = -1;
																						for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																							iend_for1_3 = is0_3;
																						}
																						if (old_def_threshold) {
																							for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1_3) && (Amat.aij[is0_3] < 0.0) && (-Amat.aij[is0_3] > maxelem_threshold_3)) maxelem_threshold_3 = -Amat.aij[is0_3];
																							}
																						}
																						else {
																							maxelem_threshold_3 = threshold_quick_all[istr_etalon1_3];
																						}

																						// ����� maxelem_threshold_3 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																						doublerealT maxelem_threshold_theta_3 = maxelem_threshold_3 * theta;
																						doublerealT maxelem_threshold_theta_strong_F_3 = maxelem_threshold_3 * magic;


																						// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																						// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																						doublerealT sumPC_3 = 0.0;
																						doublerealT sumPF_3 = 0.0;
																						doublerealT SumPall_3 = 0.0;

																						// ����� 4.1

																						for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																								if (this_is_C_node[Amat.j[is0_3]]) {

																									if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										sumPC_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																										iOPORNAQ_C++;
																									}
																								}
																								else {
																									if (this_is_F_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																											// ������������� Strong F �����.

																											sumPF_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																										}
																									}
																								}
																							}
																						}


																						if (iOPORNAQ_C >= ilimCadd) {
																							SumPall_3 = sumPC_3;

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {

																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																											Ak1 Ptmp;

																											Ptmp.j = static_cast<integer_mix_precision>(i8);
																											Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_3]]);
																											//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																											if (fabs(SumPall_3) < 1.0e-33) {
																												printf("error 1.0.5 ! division by zero. SumPall =%e\n", SumPall_3);
																												//system("pause");
																												system("PAUSE");
																												exit(1);
																											}
																											Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * (fabs(Amat.aij[is0_3]) / SumPall_3));

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
																						else {

																							SumPall_3 = sumPC_3 + sumPF_3;


																							// ����� 4.2

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																											Ak1 Ptmp;

																											Ptmp.j = static_cast<integer_mix_precision>(i8);
																											Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_3]]);

																											if (fabs(SumPall_3) < 1.0e-33) {
																												printf("error 2.0 ! division by zero. SumPall_1 =%e\n", SumPall_3);
																												//system("pause");
																												system("PAUSE");
																												exit(1);
																											}
																											Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * (fabs(Amat.aij[is0_3]) / SumPall_3));

																											if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																												P_loc[itid][isize_P_loc[itid]] = Ptmp;
																												isize_P_loc[itid]++;
																											}
																											else {
																												std::cout << "Error 2.0 Increase size P_loc memory limit\n" << std::endl;
																												system("PAUSE");
																											}
																										}
																									}
																									else {
																										if (this_is_F_node[Amat.j[is0_3]]) {

																											if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																												// ������������� Strong F �����.

																												integer i8_4 = Amat.j[is0_3];


																												doublerealT nu4 = fabs(Amat.aij[is0_3]) / SumPall_3;

																												// ����� ���������� ����� ��� ���� i8_4.

										                                                                        //********************************************************

																												doublerealT maxelem_threshold_4 = -1.0;

																												integer ii1_4 = row_startA[i8_4];
																												integer istr_etalon1_4 = Amat.i[ii1_4];
																												integer iend_for1_4 = -1;
																												for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																													iend_for1_4 = is0_4;
																												}
																												if (old_def_threshold) {
																													for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1_4) && (Amat.aij[is0_4] < 0.0) && (-Amat.aij[is0_4] > maxelem_threshold_4)) maxelem_threshold_4 = -Amat.aij[is0_4];
																													}
																												}
																												else {
																													maxelem_threshold_4 = threshold_quick_all[istr_etalon1_4];
																												}

																												// ����� maxelem_threshold_4 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																												doublerealT maxelem_threshold_theta_4 = maxelem_threshold_4 * theta;
																												doublerealT maxelem_threshold_theta_strong_F_4 = maxelem_threshold_4 * magic;


																												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																												// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																												doublerealT sumPC_4 = 0.0;
																												doublerealT sumPF_4 = 0.0;
																												doublerealT SumPall_4 = 0.0;

																												// ����� 5.1

																												for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																													if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_4)) {
																														if (this_is_C_node[Amat.j[is0_4]]) {

																															if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																sumPC_4 += fabs(Amat.aij[is0_4]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																																iOPORNAQ_C++;
																															}
																														}
																														else {
																															if (this_is_F_node[Amat.j[is0_4]]) {

																																if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_strong_F_4) {

																																	// ������������� Strong F �����.

																																	sumPF_4 += fabs(Amat.aij[is0_4]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																																}
																															}
																														}
																													}
																												}


																												if (iOPORNAQ_C >= ilimCadd) {
																													SumPall_4 = sumPC_4;

																													for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																															if (this_is_C_node[Amat.j[is0_4]]) {

																																if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																	Ak1 Ptmp;

																																	Ptmp.j = static_cast<integer_mix_precision>(i8);
																																	Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_4]]);
																																	//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																																	if (fabs(SumPall_4) < 1.0e-33) {
																																		printf("error 1.0.3 ! division by zero. SumPall =%e\n", SumPall_4);
																																		//system("pause");
																																		system("PAUSE");
																																		exit(1);
																																	}
																																	Ptmp.aij = static_cast<real_mix_precision>(nu1*nu2*nu3*nu4 * (fabs(Amat.aij[is0_4]) / SumPall_4));

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
																												else {

																													SumPall_4 = sumPC_4 + sumPF_4;


																													// ����� 5.2

																													for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																														if (((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4))) {
																															if (this_is_C_node[Amat.j[is0_4]]) {

																																if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																	Ak1 Ptmp;

																																	Ptmp.j = static_cast<integer_mix_precision>(i8);
																																	Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_4]]);

																																	if (fabs(SumPall_4) < 1.0e-33) {
																																		printf("error 2.0 ! division by zero. SumPall_1 =%e\n", SumPall_4);
																																		//system("pause");
																																		system("PAUSE");
																																		exit(1);
																																	}
																																	Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * nu4 * (fabs(Amat.aij[is0_4]) / SumPall_4));

																																	if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																																		P_loc[itid][isize_P_loc[itid]] = Ptmp;
																																		isize_P_loc[itid]++;
																																	}
																																	else {
																																		std::cout << "Error 2.0 Increase size P_loc memory limit\n" << std::endl;
																																		system("PAUSE");
																																	}
																																}
																															}
																															else {
																																if (this_is_F_node[Amat.j[is0_4]]) {

																																	if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_strong_F_4) {

																																		// ������������� Strong F �����.

																																		integer i8_5 = Amat.j[is0_4];


																																		doublerealT nu5 = fabs(Amat.aij[is0_4]) / SumPall_4;

																																		// ����� ���������� ����� ��� ���� i8_5.




																																		 //********************************************************

																																		doublerealT maxelem_threshold_5 = -1.0;

																																		integer ii1_5 = row_startA[i8_5];
																																		integer istr_etalon1_5 = Amat.i[ii1_5];
																																		integer iend_for1_5 = -1;
																																		for (integer is0_5 = ii1_5; (is0_5 <= row_startA[istr_etalon1_5 + 1] - 1); ++is0_5) {
																																			iend_for1_5 = is0_5;
																																		}
																																		if (old_def_threshold) {
																																			for (integer is0_5 = ii1_5; (is0_5 <= row_startA[istr_etalon1_5 + 1] - 1); ++is0_5) {
																																			if ((Amat.j[is0_5] != istr_etalon1_5) && (Amat.aij[is0_5] < 0.0) && (-Amat.aij[is0_5] > maxelem_threshold_5)) maxelem_threshold_5 = -Amat.aij[is0_5];
																																		}
																																	}
																																	else {
																																		maxelem_threshold_5 = threshold_quick_all[istr_etalon1_5];
																																	}

																																	// ����� maxelem_threshold_5 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																																	doublerealT maxelem_threshold_theta_5 = maxelem_threshold_5 * theta;
																																	//doublerealT maxelem_threshold_theta_strong_F_5 = maxelem_threshold_5 * magic_FFF;


																																	// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																																	// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																																	doublerealT sumPC_5 = 0.0;
																																	//doublerealT sumPF_5 = 0.0;
																																	doublerealT SumPall_5 = 0.0;

																																	// ����� 6.1

																																	for (integer is0_5 = ii1_5; is0_5 <= iend_for1_5; ++is0_5) {
																																		if ((Amat.j[is0_5] != istr_etalon1) && (Amat.j[is0_5] != istr_etalon1_1) && (Amat.j[is0_5] != istr_etalon1_2) && (Amat.j[is0_5] != istr_etalon1_3) && (Amat.j[is0_5] != istr_etalon1_4) && (Amat.j[is0_5] != istr_etalon1_5)) {
																																			if (this_is_C_node[Amat.j[is0_5]]) {

																																				if (fabs(Amat.aij[is0_5]) > maxelem_threshold_theta_5) {
																																					sumPC_5 += fabs(Amat.aij[is0_5]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																																					iOPORNAQ_C++;
																																				}
																																			}

																																		}
																																	}


																																	SumPall_5 = sumPC_5;


																																	// ����� 6.2

																																	for (integer is0_5 = ii1_5; is0_5 <= iend_for1_5; ++is0_5) {
																																	if ((Amat.j[is0_5] != istr_etalon1) && (Amat.j[is0_5] != istr_etalon1_1) && (Amat.j[is0_5] != istr_etalon1_2) && (Amat.j[is0_5] != istr_etalon1_3) && (Amat.j[is0_5] != istr_etalon1_4) && (Amat.j[is0_5] != istr_etalon1_5)) {
																																		if (this_is_C_node[Amat.j[is0_5]]) {

																																			if (fabs(Amat.aij[is0_5]) > maxelem_threshold_theta_5) {
																																				Ak1 Ptmp;

																																				Ptmp.j = static_cast<integer_mix_precision>(i8);
																																				Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_5]]);

																																				if (fabs(SumPall_5) < 1.0e-33) {
																																					printf("error 4.0 ! division by zero. SumPall_4 =%e\n", SumPall_5);
																																					//system("pause");
																																					system("PAUSE");
																																					exit(1);
																																				}
																																				Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * nu4 * nu5 * (fabs(Amat.aij[is0_5]) / SumPall_5));

																																				if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																																					P_loc[itid][isize_P_loc[itid]] = Ptmp;
																																					isize_P_loc[itid]++;
																																				}
																																				else {
																																					std::cout << "Error 4.0 Increase size P_loc memory limit\n" << std::endl;
																																					system("PAUSE");
																																				}
																																			}
																																		}

																																	}
																																}

																																//********************************************************

																																// ��������� ������ ���������� i8_5.
																																//********************************************************

																																}
																															}
																														}
																													}
																												}
																											}
																											//********************************************************

																											// ��������� ������ ���������� i8_4.
																											//********************************************************

																											}
																										}
																									}
																								}
																							}
																						}

																						//*******************************************18.3

																					// ��������� ������ ���������� i8_3.

																					}
																				}
																			}
																		}
																	}


																}
																//********************************************************

																// ��������� ������ ���������� i8_2.



															}
														}
													}
												}
											}


										}
										//********************************************************
										// ��������� ������ ���������� i8_1.

									}
								}
							}
						}
					}

				}

				if (iOPORNAQ_C < 4) {
					switch (iOPORNAQ_C) {
					case 0: iOPORNAQ_C0++;
						break;
					case 1: iOPORNAQ_C1++;
						break;
					case 2: iOPORNAQ_C2++;
						break;
					case 3: iOPORNAQ_C3++;
						break;
					case 4: iOPORNAQ_C4++;
						break;
					}
				}

				if (iOPORNAQ_C == 1) {
					//std::cout << "one uzel : " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
					P_loc[itid][isize_P_loc[itid] - 1].aij = 1.0;
				}
				if (iOPORNAQ_C == 2) {
					//std::cout << "two uzla : " << P_loc[itid][isize_P_loc[itid] - 2].aij << " " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
					doublerealT s1 = P_loc[itid][isize_P_loc[itid] - 2].aij + P_loc[itid][isize_P_loc[itid] - 1].aij;
					P_loc[itid][isize_P_loc[itid] - 2].aij /= s1;
					P_loc[itid][isize_P_loc[itid] - 1].aij /= s1;
				}
				if (iOPORNAQ_C == 3) {
					//std::cout << "three uzla : " << P_loc[itid][isize_P_loc[itid] - 3].aij << " " << P_loc[itid][isize_P_loc[itid] - 2].aij << " " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
					doublerealT s1 = P_loc[itid][isize_P_loc[itid] - 3].aij + P_loc[itid][isize_P_loc[itid] - 2].aij + P_loc[itid][isize_P_loc[itid] - 1].aij;
					P_loc[itid][isize_P_loc[itid] - 3].aij /= s1;
					P_loc[itid][isize_P_loc[itid] - 2].aij /= s1;
					P_loc[itid][isize_P_loc[itid] - 1].aij /= s1;
				}
				if (iOPORNAQ_C >= 4) {
					doublerealT s1 = 0.0;
					for (int i_5 = 1; i_5 <= iOPORNAQ_C; ++i_5) s1 += P_loc[itid][isize_P_loc[itid] - i_5].aij;
					for (int i_5 = 1; i_5 <= iOPORNAQ_C; ++i_5) P_loc[itid][isize_P_loc[itid] - i_5].aij /= s1;
				}

				//system("pause");
			}
		}

		if (my_amg_manager.iprint_log) {
			std::cout << "interpolatioon statistics: \n";
			std::cout << iOPORNAQ_C0 << " " << iOPORNAQ_C1 << " " << iOPORNAQ_C2 << " " << iOPORNAQ_C3 << " " << iOPORNAQ_C4 << " " << iadds0 << " " << std::endl;
			//system("pause");
		}

		// �������� ����������� ����������.
		for (int itid = 0; itid < inth; ++itid)
		{

			if (icount1 + isize_P_loc[itid] >= nsizePR * n) {
				printf("memory error!!!\n");
				printf("not enough memory for the interpolation operator.\n");
				//system("PAUSE");
				//exit(1);
				deallocate_prolongation(nsizePR, n, P);
			}

#pragma omp parallel for
			for (integer iscan_loc = 0; iscan_loc < isize_P_loc[itid]; ++iscan_loc) {
				P[icount1 + iscan_loc] = P_loc[itid][iscan_loc];
			}
			icount1 += isize_P_loc[itid];

		}


		for (int itid = 0; itid < inth; ++itid)
		{
			isize_P_loc[itid] = 0;

			delete[] P_loc[itid];
		}

		delete[] isize_P_loc;
		delete[] P_loc;


	}




















}


// ���� F ����� ������.
// 25.09.2021 FFFFC ������
// CFFFFFC
// ���������� ������� ���� �������� ��� PMIS2 (PMIS ������������ � �������� �������) �� ���� ������� ��������� �����.
template <typename doublerealT>
void my_long_rang_interpolation_procedure_FFFC_number3A_parallelx8(integer& the_number_of_neighbors_that_are_not_C_nodes_gl,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor_gl,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic_gl,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF_gl,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative, bool bsuffix_work, bool baggressive_interpol)
{

#ifdef _OPENMP
	int inth = number_cores();
	if (inth > 1) {
		omp_set_num_threads(inth);
	}
	//omp_set_num_threads(8);
#else 
	// ���� �����.
	int inth = 1;
#endif

	bool old_def_threshold = true;

	static int ivprjsk = 0;
	++ivprjsk;


	const int ilimCadd = 4;

	// ������
	int iadditionalCstatistic = 0;
	int iadds0 = 0;

	int iOPORNAQ_C0 = 0;
	int iOPORNAQ_C1 = 0;
	int iOPORNAQ_C2 = 0;
	int iOPORNAQ_C3 = 0;
	int iOPORNAQ_C4 = 0;

	if (ivprjsk % 2 != 0) {


#pragma omp parallel for reduction(+: iadditionalCstatistic, iadds0, iOPORNAQ_C0, iOPORNAQ_C1, iOPORNAQ_C2, iOPORNAQ_C3, iOPORNAQ_C4)
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) if (this_is_F_node[i8]) {

			bool bproblem = false;

			int iOPORNAQ_C0_1 = 0;
			int iOPORNAQ_C1_1 = 0;
			int iOPORNAQ_C2_1 = 0;
			int iOPORNAQ_C3_1 = 0;
			int iOPORNAQ_C4_1 = 0;

			int i_na_odnoi_vetke_C2 = 0;

			integer ii1 = row_startA[i8];
			integer istr_etalon1 = Amat.i[ii1];
			integer iend_for1 = -1;
			for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
				iend_for1 = is0;
			}

			if (iend_for1 != ii1) {
				// ������ ����.
				doublerealT theta_strong_F = theta83; // 3 ���� 2016
				doublerealT theta = theta_gl;

				const doublerealT magic = magic82; // 0.4 is recomended.
				const doublerealT magic_FF = magic82; // 18.09.2020
				//const doublerealT magic_FFF = magic82; // 24.09.2021

				int iOPORNAQ_C = 0;

				// ��� ����� ����������� �������� ��������.
							// 5 ������� 2015 ���� �� ��������� ��������� �������������
							// ��������� ������������ � ������ � ��������� ��������.

				

				doublerealT maxelem_threshold = -1.0;

				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
					iend_for1 = is0;
				}
				if (old_def_threshold) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
						if ((Amat.j[is0] != istr_etalon1) && (Amat.aij[is0] < 0.0) && (-Amat.aij[is0] > maxelem_threshold)) maxelem_threshold = -Amat.aij[is0];
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}

				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumPC = 0.0;
				doublerealT sumPF = 0.0;
				doublerealT SumPall = 0.0;

				// ����� 1.1

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumPC += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
								iOPORNAQ_C++;
								iOPORNAQ_C0++;
								iOPORNAQ_C0_1++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

									// ������������� Strong F �����.

									sumPF += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
								}
							}
						}
					}
				}

				if (iOPORNAQ_C >= ilimCadd) {
					SumPall = sumPC;

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									
								}
							}
						}
					}

				}
				else {

					SumPall = sumPC + sumPF;


					// ����� 1.2

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {

									if (i_na_odnoi_vetke_C2 != 2) {
										i_na_odnoi_vetke_C2 = 0;
									}

									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

										// ������������� Strong F �����.

										integer i8_1 = Amat.j[is0];


										doublerealT nu1 = fabs(Amat.aij[is0]) / SumPall;

										// ����� ���������� ����� ��� ���� i8_1.

										//********************************************************

										doublerealT maxelem_threshold_1 = -1.0;

										integer ii1_1 = row_startA[i8_1];
										integer istr_etalon1_1 = Amat.i[ii1_1];
										integer iend_for1_1 = -1;
										for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
											iend_for1_1 = is0_1;
										}
										if (old_def_threshold) {
											for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1_1) && (Amat.aij[is0_1] < 0.0) && (-Amat.aij[is0_1] > maxelem_threshold_1)) maxelem_threshold_1 = -Amat.aij[is0_1];
											}
										}
										else {
											maxelem_threshold_1 = threshold_quick_all[istr_etalon1_1];
										}

										// ����� maxelem_threshold_1 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



										doublerealT maxelem_threshold_theta_1 = maxelem_threshold_1 * theta;
										doublerealT maxelem_threshold_theta_strong_F_1 = maxelem_threshold_1 * magic;


										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
										doublerealT sumPC_1 = 0.0;
										doublerealT sumPF_1 = 0.0;
										doublerealT SumPall_1 = 0.0;

										// ����� 2.1

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														sumPC_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
														iOPORNAQ_C++;
														iOPORNAQ_C1++;
														iOPORNAQ_C1_1++;
														i_na_odnoi_vetke_C2++;
													}
												}
												else {
													if (this_is_F_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

															// ������������� Strong F �����.

															sumPF_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
														}
													}
												}
											}
										}


										if (iOPORNAQ_C >= ilimCadd) {
											SumPall_1 = sumPC_1;

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															
														}
													}
												}
											}

										}
										else {

											SumPall_1 = sumPC_1 + sumPF_1;


											// ����� 2.2

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														
														}
													}
													else {
														if (this_is_F_node[Amat.j[is0_1]]) {

															if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

																// ������������� Strong F �����.

																integer i8_2 = Amat.j[is0_1];


																doublerealT nu2 = fabs(Amat.aij[is0_1]) / SumPall_1;

																// ����� ���������� ����� ��� ���� i8_2.


																//********************************************************

																doublerealT maxelem_threshold_2 = -1.0;

																integer ii1_2 = row_startA[i8_2];
																integer istr_etalon1_2 = Amat.i[ii1_2];
																integer iend_for1_2 = -1;
																for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																	iend_for1_2 = is0_2;
																}
																if (old_def_threshold) {
																	for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1_2) && (Amat.aij[is0_2] < 0.0) && (-Amat.aij[is0_2] > maxelem_threshold_2)) maxelem_threshold_2 = -Amat.aij[is0_2];
																	}
																}
																else {
																	maxelem_threshold_2 = threshold_quick_all[istr_etalon1_2];
																}

																// ����� maxelem_threshold_2 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																doublerealT maxelem_threshold_theta_2 = maxelem_threshold_2 * theta;
																doublerealT maxelem_threshold_theta_strong_F_2 = maxelem_threshold_2 * magic_FF;


																// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																doublerealT sumPC_2 = 0.0;
																doublerealT sumPF_2 = 0.0;
																doublerealT SumPall_2 = 0.0;

																// ����� 3.1

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				sumPC_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
																				iOPORNAQ_C++;
																				iOPORNAQ_C2++;
																				iOPORNAQ_C2_1++;
																				i_na_odnoi_vetke_C2++;
																			}
																		}
																		else {
																			if (this_is_F_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																					// ������������� Strong F �����.

																					sumPF_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																				}
																			}
																		}
																	}
																}

																if (iOPORNAQ_C >= ilimCadd) {
																	SumPall_2 = sumPC_2;

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					
																				}
																			}
																		}
																	}

																}
																else {

																	SumPall_2 = sumPC_2 + sumPF_2;


																	// ����� 3.2

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					
																				}
																			}
																			else {
																				if (this_is_F_node[Amat.j[is0_2]]) {

																					if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																						// ������������� Strong F �����.

																						integer i8_3 = Amat.j[is0_2];


																						doublerealT nu3 = fabs(Amat.aij[is0_2]) / SumPall_2;


																						///******************************************************************
																						// ����� ���������� ����� ��� ���� i8_3.

										                                                //********************************************************

																						doublerealT maxelem_threshold_3 = -1.0;

																						integer ii1_3 = row_startA[i8_3];
																						integer istr_etalon1_3 = Amat.i[ii1_3];
																						integer iend_for1_3 = -1;
																						for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																							iend_for1_3 = is0_3;
																						}
																						if (old_def_threshold) {
																							for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1_3) && (Amat.aij[is0_3] < 0.0) && (-Amat.aij[is0_3] > maxelem_threshold_3)) maxelem_threshold_3 = -Amat.aij[is0_3];
																							}
																						}
																						else {
																							maxelem_threshold_3 = threshold_quick_all[istr_etalon1_3];
																						}

																						// ����� maxelem_threshold_3 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																						doublerealT maxelem_threshold_theta_3 = maxelem_threshold_3 * theta;
																						doublerealT maxelem_threshold_theta_strong_F_3 = maxelem_threshold_3 * magic;


																						// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																						// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																						doublerealT sumPC_3 = 0.0;
																						doublerealT sumPF_3 = 0.0;
																						doublerealT SumPall_3 = 0.0;

																						// ����� 4.1

																						for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																								if (this_is_C_node[Amat.j[is0_3]]) {

																									if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										sumPC_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																										iOPORNAQ_C++;
																										iOPORNAQ_C3++;
																										iOPORNAQ_C3_1++;
																										i_na_odnoi_vetke_C2++;
																									}
																								}
																								else {
																									if (this_is_F_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																											// ������������� Strong F �����.

																											sumPF_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																										}
																									}
																								}
																							}
																						}


																						if (iOPORNAQ_C >= ilimCadd) {
																							SumPall_3 = sumPC_3;

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {

																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										
																										}
																									}
																								}
																							}

																						}
																						else {

																							SumPall_3 = sumPC_3 + sumPF_3;


																							// ����� 4.2

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																											
																										}
																									}
																									else {
																										if (this_is_F_node[Amat.j[is0_3]]) {

																											if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																												// ������������� Strong F �����.

																												integer i8_4 = Amat.j[is0_3];


																												doublerealT nu4 = fabs(Amat.aij[is0_3]) / SumPall_3;

																												// ����� ���������� ����� ��� ���� i8_4.




																						                        //********************************************************

																												doublerealT maxelem_threshold_4 = -1.0;

																												integer ii1_4 = row_startA[i8_4];
																												integer istr_etalon1_4 = Amat.i[ii1_4];
																												integer iend_for1_4 = -1;
																												for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																													iend_for1_4 = is0_4;
																												}
																												if (old_def_threshold) {
																													for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1_4) && (Amat.aij[is0_4] < 0.0) && (-Amat.aij[is0_4] > maxelem_threshold_4)) maxelem_threshold_4 = -Amat.aij[is0_4];
																													}
																												}
																												else {
																													maxelem_threshold_4 = threshold_quick_all[istr_etalon1_4];
																												}

																												// ����� maxelem_threshold_4 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																												doublerealT maxelem_threshold_theta_4 = maxelem_threshold_4 * theta;
																												//doublerealT maxelem_threshold_theta_strong_F_4 = maxelem_threshold_4 * magic_FFF;


																												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																												// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																												doublerealT sumPC_4 = 0.0;
																												//doublerealT sumPF_4 = 0.0;
																												//doublerealT SumPall_4 = 0.0;

																												// ����� 5.1

																												for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																													if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																														if (this_is_C_node[Amat.j[is0_4]]) {

																															if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																sumPC_4 += fabs(Amat.aij[is0_4]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																																iOPORNAQ_C++;
																																iOPORNAQ_C4++;
																																iOPORNAQ_C4_1++;
																																i_na_odnoi_vetke_C2++;
																															}
																														}

																													}
																												}


																												//SumPall_4 = sumPC_4;


																												// ����� 5.2

																												for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																													if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																														if (this_is_C_node[Amat.j[is0_4]]) {

																															if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																
																															}
																														}

																													}
																												}

																												//********************************************************

																												// ��������� ������ ���������� i8_4.
																												//********************************************************

																											}
																										}
																									}
																								}
																							}
																						}

																						//*******************************************18.3

																					// ��������� ������ ���������� i8_3.

																					}
																				}
																			}
																		}
																	}


																}
																//********************************************************

																// ��������� ������ ���������� i8_2.



															}
														}
													}
												}
											}


										}
										//********************************************************
										// ��������� ������ ���������� i8_1.



									}
								}
							}
						}
					}

				}
			
				// ��� ����� ��������� �������� ����� �������� ������� �����.
				// ����������� ����� �������� ����� � ������������ � ����� ����������.
				if ((iOPORNAQ_C == 0) ||
					((iOPORNAQ_C == 1) && (iOPORNAQ_C0_1 == 0) && (iOPORNAQ_C1_1 == 0)) ||
					((iOPORNAQ_C == 2) && (i_na_odnoi_vetke_C2==2) && (iOPORNAQ_C0_1 == 0) && (iOPORNAQ_C1_1 == 0))
					) {

					    bproblem = true;
					
				    }
					// ���� ��� � ����� ������, ���� ���� ������ ���� � �����, �� �� ���������� ������ �� ���������� 2 � ����.
				    // ���� ���� ������ ��� ���� ������ �� ����� ����� � ������ �������.
                }
			
				if (ii1 == (row_startA[istr_etalon1 + 1] - 1)) bproblem = false; // ������� ������� �� �� ��������� � � ����. 01.10.2021
				if (ii1 + 1 == (row_startA[istr_etalon1 + 1] - 1)) bproblem = false; // ���������� ������� ������� 01.10.2021
				//if (ii1 + 2 <= (row_startA[istr_etalon1 + 1] - 1)) {
					// ���� ���� ��� ������
					//if (iOPORNAQ_C == 0) {
						//� ����� ������� ���� �������� � ����.
					//    bproblem = false;
					//}
			//	}


			if (bproblem)
			{

				//int ir = rand() % 1000;
				//if (ir < 1000) {

					// ������ 10%.

					this_is_C_node[istr_etalon1] = true;
					this_is_F_node[istr_etalon1] = false;

					iadditionalCstatistic++;
				//}
			}

		}

		iadditionalCstatistic_gl += iadditionalCstatistic;
		iadds0 = iadditionalCstatistic;
		iadditionalCstatistic = 0;

		if (my_amg_manager.iprint_log) {

			std::cout << "iadditional C=" << iadds0 << std::endl;

			std::cout << iOPORNAQ_C0 << " " << iOPORNAQ_C1 << " " << iOPORNAQ_C2 << " " << iOPORNAQ_C3 << " " << iOPORNAQ_C4 << " " << std::endl;
		}

		bweSholdbeContinue = true;

	}
	else {
		bweSholdbeContinue = false;
	}

	//if (iadditionalCstatistic == 0) 
	{

		doublerealT theta_strong_F = theta83; // 3 ���� 2016
		doublerealT theta = theta_gl;

		const doublerealT magic = magic82; // 0.4 is recomended.
		const doublerealT magic_FF = magic82; // 18.09.2020
		//const doublerealT magic_FFF = magic82; // 24.09.2021


		//int itid = 0;

		integer LOC_P_LIM_SIZE = -1;
		double dsize =  2.2 *( 1.0/ (1.0*inth)) * nsizePR * n_a[0];
		LOC_P_LIM_SIZE = static_cast<integer>(dsize);

		Ak1** P_loc = new Ak1 * [inth];
		integer* isize_P_loc = new integer[inth];
		for (int itid = 0; itid < inth; ++itid)
		{
			isize_P_loc[itid] = 0;

			P_loc[itid] = new Ak1[LOC_P_LIM_SIZE];
		}

		
		//int istat[5] = { 0,0,0,0,0 };
		//istat[4] = iadds0;

		iOPORNAQ_C0 = 0;
		iOPORNAQ_C1 = 0;
		iOPORNAQ_C2 = 0;
		iOPORNAQ_C3 = 0;
		iOPORNAQ_C4 = 0;

#pragma omp parallel for reduction(+:  iOPORNAQ_C0, iOPORNAQ_C1, iOPORNAQ_C2, iOPORNAQ_C3, iOPORNAQ_C4)
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8)
		{
#ifdef _OPENMP
				int itid = omp_get_thread_num();
#else
			int itid = 0;
#endif

			if (this_is_F_node[i8]) {


				int iOPORNAQ_C = 0;

				// ��� ����� ����������� �������� ��������.
							// 5 ������� 2015 ���� �� ��������� ��������� �������������
							// ��������� ������������ � ������ � ��������� ��������.



				doublerealT maxelem_threshold = -1.0;

				integer ii1 = row_startA[i8];
				integer istr_etalon1 = Amat.i[ii1];
				integer iend_for1 = -1;
				for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
					iend_for1 = is0;
				}
				if (old_def_threshold) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
						if ((Amat.j[is0] != istr_etalon1) && (Amat.aij[is0] < 0.0) && (-Amat.aij[is0] > maxelem_threshold)) maxelem_threshold = -Amat.aij[is0];
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}

				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumPC = 0.0;
				doublerealT sumPF = 0.0;
				doublerealT SumPall = 0.0;

				// ����� 1.1

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumPC += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
								iOPORNAQ_C++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

									// ������������� Strong F �����.

									sumPF += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
								}
							}
						}
					}
				}

				if (iOPORNAQ_C >= ilimCadd) {
					SumPall = sumPC;

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									Ak1 Ptmp;

									Ptmp.j = static_cast<integer_mix_precision>(i8);
									Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
									//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
									if (fabs(SumPall) < 1.0e-33) {
										printf("error 1.0.1 ! division by zero. SumPall =%e\n", SumPall);
										//system("pause");
										system("PAUSE");
										exit(1);
									}
									Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
				else {

					SumPall = sumPC + sumPF;


					// ����� 1.2

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									Ak1 Ptmp;

									Ptmp.j = static_cast<integer_mix_precision>(i8);
									Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
									//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
									if (fabs(SumPall) < 1.0e-33) {
										printf("error 1.0.2 ! division by zero. SumPall =%e\n", SumPall);
										//system("pause");
										system("PAUSE");
										exit(1);
									}
									Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
							else {
								if (this_is_F_node[Amat.j[is0]]) {

									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

										// ������������� Strong F �����.

										integer i8_1 = Amat.j[is0];


										doublerealT nu1 = fabs(Amat.aij[is0]) / SumPall;

										// ����� ���������� ����� ��� ���� i8_1.

										//********************************************************

										doublerealT maxelem_threshold_1 = -1.0;

										integer ii1_1 = row_startA[i8_1];
										integer istr_etalon1_1 = Amat.i[ii1_1];
										integer iend_for1_1 = -1;
										for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
											iend_for1_1 = is0_1;
										}
										if (old_def_threshold) {
											for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1_1) && (Amat.aij[is0_1] < 0.0) && (-Amat.aij[is0_1] > maxelem_threshold_1)) maxelem_threshold_1 = -Amat.aij[is0_1];
											}
										}
										else {
											maxelem_threshold_1 = threshold_quick_all[istr_etalon1_1];
										}

										// ����� maxelem_threshold_1 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



										doublerealT maxelem_threshold_theta_1 = maxelem_threshold_1 * theta;
										doublerealT maxelem_threshold_theta_strong_F_1 = maxelem_threshold_1 * magic;


										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
										doublerealT sumPC_1 = 0.0;
										doublerealT sumPF_1 = 0.0;
										doublerealT SumPall_1 = 0.0;

										// ����� 2.1

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														sumPC_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
														iOPORNAQ_C++;
													}
												}
												else {
													if (this_is_F_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

															// ������������� Strong F �����.

															sumPF_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
														}
													}
												}
											}
										}


										if (iOPORNAQ_C >= ilimCadd) {
											SumPall_1 = sumPC_1;

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															Ak1 Ptmp;

															Ptmp.j = static_cast<integer_mix_precision>(i8);
															Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_1]]);
															//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
															if (fabs(SumPall_1) < 1.0e-33) {
																printf("error 1.0.3 ! division by zero. SumPall =%e\n", SumPall_1);
																//system("pause");
																system("PAUSE");
																exit(1);
															}
															Ptmp.aij = static_cast<real_mix_precision>(nu1 * (fabs(Amat.aij[is0_1]) / SumPall_1));

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
										else {

											SumPall_1 = sumPC_1 + sumPF_1;


											// ����� 2.2

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															Ak1 Ptmp;

															Ptmp.j = static_cast<integer_mix_precision>(i8);
															Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_1]]);

															if (fabs(SumPall_1) < 1.0e-33) {
																printf("error 2.0 ! division by zero. SumPall_1 =%e\n", SumPall_1);
																//system("pause");
																system("PAUSE");
																exit(1);
															}
															Ptmp.aij = static_cast<real_mix_precision>(nu1 * (fabs(Amat.aij[is0_1]) / SumPall_1));

															if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																P_loc[itid][isize_P_loc[itid]] = Ptmp;
																isize_P_loc[itid]++;
															}
															else {
																std::cout << "Error 2.0 Increase size P_loc memory limit\n" << std::endl;
																system("PAUSE");
															}
														}
													}
													else {
														if (this_is_F_node[Amat.j[is0_1]]) {

															if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

																// ������������� Strong F �����.

																integer i8_2 = Amat.j[is0_1];


																doublerealT nu2 = fabs(Amat.aij[is0_1]) / SumPall_1;

																// ����� ���������� ����� ��� ���� i8_2.


																//********************************************************

																doublerealT maxelem_threshold_2 = -1.0;

																integer ii1_2 = row_startA[i8_2];
																integer istr_etalon1_2 = Amat.i[ii1_2];
																integer iend_for1_2 = -1;
																for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																	iend_for1_2 = is0_2;
																}
																if (old_def_threshold) {
																	for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1_2) && (Amat.aij[is0_2] < 0.0) && (-Amat.aij[is0_2] > maxelem_threshold_2)) maxelem_threshold_2 = -Amat.aij[is0_2];
																	}
																}
																else {
																	maxelem_threshold_2 = threshold_quick_all[istr_etalon1_2];
																}

																// ����� maxelem_threshold_2 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																doublerealT maxelem_threshold_theta_2 = maxelem_threshold_2 * theta;
																doublerealT maxelem_threshold_theta_strong_F_2 = maxelem_threshold_2 * magic_FF;


																// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																doublerealT sumPC_2 = 0.0;
																doublerealT sumPF_2 = 0.0;
																doublerealT SumPall_2 = 0.0;

																// ����� 3.1

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				sumPC_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
																				iOPORNAQ_C++;
																			}
																		}
																		else {
																			if (this_is_F_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																					// ������������� Strong F �����.

																					sumPF_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																				}
																			}
																		}
																	}
																}

																if (iOPORNAQ_C >= ilimCadd) {
																	SumPall_2 = sumPC_2;

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					Ak1 Ptmp;

																					Ptmp.j = static_cast<integer_mix_precision>(i8);
																					Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_2]]);
																					//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																					if (fabs(SumPall_2) < 1.0e-33) {
																						printf("error 1.0.4 ! division by zero. SumPall =%e\n", SumPall_2);
																						//system("pause");
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * (fabs(Amat.aij[is0_2]) / SumPall_2));

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
																else {

																	SumPall_2 = sumPC_2 + sumPF_2;


																	// ����� 3.2

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					Ak1 Ptmp;

																					Ptmp.j = static_cast<integer_mix_precision>(i8);
																					Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_2]]);

																					if (fabs(SumPall_2) < 1.0e-33) {
																						printf("error 3.0 ! division by zero. SumPall_2 =%e\n", SumPall_2);
																						//system("pause");
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * (fabs(Amat.aij[is0_2]) / SumPall_2));

																					if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																						P_loc[itid][isize_P_loc[itid]] = Ptmp;
																						isize_P_loc[itid]++;
																					}
																					else {
																						std::cout << "Error 3.0 Increase size P_loc memory limit\n" << std::endl;
																						system("PAUSE");
																					}
																				}
																			}
																			else {
																				if (this_is_F_node[Amat.j[is0_2]]) {

																					if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																						// ������������� Strong F �����.

																						integer i8_3 = Amat.j[is0_2];


																						doublerealT nu3 = fabs(Amat.aij[is0_2]) / SumPall_2;


																						///******************************************************************
																						// ����� ���������� ����� ��� ���� i8_3.

										//********************************************************

																						doublerealT maxelem_threshold_3 = -1.0;

																						integer ii1_3 = row_startA[i8_3];
																						integer istr_etalon1_3 = Amat.i[ii1_3];
																						integer iend_for1_3 = -1;
																						for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																							iend_for1_3 = is0_3;
																						}
																						if (old_def_threshold) {
																							for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1_3) && (Amat.aij[is0_3] < 0.0) && (-Amat.aij[is0_3] > maxelem_threshold_3)) maxelem_threshold_3 = -Amat.aij[is0_3];
																							}
																						}
																						else {
																							maxelem_threshold_3 = threshold_quick_all[istr_etalon1_3];
																						}

																						// ����� maxelem_threshold_3 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																						doublerealT maxelem_threshold_theta_3 = maxelem_threshold_3 * theta;
																						doublerealT maxelem_threshold_theta_strong_F_3 = maxelem_threshold_3 * magic;


																						// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																						// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																						doublerealT sumPC_3 = 0.0;
																						doublerealT sumPF_3 = 0.0;
																						doublerealT SumPall_3 = 0.0;

																						// ����� 4.1

																						for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																								if (this_is_C_node[Amat.j[is0_3]]) {

																									if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										sumPC_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																										iOPORNAQ_C++;
																									}
																								}
																								else {
																									if (this_is_F_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																											// ������������� Strong F �����.

																											sumPF_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																										}
																									}
																								}
																							}
																						}


																						if (iOPORNAQ_C >= ilimCadd) {
																							SumPall_3 = sumPC_3;

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {

																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																											Ak1 Ptmp;

																											Ptmp.j = static_cast<integer_mix_precision>(i8);
																											Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_3]]);
																											//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																											if (fabs(SumPall_3) < 1.0e-33) {
																												printf("error 1.0.5 ! division by zero. SumPall =%e\n", SumPall_3);
																												//system("pause");
																												system("PAUSE");
																												exit(1);
																											}
																											Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * (fabs(Amat.aij[is0_3]) / SumPall_3));

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
																						else {

																							SumPall_3 = sumPC_3 + sumPF_3;


																							// ����� 4.2

																							for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																									if (this_is_C_node[Amat.j[is0_3]]) {

																										if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																											Ak1 Ptmp;

																											Ptmp.j = static_cast<integer_mix_precision>(i8);
																											Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_3]]);

																											if (fabs(SumPall_3) < 1.0e-33) {
																												printf("error 2.0 ! division by zero. SumPall_1 =%e\n", SumPall_3);
																												//system("pause");
																												system("PAUSE");
																												exit(1);
																											}
																											Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * (fabs(Amat.aij[is0_3]) / SumPall_3));

																											if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																												P_loc[itid][isize_P_loc[itid]] = Ptmp;
																												isize_P_loc[itid]++;
																											}
																											else {
																												std::cout << "Error 2.0 Increase size P_loc memory limit\n" << std::endl;
																												system("PAUSE");
																											}
																										}
																									}
																									else {
																										if (this_is_F_node[Amat.j[is0_3]]) {

																											if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_strong_F_3) {

																												// ������������� Strong F �����.

																												integer i8_4 = Amat.j[is0_3];


																												doublerealT nu4 = fabs(Amat.aij[is0_3]) / SumPall_3;

																												// ����� ���������� ����� ��� ���� i8_4.




																						//********************************************************

																												doublerealT maxelem_threshold_4 = -1.0;

																												integer ii1_4 = row_startA[i8_4];
																												integer istr_etalon1_4 = Amat.i[ii1_4];
																												integer iend_for1_4 = -1;
																												for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																													iend_for1_4 = is0_4;
																												}
																												if (old_def_threshold) {
																													for (integer is0_4 = ii1_4; (is0_4 <= row_startA[istr_etalon1_4 + 1] - 1); ++is0_4) {
																														if ((Amat.j[is0_4] != istr_etalon1_4) && (Amat.aij[is0_4] < 0.0) && (-Amat.aij[is0_4] > maxelem_threshold_4)) maxelem_threshold_4 = -Amat.aij[is0_4];
																													}
																												}
																												else {
																													maxelem_threshold_4 = threshold_quick_all[istr_etalon1_4];
																												}

																												// ����� maxelem_threshold_4 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																												doublerealT maxelem_threshold_theta_4 = maxelem_threshold_4 * theta;
																												//doublerealT maxelem_threshold_theta_strong_F_4 = maxelem_threshold_4 * magic_FFF;


																												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																												// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																												doublerealT sumPC_4 = 0.0;
																												//doublerealT sumPF_4 = 0.0;
																												doublerealT SumPall_4 = 0.0;

																												// ����� 5.1

																												for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																													if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																														if (this_is_C_node[Amat.j[is0_4]]) {

																															if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																sumPC_4 += fabs(Amat.aij[is0_4]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
																																iOPORNAQ_C++;
																															}
																														}

																													}
																												}


																												SumPall_4 = sumPC_4;


																												// ����� 5.2

																												for (integer is0_4 = ii1_4; is0_4 <= iend_for1_4; ++is0_4) {
																													if ((Amat.j[is0_4] != istr_etalon1) && (Amat.j[is0_4] != istr_etalon1_1) && (Amat.j[is0_4] != istr_etalon1_2) && (Amat.j[is0_4] != istr_etalon1_3) && (Amat.j[is0_4] != istr_etalon1_4)) {
																														if (this_is_C_node[Amat.j[is0_4]]) {

																															if (fabs(Amat.aij[is0_4]) > maxelem_threshold_theta_4) {
																																Ak1 Ptmp;

																																Ptmp.j = static_cast<integer_mix_precision>(i8);
																																Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_4]]);

																																if (fabs(SumPall_4) < 1.0e-33) {
																																	printf("error 4.0 ! division by zero. SumPall_4 =%e\n", SumPall_4);
																																	//system("pause");
																																	system("PAUSE");
																																	exit(1);
																																}
																																Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * nu4 * (fabs(Amat.aij[is0_4]) / SumPall_4));

																																if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																																	P_loc[itid][isize_P_loc[itid]] = Ptmp;
																																	isize_P_loc[itid]++;
																																}
																																else {
																																	std::cout << "Error 4.0 Increase size P_loc memory limit\n" << std::endl;
																																	system("PAUSE");
																																}
																															}
																														}

																													}
																												}

																												//********************************************************

																												// ��������� ������ ���������� i8_4.
																												//********************************************************

																											}
																										}
																									}
																								}
																							}
																						}

																						//*******************************************18.3

																					// ��������� ������ ���������� i8_3.

																					}
																				}
																			}
																		}
																	}


																}
																//********************************************************

																// ��������� ������ ���������� i8_2.



															}
														}
													}
												}
											}


										}
										//********************************************************
										// ��������� ������ ���������� i8_1.



									}
								}
							}
						}
					}

				}

				if (iOPORNAQ_C < 4) {
					switch (iOPORNAQ_C) {
					case 0: iOPORNAQ_C0++;
						break;
					case 1: iOPORNAQ_C1++;
						break;
					case 2: iOPORNAQ_C2++;
						break;
					case 3: iOPORNAQ_C3++;
						break;
					}
				}

				if (iOPORNAQ_C == 1) {
					//std::cout << "one uzel : " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
					P_loc[itid][isize_P_loc[itid] - 1].aij = 1.0;
				}
				if (iOPORNAQ_C == 2) {
					//std::cout << "two uzla : " << P_loc[itid][isize_P_loc[itid] - 2].aij << " " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
					doublerealT s1 = P_loc[itid][isize_P_loc[itid] - 2].aij + P_loc[itid][isize_P_loc[itid] - 1].aij;
					P_loc[itid][isize_P_loc[itid] - 2].aij /= s1;
					P_loc[itid][isize_P_loc[itid] - 1].aij /= s1;
				}
				if (iOPORNAQ_C == 3) {
					//std::cout << "three uzla : " << P_loc[itid][isize_P_loc[itid] - 3].aij << " " << P_loc[itid][isize_P_loc[itid] - 2].aij << " " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
					doublerealT s1 = P_loc[itid][isize_P_loc[itid] - 3].aij + P_loc[itid][isize_P_loc[itid] - 2].aij + P_loc[itid][isize_P_loc[itid] - 1].aij;
					P_loc[itid][isize_P_loc[itid] - 3].aij /= s1;
					P_loc[itid][isize_P_loc[itid] - 2].aij /= s1;
					P_loc[itid][isize_P_loc[itid] - 1].aij /= s1;
				}
				if (iOPORNAQ_C >= 4) {
					doublerealT s1 = 0.0;
					for (int i_5 = 1; i_5 <= iOPORNAQ_C; ++i_5) s1 += P_loc[itid][isize_P_loc[itid] - i_5].aij;
					for (int i_5 = 1; i_5 <= iOPORNAQ_C; ++i_5) P_loc[itid][isize_P_loc[itid] - i_5].aij /= s1;
				}

				//system("pause");
			}
		}

		if (my_amg_manager.iprint_log) {
			std::cout << "interpolatioon statistics: \n";
			std::cout << iOPORNAQ_C0 << " " << iOPORNAQ_C1 << " " << iOPORNAQ_C2 << " " << iOPORNAQ_C3 << " " << iadds0 << " " << std::endl;
			//system("pause");
		}

		// �������� ����������� ����������.
		for (int itid = 0; itid < inth; ++itid)
		{

			if (icount1 + isize_P_loc[itid] >= nsizePR * n) {
				printf("memory error!!!\n");
				printf("not enough memory for the interpolation operator.\n");
				//system("PAUSE");
				//exit(1);
				deallocate_prolongation(nsizePR, n, P);
			}

#pragma omp parallel for
			for (integer iscan_loc = 0; iscan_loc < isize_P_loc[itid]; ++iscan_loc) {
				P[icount1 + iscan_loc] = P_loc[itid][iscan_loc];
			}
			icount1 += isize_P_loc[itid];

		}


		for (int itid = 0; itid < inth; ++itid)
		{
			isize_P_loc[itid] = 0;

			delete[] P_loc[itid];
		}

		delete[] isize_P_loc;
		delete[] P_loc;


	}




















}


// ������ F ���� ������.
// 25.09.2021 FFFC ������
// 03.10.2021 ok
// CFFFFC
// ������������ � �������� long range interpolation.
// � �������� long range interpolation ����������� ����� ���� ��������� ������  CFFFFFC.
template <typename doublerealT>
void my_long_rang_interpolation_procedure_FFFC4_number3A_parallelx8(integer& the_number_of_neighbors_that_are_not_C_nodes_gl,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor_gl,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic_gl,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF_gl,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative, bool bsuffix_work, bool baggressive_interpol)
{

#ifdef _OPENMP
	int inth = number_cores();
	if (inth > 1) {
		omp_set_num_threads(inth);
	}
	//omp_set_num_threads(8);
#else 
	// ���� �����.
	int inth = 1;
#endif

	bool old_def_threshold = true;

	static int ivprjsk = 0;
	++ivprjsk;


	const int ilimCadd = 4;

	// ������
	int iadditionalCstatistic = 0;
	int iadds0 = 0;

	int iOPORNAQ_C0 = 0;
	int iOPORNAQ_C1 = 0;
	int iOPORNAQ_C2 = 0;
	int iOPORNAQ_C3 = 0;
	

	if (ivprjsk % 2 != 0) {

#pragma omp parallel for reduction(+: iadditionalCstatistic, iadds0, iOPORNAQ_C0, iOPORNAQ_C1, iOPORNAQ_C2, iOPORNAQ_C3)
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) if (this_is_F_node[i8]) {

			bool bproblem = false;

			int iOPORNAQ_C0_1 = 0;
			int iOPORNAQ_C1_1 = 0;
			int iOPORNAQ_C2_1 = 0;
			int iOPORNAQ_C3_1 = 0;
			

			int i_na_odnoi_vetke_C2 = 0;


			// ��� ����� ����������� �������� ��������.
						// 5 ������� 2015 ���� �� ��������� ��������� �������������
						// ��������� ������������ � ������ � ��������� ��������.
			doublerealT maxelem_threshold = -1.0;

			integer ii1 = row_startA[i8];
			integer istr_etalon1 = Amat.i[ii1];
			integer iend_for1 = -1;
			for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
				iend_for1 = is0;
			}

			if (iend_for1 != ii1) {
				// ������ ����.

				doublerealT theta_strong_F = theta83; // 3 ���� 2016
				doublerealT theta = theta_gl;

				const doublerealT magic = magic82; // 0.4 is recomended.
				const doublerealT magic_FF = magic82; // 18.09.2020
				//const doublerealT magic_FFF = magic82; // 24.09.2021

				int iOPORNAQ_C = 0;

				if (old_def_threshold) {
					for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
						if ((Amat.j[is0] != istr_etalon1) && (Amat.aij[is0] < 0.0) && (-Amat.aij[is0] > maxelem_threshold)) maxelem_threshold = -Amat.aij[is0];
					}
				}
				else {
					maxelem_threshold = threshold_quick_all[istr_etalon1];
				}

				// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



				doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
				doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


				// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
				// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
				doublerealT sumPC = 0.0;
				doublerealT sumPF = 0.0;
				doublerealT SumPall = 0.0;

				// ����� 1.1

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								sumPC += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
								iOPORNAQ_C++;
								iOPORNAQ_C0++;
								iOPORNAQ_C0_1++;
							}
						}
						else {
							if (this_is_F_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

									// ������������� Strong F �����.

									sumPF += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
								}
							}
						}
					}
				}

				if (iOPORNAQ_C >= ilimCadd) {
					SumPall = sumPC;

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									
								}
							}
						}
					}

				}
				else {

					SumPall = sumPC + sumPF;


					// ����� 1.2

					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {

									if (i_na_odnoi_vetke_C2 != 2) {
										i_na_odnoi_vetke_C2 = 0;
									}


									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

										// ������������� Strong F �����.

										integer i8_1 = Amat.j[is0];


										doublerealT nu1 = fabs(Amat.aij[is0]) / SumPall;

										// ����� ���������� ����� ��� ���� i8_1.

										//********************************************************

										doublerealT maxelem_threshold_1 = -1.0;

										integer ii1_1 = row_startA[i8_1];
										integer istr_etalon1_1 = Amat.i[ii1_1];
										integer iend_for1_1 = -1;
										for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
											iend_for1_1 = is0_1;
										}
										if (old_def_threshold) {
											for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1_1) && (Amat.aij[is0_1] < 0.0) && (-Amat.aij[is0_1] > maxelem_threshold_1)) maxelem_threshold_1 = -Amat.aij[is0_1];
											}
										}
										else {
											maxelem_threshold_1 = threshold_quick_all[istr_etalon1_1];
										}

										// ����� maxelem_threshold_1 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



										doublerealT maxelem_threshold_theta_1 = maxelem_threshold_1 * theta;
										doublerealT maxelem_threshold_theta_strong_F_1 = maxelem_threshold_1 * magic;


										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
										doublerealT sumPC_1 = 0.0;
										doublerealT sumPF_1 = 0.0;
										doublerealT SumPall_1 = 0.0;

										// ����� 2.1

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														sumPC_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
														iOPORNAQ_C++;
														iOPORNAQ_C1++;
														iOPORNAQ_C1_1++;
														i_na_odnoi_vetke_C2++;
													}
												}
												else {
													if (this_is_F_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

															// ������������� Strong F �����.

															sumPF_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
														}
													}
												}
											}
										}


										if (iOPORNAQ_C >= ilimCadd) {
											SumPall_1 = sumPC_1;

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if (Amat.j[is0_1] != istr_etalon1_1) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															
														}
													}
												}
											}

										}
										else {

											SumPall_1 = sumPC_1 + sumPF_1;


											// ����� 2.2

											for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
												if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
													if (this_is_C_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
															
														}
													}
													else {
														if (this_is_F_node[Amat.j[is0_1]]) {

															if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

																// ������������� Strong F �����.

																integer i8_2 = Amat.j[is0_1];


																doublerealT nu2 = fabs(Amat.aij[is0_1]) / SumPall_1;

																// ����� ���������� ����� ��� ���� i8_2.


																//********************************************************

																doublerealT maxelem_threshold_2 = -1.0;

																integer ii1_2 = row_startA[i8_2];
																integer istr_etalon1_2 = Amat.i[ii1_2];
																integer iend_for1_2 = -1;
																for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																	iend_for1_2 = is0_2;
																}
																if (old_def_threshold) {
																	for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1_2) && (Amat.aij[is0_2] < 0.0) && (-Amat.aij[is0_2] > maxelem_threshold_2)) maxelem_threshold_2 = -Amat.aij[is0_2];
																	}
																}
																else {
																	maxelem_threshold_2 = threshold_quick_all[istr_etalon1_2];
																}

																// ����� maxelem_threshold_2 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																doublerealT maxelem_threshold_theta_2 = maxelem_threshold_2 * theta;
																doublerealT maxelem_threshold_theta_strong_F_2 = maxelem_threshold_2 * magic_FF;


																// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																doublerealT sumPC_2 = 0.0;
																doublerealT sumPF_2 = 0.0;
																doublerealT SumPall_2 = 0.0;

																// ����� 3.1

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				sumPC_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
																				iOPORNAQ_C++;
																				iOPORNAQ_C2++;
																				iOPORNAQ_C2_1++;
																				i_na_odnoi_vetke_C2++;
																			}
																		}
																		else {
																			if (this_is_F_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																					// ������������� Strong F �����.

																					sumPF_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																				}
																			}
																		}
																	}
																}

																if (iOPORNAQ_C >= ilimCadd) {
																	SumPall_2 = sumPC_2;

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if (Amat.j[is0_2] != istr_etalon1_2) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					
																				}
																			}
																		}
																	}

																}
																else {

																	SumPall_2 = sumPC_2 + sumPF_2;


																	// ����� 3.2

																	for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																		if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																			if (this_is_C_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																					
																				}
																			}
																			else {
																				if (this_is_F_node[Amat.j[is0_2]]) {

																					if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																						// ������������� Strong F �����.

																						integer i8_3 = Amat.j[is0_2];


																						doublerealT nu3 = fabs(Amat.aij[is0_2]) / SumPall_2;

																						// ����� ���������� ����� ��� ���� i8_3.


																						//********************************************************

																						doublerealT maxelem_threshold_3 = -1.0;

																						integer ii1_3 = row_startA[i8_3];
																						integer istr_etalon1_3 = Amat.i[ii1_3];
																						integer iend_for1_3 = -1;
																						for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																							iend_for1_3 = is0_3;
																						}
																						if (old_def_threshold) {
																							for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																								if ((Amat.j[is0_3] != istr_etalon1_3) && (Amat.aij[is0_3] < 0.0) && (-Amat.aij[is0_3] > maxelem_threshold_3)) maxelem_threshold_3 = -Amat.aij[is0_3];
																							}
																						}
																						else {
																							maxelem_threshold_3 = threshold_quick_all[istr_etalon1_3];
																						}

																						// ����� maxelem_threshold_3 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																						doublerealT maxelem_threshold_theta_3 = maxelem_threshold_3 * theta;
																						//doublerealT maxelem_threshold_theta_strong_F_3 = maxelem_threshold_3 * magic_FFF;


																						// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																						// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																						doublerealT sumPC_3 = 0.0;
																						//doublerealT sumPF_3 = 0.0;
																						//doublerealT SumPall_3 = 0.0;

																						// ����� 4.1

																						for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																								if (this_is_C_node[Amat.j[is0_3]]) {

																									if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										sumPC_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.		
																										iOPORNAQ_C++;
																										iOPORNAQ_C3++;
																										iOPORNAQ_C3_1++;
																										i_na_odnoi_vetke_C2++;
																									}
																								}

																							}
																						}


																						//SumPall_3 = sumPC_3;


																						// ����� 4.2

																						for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																								if (this_is_C_node[Amat.j[is0_3]]) {

																									if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																										
																									}
																								}

																							}
																						}

																						//********************************************************

																						// ��������� ������ ���������� i8_3.



																					}
																				}
																			}
																		}
																	}


																}
																//********************************************************

																// ��������� ������ ���������� i8_2.



															}
														}
													}
												}
											}


										}
										//********************************************************
										// ��������� ������ ���������� i8_1.



									}
								}
							}
						}
					}

				}


				// ��� ����� ��������� �������� ����� �������� ������� �����.
				// ����������� ����� �������� ����� � ������������ � ����� ����������.
				if ((iOPORNAQ_C == 0) ||
					((iOPORNAQ_C == 1) && (iOPORNAQ_C0_1 == 0) && (iOPORNAQ_C1_1 == 0)) ||
					((iOPORNAQ_C == 2) && (i_na_odnoi_vetke_C2 == 2) && (iOPORNAQ_C0_1 == 0) && (iOPORNAQ_C1_1 == 0))) {
					bproblem = true;
				}
				// ���� ��� � ����� ������, ���� ���� ������ ���� � �����, �� �� ���������� ������ �� ���������� 2 � ����.
				// ���� ���� ������ ��� ���� ������ �� ����� ����� � ������ �������.

			}

			if (ii1 == (row_startA[istr_etalon1 + 1] - 1)) bproblem = false; // ������� ������� �� �� ��������� � � ����. 01.10.2021
			if (ii1 + 1 == (row_startA[istr_etalon1 + 1] - 1)) bproblem = false; // ���������� ������� ������� 01.10.2021

			if (bproblem)
			{

				//int ir = rand() % 1000;
				//if (ir < 1000) {

					// ������ 10%.

					this_is_C_node[istr_etalon1] = true;
					this_is_F_node[istr_etalon1] = false;

					iadditionalCstatistic++;
				//}
			}
		}

		iadditionalCstatistic_gl += iadditionalCstatistic;
		iadds0 = iadditionalCstatistic;
		iadditionalCstatistic = 0;

		if (my_amg_manager.iprint_log) {

			std::cout << "iadditional C=" << iadds0 << std::endl;

			std::cout << iOPORNAQ_C0 << " " << iOPORNAQ_C1 << " " << iOPORNAQ_C2 << " " << iOPORNAQ_C3 << " " << std::endl;
		}
		bweSholdbeContinue = true;

	}
	else {
		bweSholdbeContinue = false;
	}

	//if (iadditionalCstatistic == 0) 
	{


		doublerealT theta_strong_F = theta83; // 3 ���� 2016
		doublerealT theta = theta_gl;

		const doublerealT magic = magic82; // 0.4 is recomended.
		const doublerealT magic_FF = magic82; // 18.09.2020
		//const doublerealT magic_FFF = magic82; // 24.09.2021

		const int ilimCadd = 100;

		//int itid = 0;

		integer LOC_P_LIM_SIZE = -1;
		double dsize = 2.2 * (1.0 / (1.0 * inth)) * nsizePR * n_a[0];
		LOC_P_LIM_SIZE = static_cast<integer>(dsize);

		Ak1** P_loc = new Ak1 * [inth];
		integer* isize_P_loc = new integer[inth];
		for (int itid = 0; itid < inth; ++itid)
		{
			isize_P_loc[itid] = 0;

			P_loc[itid] = new Ak1[LOC_P_LIM_SIZE];
		}

		bool old_def_threshold = true;


		iOPORNAQ_C0 = 0;
		iOPORNAQ_C1 = 0;
		iOPORNAQ_C2 = 0;
		iOPORNAQ_C3 = 0;

#pragma omp parallel for reduction(+:  iOPORNAQ_C0, iOPORNAQ_C1, iOPORNAQ_C2, iOPORNAQ_C3)
		for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) 
		{

#ifdef _OPENMP
			int itid = omp_get_thread_num();
#else
			int itid = 0;
#endif

			if (this_is_F_node[i8]) {


			int iOPORNAQ_C = 0;

			// ��� ����� ����������� �������� ��������.
						// 5 ������� 2015 ���� �� ��������� ��������� �������������
						// ��������� ������������ � ������ � ��������� ��������.
			doublerealT maxelem_threshold = -1.0;

			integer ii1 = row_startA[i8];
			integer istr_etalon1 = Amat.i[ii1];
			integer iend_for1 = -1;
			for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
				iend_for1 = is0;
			}
			if (old_def_threshold) {
				for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
					if ((Amat.j[is0] != istr_etalon1) && (Amat.aij[is0] < 0.0) && (-Amat.aij[is0] > maxelem_threshold)) maxelem_threshold = -Amat.aij[is0];
				}
			}
			else {
				maxelem_threshold = threshold_quick_all[istr_etalon1];
			}

			// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



			doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
			doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;


			// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
			// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
			doublerealT sumPC = 0.0;
			doublerealT sumPF = 0.0;
			doublerealT SumPall = 0.0;

			// ����� 1.1

			for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
				if (Amat.j[is0] != istr_etalon1) {
					if (this_is_C_node[Amat.j[is0]]) {

						if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
							sumPC += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
							iOPORNAQ_C++;
						}
					}
					else {
						if (this_is_F_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

								// ������������� Strong F �����.

								sumPF += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
							}
						}
					}
				}
			}

			if (iOPORNAQ_C >= ilimCadd) {
				SumPall = sumPC;

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								Ak1 Ptmp;

								Ptmp.j = static_cast<integer_mix_precision>(i8);
								Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
								//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
								if (fabs(SumPall) < 1.0e-33) {
									printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
									//system("pause");
									system("PAUSE");
									exit(1);
								}
								Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
			else {

				SumPall = sumPC + sumPF;


				// ����� 1.2

				for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
					if (Amat.j[is0] != istr_etalon1) {
						if (this_is_C_node[Amat.j[is0]]) {

							if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
								Ak1 Ptmp;

								Ptmp.j = static_cast<integer_mix_precision>(i8);
								Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
								//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
								if (fabs(SumPall) < 1.0e-33) {
									printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
									//system("pause");
									system("PAUSE");
									exit(1);
								}
								Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
						else {
							if (this_is_F_node[Amat.j[is0]]) {

								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {

									// ������������� Strong F �����.

									integer i8_1 = Amat.j[is0];


									doublerealT nu1 = fabs(Amat.aij[is0]) / SumPall;

									// ����� ���������� ����� ��� ���� i8_1.

									//********************************************************

									doublerealT maxelem_threshold_1 = -1.0;

									integer ii1_1 = row_startA[i8_1];
									integer istr_etalon1_1 = Amat.i[ii1_1];
									integer iend_for1_1 = -1;
									for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
										iend_for1_1 = is0_1;
									}
									if (old_def_threshold) {
										for (integer is0_1 = ii1_1; (is0_1 <= row_startA[istr_etalon1_1 + 1] - 1); ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1_1) && (Amat.aij[is0_1] < 0.0) && (-Amat.aij[is0_1] > maxelem_threshold_1)) maxelem_threshold_1 = -Amat.aij[is0_1];
										}
									}
									else {
										maxelem_threshold_1 = threshold_quick_all[istr_etalon1_1];
									}

									// ����� maxelem_threshold_1 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



									doublerealT maxelem_threshold_theta_1 = maxelem_threshold_1 * theta;
									doublerealT maxelem_threshold_theta_strong_F_1 = maxelem_threshold_1 * magic;


									// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
									// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
									doublerealT sumPC_1 = 0.0;
									doublerealT sumPF_1 = 0.0;
									doublerealT SumPall_1 = 0.0;

									// ����� 2.1

									for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
										if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
											if (this_is_C_node[Amat.j[is0_1]]) {

												if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
													sumPC_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.	
													iOPORNAQ_C++;
												}
											}
											else {
												if (this_is_F_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

														// ������������� Strong F �����.

														sumPF_1 += fabs(Amat.aij[is0_1]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
													}
												}
											}
										}
									}


									if (iOPORNAQ_C >= ilimCadd) {
										SumPall_1 = sumPC_1;

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if (Amat.j[is0_1] != istr_etalon1_1) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														Ak1 Ptmp;

														Ptmp.j = static_cast<integer_mix_precision>(i8);
														Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_1]]);
														//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
														if (fabs(SumPall_1) < 1.0e-33) {
															printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall_1);
															//system("pause");
															system("PAUSE");
															exit(1);
														}
														Ptmp.aij = static_cast<real_mix_precision>(nu1 * (fabs(Amat.aij[is0_1]) / SumPall_1));

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
									else {

										SumPall_1 = sumPC_1 + sumPF_1;


										// ����� 2.2

										for (integer is0_1 = ii1_1; is0_1 <= iend_for1_1; ++is0_1) {
											if ((Amat.j[is0_1] != istr_etalon1) && (Amat.j[is0_1] != istr_etalon1_1)) {
												if (this_is_C_node[Amat.j[is0_1]]) {

													if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_1) {
														Ak1 Ptmp;

														Ptmp.j = static_cast<integer_mix_precision>(i8);
														Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_1]]);

														if (fabs(SumPall_1) < 1.0e-33) {
															printf("error 2.0 ! division by zero. SumPall_1 =%e\n", SumPall_1);
															//system("pause");
															system("PAUSE");
															exit(1);
														}
														Ptmp.aij = static_cast<real_mix_precision>(nu1 * (fabs(Amat.aij[is0_1]) / SumPall_1));

														if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

															P_loc[itid][isize_P_loc[itid]] = Ptmp;
															isize_P_loc[itid]++;
														}
														else {
															std::cout << "Error 2.0 Increase size P_loc memory limit\n" << std::endl;
															system("PAUSE");
														}
													}
												}
												else {
													if (this_is_F_node[Amat.j[is0_1]]) {

														if (fabs(Amat.aij[is0_1]) > maxelem_threshold_theta_strong_F_1) {

															// ������������� Strong F �����.

															integer i8_2 = Amat.j[is0_1];


															doublerealT nu2 = fabs(Amat.aij[is0_1]) / SumPall_1;

															// ����� ���������� ����� ��� ���� i8_2.


															//********************************************************

															doublerealT maxelem_threshold_2 = -1.0;

															integer ii1_2 = row_startA[i8_2];
															integer istr_etalon1_2 = Amat.i[ii1_2];
															integer iend_for1_2 = -1;
															for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																iend_for1_2 = is0_2;
															}
															if (old_def_threshold) {
																for (integer is0_2 = ii1_2; (is0_2 <= row_startA[istr_etalon1_2 + 1] - 1); ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1_2) && (Amat.aij[is0_2] < 0.0) && (-Amat.aij[is0_2] > maxelem_threshold_2)) maxelem_threshold_2 = -Amat.aij[is0_2];
																}
															}
															else {
																maxelem_threshold_2 = threshold_quick_all[istr_etalon1_2];
															}

															// ����� maxelem_threshold_2 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



															doublerealT maxelem_threshold_theta_2 = maxelem_threshold_2 * theta;
															doublerealT maxelem_threshold_theta_strong_F_2 = maxelem_threshold_2 * magic_FF;


															// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
															// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
															doublerealT sumPC_2 = 0.0;
															doublerealT sumPF_2 = 0.0;
															doublerealT SumPall_2 = 0.0;

															// ����� 3.1

															for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																	if (this_is_C_node[Amat.j[is0_2]]) {

																		if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																			sumPC_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
																			iOPORNAQ_C++;
																		}
																	}
																	else {
																		if (this_is_F_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																				// ������������� Strong F �����.

																				sumPF_2 += fabs(Amat.aij[is0_2]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
																			}
																		}
																	}
																}
															}

															if (iOPORNAQ_C >= ilimCadd) {
																SumPall_2 = sumPC_2;

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if (Amat.j[is0_2] != istr_etalon1_2) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				Ak1 Ptmp;

																				Ptmp.j = static_cast<integer_mix_precision>(i8);
																				Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_2]]);
																				//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																				if (fabs(SumPall_2) < 1.0e-33) {
																					printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall_2);
																					//system("pause");
																					system("PAUSE");
																					exit(1);
																				}
																				Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * (fabs(Amat.aij[is0_2]) / SumPall_2));

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
															else {

																SumPall_2 = sumPC_2 + sumPF_2;


																// ����� 3.2

																for (integer is0_2 = ii1_2; is0_2 <= iend_for1_2; ++is0_2) {
																	if ((Amat.j[is0_2] != istr_etalon1) && (Amat.j[is0_2] != istr_etalon1_1) && (Amat.j[is0_2] != istr_etalon1_2)) {
																		if (this_is_C_node[Amat.j[is0_2]]) {

																			if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_2) {
																				Ak1 Ptmp;

																				Ptmp.j = static_cast<integer_mix_precision>(i8);
																				Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_2]]);

																				if (fabs(SumPall_2) < 1.0e-33) {
																					printf("error 3.0 ! division by zero. SumPall_2 =%e\n", SumPall_2);
																					//system("pause");
																					system("PAUSE");
																					exit(1);
																				}
																				Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * (fabs(Amat.aij[is0_2]) / SumPall_2));

																				if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																					P_loc[itid][isize_P_loc[itid]] = Ptmp;
																					isize_P_loc[itid]++;
																				}
																				else {
																					std::cout << "Error 3.0 Increase size P_loc memory limit\n" << std::endl;
																					system("PAUSE");
																				}
																			}
																		}
																		else {
																			if (this_is_F_node[Amat.j[is0_2]]) {

																				if (fabs(Amat.aij[is0_2]) > maxelem_threshold_theta_strong_F_2) {

																					// ������������� Strong F �����.

																					integer i8_3 = Amat.j[is0_2];


																					doublerealT nu3 = fabs(Amat.aij[is0_2]) / SumPall_2;

																					// ����� ���������� ����� ��� ���� i8_3.


																					//********************************************************

																					doublerealT maxelem_threshold_3 = -1.0;

																					integer ii1_3 = row_startA[i8_3];
																					integer istr_etalon1_3 = Amat.i[ii1_3];
																					integer iend_for1_3 = -1;
																					for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																						iend_for1_3 = is0_3;
																					}
																					if (old_def_threshold) {
																						for (integer is0_3 = ii1_3; (is0_3 <= row_startA[istr_etalon1_3 + 1] - 1); ++is0_3) {
																							if ((Amat.j[is0_3] != istr_etalon1_3) && (Amat.aij[is0_3] < 0.0) && (-Amat.aij[is0_3] > maxelem_threshold_3)) maxelem_threshold_3 = -Amat.aij[is0_3];
																						}
																					}
																					else {
																						maxelem_threshold_3 = threshold_quick_all[istr_etalon1_3];
																					}

																					// ����� maxelem_threshold_3 ��� ������ ������������� ���������������� �������� � ������ ����� � �������.



																					doublerealT maxelem_threshold_theta_3 = maxelem_threshold_3 * theta;
																					//doublerealT maxelem_threshold_theta_strong_F_3 = maxelem_threshold_3 * magic_FFF;


																					// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																					// ��� ������� ������ ����� ����� ����������� ������� iscos. �� ���� iscos ������ ���� 2 � �����.
																					doublerealT sumPC_3 = 0.0;
																					//doublerealT sumPF_3 = 0.0;
																					doublerealT SumPall_3 = 0.0;

																					// ����� 4.1

																					for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																						if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																							if (this_is_C_node[Amat.j[is0_3]]) {

																								if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																									sumPC_3 += fabs(Amat.aij[is0_3]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.						
																								}
																							}

																						}
																					}


																					SumPall_3 = sumPC_3;


																					// ����� 4.2

																					for (integer is0_3 = ii1_3; is0_3 <= iend_for1_3; ++is0_3) {
																						if ((Amat.j[is0_3] != istr_etalon1) && (Amat.j[is0_3] != istr_etalon1_1) && (Amat.j[is0_3] != istr_etalon1_2) && (Amat.j[is0_3] != istr_etalon1_3)) {
																							if (this_is_C_node[Amat.j[is0_3]]) {

																								if (fabs(Amat.aij[is0_3]) > maxelem_threshold_theta_3) {
																									Ak1 Ptmp;

																									Ptmp.j = static_cast<integer_mix_precision>(i8);
																									Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_3]]);

																									if (fabs(SumPall_3) < 1.0e-33) {
																										printf("error 4.0 ! division by zero. SumPall_3 =%e\n", SumPall_3);
																										//system("pause");
																										system("PAUSE");
																										exit(1);
																									}
																									Ptmp.aij = static_cast<real_mix_precision>(nu1 * nu2 * nu3 * (fabs(Amat.aij[is0_3]) / SumPall_3));

																									if (isize_P_loc[itid] < LOC_P_LIM_SIZE) {

																										P_loc[itid][isize_P_loc[itid]] = Ptmp;
																										isize_P_loc[itid]++;
																									}
																									else {
																										std::cout << "Error 4.0 Increase size P_loc memory limit\n" << std::endl;
																										system("PAUSE");
																									}
																								}
																							}

																						}
																					}

																					//********************************************************

																					// ��������� ������ ���������� i8_3.



																				}
																			}
																		}
																	}
																}


															}
															//********************************************************

															// ��������� ������ ���������� i8_2.



														}
													}
												}
											}
										}


									}
									//********************************************************
									// ��������� ������ ���������� i8_1.



								}
							}
						}
					}
				}

			}


			if (iOPORNAQ_C < 4) {
				switch (iOPORNAQ_C) {
				case 0: iOPORNAQ_C0++;
					break;
				case 1: iOPORNAQ_C1++;
					break;
				case 2: iOPORNAQ_C2++;
					break;
				case 3: iOPORNAQ_C3++;
					break;
				}
			}

			if (iOPORNAQ_C == 1) {
				//std::cout << "one uzel : " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
				P_loc[itid][isize_P_loc[itid] - 1].aij = 1.0;
			}
			if (iOPORNAQ_C == 2) {
				//std::cout << "two uzla : " << P_loc[itid][isize_P_loc[itid] - 2].aij << " " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
				doublerealT s1 = P_loc[itid][isize_P_loc[itid] - 2].aij + P_loc[itid][isize_P_loc[itid] - 1].aij;
				P_loc[itid][isize_P_loc[itid] - 2].aij /= s1;
				P_loc[itid][isize_P_loc[itid] - 1].aij /= s1;
			}
			if (iOPORNAQ_C == 3) {
				//std::cout << "three uzla : " << P_loc[itid][isize_P_loc[itid] - 3].aij << " " << P_loc[itid][isize_P_loc[itid] - 2].aij << " " << P_loc[itid][isize_P_loc[itid] - 1].aij << std::endl;
				doublerealT s1 = P_loc[itid][isize_P_loc[itid] - 3].aij + P_loc[itid][isize_P_loc[itid] - 2].aij + P_loc[itid][isize_P_loc[itid] - 1].aij;
				P_loc[itid][isize_P_loc[itid] - 3].aij /= s1;
				P_loc[itid][isize_P_loc[itid] - 2].aij /= s1;
				P_loc[itid][isize_P_loc[itid] - 1].aij /= s1;
			}
			if (iOPORNAQ_C >= 4) {
				doublerealT s1 = 0.0;
				for (int i_5 = 1; i_5 <= iOPORNAQ_C; ++i_5) s1 += P_loc[itid][isize_P_loc[itid] - i_5].aij;
				for (int i_5 = 1; i_5 <= iOPORNAQ_C; ++i_5) P_loc[itid][isize_P_loc[itid] - i_5].aij /= s1;
			}

			//system("pause");

		}

	    }

		if (my_amg_manager.iprint_log) {
			std::cout << "interpolatioon statistics: \n";
			std::cout << iOPORNAQ_C0 << " " << iOPORNAQ_C1 << " " << iOPORNAQ_C2 << " " << iOPORNAQ_C3 << " " << iadds0 << " " << std::endl;
			//system("pause");
		}

		// �������� ����������� ����������.
		for (int itid = 0; itid < inth; ++itid) 
		{

			if (icount1 + isize_P_loc[itid] >= nsizePR * n) {
				printf("memory error!!!\n");
				printf("not enough memory for the interpolation operator.\n");
				//system("PAUSE");
				//exit(1);
				deallocate_prolongation(nsizePR, n, P);
			}

#pragma omp parallel for
			for (integer iscan_loc = 0; iscan_loc < isize_P_loc[itid]; ++iscan_loc) {
				P[icount1 + iscan_loc] = P_loc[itid][iscan_loc];
			}
			icount1 += isize_P_loc[itid];

		}


		for (int itid = 0; itid < inth; ++itid)
		{
			isize_P_loc[itid] = 0;

			delete[] P_loc[itid];
		}

		delete[] isize_P_loc;
		delete[] P_loc;


	}





















}

// ������������ ������������ �� ������ ��������������� ��������� �3.
  // ���������������� ������� ������������� �� ������� ����.
  // �� ���� ������������ �� ������� ���� �������� ������� ����.
  // 12.05.2018 - ������ ���������� � ���������� ������������.
  // 18.09.2020. 
  // 12.10.2020 ������ ��������� ������������ ������ ������������ � ����������� 3.
// 25.09.2021 FFFC ������
template <typename doublerealT>
void my_long_rang_interpolation_procedure_FFFC_number3A_parallelx8_old(integer& the_number_of_neighbors_that_are_not_C_nodes_gl,
	integer& number_of_F_nodes_with_one_single_strong_C_neighbor_gl,
	integer*& n_a, bool*& this_is_F_node, integer*& row_startA,
	integer*& nnz_a, bool& bpositive_connections, Ak2& Amat,
	bool& bweSholdbeContinue, bool*& this_is_C_node, integer& iadditionalCstatistic_gl,
	const doublerealT RealZERO, integer& icount1, Ak1*& P, integer& nsizePR, integer& ilevel,
	integer& iadd, doublerealT& theta_gl, integer& n, integer*& C_numerate,
	integer& number_of_F_nodes_with_one_single_strong_C_neighborF_gl,
	doublerealT& theta83, bool& btreshold_on_new_vetv, integer& ifrom_re_operation_protection,
	bool& from_re_operation_protection0, doublerealT& magic82, doublerealT*& threshold_quick_all,
	doublerealT*& threshold_quick_only_negative, bool bsuffix_work, bool baggressive_interpol)
{

	// ���������� ���������� ����� � ����� 01.10.2020
	bool weakening_the_addition_of_new_C_nodes = true;
	if (!b_on_adaptive_local_refinement_mesh) {
		// ������ ����������� ������������� ������ ��� ���� ����� �
		// �� ����������������� ����� �� ��������.
		weakening_the_addition_of_new_C_nodes = false;
	}


	//bool bprint_mesage_diagnostic = true;
	//if (my_amg_manager.iprint_log == 0) {
		//bprint_mesage_diagnostic = false;
	//}


#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(8); // ���������� 8 �������, 10 ������� ��� �������� �� �������.
#endif

	//printf("incomming long range%e\n", free_debug_parametr1);
	//system("pause");

	// do_not_let_one_neighbor = false; �������� ���� � ����� � ����� ������ ���������� � �������.
	// do_not_let_one_neighbor =true; ����������� ����� � ���� � �������� �������.
	bool do_not_let_one_neighbor = true;// bsuffix_work; // true; // ����� � ���� ����������� ����� ��������� ����� ���������� �����������.

	//integer icount1_memo = icount1;

	//theta = 0.24;
	// theta_strong_F iter_number time,s
	// 0.21 56 22.63
	// 0.22 55 21.769
	// 0.23 52 21.488
	// 0.24 52 21.741 == theta // optimum
	// 0.26 69 24.623
	//doublerealT theta_strong_F = 0.23; // ����������� �����.
	doublerealT theta_strong_F = theta83; // 3 ���� 2016
	doublerealT theta = theta_gl;


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

	//integer icount_93 = 0;

	// the begining

	integer iadditionalCstatistic = 0;
	integer number_of_F_nodes_with_one_single_strong_C_neighborF = 0;

	bool byes_add = false;
	if (do_not_let_one_neighbor) {
		// ������� ���������� ����������� � �����.
		if (1) {
			// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
			//doublerealTT magic = 0.4; // 0.4 optimum


			integer the_number_of_neighbors_that_are_not_C_nodes = 0;
			integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;


			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
#pragma omp parallel for reduction(+: number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes, iadditionalCstatistic, number_of_F_nodes_with_one_single_strong_C_neighborF)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8)
			{
				if (this_is_F_node[i8]) {


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

					//system("pause");


					// ��� ����� ����������� �������� ��������.
					// 5 ������� 2015 ���� �� ��������� ��������� �������������
					// ��������� ������������ � ������ � ��������� ��������.
						doublerealT maxelem_threshold = -1.0;

						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {

						if (!btreshold_on_new_vetv) {
							for (integer is0 = ii1; (is0 <= iend_marker_position); ++is0) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
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
						// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


						doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
						doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); ++is0) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (this_is_C_node[Amat.j[is0]]) {
									//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
										sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
										icsos++;
									}
								}
								else {
									if (this_is_F_node[Amat.j[is0]]) {
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

						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {

						if (!btreshold_on_new_vetv) {
							for (integer is0 = ii1; (is0 <= iend_marker_position); ++is0) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
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
						// ����� maxelem_threshold ��� ������ ������������� ���������������� �������� � ������ ����� � �������.


						doublerealT maxelem_threshold_theta = maxelem_threshold * theta;
						doublerealT maxelem_threshold_theta_strong_F = maxelem_threshold * theta_strong_F;
						for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); ++is0) {
							if (Amat.j[is0] != Amat.i[ii1]) {
								if (this_is_C_node[Amat.j[is0]]) {
									//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
										sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
										icsos++;
									}
								}
								else {
									if (this_is_F_node[Amat.j[is0]]) {
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

					bool baddC = true;
					for (integer is0 = ii1; (is0 <= iend_marker_position); ++is0) {
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

					if (baggressive_interpol) {
						baddC = false;
					}

					if (baddC) {


						if (do_not_let_one_neighbor) {
							// �� ������ �� F ���� C ���� ������ � ��� ������ ���� � ��� ���� � ������.
							//integer iend_marker_position = row_startA[Amat.i[ii1] + 1] - 1;
							for (integer is0 = ii1; (is0 <= iend_marker_position); ++is0) {
								if (this_is_F_node[Amat.j[is0]]) {
									if (Amat.j[is0] != Amat.i[ii1]) {


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

												// � ������ �������� ��� ����������� ��������� PMIS ������ ��������� ������
												// ��� � ���� � ���������. ���������� ����������� ��������� ��� ������� 
												// ����� � ������ ����� ���� �������� � ������ �������� ����������.
												// �������� ���������� ����������� ����� ������.
												//if ((icount_93 % (static_cast<integer>(free_debug_parametr1)) == 0))
												{
													// ���������� ����� ����� ������ �����������, ����� ���������� 
													// �����������. 21.09.2020
													if (this_is_C_node[i8] == false) iadditionalCstatistic++;
													this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
													this_is_C_node[i8] = true;
													//bweSholdbeContinue = true;

#pragma omp critical
													{// ����� ������ ������ ���, ��� ��� ������ Acomp
														if (!byes_add) {
															byes_add = true;  // ���� ���������� �����.
														}
													}


													//exit(1);
													// ����� ����� �������� ������������.
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
	1. ���������� threshold�.
	2. ������� �����. SumP - ��� strong C, SumPall - ��� strong F.
	3. icsos - ���������� strong C. icsosF - ���������� strong F.
	4. ���� SumP ==0 �� ���� ���������� � �����. iadditionalCstatistic++;
	*/

	// ���� do_not_let_one_neighbor =true, �� byes_add==false.

	if (1 && (!byes_add)) {

		// � ���������� 0.4 ��������� ������������ �������� ����� ������� ������.
		//--->doublerealT magic = 0.4; // 0.4 optimum
		//magic = 0.3; // 3 ���� 2016 ��� ������������ �����
		// �������� ������� �� ���� ���������
		// �� �� �������������� �� �� ����� V ������.
		//magic = 0.5 - 0.2*ilevel / 12.0;
		const doublerealT magic = magic82; // 0.4 is recomended.

		const doublerealT magic_FF = magic82; // 18.09.2020
										//magic_FF =  0.45; // 20.09.2020 //free_debug_parametr1;//

		const doublerealT magic_FFF = magic82; // 24.09.2021

		integer LOC_P_LIM_SIZE = -1;
		double dsize = 2.2 * 0.125 * nsizePR * n_a[0];
		LOC_P_LIM_SIZE = static_cast<integer>(dsize);

		Ak1** P_loc = new Ak1 * [8];
		integer* isize_P_loc = new integer[8];
		for (int itid = 0; itid < 8; ++itid) {
			isize_P_loc[itid] = 0;

			P_loc[itid] = new Ak1[LOC_P_LIM_SIZE];
		}


		bool** permission_to_use_F = new bool* [8];// ���������� �� ������������� F.
		bool** bactiveF = new bool* [8];


		for (int itid = 0; itid < 8; ++itid) {

			permission_to_use_F[itid] = new bool[n_a[ilevel - 1] + 1];
#pragma omp parallel for
			for (integer i8 = 0; i8 <= n_a[ilevel - 1]; ++i8) {
				permission_to_use_F[itid][i8] = true; // �������������, �� ��������� �� ���������.
			}
			bactiveF[itid] = new bool[n_a[ilevel - 1] + 1];
#pragma omp parallel for
			for (integer i8 = 0; i8 <= n_a[ilevel - 1]; ++i8) {
				bactiveF[itid][i8] = true; // �������������, �� ��������� �� ���������.
			}
		}

		integer the_number_of_neighbors_that_are_not_C_nodes = 0;
		integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;

		if (bpositive_connections) {

			// ���������� ����������� ��� ����� ������� ���������� F-nodes.
			// ������ F-nodes ������ C-nodes.
#pragma omp parallel for reduction(+ : number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes,  number_of_F_nodes_with_one_single_strong_C_neighborF, iadditionalCstatistic) schedule (guided)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) {

#ifdef _OPENMP 
				int itid = omp_get_thread_num();
#else
				int itid = 0;
#endif		

				if (this_is_F_node[i8]) {

					permission_to_use_F[itid][i8] = false;

					// ��� ����� ����������� �������� ��������.
					// 5 ������� 2015 ���� �� ��������� ��������� �������������
					// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;
					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer istr_etalon1 = Amat.i[ii1];
					integer iend_for1 = -1;
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
							iend_for1 = is0;
							if (Amat.j[is0] != istr_etalon1) {
								// ���� ���������������� �� ��������� ������������ ��������������� ������� � ������.
								//if (this_is_C_node[Amat.j[is0]] ) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold) {
									maxelem_threshold = fabs(Amat.aij[is0]);
								}
								//}
							}
						}
					}
					else {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
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

					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {
					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
								if (fabs(Amat.aij[is0]) > maxelem_threshold_theta) {
									sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


										// ����� SumPall ����������� ������ ���� � ��� F ���� ������� � ������ ����� magic.

										// ������������� Strong F �����.

										// �����:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
										// � ������� F ������.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// ��������� ���� ������� iFpointeger 
										// ����� ����� ����� ��� � ����.

										// �������������� ��������� �����.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); ++is0_loc) {
											iend_for = is0_loc;
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
												// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
												// ����� ������ ����������.
												if (this_is_C_node[Amat.j[is0_loc]]) {
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
										integer icsos_locall = 0;

										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

											// ����� �������� ����� ���������� ����� ���������.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

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

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;

										if (fabs(sumP_loc) > 1.0e-33) {
											// ��������� ���� �������������� F ���� ������ � ��� ������ ���� � ���� ���� ������� C ������.
											bactiveF[itid][Amat.j[is0]] = true;

											// ���� ��� ������ ��������������� �� ����������
											// ����� ���������. 247 �������� ������ 20 ��������.
											permission_to_use_F[itid][Amat.j[is0]] = false;

											SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
											icount_StronglyF++;
											icsosF++;
										}
										else {
											bactiveF[itid][Amat.j[is0]] = false;
										}
									}
								}
								// ������������ ���������� ������� ������� �� �������� � ������.
								the_number_of_neighbors_that_are_not_C_nodes++; // ������������ �������� ������������ 
							}
						}
					}


					if (weakening_the_addition_of_new_C_nodes && (icsos == 0)) {
						bool baddC = true;
						doublerealT Cmax = -1.0; // �������� ������� � ����� �� ������ � ������.
						for (integer is0 = ii1; (is0 <= iend_for1); ++is0) {
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
							maxelem_threshold_theta = (doublerealT)(0.98 * Cmax);
							for (integer is0 = ii1; (is0 <= iend_for1); ++is0) {
								if (this_is_C_node[Amat.j[is0]]) {
									if (Amat.j[is0] != Amat.i[ii1]) {
										doublerealT aij_abs = Amat.abs_aij[is0];
										if (aij_abs > 0.98 * Cmax) {
											sumP += aij_abs; // ����� ������� ��������������� ��������� ������� ����������� Strongly C �����.
											icsos++;
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
						if (icsosF == 0) number_of_F_nodes_with_one_single_strong_C_neighborF++; // ���������� F ����� � ����� ������������ �������  C �������
																								 // � � ����-�� �� ������� ������� F �������.
					}



					{

						//if (((icsos == 1) || (icsos == 2) || (icsos == 3)) && (icsosF != 0)) {
						if (1) {
							//if (((icsos == 1) || (icsos == 2) || (icsos == 3) || (icsos >= 4)) && (icsosF != 0)) {
							// ������ ������ Strong C ������ � ������� � �������� ���� �� ������� ���� ���� Strong F �����.
							//
							SumPall += sumP;

							for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); ++is0) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// ��� ���������� ������ Strong �����.


									if (this_is_C_node[Amat.j[is0]]) {
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
#pragma omp critical
												{
													if (!bweSholdbeContinue) {
														bweSholdbeContinue = true;
													}
												}
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

													Ak1 Ptmp;

													Ptmp.j = static_cast<integer_mix_precision>(i8);
													Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
													//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
													if (fabs(SumPall) < 1.0e-33) {
														printf("error 1.0 ! division by zero. SumPall =%e\n", SumPall);
														//system("pause");
														system("PAUSE");
														exit(1);
													}
													Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
												// ������������� Strong F �����.

												// �����:
												// 



												//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
												// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
												// � ������� F ������.
												//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

												integer iFpoint = Amat.j[is0];
												if (fabs(SumPall) < 1.0e-33) {
													printf("error 2.0 ! division by zero. SumPall =%e\n", SumPall);
													//system("pause");
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
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); ++is0_loc) {
													iend_for = is0_loc;
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
														// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
														// ����� ������ ����������.
														if (this_is_C_node[Amat.j[is0_loc]]) {
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
												integer icsos_locall = 0;

												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
												// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
												doublerealT sumP_loc = 0.0;
												doublerealT sumP_locall = 0.0;
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

													// ����� �������� ����� ���������� ����� ���������.
													if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

														if (this_is_C_node[Amat.j[is0_loc]]) {

															if (Amat.j[is0_loc] != istr_etalon) {

																sumP_loc += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																icsos_loc++;
															}

														}
														else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {

															//if (Amat.j[is0_loc] != istr_etalon) {
															// ������������ ���������� ������� ������� �� �������� � ������.
															//the_number_of_neighbors_that_are_not_C_nodes_loc++; // ������������ �������� ������������ 
															//}

															if (Amat.j[is0_loc] != istr_etalon) {



																// F node

																integer iFpoint_loc = Amat.j[is0_loc];

																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// ��������� ���� ������� iFpointeger_loc 
																// ����� ����� ����� ��� � ����.

																// �������������� ��������� �����.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); ++is0_loc2) {
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2) {
																		// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
																		// ����� ������ ����������.
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
																	// ���� F ��� � ����� ������ �� �� �����������.

																	sumP_locall += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																	icsos_locall++;
																}
															}

														}

													}
												}

												//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												sumP_locall += sumP_loc;

												if (fabs(sumP_locall) > 1.0e-33)
												{
													// ����������� ������ � ��� ������ ���� ���� ������� � ������.

													// � ����� ��� ������� ���������������� ����� 
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

														// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


														// ��������������� ������� �� ��������� � �����.

														// ������ ������� ������ ����������� ��������� 
														// ������������� ��������� �������� �� �������� 
														// �������� �������.
														// ����������� 5 ������� 2015.
														if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {
																	if (fabs(sumP_loc) > 1.0e-33)
																	{

																		Ak1 Ptmp;

																		Ptmp.j = static_cast<integer_mix_precision>(i8);
																		Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_loc]]);
																		//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																		if (fabs(sumP_loc) < 1.0e-33) {
																			printf("error 3.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																			//system("pause");
																			system("PAUSE");
																			exit(1);
																		}
																		Ptmp.aij = static_cast<real_mix_precision>(multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall);

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
															else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {
																// F node

																integer iFpoint_loc = Amat.j[is0_loc];
																if (fabs(sumP_locall) < 1.0e-33) {
																	printf("error 2.a ! division by zero. SumP_locall =%e\n", sumP_locall);
																	//system("pause");
																	system("PAUSE");
																	exit(1);
																}
																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// ��������� ���� ������� iFpointeger_loc 
																// ����� ����� ����� ��� � ����.

																// �������������� ��������� �����.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); ++is0_loc2) {
																	if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2) {
																		// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
																		// ����� ������ ����������.
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
																	// ����� maxelem_threshold_loc2 ��� ������ ������������� ���������������� 
																	// �������� � ������ ����� � ������� ��������.

																	// ����� ������� ������� F-node ������� C-node.
																	integer icsos_loc2 = 0;

																	// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																	// ��� ������� ������ ����� ����� ����������� ������� iscos_loc2. 
																	// �� ���� iscos_loc2 ������ ���� 2 � �����.
																	doublerealT sumP_loc2 = 0.0;
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); ++is0_loc2) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; ++is0_loc2) {

																		// ����� �������� ����� ���������� ����� ���������.
																		if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {

																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					sumP_loc2 += fabs(Amat.aij[is0_loc2]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																					icsos_loc2++;
																				}

																			}
																			else {

																				//if (Amat.j[is0_loc2] != istr_etalon2) {

																			}

																		}
																	}

																	//doublerealT maxelem_threshold_loc_magic_minus2 = -maxelem_threshold_loc_magic2;

																	// � ����� ��� ������� ���������������� ����� 
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); ++is0_loc2) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; ++is0_loc2) {

																		// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


																		// ��������������� ������� �� ��������� � �����.

																		// ������ ������� ������ ����������� ��������� 
																		// ������������� ��������� �������� �� �������� 
																		// �������� �������.
																		// ����������� 5 ������� 2015.
																		if (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {
																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					

																					Ak1 Ptmp;

																					Ptmp.j = static_cast<integer_mix_precision>(i8);
																					Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_loc2]]);
																					//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																					if (fabs(sumP_loc2) < 1.0e-33) {
																						printf("error 3.a ! division by zero. sumP_loc2 =%e\n", sumP_loc2);
																						//system("pause");
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = static_cast<real_mix_precision>(multiplyer_nu_loc * fabs(Amat.aij[is0_loc2]) / sumP_loc2);

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

					// ������� ���.
					permission_to_use_F[itid][i8] = true;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {
					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							{
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


										// ����� SumPall ����������� ������ ���� � ��� F ���� ������� � ������ ����� magic.

										// ������������� Strong F �����.

										// �����:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
										// � ������� F ������.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// ��������� ���� ������� iFpointeger 
										// ����� ����� ����� ��� � ����.

										// �������������� ��������� �����.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); ++is0_loc) {
											iend_for = is0_loc;
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc) {
												// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
												// ����� ������ ����������.
												if (this_is_C_node[Amat.j[is0_loc]]) {
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
										integer icsos_locall = 0;

										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

											// ����� �������� ����� ���������� ����� ���������.
											if (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

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

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;
										if (fabs(sumP_loc) > 1.0e-33) {
											// ��������� ���� �������������� F ���� ������ � ��� ������ ���� � ���� ���� ������� C ������.
											permission_to_use_F[itid][Amat.j[is0]] = true;
										}
										bactiveF[itid][Amat.j[is0]] = true;// �������������.
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

			std::cout << "only negative connections" << std::endl;
			system("pause");

			// ���������� ����������� ��� ����� ������� ���������� F nodes.
			// ������ F-nodes ������ C-nodes.
#pragma omp parallel for reduction(+ : number_of_F_nodes_with_one_single_strong_C_neighbor,	the_number_of_neighbors_that_are_not_C_nodes,  number_of_F_nodes_with_one_single_strong_C_neighborF, iadditionalCstatistic)
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) {

#ifdef _OPENMP 
				int itid = omp_get_thread_num();
#else
				int itid = 0;
#endif

				if (this_is_F_node[i8]) {

					permission_to_use_F[itid][i8] = false;

					// ��� ����� ����������� �������� ��������.
					// 5 ������� 2015 ���� �� ��������� ��������� �������������
					// ��������� ������������ � ������ � ��������� ��������.
					doublerealT maxelem_threshold = -1.0;
					//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					integer ii1 = row_startA[i8];
					integer istr_etalon1 = Amat.i[ii1];
					integer iend_for1 = -1;
					if (!btreshold_on_new_vetv) {
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
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
						for (integer is0 = ii1; (is0 <= row_startA[istr_etalon1 + 1] - 1); ++is0) {
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
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {
					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							if (this_is_C_node[Amat.j[is0]]) {
								//	if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
								//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
								if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
									sumP += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly � �����.
									icsos++;
								}
							}
							else {
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F)) {


										// ������������� Strong F �����.

										// �����:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
										// � ������� F ������.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// ��������� ���� ������� iFpointeger 
										// ����� ����� ����� ��� � ����.

										// �������������� ��������� �����.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); ++is0_loc) {
											iend_for = is0_loc;
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
												// ����� ������ ����������.
												if (this_is_C_node[Amat.j[is0_loc]]) {
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
										integer icsos_locall = 0;

										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

											// ����� �������� ����� ���������� ����� ���������.
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

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

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;

										if (fabs(sumP_loc) > 1.0e-33) {
											// ��������� ���� �������������� F ���� ������ � ��� ������ ���� � ���� ���� ������� C ������.

											bactiveF[itid][Amat.j[is0]] = true;

											permission_to_use_F[itid][Amat.j[is0]] = false;

											SumPall += fabs(Amat.aij[is0]); // ����� ������� ��������������� ��������� ������� ����������� Strongly F �����.
											icount_StronglyF++;
											icsosF++;
										}
										else {
											bactiveF[itid][Amat.j[is0]] = false;
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

							for (integer is0 = ii1; (is0 <= row_startA[Amat.i[ii1] + 1] - 1); ++is0) {
								if (Amat.j[is0] != Amat.i[ii1]) {
									// ��� ���������� ������ Strong �����.


									if (this_is_C_node[Amat.j[is0]]) {
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta) {
										if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {
											if ((!baggressive_interpol) && (fabs(sumP) < RealZERO)) {
												//printf("error interpolation zero diagonal sumP.\n");
												//printf("Fnode all sosed is F");
												//system("pause");
												//	printf("i8 is Dirichlet node\n");
												if (this_is_C_node[i8] == false) iadditionalCstatistic++;
												this_is_F_node[i8] = false; // ���� ���� ������� ������ � �����.
												this_is_C_node[i8] = true;

#pragma omp critical
												{
													if (!bweSholdbeContinue) {
														bweSholdbeContinue = true;
													}
												}

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
												if ((Amat.aij[is0] < 0.0) && (fabs(Amat.aij[is0]) > maxelem_threshold_theta)) {

													Ak1 Ptmp;

													Ptmp.j = static_cast<integer_mix_precision>(i8);
													Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0]]);
													//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
													if (fabs(SumPall) < 1.0e-23) {
														printf("error 5.0 ! division by zero. SumPall =%e\n", SumPall);
														//system("pause");
														system("PAUSE");
														exit(1);
													}
													Ptmp.aij = static_cast<real_mix_precision>(fabs(Amat.aij[is0]) / SumPall);

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
												integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
												for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); ++is0_loc) {
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
												// ����� maxelem_threshold_loc ��� ������ ������������� ���������������� �������� � ������ ����� � ������� ��������.

												// ����� ������� ������� F-node ������� C-node.
												integer icsos_loc = 0;
												integer icsos_locall = 0;

												// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
												// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
												doublerealT sumP_loc = 0.0;
												doublerealT sumP_locall = 0.0;
												//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
												for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

													// ����� �������� ����� ���������� ����� ���������.
													if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

														if (this_is_C_node[Amat.j[is0_loc]]) {

															if (Amat.j[is0_loc] != istr_etalon) {

																sumP_loc += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																icsos_loc++;
															}

														}
														else if ((permission_to_use_F[itid][Amat.j[is0_loc]]) && (this_is_F_node[Amat.j[is0_loc]])) {

															if (Amat.j[is0_loc] != istr_etalon) {

																// F node

																integer iFpoint_loc = Amat.j[is0_loc];

																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// ��������� ���� ������� iFpointeger_loc 
																// ����� ����� ����� ��� � ����.

																// �������������� ��������� �����.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); ++is0_loc2) {
																	if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2)) {
																		// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
																		// ����� ������ ����������.
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
																	// ���� F ��� � ����� ������ �� �� �����������.

																	sumP_locall += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																	icsos_locall++;
																}

															}
														}

													}
												}

												//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

												sumP_locall += sumP_loc;

												if (fabs(sumP_loc) > 1.0e-33) {

													// ����������� ������ � ��� ������ ���� ���� ������� � ������.


													// � ����� ��� ������� ���������������� ����� 
													//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
													for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

														// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


														// ��������������� ������� �� ��������� � �����.

														// ������ ������� ������ ����������� ��������� 
														// ������������� ��������� �������� �� �������� 
														// �������� �������.
														// ����������� 5 ������� 2015.
														if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

															if (this_is_C_node[Amat.j[is0_loc]]) {
																if (Amat.j[is0_loc] != istr_etalon) {

																	Ak1 Ptmp;

																	Ptmp.j = static_cast<integer_mix_precision>(i8);
																	Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_loc]]);
																	//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																	if (fabs(sumP_loc) < 1.0e-23) {
																		printf("error 6.0 ! division by zero. sumP_loc =%e\n", sumP_loc);
																		//system("pause");
																		system("PAUSE");
																		exit(1);
																	}
																	Ptmp.aij = static_cast<real_mix_precision>(multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_loc);


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
																	//system("pause");
																	system("PAUSE");
																	exit(1);
																}
																doublerealT multiplyer_nu_loc = multiplyer_nu * fabs(Amat.aij[is0_loc]) / sumP_locall;

																// ��������� ���� ������� iFpointeger_loc 
																// ����� ����� ����� ��� � ����.

																// �������������� ��������� �����.
																doublerealT maxelem_threshold_loc2 = -1.0;
																//integer ii1_loc2 = BinarySearchAi(Amat, iFpoint_loc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
																integer ii1_loc2 = row_startA[iFpoint_loc];
																integer istr_etalon2 = Amat.i[ii1_loc2];
																integer iend_marker_position_loc2 = row_startA[istr_etalon2 + 1] - 1;
																for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= iend_marker_position_loc2); ++is0_loc2) {
																	if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc2)) {
																		// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
																		// ����� ������ ����������.
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
																	// ����� maxelem_threshold_loc2 ��� ������ ������������� ���������������� 
																	// �������� � ������ ����� � ������� ��������.


																	// ����� ������� ������� F-node ������� C-node.
																	integer icsos_loc2 = 0;

																	// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
																	// ��� ������� ������ ����� ����� ����������� ������� iscos_loc2. 
																	// �� ���� iscos_loc2 ������ ���� 2 � �����.
																	doublerealT sumP_loc2 = 0.0;
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); ++is0_loc2) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; ++is0_loc2) {

																		// ����� �������� ����� ���������� ����� ���������.
																		if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2)) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {

																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					sumP_loc2 += fabs(Amat.aij[is0_loc2]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
																					icsos_loc2++;
																				}

																			}
																			else {

																				//if (Amat.j[is0_loc2] != istr_etalon2) {

																			}

																		}
																	}

																	//doublerealT maxelem_threshold_loc_magic_minus2 = -maxelem_threshold_loc_magic2;

																	// � ����� ��� ������� ���������������� ����� 
																	//for (integer is0_loc2 = ii1_loc2; (is0_loc2 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc2] == Amat.i[ii1_loc2]); ++is0_loc2) {
																	for (integer is0_loc2 = ii1_loc2; is0_loc2 <= iend_marker_position_loc2; ++is0_loc2) {

																		// ��� �������������� ��������� � ������ ������� ���� �������� ���������� ����� ���������.


																		// ��������������� ������� �� ��������� � �����.

																		// ������ ������� ������ ����������� ��������� 
																		// ������������� ��������� �������� �� �������� 
																		// �������� �������.
																		// ����������� 5 ������� 2015.
																		if ((Amat.aij[is0_loc2] < 0.0) && (fabs(Amat.aij[is0_loc2]) > maxelem_threshold_loc_magic2)) {

																			if (this_is_C_node[Amat.j[is0_loc2]]) {
																				if (Amat.j[is0_loc2] != istr_etalon2) {

																					Ak1 Ptmp;

																					Ptmp.j = static_cast<integer_mix_precision>(i8);
																					Ptmp.i = static_cast<integer_mix_precision>(C_numerate[Amat.j[is0_loc2]]);
																					//Ptmp.aij = fabs(Amat.aij[is0]) / sumP;
																					if (fabs(sumP_loc2) < 1.0e-33) {
																						printf("error 3.a ! division by zero. sumP_loc2 =%e\n", sumP_loc2);
																						//system("pause");
																						system("PAUSE");
																						exit(1);
																					}
																					Ptmp.aij = static_cast<real_mix_precision>(multiplyer_nu_loc * fabs(Amat.aij[is0_loc2]) / sumP_loc2);

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

					// ������� ���.
					permission_to_use_F[itid][i8] = true;
					//	for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {
					for (integer is0 = ii1; is0 <= iend_for1; ++is0) {
						if (Amat.j[is0] != istr_etalon1) {
							{
								if (this_is_F_node[Amat.j[is0]]) {
									//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {
									if (fabs(Amat.aij[is0]) > maxelem_threshold_theta_strong_F) {


										// ����� SumPall ����������� ������ ���� � ��� F ���� ������� � ������ ����� magic.

										// ������������� Strong F �����.

										// �����:
										// 



										//if (fabs(Amat.aij[is0]) > maxelem_threshold*barjer) {
										// ��� ������ �������, ����� ��� ���� ��������� ��� �� ����� ����
										// � ������� F ������.
										//if (fabs(Amat.aij[is0]) > maxelem_threshold*theta_strong_F) {

										integer iFpoint = Amat.j[is0];


										// ��������� ���� ������� iFpointeger 
										// ����� ����� ����� ��� � ����.

										// �������������� ��������� �����.
										doublerealT maxelem_threshold_loc = -1.0;
										//integer ii1_loc = BinarySearchAi(Amat, iFpoint, 1 + iadd, nnz_a[ilevel - 1] + iadd);
										integer ii1_loc = row_startA[iFpoint];
										integer istr_etalon = Amat.i[ii1_loc];
										integer iend_for = -1;
										integer iend_marker_position_loc = row_startA[istr_etalon + 1] - 1;
										for (integer is0_loc = ii1_loc; (is0_loc <= iend_marker_position_loc); ++is0_loc) {
											iend_for = is0_loc;
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc)) {
												// ����� ��� threshold ������ ���� ��������� ������ �� � �����.
												// ����� ������ ����������.
												if (this_is_C_node[Amat.j[is0_loc]]) {
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
										integer icsos_locall = 0;

										// ��������� ������ ��������������� ��������� �� � ����� ������� ������ ������.
										// ��� ������� ������ ����� ����� ����������� ������� iscos_loc. �� ���� iscos_loc ������ ���� 2 � �����.
										doublerealT sumP_loc = 0.0;
										doublerealT sumP_locall = 0.0;
										//for (integer is0_loc = ii1_loc; (is0_loc <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0_loc] == Amat.i[ii1_loc]); ++is0_loc) {
										for (integer is0_loc = ii1_loc; is0_loc <= iend_for; ++is0_loc) {

											// ����� �������� ����� ���������� ����� ���������.
											if ((Amat.aij[is0_loc] < 0.0) && (fabs(Amat.aij[is0_loc]) > maxelem_threshold_loc_magic)) {

												if (this_is_C_node[Amat.j[is0_loc]]) {

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

													if (Amat.j[is0_loc] != istr_etalon) {

														sumP_locall += fabs(Amat.aij[is0_loc]); // ����� ������� ��������������� ��������� ������� ����������� � �����.
														icsos_locall++;
													}

												}

											}
										}

										//doublerealT maxelem_threshold_loc_magic_minus = -maxelem_threshold_loc_magic;

										sumP_locall += sumP_loc;
										if (fabs(sumP_loc) > 1.0e-33) {
											// ��������� ���� �������������� F ���� ������ � ��� ������ ���� � ���� ���� ������� C ������.
											permission_to_use_F[itid][Amat.j[is0]] = true;
										}
										bactiveF[itid][Amat.j[is0]] = true;// �������������.
									}
								}

							}
						}
					}

				} // end only negative connections

			}
		}

		iadditionalCstatistic_gl += iadditionalCstatistic;

		// �������� ����������� ����������.
		for (int itid = 0; itid < 8; ++itid) {

			if (icount1 + isize_P_loc[itid] >= nsizePR * n) {
				printf("memory error!!!\n");
				printf("not enough memory for the interpolation operator.\n");
				//system("PAUSE");
				//exit(1);
				deallocate_prolongation(nsizePR, n, P);
			}

#pragma omp parallel for
			for (integer iscan_loc = 0; iscan_loc < isize_P_loc[itid]; ++iscan_loc) {
				P[icount1 + iscan_loc] = P_loc[itid][iscan_loc];
			}
			icount1 += isize_P_loc[itid];

		}


		for (int itid = 0; itid < 8; ++itid) {
			delete[] permission_to_use_F[itid];
			delete[] bactiveF[itid];
		}

		delete[] permission_to_use_F;
		delete[] bactiveF;


		// �� �������� ����������� ������.
		for (int itid = 0; itid < 8; ++itid)
		{
			isize_P_loc[itid] = 0;

			delete[] P_loc[itid];
		}

		number_of_F_nodes_with_one_single_strong_C_neighbor_gl = number_of_F_nodes_with_one_single_strong_C_neighbor;
		the_number_of_neighbors_that_are_not_C_nodes_gl = the_number_of_neighbors_that_are_not_C_nodes;

	}

	


	number_of_F_nodes_with_one_single_strong_C_neighborF_gl += number_of_F_nodes_with_one_single_strong_C_neighborF;

	/*
	if (0) {
	// ���� ������ ������.
	// 26.09.2020 ���������� ����������.

	//Ak1* P_copy = new Ak1[nsizePR * n_a[0] + 1];

	//icount1
	//for (integer i_1 = 1; i_1 < icount1; ++i_1) {
	//P_copy[i_1] = P[i_1];// �����.
	//}

	bool* bvisible = new bool[n_a[ilevel - 1] + 1];
	for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
	bvisible[i_1] = false;
	}

	doublerealT* xP = new doublerealT[n_a[ilevel - 1] + 1];
	doublerealT* xP2 = new doublerealT[n_a[ilevel - 1] + 1];
	for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
	xP[i_1] = 0.0;
	xP2[i_1] = 0.0;
	}

	integer k_1 = icount1_memo;

	//for (integer i_1 = icount1_memo; i_1 < icount1; ++i_1)
	for (integer i_1 = 1; i_1 < n_a[ilevel - 1]; ++i_1) {
	if (this_is_F_node[i_1])
	{
	//if (!bvisible[P_copy[i_1].i])
	{

	// �������������.
	//for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); ++j_1) {
	//xP[P_copy[i_1].j] = P_copy[i_1].aij;
	//}


	// smoother
	doublerealT omega = 0.6667;

	// ���� ��������.
	// ������������ ������� ���������.
	//for (integer ii = 1; ii <= n_a[ilevel-1]; ++ii)
	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); ++j_1)
	{
	//integer ii = P_copy[i_1].i;
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; ++s_1)
	{
	//integer	ii = P_copy[i_1].j;
	integer ii = Amat.j[s_1];

	integer istr = ii;
	doublerealT rold = xP[istr];

	// ������� ������ �����.

	doublerealT rsum = 0.0;
	doublerealT rsum_ap = 0.0;
	integer is1 = row_startA[ii];
	integer is2 = row_startA[ii + 1] - 1;
	doublerealT ap = 0.0;
	// ����������������� �������� �������� ����� ������.
	//#pragma omp parallel for reduction(+:rsum)
	for (integer ii1 = is1; ii1 <= is2; ++ii1)
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
	//xP2[istr] *= ap;// ����� �� ������������ �������, �.�. �
	//��������� �������� �������� �������� ��� ������������� ������������.
	xP2[istr] /= rsum_ap;

	// ����������� ������ ����� ������� ���� ����� ��������� ������
	// �.�. ��������� �������� �� �������� ����������.
	xP2[istr] = omega * xP2[istr] + (1.0 - omega) * rold; // this is SOR
	}
	}

	// find_max
	doublerealT xPmax = -1.0;
	doublerealT sumxP = 0.0;
	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); ++j_1)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; ++s_1)
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
	// �����������.
	doublerealT sumxPnew = 0.0;
	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); ++j_1)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; ++s_1)
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

	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); ++j_1)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; ++s_1)
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

	for (integer j_1 = i_1; (j_1 < icount1) && (P_copy[i_1].i == P_copy[j_1].i); ++j_1)
	{
	for (integer s_1 = row_startA[P_copy[i_1].j]; s_1 < row_startA[P_copy[i_1].j + 1]; ++s_1)
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


} // my_long_rang_interpolation_procedure_FFFC_number3A_parallelx8
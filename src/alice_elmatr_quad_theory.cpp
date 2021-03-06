// � ������ alice_elmatr_quad_theory.cpp 
// ���������� ������� ������������� ������ 
// (��� ������������� ������������ � �������������� ���������������� ������� ��� �������� ����� E).

// �� ������ ������ � ��������� �.�. � ��� ����������� �������.

/*
// 28.09.2016. ����� �� ������������� ������� ����.





// ������ ������� �� ���� �����:
if (iE > -1) {
	if (bE) {
		// ��������� ����.
		sl[iP].ae += Ge*border_neighbor[iE - maxelm].dS / dxe;
	}
	else {
		// iE ���������� ����.
		// �������� 3 ������:
		// 1. ������ ���������.
		// 2. ������� ������ ���� (� ����� ��� ��������� �� 2, 3 ��� 4 ������).
		// 3. ������� ������ ���� (����� ������� ������ ��� ����������� ������������ ������������).
		if (ilevel_alice[iP] == ilevel_alice[iE]) {
			// ������ �����:
			sl[iP].ae += Ge*dy*dz / dxe;
		}
		else if (ilevel_alice[iE] > ilevel_alice[iP]) {
			// ��� ������ ���� �� �������.
			// ����� ���������� ������ �� ae2, ae3, ae4. ��� ����� ���������� ���� �� ������.

			// ����� ����� �� ������ ����:
			// ���������.
			// ���������� �������� ��������� ������������ ������:
			doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
			volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);
			// ������� ������ �� ����� ���� ������ ��� ������� ����� ���� ������ iP.
			sl[iP].ae += Ge*dy_loc*dz_loc / dxe;
		}
		else {
			// ������� iP ���� ������ ������.
			// ���������� ������������� ������������.
			// ����� ������� ������. (24 ������ ��� ����� iE).
			//sl[iP].ae += Ge*dy*dz / dxe;
			// ���������� �������� ��������� ������������ ������:
			doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
			volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);
			sl[iP].ae += Ge*dy_loc*dz_loc / dxe;
			// ����������� �������� iN, iB, iS � iT.
			// �����������: �������� iN2 ��� iN2>-1 ��� iN==iN3==iN4==-1.
			// �.�. ���������� ���������� �������� iN �� (iN,iN2,iN3,iN4).
			// ���������� � iB, iT, iS.
			int iN_now = -1, iS_now = -1, iB_now = -1, iT_now = -1;
			if ((iN > -1) && (iN < maxelm) && (iN2 == -1) && (iN3 == -1) && (iN4 == -1)) {
				iN_now = iN;
			}
			if ((iN2 > -1) && (iN2 < maxelm) && (iN == -1) && (iN3 == -1) && (iN4 == -1)) {
				iN_now = iN2;
			}
			if ((iN3 > -1) && (iN3 < maxelm) && (iN2 == -1) && (iN == -1) && (iN4 == -1)) {
				iN_now = iN3;
			}
			if ((iN4 > -1) && (iN4 < maxelm) && (iN2 == -1) && (iN3 == -1) && (iN == -1)) {
				iN_now = iN4;
			}
			// iS
			if ((iS > -1) && (iS < maxelm) && (iS2 == -1) && (iS3 == -1) && (iS4 == -1)) {
				iS_now = iS;
			}
			if ((iS2 > -1) && (iS2 < maxelm) && (iS == -1) && (iS3 == -1) && (iS4 == -1)) {
				iS_now = iS2;
			}
			if ((iS3 > -1) && (iS3 < maxelm) && (iS2 == -1) && (iS == -1) && (iS4 == -1)) {
				iS_now = iS3;
			}
			if ((iS4 > -1) && (iS4 < maxelm) && (iS2 == -1) && (iS3 == -1) && (iS == -1)) {
				iS_now = iS4;
			}
			// iB
			if ((iB > -1) && (iB < maxelm) && (iB2 == -1) && (iB3 == -1) && (iB4 == -1)) {
				iB_now = iB;
			}
			if ((iB2 > -1) && (iB2 < maxelm) && (iB == -1) && (iB3 == -1) && (iB4 == -1)) {
				iB_now = iB2;
			}
			if ((iB3 > -1) && (iB3 < maxelm) && (iB2 == -1) && (iB == -1) && (iB4 == -1)) {
				iB_now = iB3;
			}
			if ((iB4 > -1) && (iB4 < maxelm) && (iB2 == -1) && (iB3 == -1) && (iB == -1)) {
				iB_now = iB4;
			}
			// iT
			if ((iT > -1) && (iT < maxelm) && (iT2 == -1) && (iT3 == -1) && (iT4 == -1)) {
				iT_now = iT;
			}
			if ((iT2 > -1) && (iT2 < maxelm) && (iT == -1) && (iT3 == -1) && (iT4 == -1)) {
				iT_now = iT2;
			}
			if ((iT3 > -1) && (iT3 < maxelm) && (iT2 == -1) && (iT == -1) && (iT4 == -1)) {
				iT_now = iT3;
			}
			if ((iT4 > -1) && (iT4 < maxelm) && (iT2 == -1) && (iT3 == -1) && (iT == -1)) {
				iT_now = iT4;
			}
			// ����� ������������� ������ ����������� -1.

			// ����� ����� ��� �� ����� ���� ���.
			if (neighbors_for_the_internal_node[W][iE].iNODE1 == iP) {
				// � iN � iS ������ 6 ���������.
				if ((neighbors_for_the_internal_node[W][iE].iNODE2>-1) && (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 > -1)) {
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iN_now)) {
						// �������� ���� iNODE3 && iNODE4

						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iN_now)) {

						// �������� ���� iNODE2 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iN_now)) {

						// �������� ���� iNODE2 && iNODE3
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iS_now)) {

						// �������� ���� iNODE3 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iS_now)) {


						// �������� ���� iNODE2 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iS_now)) {

						// �������� ���� iNODE2 && iNODE3
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
				} // NODE1==iP
				else if ((neighbors_for_the_internal_node[W][iE].iNODE3 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1 && NODE2. // ����������.
					// NODE3 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1 && NODE3. // ����������.
					// NODE2 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}


				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == -1)) {
					// ������������� ������ ���� NODE1 && NODE4. // ����������.
					// NODE2 && NODE3 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1  // ����������.
					// NODE2 && NODE3 && NODE4 �� ����������.

					// sl[iP].ae ��� ������.
				}
				else {
					printf("matrix assembles gran E in my_elmatr_quad_T3D error!!! nepredusmotrennaq situaciq.");
					getchar();
					exit(1);
				}
			} // NODE1 == iP
			else if (neighbors_for_the_internal_node[W][iE].iNODE2 == iP) {
				// � iN � iS ������ 6 ���������.
				if ((neighbors_for_the_internal_node[W][iE].iNODE1>-1) && (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 > -1)) {
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iN_now)) {
						// �������� ���� iNODE3 && iNODE4

						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iN_now)) {

						// �������� ���� iNODE1 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iN_now)) {

						// �������� ���� iNODE1 && iNODE3
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iS_now)) {

						// �������� ���� iNODE3 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iS_now)) {


						// �������� ���� iNODE1 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iS_now)) {

						// �������� ���� iNODE1 && iNODE3
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
				} // NODE1==iP
				else if ((neighbors_for_the_internal_node[W][iE].iNODE3 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1 && NODE2. // ����������.
					// NODE3 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE1 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1 && NODE3. // ����������.
					// NODE2 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}


				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE1 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == -1)) {
					// ������������� ������ ���� NODE1 && NODE4. // ����������.
					// NODE2 && NODE3 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE1 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1  // ����������.
					// NODE2 && NODE3 && NODE4 �� ����������.

					// sl[iP].ae ��� ������.
				}
				else {
					printf("matrix assembles gran E in my_elmatr_quad_T3D error!!! nepredusmotrennaq situaciq.");
					getchar();
					exit(1);
				}
			} // NODE1 == iP
			else if (neighbors_for_the_internal_node[W][iE].iNODE3 == iP) {
				// � iN � iS ������ 6 ���������.
				if ((neighbors_for_the_internal_node[W][iE].iNODE2>-1) && (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 > -1)) {
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iN_now)) {
						// �������� ���� iNODE1 && iNODE4

						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iN_now)) {

						// �������� ���� iNODE2 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iN_now)) {

						// �������� ���� iNODE2 && iNODE1
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iS_now)) {

						// �������� ���� iNODE1 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iS_now)) {


						// �������� ���� iNODE2 && iNODE4
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE4
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE4]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE4 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE4, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE4;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iS_now)) {

						// �������� ���� iNODE2 && iNODE1
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
				} // NODE1==iP
				else if ((neighbors_for_the_internal_node[W][iE].iNODE1 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1 && NODE2. // ����������.
					// NODE3 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1 && NODE3. // ����������.
					// NODE2 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}


				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == -1)) {
					// ������������� ������ ���� NODE1 && NODE4. // ����������.
					// NODE2 && NODE3 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE4 == -1)) {
					// ������������� ������ ���� NODE1  // ����������.
					// NODE2 && NODE3 && NODE4 �� ����������.

					// sl[iP].ae ��� ������.
				}
				else {
					printf("matrix assembles gran E in my_elmatr_quad_T3D error!!! nepredusmotrennaq situaciq.");
					getchar();
					exit(1);
				}
			} // NODE1 == iP
			else if (neighbors_for_the_internal_node[W][iE].iNODE4 == iP) {
				// � iN � iS ������ 6 ���������.
				if ((neighbors_for_the_internal_node[W][iE].iNODE2>-1) && (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 > -1)) {
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iN_now)) {
						// �������� ���� iNODE3 && iNODE1

						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iN_now)) {

						// �������� ���� iNODE2 && iNODE1
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iN_now)) {

						// �������� ���� iNODE2 && iNODE3
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iN_now == iN) {
								sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN2) {
								sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN3) {
								sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iN_now == iN4) {
								sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iS_now)) {

						// �������� ���� iNODE3 && iNODE1
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.

						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iS_now)) {


						// �������� ���� iNODE2 && iNODE1
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}


							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,B NODE4==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE1
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE1]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE1 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE1, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE1;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iS_now)) {

						// �������� ���� iNODE2 && iNODE3
						if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}


							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,B NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
							// ����� ���� iNODE3
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE3]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE3 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE3, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE3;
							}
							else {
								printf("E: N,T NODE3==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iB_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � B:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iB_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iB_now == iB) {
								sl[iP].ab -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB2) {
								sl[iP].ab2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB3) {
								sl[iP].ab3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iB_now == iB4) {
								sl[iP].ab4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,B NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}
						else if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
							// ����� ���� iNODE2
							// ����� ����� �� ����� ������ ���������������� ���� ����� � ���-�� ��� ��������� ������ ��.
							doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now] + prop[LAM][iT_now] + prop[LAM][neighbors_for_the_internal_node[W][iE].iNODE2]) / 5.0;
							// ������ N � T:
							doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
							volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
							doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
							// ���������� Ge � dxe.
							if (iS_now == iS) {
								sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS2) {
								sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS3) {
								sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}
							else if (iS_now == iS4) {
								sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
							}

							doublereal dx_loc2 = 0.0, dy_loc2 = 0.0, dz_loc2 = 0.0; // ����� �������� ������������ ������
							volume3D(iT_now, nvtx, pa, dx_loc2, dy_loc2, dz_loc2);
							doublereal dxe_loc2 = 0.5*(dx_loc2 + dx_loc);
							// ���������� Ge � dxe.
							if (iT_now == iT) {
								sl[iP].at -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT2) {
								sl[iP].at2 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT3) {
								sl[iP].at3 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}
							else if (iT_now == iT4) {
								sl[iP].at4 -= Ge_loc*dy_loc2*dz_loc2 / dxe_loc2;
							}

							if (neighbors_for_the_internal_node[W][iE].iNODE2 > -1) {
								doublereal dx_loc3 = 0.0, dy_loc3 = 0.0, dz_loc3 = 0.0; // ����� �������� ������������ ������
								volume3D(neighbors_for_the_internal_node[W][iE].iNODE2, nvtx, pa, dx_loc3, dy_loc3, dz_loc3);
								doublereal dxe_loc3 = 0.5*(dx_loc3 + dx_loc);
								// ���������� Ge � dxe.
								sl[iP].ae_dop -= Ge_loc*dy_loc3*dz_loc3 / dxe_loc3;
								sl[iP].iE_dop = neighbors_for_the_internal_node[W][iE].iNODE2;
							}
							else {
								printf("E: N,T NODE2==-1 in my_elmatr_quad_T3D.\n");
								getchar();
								exit(1);
							}
							// 27.09.2016 ���� ���������.
						}

					}
				} // NODE1==iP
				else if ((neighbors_for_the_internal_node[W][iE].iNODE3 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == -1)) {
					// ������������� ������ ���� NODE1 && NODE2. // ����������.
					// NODE3 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE2 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == -1)) {
					// ������������� ������ ���� NODE1 && NODE3. // ����������.
					// NODE2 && NODE4 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}


				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == -1)) {
					// ������������� ������ ���� NODE1 && NODE4. // ����������.
					// NODE2 && NODE3 �� ����������.
					if ((iN_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iN_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iN_now]) / 3.0;
						// ����� N:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iN_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iN_now == iN) {
							sl[iP].an -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN2) {
							sl[iP].an2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN3) {
							sl[iP].an3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iN_now == iN4) {
							sl[iP].an4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}

					}
					else if ((iS_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iS_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iS_now]) / 3.0;
						// ����� S:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iS_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iS_now == iS) {
							sl[iP].as -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS2) {
							sl[iP].as2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS3) {
							sl[iP].as3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iS_now == iS4) {
							sl[iP].as4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					if ((iT_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iT_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iT_now]) / 3.0;
						// ����� T:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iT_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iT_now == iT) {
							sl[iP].at -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT2) {
							sl[iP].at2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT3) {
							sl[iP].at3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iT_now == iT4) {
							sl[iP].at4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
					else if ((iB_now > -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == iB_now)) {
						doublereal Ge_loc = (prop[LAM][iP] + prop[LAM][iE] + prop[LAM][iB_now]) / 3.0;
						// ����� B:
						doublereal dx_loc1 = 0.0, dy_loc1 = 0.0, dz_loc1 = 0.0; // ����� �������� ������������ ������
						volume3D(iB_now, nvtx, pa, dx_loc1, dy_loc1, dz_loc1);
						doublereal dxe_loc1 = 0.5*(dx_loc1 + dx_loc);
						// ���������� Ge � dxe.
						if (iB_now == iB) {
							sl[iP].ab -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB2) {
							sl[iP].ab2 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB3) {
							sl[iP].ab3 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
						else if (iB_now == iB4) {
							sl[iP].ab4 -= Ge_loc*dy_loc1*dz_loc1 / dxe_loc1;
						}
					}
				}
				else if ((neighbors_for_the_internal_node[W][iE].iNODE2 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE3 == -1) && (neighbors_for_the_internal_node[W][iE].iNODE1 == -1)) {
					// ������������� ������ ���� NODE1  // ����������.
					// NODE2 && NODE3 && NODE4 �� ����������.

					// sl[iP].ae ��� ������.
				}
				else {
					printf("matrix assembles gran E in my_elmatr_quad_T3D error!!! nepredusmotrennaq situaciq.");
					getchar();
					exit(1);
				}
			} // NODE1 == iP 

			// END VARIANTS E
		}

	}
}
else {
	sl[iP].ae = 0.0;
}
// ������ ����� iE ���������.
*/
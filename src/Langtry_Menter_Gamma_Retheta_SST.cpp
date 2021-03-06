// ���� Langtry_Menter_Gamma_Retheta_SST.cpp ������ ��������� ������������� �������� 
// Langtry � Menter�.
// ������ 12.01.2021 10.43.
// 18.10  12.01.2021 ������ ��������� ��� ���������� ����� ��� �������������� gamma : my_elmatr_quad_turbulent_Gamma_Langtry_Menter_SST_3D.
// 13.37  13.01.2021 ������ ��������� ��� ���������� ����� ���  ������� Re_Theta: my_elmatr_quad_turbulent_Re_Theta_Langtry_MenterSST_3D.
// 14.43  13.01.2021 ������ ��������� ������� ��� �������������� gamma � ������ ������� �������: my_elmatr_quad_Gamma_Langtry_Mentor_3D_bound.
// 15.28  13.01.2021 ������ ��������� ������� ��� ������� Re_theta_T � ������ ������� �������: my_elmatr_quad_Re_Theta_T_Langtry_Mentor_3D_bound.
// 15.44  13.01.2021 ����� �� ���������� ���������� ��� gamma � Re_theta_T �.�. �� ��������� �� ������������ � ������ ������� ��������.
// 16.44  13.01.2021 ���� �������� � ������������ ����� ��������� ������������ ������� �������������� SST ��������� ��������� ������� ��������. 
// �������� � ������������ ������ � ��������� ������������ ������� �������������� ���������� ���� ���� bMentor_Langtry = true.
// 17.05 13.01.2021 ������� ������������� F1 �������������� � ������ ������ ������� ��������.
// 15.00-19.00 14.01.2021 - 11.45 15.01.2021 ������� solve() � ����� mysolverv0_03.c. 
// 12.00 15.01.2021 fluent_residual_for_x. ���� my_linalg.cpp.
// 12.47 15.01.2021  ���� my_linalg.cpp �������.
// 14.43 15.01.2021 ��������� ����������� ������� my_version_SIMPLE_Algorithm3D ��� ������ ������� ��������.
// 15.48 21.01.2021 my_unsteady_temperature.cpp
// 15.56 21.01.2021 my_amgcl_alg.cpp
// 16.34 21.01.2021 my_export_tecplot3 �� ������� �� ������ ����� � ReTheta. ������ ������� ������ ��� � ������ ������ SST, ��� �� �������.
// 16.40 21.01.2021 inputlaplas.cpp
// 16.49 21.01.2021 Mingw_input_laplas.cpp
// 18.13 21.01.2021 constr_struct.cpp
// 11.11 22.01.2021 classic_aglomerative_amg6_2018year.cpp
// 11.18 22.01.2021 my_agregat_amg.cpp

#pragma once
#ifndef LANGTRY_MENTER_GAMMA_RETHETA_SST_CPP
#define LANGTRY_MENTER_GAMMA_RETHETA_SST_CPP 1


doublereal F_length_1(doublereal &Rethetat)
{
	if (Rethetat < 0.0) {
		Rethetat = 0.0;
	}

	if (Rethetat < 400.0) {

		return 39.8189 + (-119.270 * 1e-4) * Rethetat + (-132.567 * 1e-6) * Rethetat * Rethetat;
	}
	else if ((Rethetat >= 400.0) && (Rethetat < 596)) {

		return 263.404 + (-123.939 * 1e-2) * Rethetat + (194.548 * 1e-5) * Rethetat * Rethetat + (-101.695 * 1e-8) * Rethetat * Rethetat * Rethetat;
	}
	else if ((Rethetat >= 596) && (Rethetat < 1200.0)) {

		return 0.5 - (3.0 * 1e-4) * (Rethetat - 596.0);
	}
	else {
		return 0.3188;
	}
} // F_length_1


doublereal F_lambda_theta(doublereal &lambda_theta, doublereal &Tu_val)
{
	if (lambda_theta <= 0) {
		return 1.0 + (12.986 * lambda_theta + 123.66 * lambda_theta * lambda_theta + 405.689 * lambda_theta * lambda_theta * lambda_theta) * exp(-pow(Tu_val / 1.5, 1.5));
	}
	else {
		return 1.0 + 0.275 * (1.0 - exp(-35.0 * lambda_theta)) * exp(-Tu_val / 0.5);
	}
} // F_lambda_theta

doublereal Re_theta_eq(doublereal& lambda_theta, doublereal& Tu_val)
{
	if (Tu_val < 1.0e-10) {
		printf("error Langtry_Menter_Gamma_Retheta_SST.cpp in function Re_theta_eq()\n");
		printf("negative or zero Tu\n");
		system("pause");
		//exit(1);
		Tu_val = 1.0e-10;
	}

	if (Tu_val <= 1.3)
	{
		return (1173.51 - 589.428 * Tu_val + 0.2196 / (Tu_val * Tu_val))* F_lambda_theta(lambda_theta, Tu_val);
	}
	else
	{
		return 331.50 * pow(Tu_val - 0.5658, -0.671) * F_lambda_theta( lambda_theta, Tu_val);
	}
} // Re_theta_eq

// Re_theta_t � ������� ����, �� ������� �������.
doublereal Re_theta_T_farfield(doublereal& Tu_infinity) {

	if (Tu_infinity < 1.0e-10) {
		printf("error Langtry_Menter_Gamma_Retheta_SST.cpp in function Re_theta_T_farfield()\n");
		printf("negative or zero Tu_infinity\n");
		system("pause");
		//exit(1);
		Tu_infinity = 1.0e-10;
	}

	if (Tu_infinity <= 1.3) {
		return (1173.51 - 589.428* Tu_infinity + 0.2196/(Tu_infinity* Tu_infinity));
	}
	else {
		return 331.50 * pow(Tu_infinity - 0.5658,-0.671);
	}

} // Re_theta_T_farfield


// �������� ���� ��������� ������� ���� ��� ����������� ��������� 
// ��������� - ��������, ��� ������������ ����������� ������������ ������.
// ��� ������������� ��������� ������������ Hex �����.
// ��� ������ ����������� ������ ��� �������������� (Gamma) � ������ ��������� ������������� �������� ������� ��������.
void my_elmatr_quad_turbulent_Gamma_Langtry_Menter_SST_3D(
	integer iP,
	BOUND* border_neighbor,
	integer lw,
	integer ls,
	equation3D**& sl,
	equation3D_bon**& slb,
	//doublereal** diag_coef,
	//integer iVar,
	bool btimedep,
	doublereal tauparam,
	int* ptr,
	int** nvtx,
	doublereal** potent,
	//doublereal* potent_temper,
	TOCHKA* pa,
	float** prop,
	float** prop_b,
	integer maxelm,
	int*** neighbors_for_the_internal_node,
	//doublereal* alpha,
	//doublereal dgx,
	//doublereal dgy,
	//doublereal dgz,
	//doublereal dbeta, 
	integer ishconvection,
	doublereal*& gamma_Langtry_Menter_old_time_step,
	//bool bBussineskApproach,
	//doublereal temp_ref,
	//bool bfirst_start,
	//doublereal RCh, 
	//integer iflowregime,
	//doublereal* speedoldtimestep, 
	//doublereal* toldtimestep,
	BLOCK* b,
	integer lb,
	TPROP* matlist,
	doublereal**& mf,
	//bool bVERYStable,
	//doublereal &sumanb,
	integer* ilevel_alice,
	doublereal*& distance_to_wall,
	doublereal*& SInvariantStrainRateTensor,
	bool brthdsd_ON,
	TOCHKA*& center_coord_loc,
	TOCHKA*& volume_loc,
	integer maxbound
) {

	// iP - ����� ����������� ������������ ������
	// iP ���������� �� 0 �� maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // ������ �������� ����������� �������
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];
	sl[GAMMA_LANGTRY_MENTER_SL][iP].iE = iE; sl[GAMMA_LANGTRY_MENTER_SL][iP].iN = iN; sl[GAMMA_LANGTRY_MENTER_SL][iP].iT = iT;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].iS = iS; sl[GAMMA_LANGTRY_MENTER_SL][iP].iW = iW; sl[GAMMA_LANGTRY_MENTER_SL][iP].iB = iB;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].iP = iP;

	// 26.09.2016 ������� ��� ���� �����.
	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // ������ �������� ����������� �������
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // ������ �������� ����������� �������
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // ������ �������� ����������� �������

	if (b_on_adaptive_local_refinement_mesh) {
		// NON_EXISTENT_NODE ���� �� ������������ � [0..maxelm+maxbound-1] ���� ������������.
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}


	sl[GAMMA_LANGTRY_MENTER_SL][iP].iE2 = iE2; sl[GAMMA_LANGTRY_MENTER_SL][iP].iN2 = iN2; sl[GAMMA_LANGTRY_MENTER_SL][iP].iT2 = iT2;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].iS2 = iS2; sl[GAMMA_LANGTRY_MENTER_SL][iP].iW2 = iW2; sl[GAMMA_LANGTRY_MENTER_SL][iP].iB2 = iB2;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].iE3 = iE3; sl[GAMMA_LANGTRY_MENTER_SL][iP].iN3 = iN3; sl[GAMMA_LANGTRY_MENTER_SL][iP].iT3 = iT3;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].iS3 = iS3; sl[GAMMA_LANGTRY_MENTER_SL][iP].iW3 = iW3; sl[GAMMA_LANGTRY_MENTER_SL][iP].iB3 = iB3;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].iE4 = iE4; sl[GAMMA_LANGTRY_MENTER_SL][iP].iN4 = iN4; sl[GAMMA_LANGTRY_MENTER_SL][iP].iT4 = iT4;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].iS4 = iS4; sl[GAMMA_LANGTRY_MENTER_SL][iP].iW4 = iW4; sl[GAMMA_LANGTRY_MENTER_SL][iP].iB4 = iB4;



	// ���������������� ���������.
	sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].an = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].as = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].at = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = 0.0;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = 0.0;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = 0.0;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = 0.0;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = 0.0;


	// ������� ����������� �����.
	// �� ������� ������ ����� ���������� � ����� �������� ������ ���.
	sl[GAMMA_LANGTRY_MENTER_SL][iP].bE2 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bW2 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bS2 = false;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].bN2 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bB2 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bT2 = false;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].bE3 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bW3 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bS3 = false;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].bN3 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bB3 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bT3 = false;

	sl[GAMMA_LANGTRY_MENTER_SL][iP].bE4 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bW4 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bS4 = false;
	sl[GAMMA_LANGTRY_MENTER_SL][iP].bN4 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bB4 = false; sl[GAMMA_LANGTRY_MENTER_SL][iP].bT4 = false;

	if (CHECK_NODE_FOR_EXISTENCE(iE2)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bE2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iW2)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bW2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iN2)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bN2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iS2)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bS2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iT2)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bT2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iB2)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bB2 = true;

	if (CHECK_NODE_FOR_EXISTENCE(iE3)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bE3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iW3)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bW3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iN3)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bN3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iS3)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bS3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iT3)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bT3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iB3)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bB3 = true;

	if (CHECK_NODE_FOR_EXISTENCE(iE4)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bE4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iW4)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bW4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iN4)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bN4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iS4)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bS4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iT4)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bT4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iB4)) sl[GAMMA_LANGTRY_MENTER_SL][iP].bB4 = true;

	// ���������� ��.	

	// ���� � ����� �� ������ ����� ������� ��������� �������
	// �� ��������������� ���������� ����� true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if ((iE >= maxelm) && (iE < maxelm + maxbound)) bE = true;
	if ((iN >= maxelm) && (iN < maxelm + maxbound)) bN = true;
	if ((iT >= maxelm) && (iT < maxelm + maxbound)) bT = true;
	if ((iW >= maxelm) && (iW < maxelm + maxbound)) bW = true;
	if ((iS >= maxelm) && (iS < maxelm + maxbound)) bS = true;
	if ((iB >= maxelm) && (iB < maxelm + maxbound)) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if ((iE2 >= maxelm) && (iE2 < maxelm + maxbound)) bE2 = true;
	if ((iN2 >= maxelm) && (iN2 < maxelm + maxbound)) bN2 = true;
	if ((iT2 >= maxelm) && (iT2 < maxelm + maxbound)) bT2 = true;
	if ((iW2 >= maxelm) && (iW2 < maxelm + maxbound)) bW2 = true;
	if ((iS2 >= maxelm) && (iS2 < maxelm + maxbound)) bS2 = true;
	if ((iB2 >= maxelm) && (iB2 < maxelm + maxbound)) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if ((iE3 >= maxelm) && (iE3 < maxelm + maxbound)) bE3 = true;
	if ((iN3 >= maxelm) && (iN3 < maxelm + maxbound)) bN3 = true;
	if ((iT3 >= maxelm) && (iT3 < maxelm + maxbound)) bT3 = true;
	if ((iW3 >= maxelm) && (iW3 < maxelm + maxbound)) bW3 = true;
	if ((iS3 >= maxelm) && (iS3 < maxelm + maxbound)) bS3 = true;
	if ((iB3 >= maxelm) && (iB3 < maxelm + maxbound)) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if ((iE4 >= maxelm) && (iE4 < maxelm + maxbound)) bE4 = true;
	if ((iN4 >= maxelm) && (iN4 < maxelm + maxbound)) bN4 = true;
	if ((iT4 >= maxelm) && (iT4 < maxelm + maxbound)) bT4 = true;
	if ((iW4 >= maxelm) && (iW4 < maxelm + maxbound)) bW4 = true;
	if ((iS4 >= maxelm) && (iS4 < maxelm + maxbound)) bS4 = true;
	if ((iB4 >= maxelm) && (iB4 < maxelm + maxbound)) bB4 = true;

	// ���������� �������� �������� ������������ ������:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
	//volume3D(iP, nvtx, pa, dx, dy, dz);

	TOCHKA pvol = volume_loc[iP];
	dx = pvol.x;
	dy = pvol.y;
	dz = pvol.z;

	//printf("%.2f %.2f\n",dx,dy); // debug GOOD
	//getchar();

	doublereal dxe = 0.5 * dx, dxw = 0.5 * dx, dyn = 0.5 * dy, dys = 0.5 * dy, dzt = 0.5 * dz, dzb = 0.5 * dz;
	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE)) {
		if (!bE) dxe = 0.5 * (pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW)) {
		if (!bW) dxw = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5 * (pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN)) {
		if (!bN) dyn = 0.5 * (pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS)) {
		if (!bS) dys = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5 * (pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT)) {
		if (!bT) dzt = 0.5 * (pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB)) {
		if (!bB) dzb = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5 * (pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5 * dx, dxw2 = 0.5 * dx, dyn2 = 0.5 * dy, dys2 = 0.5 * dy, dzt2 = 0.5 * dz, dzb2 = 0.5 * dz;
	doublereal dxe3 = 0.5 * dx, dxw3 = 0.5 * dx, dyn3 = 0.5 * dy, dys3 = 0.5 * dy, dzt3 = 0.5 * dz, dzb3 = 0.5 * dz;
	doublereal dxe4 = 0.5 * dx, dxw4 = 0.5 * dx, dyn4 = 0.5 * dy, dys4 = 0.5 * dy, dzt4 = 0.5 * dz, dzb4 = 0.5 * dz;

	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE2)) {
		if (!bE2) dxe2 = 0.5 * (pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW2)) {
		if (!bW2) dxw2 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5 * (pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN2)) {
		if (!bN2) dyn2 = 0.5 * (pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS2)) {
		if (!bS2) dys2 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5 * (pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT2)) {
		if (!bT2) dzt2 = 0.5 * (pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB2)) {
		if (!bB2) dzb2 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5 * (pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE3)) {
		if (!bE3) dxe3 = 0.5 * (pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW3)) {
		if (!bW3) dxw3 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5 * (pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN3)) {
		if (!bN3) dyn3 = 0.5 * (pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS3)) {
		if (!bS3) dys3 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5 * (pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT3)) {
		if (!bT3) dzt3 = 0.5 * (pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB3)) {
		if (!bB3) dzb3 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5 * (pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE4)) {
		if (!bE4) dxe4 = 0.5 * (pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW4)) {
		if (!bW4) dxw4 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5 * (pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN4)) {
		if (!bN4) dyn4 = 0.5 * (pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS4)) {
		if (!bS4) dys4 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5 * (pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT4)) {
		if (!bT4) dzt4 = 0.5 * (pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB4)) {
		if (!bB4) dzb4 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5 * (pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}


	// ���� ��������������� ��������� �����:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5 * dx / dxe;
	fwplus = 0.5 * dx / dxw;
	// y-direction
	fnplus = 0.5 * dy / dyn;
	fsplus = 0.5 * dy / dys;
	// z-direction
	ftplus = 0.5 * dz / dzt;
	fbplus = 0.5 * dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5 * dx / dxe2;
	fwplus2 = 0.5 * dx / dxw2;
	// y-direction
	fnplus2 = 0.5 * dy / dyn2;
	fsplus2 = 0.5 * dy / dys2;
	// z-direction
	ftplus2 = 0.5 * dz / dzt2;
	fbplus2 = 0.5 * dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5 * dx / dxe3;
	fwplus3 = 0.5 * dx / dxw3;
	// y-direction
	fnplus3 = 0.5 * dy / dyn3;
	fsplus3 = 0.5 * dy / dys3;
	// z-direction
	ftplus3 = 0.5 * dz / dzt3;
	fbplus3 = 0.5 * dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5 * dx / dxe4;
	fwplus4 = 0.5 * dx / dxw4;
	// y-direction
	fnplus4 = 0.5 * dy / dyn4;
	fsplus4 = 0.5 * dy / dys4;
	// z-direction
	ftplus4 = 0.5 * dz / dzt4;
	fbplus4 = 0.5 * dz / dzb4;

	//printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	//getchar();

	// ��������� �� ����� �� ���������������� ������� �������������	
	doublereal rP, rE = 0.0, rN = 0.0, rT = 0.0, rW = 0.0, rS = 0.0, rB = 0.0;
	rP = prop[RHO][iP];
	if (iE > -1) {
		if (!bE) rE = prop[RHO][iE]; else rE = prop_b[RHO][iE - maxelm];
	}
	if (iN > -1) {
		if (!bN) rN = prop[RHO][iN]; else rN = prop_b[RHO][iN - maxelm];
	}
	if (iT > -1) {
		if (!bT) rT = prop[RHO][iT]; else rT = prop_b[RHO][iT - maxelm];
	}
	if (iW > -1) {
		if (!bW) rW = prop[RHO][iW]; else rW = prop_b[RHO][iW - maxelm];
	}
	if (iS > -1) {
		if (!bS) rS = prop[RHO][iS]; else rS = prop_b[RHO][iS - maxelm];
	}
	if (iB > -1) {
		if (!bB) rB = prop[RHO][iB]; else rB = prop_b[RHO][iB - maxelm];
	}

	doublereal  rE2 = 0.0, rN2 = 0.0, rT2 = 0.0, rW2 = 0.0, rS2 = 0.0, rB2 = 0.0;
	if (iE2 > -1) {
		if (!bE2) rE2 = prop[RHO][iE2]; else rE2 = prop_b[RHO][iE2 - maxelm];
	}
	if (iN2 > -1) {
		if (!bN2) rN2 = prop[RHO][iN2]; else rN2 = prop_b[RHO][iN2 - maxelm];
	}
	if (iT2 > -1) {
		if (!bT2) rT2 = prop[RHO][iT2]; else rT2 = prop_b[RHO][iT2 - maxelm];
	}
	if (iW2 > -1) {
		if (!bW2) rW2 = prop[RHO][iW2]; else rW2 = prop_b[RHO][iW2 - maxelm];
	}
	if (iS2 > -1) {
		if (!bS2) rS2 = prop[RHO][iS2]; else rS2 = prop_b[RHO][iS2 - maxelm];
	}
	if (iB2 > -1) {
		if (!bB2) rB2 = prop[RHO][iB2]; else rB2 = prop_b[RHO][iB2 - maxelm];
	}

	doublereal  rE3 = 0.0, rN3 = 0.0, rT3 = 0.0, rW3 = 0.0, rS3 = 0.0, rB3 = 0.0;
	if (iE3 > -1) {
		if (!bE3) rE3 = prop[RHO][iE3]; else rE3 = prop_b[RHO][iE3 - maxelm];
	}
	if (iN3 > -1) {
		if (!bN3) rN3 = prop[RHO][iN3]; else rN3 = prop_b[RHO][iN3 - maxelm];
	}
	if (iT3 > -1) {
		if (!bT3) rT3 = prop[RHO][iT3]; else rT3 = prop_b[RHO][iT3 - maxelm];
	}
	if (iW3 > -1) {
		if (!bW3) rW3 = prop[RHO][iW3]; else rW3 = prop_b[RHO][iW3 - maxelm];
	}
	if (iS3 > -1) {
		if (!bS3) rS3 = prop[RHO][iS3]; else rS3 = prop_b[RHO][iS3 - maxelm];
	}
	if (iB3 > -1) {
		if (!bB3) rB3 = prop[RHO][iB3]; else rB3 = prop_b[RHO][iB3 - maxelm];
	}

	doublereal  rE4 = 0.0, rN4 = 0.0, rT4 = 0.0, rW4 = 0.0, rS4 = 0.0, rB4 = 0.0;
	if (iE4 > -1) {
		if (!bE4) rE4 = prop[RHO][iE4]; else rE4 = prop_b[RHO][iE4 - maxelm];
	}
	if (iN4 > -1) {
		if (!bN4) rN4 = prop[RHO][iN4]; else rN4 = prop_b[RHO][iN4 - maxelm];
	}
	if (iT4 > -1) {
		if (!bT4) rT4 = prop[RHO][iT4]; else rT4 = prop_b[RHO][iT4 - maxelm];
	}
	if (iW4 > -1) {
		if (!bW4) rW4 = prop[RHO][iW4]; else rW4 = prop_b[RHO][iW4 - maxelm];
	}
	if (iS4 > -1) {
		if (!bS4) rS4 = prop[RHO][iS4]; else rS4 = prop_b[RHO][iS4 - maxelm];
	}
	if (iB4 > -1) {
		if (!bB4) rB4 = prop[RHO][iB4]; else rB4 = prop_b[RHO][iB4 - maxelm];
	}

	/*
	doublereal rhoe, rhow, rhon, rhos, rhot, rhob;
	// �������� ���������  �� ����� ��:
	rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP);  // ���������.
	rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP);
	rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP);
	rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP);
	rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP);
	rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP);


	doublereal rhoe2, rhow2, rhon2, rhos2, rhot2, rhob2;
	doublereal rhoe3, rhow3, rhon3, rhos3, rhot3, rhob3;
	doublereal rhoe4, rhow4, rhon4, rhos4, rhot4, rhob4;

	rhoe2 = rE2 * rP / (feplus2*rE2 + (1.0 - feplus2)*rP);
	rhow2 = rW2 * rP / (fwplus2*rW2 + (1.0 - fwplus2)*rP);
	rhon2 = rN2 * rP / (fnplus2*rN2 + (1.0 - fnplus2)*rP);
	rhos2 = rS2 * rP / (fsplus2*rS2 + (1.0 - fsplus2)*rP);
	rhot2 = rT2 * rP / (ftplus2*rT2 + (1.0 - ftplus2)*rP);
	rhob2 = rB2 * rP / (fbplus2*rB2 + (1.0 - fbplus2)*rP);

	rhoe3 = rE3 * rP / (feplus3*rE3 + (1.0 - feplus3)*rP);
	rhow3 = rW3 * rP / (fwplus3*rW3 + (1.0 - fwplus3)*rP);
	rhon3 = rN3 * rP / (fnplus3*rN3 + (1.0 - fnplus3)*rP);
	rhos3 = rS3 * rP / (fsplus3*rS3 + (1.0 - fsplus3)*rP);
	rhot3 = rT3 * rP / (ftplus3*rT3 + (1.0 - ftplus3)*rP);
	rhob3 = rB3 * rP / (fbplus3*rB3 + (1.0 - fbplus3)*rP);

	rhoe4 = rE4 * rP / (feplus4*rE4 + (1.0 - feplus4)*rP);
	rhow4 = rW4 * rP / (fwplus4*rW4 + (1.0 - fwplus4)*rP);
	rhon4 = rN4 * rP / (fnplus4*rN4 + (1.0 - fnplus4)*rP);
	rhos4 = rS4 * rP / (fsplus4*rS4 + (1.0 - fsplus4)*rP);
	rhot4 = rT4 * rP / (ftplus4*rT4 + (1.0 - ftplus4)*rP);
	rhob4 = rB4 * rP / (fbplus4*rB4 + (1.0 - fbplus4)*rP);
	*/
	doublereal rhoe = 0.0, rhow = 0.0, rhon = 0.0, rhos = 0.0, rhot = 0.0, rhob = 0.0;
	// ������������ ��������� ������� ���, ����� ����������� 
	// ���������� �����������.
	if (iE > -1) {
		if (!bE) rhoe = rE * rP / (feplus * rE + (1.0 - feplus) * rP); else rhoe = rE; // ��������� !
	}
	if (iW > -1) {
		if (!bW) rhow = rW * rP / (fwplus * rW + (1.0 - fwplus) * rP); else rhow = rW;
	}
	if (iN > -1) {
		if (!bN) rhon = rN * rP / (fnplus * rN + (1.0 - fnplus) * rP); else rhon = rN;
	}
	if (iS > -1) {
		if (!bS) rhos = rS * rP / (fsplus * rS + (1.0 - fsplus) * rP); else rhos = rS;
	}
	if (iT > -1) {
		if (!bT) rhot = rT * rP / (ftplus * rT + (1.0 - ftplus) * rP); else rhot = rT;
	}
	if (iB > -1) {
		if (!bB) rhob = rB * rP / (fbplus * rB + (1.0 - fbplus) * rP); else rhob = rB;
	}

	doublereal rhoe2 = 0.0, rhow2 = 0.0, rhon2 = 0.0, rhos2 = 0.0, rhot2 = 0.0, rhob2 = 0.0;
	doublereal rhoe3 = 0.0, rhow3 = 0.0, rhon3 = 0.0, rhos3 = 0.0, rhot3 = 0.0, rhob3 = 0.0;
	doublereal rhoe4 = 0.0, rhow4 = 0.0, rhon4 = 0.0, rhos4 = 0.0, rhot4 = 0.0, rhob4 = 0.0;

	if (iE2 > -1) {
		if (!bE2)  rhoe2 = rE2 * rP / (feplus2 * rE2 + (1.0 - feplus2) * rP); else rhoe2 = rE2; // ��������� !
	}
	if (iW2 > -1) {
		if (!bW2)  rhow2 = rW2 * rP / (fwplus2 * rW2 + (1.0 - fwplus2) * rP); else rhow2 = rW2;
	}
	if (iN2 > -1) {
		if (!bN2) rhon2 = rN2 * rP / (fnplus2 * rN2 + (1.0 - fnplus2) * rP); else rhon2 = rN2;
	}
	if (iS2 > -1) {
		if (!bS2)  rhos2 = rS2 * rP / (fsplus2 * rS2 + (1.0 - fsplus2) * rP); else rhos2 = rS2;
	}
	if (iT2 > -1) {
		if (!bT2)  rhot2 = rT2 * rP / (ftplus2 * rT2 + (1.0 - ftplus2) * rP); else rhot2 = rT2;
	}
	if (iB2 > -1) {
		if (!bB2) rhob2 = rB2 * rP / (fbplus2 * rB2 + (1.0 - fbplus2) * rP); else rhob2 = rB2;
	}

	if (iE3 > -1) {
		if (!bE3) rhoe3 = rE3 * rP / (feplus3 * rE3 + (1.0 - feplus3) * rP); else rhoe3 = rE3;
	}
	if (iW3 > -1) {
		if (!bW3) rhow3 = rW3 * rP / (fwplus3 * rW3 + (1.0 - fwplus3) * rP); else rhow3 = rW3;
	}
	if (iN3 > -1) {
		if (!bN3) rhon3 = rN3 * rP / (fnplus3 * rN3 + (1.0 - fnplus3) * rP); else rhon3 = rN3;
	}
	if (iS3 > -1) {
		if (!bS3) rhos3 = rS3 * rP / (fsplus3 * rS3 + (1.0 - fsplus3) * rP); else rhos3 = rS3;
	}
	if (iT3 > -1) {
		if (!bT3) rhot3 = rT3 * rP / (ftplus3 * rT3 + (1.0 - ftplus3) * rP); else rhot3 = rT3;
	}
	if (iB3 > -1) {
		if (!bB3) rhob3 = rB3 * rP / (fbplus3 * rB3 + (1.0 - fbplus3) * rP); else rhob3 = rB3;
	}

	if (iE4 > -1) {
		if (!bE4) rhoe4 = rE4 * rP / (feplus4 * rE4 + (1.0 - feplus4) * rP); else rhoe4 = rE4;
	}
	if (iW4 > -1) {
		if (!bW4) rhow4 = rW4 * rP / (fwplus4 * rW4 + (1.0 - fwplus4) * rP); else rhow4 = rW4;
	}
	if (iN4 > -1) {
		if (!bN4) rhon4 = rN4 * rP / (fnplus4 * rN4 + (1.0 - fnplus4) * rP); else rhon4 = rN4;
	}
	if (iS4 > -1) {
		if (!bS4) rhos4 = rS4 * rP / (fsplus4 * rS4 + (1.0 - fsplus4) * rP); else rhos4 = rS4;
	}
	if (iT4 > -1) {
		if (!bT4) rhot4 = rT4 * rP / (ftplus4 * rT4 + (1.0 - ftplus4) * rP); else rhot4 = rT4;
	}
	if (iB4 > -1) {
		if (!bB4) rhob4 = rB4 * rP / (fbplus4 * rB4 + (1.0 - fbplus4) * rP); else rhob4 = rB4;
	}



	/*
	   ����������� ����������:
	   �� �������� ������������ ���-��� �������� �� ����� ��
	   ������ ���� ��������� � � ���������� ��� ��������� ��������.
	   ��� ������������� ������ ���� ��������� ��� ���������� ������� �� ����� ������������ ������
	   � ���������� ���������� ��������, ���������������� � ������������ ������������� �.�. � ����������
	   ����������� ��� ������ ����� ������������� ��������� �������������. ���� ��� ���������� ������� ��
	   ��������� ������������ ���-��� �� ���� ��������� ���������� � �� ��������� (���� ��� � ������������
	   ������ � ��������� ���������� �������� � ����������������. ��� �������� �������� ��������� �����������
	   ������� ����� ��������� ��������� ����������) �� ������ ����� �� ����� ������������� ��������� �������������
	   � ������������� ������� ����� ��������.
	   ������������ ���������� ������������ �������� �� ��� ��� ������������, �
	   �� ����������� ������ ���. ��� ����������� ����� ��������� ������������� �������� ���� ��� �� ������
	   ����������������� �������� � ��������. ��� ���������� ��������� ������������ ������������ �
	   ���������� ��� ��������� ��������. ��� ������� � ������� �������� ���������
	   SIMPLE (�� ������ ����������������� ���������� ������� ������������ ������� � ������� ��������, ��
	   ���� � ��� ��� ����� �� �� ��������� �������� ���������� ����������� ����� ������ ����� (������� ���� ��������� � ������)
	   � ������� ������� ��� ������������ ������������ �������� � ���������� ��������).
	   ��������� ������������� ��������...
	   �������� ������ �������������� ������ ��������, �.�. �� ��� ������������ ������������
	   ��� ���� ����� ��� �� ���������. ������� ������������ �������� ������������ ���-��� ������
	   �� ������ �������� ��������� SIMPLE. �� ������ �������� ������������� ������� ������������
	   �������� ����� ���������� ������� �������� ������������� ��������.
	*/


	// ������������ ����� ����� ����� ��.
	// � ���������� �������� � ������ ���������������� ��������������� 
	// ���������� ���������������� �������� ���-���.
	doublereal Fe = 0.0, Fw = 0.0, Fn = 0.0, Fs = 0.0, Ft = 0.0, Fb = 0.0;


	// ��� ���� �����.
	doublereal Fe1 = 0.0, Fe2 = 0.0, Fe3 = 0.0, Fe4 = 0.0;
	doublereal Fw1 = 0.0, Fw2 = 0.0, Fw3 = 0.0, Fw4 = 0.0;
	doublereal Fn1 = 0.0, Fn2 = 0.0, Fn3 = 0.0, Fn4 = 0.0;
	doublereal Fs1 = 0.0, Fs2 = 0.0, Fs3 = 0.0, Fs4 = 0.0;
	doublereal Ft1 = 0.0, Ft2 = 0.0, Ft3 = 0.0, Ft4 = 0.0;
	doublereal Fb1 = 0.0, Fb2 = 0.0, Fb3 = 0.0, Fb4 = 0.0;

	// �������� ����� ������������ � �� �����������, 
	// ��� ������ ����������� ��������� �� �������� ������ ������� ����.
	if (!b_on_adaptive_local_refinement_mesh) {
		if (Fe != Fe) {
			printf("Fe=%e\n", Fe);
		}
		Fe = mf[iP][E_SIDE];
		if (Fe != Fe) {
			printf("Fe=%e\n", Fe);
			system("pause");
		}
		Fn = mf[iP][N_SIDE];
		Ft = mf[iP][T_SIDE];
		Fw = mf[iP][W_SIDE];
		Fs = mf[iP][S_SIDE];
		Fb = mf[iP][B_SIDE];
	}
	else {
		// ����� �� ����. 24.11.2018

		if (iE > -1) {
			if (bE) {
				// ��������� ����.
				Fe = mf[iP][E_SIDE] * (border_neighbor[iE - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
					Fe = mf[iP][E_SIDE];
				}
				else {

					Fe = mf[iE][W_SIDE];

				}
			}
		}

		if (iW > -1) {
			if (bW) {
				// ��������� ����.
				Fw = mf[iP][W_SIDE] * (border_neighbor[iW - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
					Fw = mf[iP][W_SIDE];
				}
				else {

					Fw = mf[iW][E_SIDE];

				}
			}
		}

		if (iN > -1) {
			if (bN) {
				// ��������� ����.
				Fn = mf[iP][N_SIDE] * (border_neighbor[iN - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
					Fn = mf[iP][N_SIDE];
				}
				else {

					Fn = mf[iN][S_SIDE];

				}
			}
		}

		if (iS > -1) {
			if (bS) {
				// ��������� ����.
				Fs = mf[iP][S_SIDE] * (border_neighbor[iS - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
					Fs = mf[iP][S_SIDE];
				}
				else {

					Fs = mf[iS][N_SIDE];

				}
			}
		}

		if (iT > -1) {
			if (bT) {
				// ��������� ����.
				Ft = mf[iP][T_SIDE] * (border_neighbor[iT - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
					Ft = mf[iP][T_SIDE];
				}
				else {

					Ft = mf[iT][B_SIDE];

				}
			}
		}

		if (iB > -1) {
			if (bB) {
				// ��������� ����.
				Fb = mf[iP][B_SIDE] * (border_neighbor[iB - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
					Fb = mf[iP][B_SIDE];
				}
				else {

					Fb = mf[iB][T_SIDE];

				}
			}
		}

		if (iE2 > -1) {
			if (bE2) {
				// ��������� ����.
				Fe2 = mf[iP][E_SIDE] * (border_neighbor[iE2 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
					Fe2 = mf[iP][E_SIDE];
				}
				else {

					Fe2 = mf[iE2][W_SIDE];

				}
			}
		}

		if (iW2 > -1) {
			if (bW2) {
				// ��������� ����.
				Fw2 = mf[iP][W_SIDE] * (border_neighbor[iW2 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW2]]) {
					Fw = mf[iP][W_SIDE];
				}
				else {

					Fw = mf[iW2][E_SIDE];

				}
			}
		}

		if (iN2 > -1) {
			if (bN2) {
				// ��������� ����.
				Fn2 = mf[iP][N_SIDE] * (border_neighbor[iN2 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
					Fn2 = mf[iP][N_SIDE];
				}
				else {

					Fn2 = mf[iN2][S_SIDE];

				}
			}
		}

		if (iS2 > -1) {
			if (bS2) {
				// ��������� ����.
				Fs2 = mf[iP][S_SIDE] * (border_neighbor[iS2 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
					Fs2 = mf[iP][S_SIDE];
				}
				else {

					Fs2 = mf[iS2][N_SIDE];

				}
			}
		}

		if (iT2 > -1) {
			if (bT2) {
				// ��������� ����.
				Ft2 = mf[iP][T_SIDE] * (border_neighbor[iT2 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
					Ft2 = mf[iP][T_SIDE];
				}
				else {

					Ft2 = mf[iT2][B_SIDE];

				}
			}
		}

		if (iB2 > -1) {
			if (bB2) {
				// ��������� ����.
				Fb2 = mf[iP][B_SIDE] * (border_neighbor[iB2 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
					Fb2 = mf[iP][B_SIDE];
				}
				else {

					Fb2 = mf[iB2][T_SIDE];

				}
			}
		}


		if (iE3 > -1) {
			if (bE3) {
				// ��������� ����.
				Fe3 = mf[iP][E_SIDE] * (border_neighbor[iE3 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
					Fe3 = mf[iP][E_SIDE];
				}
				else {

					Fe3 = mf[iE3][W_SIDE];

				}
			}
		}

		if (iW3 > -1) {
			if (bW3) {
				// ��������� ����.
				Fw3 = mf[iP][W_SIDE] * (border_neighbor[iW3 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
					Fw3 = mf[iP][W_SIDE];
				}
				else {

					Fw3 = mf[iW3][E_SIDE];

				}
			}
		}

		if (iN3 > -1) {
			if (bN3) {
				// ��������� ����.
				Fn3 = mf[iP][N_SIDE] * (border_neighbor[iN3 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
					Fn3 = mf[iP][N_SIDE];
				}
				else {

					Fn3 = mf[iN3][S_SIDE];

				}
			}
		}

		if (iS3 > -1) {
			if (bS3) {
				// ��������� ����.
				Fs3 = mf[iP][S_SIDE] * (border_neighbor[iS3 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
					Fs3 = mf[iP][S_SIDE];
				}
				else {

					Fs3 = mf[iS3][N_SIDE];

				}
			}
		}

		if (iT3 > -1) {
			if (bT3) {
				// ��������� ����.
				Ft3 = mf[iP][T_SIDE] * (border_neighbor[iT3 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
					Ft3 = mf[iP][T_SIDE];
				}
				else {

					Ft3 = mf[iT3][B_SIDE];

				}
			}
		}

		if (iB3 > -1) {
			if (bB3) {
				// ��������� ����.
				Fb3 = mf[iP][B_SIDE] * (border_neighbor[iB3 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
					Fb3 = mf[iP][B_SIDE];
				}
				else {

					Fb3 = mf[iB3][T_SIDE];

				}
			}
		}

		if (iE4 > -1) {
			if (bE4) {
				// ��������� ����.
				Fe4 = mf[iP][E_SIDE] * (border_neighbor[iE4 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
					Fe4 = mf[iP][E_SIDE];
				}
				else {

					Fe4 = mf[iE4][W_SIDE];

				}
			}
		}

		if (iW4 > -1) {
			if (bW4) {
				// ��������� ����.
				Fw4 = mf[iP][W_SIDE] * (border_neighbor[iW4 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
					Fw4 = mf[iP][W_SIDE];
				}
				else {

					Fw4 = mf[iW4][E_SIDE];

				}
			}
		}

		if (iN4 > -1) {
			if (bN4) {
				// ��������� ����.
				Fn4 = mf[iP][N_SIDE] * (border_neighbor[iN4 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
					Fn4 = mf[iP][N_SIDE];
				}
				else {

					Fn4 = mf[iN4][S_SIDE];

				}
			}
		}

		if (iS4 > -1) {
			if (bS4) {
				// ��������� ����.
				Fs4 = mf[iP][S_SIDE] * (border_neighbor[iS4 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
					Fs4 = mf[iP][S_SIDE];
				}
				else {

					Fs4 = mf[iS4][N_SIDE];

				}
			}
		}

		if (iT4 > -1) {
			if (bT4) {
				// ��������� ����.
				Ft4 = mf[iP][T_SIDE] * (border_neighbor[iT4 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
					Ft4 = mf[iP][T_SIDE];
				}
				else {

					Ft4 = mf[iT4][B_SIDE];

				}
			}
		}

		if (iB4 > -1) {
			if (bB4) {
				// ��������� ����.
				Fb4 = mf[iP][B_SIDE] * (border_neighbor[iB4 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
					Fb4 = mf[iP][B_SIDE];
				}
				else {

					Fb4 = mf[iB4][T_SIDE];

				}
			}
		}



	}



	doublereal eps = 1e-37; // ��� ��������� ������������� ����.

	//eqin.fluidinfo[0].sigma_nu
	// ������������ ��������:
	doublereal  GP, GE, GW, GN, GS, GT, GB;
	doublereal  GE2, GW2, GN2, GS2, GT2, GB2;
	doublereal  GE3, GW3, GN3, GS3, GT3, GB3;
	doublereal  GE4, GW4, GN4, GS4, GT4, GB4;

	
	doublereal sigmaf = eqin.fluidinfo[0].sigmaf;

	doublereal sigmafE = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafW = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafN = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafS = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafT = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafB = eqin.fluidinfo[0].sigmaf;

	doublereal sigmafE2 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafW2 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafN2 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafS2 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafT2 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafB2 = eqin.fluidinfo[0].sigmaf;

	doublereal sigmafE3 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafW3 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafN3 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafS3 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafT3 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafB3 = eqin.fluidinfo[0].sigmaf;

	doublereal sigmafE4 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafW4 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafN4 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafS4 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafT4 = eqin.fluidinfo[0].sigmaf;
	doublereal sigmafB4 = eqin.fluidinfo[0].sigmaf;



	// ���������� ������������ ��������:
	GP = ((prop[MU_DYNAMIC_VISCOSITY][iP]) + fmax(0.0,   potent[MUT][iP] / sigmaf)); // � ������ ����������� ��.
	if (iE > -1) {
		if (!bE) GE = ((prop[MU_DYNAMIC_VISCOSITY][iE]) + fmax(0.0, potent[MUT][iE] / sigmafE)); else GE = ((prop_b[MU_DYNAMIC_VISCOSITY][iE - maxelm]) + fmax(0.0, potent[MUT][iE] / sigmaf));
	}
	if (iN > -1) {
		if (!bN) GN = ((prop[MU_DYNAMIC_VISCOSITY][iN]) + fmax(0.0, potent[MUT][iN] / sigmafN)); else GN = ((prop_b[MU_DYNAMIC_VISCOSITY][iN - maxelm]) + fmax(0.0,  potent[MUT][iN] / sigmaf));
	}
	if (iT > -1) {
		if (!bT) GT = ((prop[MU_DYNAMIC_VISCOSITY][iT]) + fmax(0.0,  potent[MUT][iT] / sigmafT)); else GT = ((prop_b[MU_DYNAMIC_VISCOSITY][iT - maxelm]) + fmax(0.0, potent[MUT][iT] / sigmaf));
	}
	if (iW > -1) {
		if (!bW) GW = ((prop[MU_DYNAMIC_VISCOSITY][iW]) + fmax(0.0,  potent[MUT][iW] / sigmafW)); else GW = ((prop_b[MU_DYNAMIC_VISCOSITY][iW - maxelm]) + fmax(0.0, potent[MUT][iW] / sigmaf));
	}
	if (iS > -1) {
		if (!bS) GS = ((prop[MU_DYNAMIC_VISCOSITY][iS]) + fmax(0.0,  potent[MUT][iS] / sigmafS)); else GS = ((prop_b[MU_DYNAMIC_VISCOSITY][iS - maxelm]) + fmax(0.0, potent[MUT][iS] / sigmaf));
	}
	if (iB > -1) {
		if (!bB) GB = ((prop[MU_DYNAMIC_VISCOSITY][iB]) + fmax(0.0,  potent[MUT][iB] / sigmafB)); else GB = ((prop_b[MU_DYNAMIC_VISCOSITY][iB - maxelm]) + fmax(0.0, potent[MUT][iB] / sigmaf));
	}

	if (iE2 > -1) {
		if (!bE2) GE2 = ((prop[MU_DYNAMIC_VISCOSITY][iE2]) + fmax(0.0,  potent[MUT][iE2] / sigmafE2)); else GE2 = ((prop_b[MU_DYNAMIC_VISCOSITY][iE2 - maxelm]) + fmax(0.0,  potent[MUT][iE2] / sigmaf));
	}
	if (iN2 > -1) {
		if (!bN2) GN2 = ((prop[MU_DYNAMIC_VISCOSITY][iN2]) + fmax(0.0, potent[MUT][iN2] / sigmafN2)); else GN2 = ((prop_b[MU_DYNAMIC_VISCOSITY][iN2 - maxelm]) + fmax(0.0,  potent[MUT][iN2] / sigmaf));
	}
	if (iT2 > -1) {
		if (!bT2) GT2 = ((prop[MU_DYNAMIC_VISCOSITY][iT2]) + fmax(0.0,  potent[MUT][iT2] / sigmafT2)); else GT2 = ((prop_b[MU_DYNAMIC_VISCOSITY][iT2 - maxelm]) + fmax(0.0,  potent[MUT][iT2] / sigmaf));
	}
	if (iW2 > -1) {
		if (!bW2) GW2 = ((prop[MU_DYNAMIC_VISCOSITY][iW2]) + fmax(0.0,  potent[MUT][iW2] / sigmafW2)); else GW2 = ((prop_b[MU_DYNAMIC_VISCOSITY][iW2 - maxelm]) + fmax(0.0,  potent[MUT][iW2] / sigmaf));
	}
	if (iS2 > -1) {
		if (!bS2) GS2 = ((prop[MU_DYNAMIC_VISCOSITY][iS2]) + fmax(0.0, potent[MUT][iS2] / sigmafS2)); else GS2 = ((prop_b[MU_DYNAMIC_VISCOSITY][iS2 - maxelm]) + fmax(0.0,  potent[MUT][iS2] / sigmaf));
	}
	if (iB2 > -1) {
		if (!bB2) GB2 = ((prop[MU_DYNAMIC_VISCOSITY][iB2]) + fmax(0.0,  potent[MUT][iB2] / sigmafB2)); else GB2 = ((prop_b[MU_DYNAMIC_VISCOSITY][iB2 - maxelm]) + fmax(0.0,  potent[MUT][iB2] / sigmaf));
	}

	if (iE3 > -1) {
		if (!bE3) GE3 = ((prop[MU_DYNAMIC_VISCOSITY][iE3]) + fmax(0.0,  potent[MUT][iE3] / sigmafE3)); else GE3 = ((prop_b[MU_DYNAMIC_VISCOSITY][iE3 - maxelm]) + fmax(0.0,  potent[MUT][iE3] / sigmaf));
	}
	if (iN3 > -1) {
		if (!bN3) GN3 = ((prop[MU_DYNAMIC_VISCOSITY][iN3]) + fmax(0.0,  potent[MUT][iN3] / sigmafN3)); else GN3 = ((prop_b[MU_DYNAMIC_VISCOSITY][iN3 - maxelm]) + fmax(0.0,  potent[MUT][iN3] / sigmaf));
	}
	if (iT3 > -1) {
		if (!bT3) GT3 = ((prop[MU_DYNAMIC_VISCOSITY][iT3]) + fmax(0.0,  potent[MUT][iT3] / sigmafT3)); else GT3 = ((prop_b[MU_DYNAMIC_VISCOSITY][iT3 - maxelm]) + fmax(0.0,  potent[MUT][iT3] / sigmaf));
	}
	if (iW3 > -1) {
		if (!bW3) GW3 = ((prop[MU_DYNAMIC_VISCOSITY][iW3]) + fmax(0.0,  potent[MUT][iW3] / sigmafW3)); else GW3 = ((prop_b[MU_DYNAMIC_VISCOSITY][iW3 - maxelm]) + fmax(0.0,  potent[MUT][iW3] / sigmaf));
	}
	if (iS3 > -1) {
		if (!bS3) GS3 = ((prop[MU_DYNAMIC_VISCOSITY][iS3]) + fmax(0.0,  potent[MUT][iS3] / sigmafS3)); else GS3 = ((prop_b[MU_DYNAMIC_VISCOSITY][iS3 - maxelm]) + fmax(0.0, potent[MUT][iS3] / sigmaf));
	}
	if (iB3 > -1) {
		if (!bB3) GB3 = ((prop[MU_DYNAMIC_VISCOSITY][iB3]) + fmax(0.0,  potent[MUT][iB3] / sigmafB3)); else GB3 = ((prop_b[MU_DYNAMIC_VISCOSITY][iB3 - maxelm]) + fmax(0.0, potent[MUT][iB3] / sigmaf));
	}

	if (iE4 > -1) {
		if (!bE4) GE4 = ((prop[MU_DYNAMIC_VISCOSITY][iE4]) + fmax(0.0,  potent[MUT][iE4] / sigmafE4)); else GE4 = ((prop_b[MU_DYNAMIC_VISCOSITY][iE4 - maxelm]) + fmax(0.0,  potent[MUT][iE4] / sigmaf));
	}
	if (iN4 > -1) {
		if (!bN4) GN4 = ((prop[MU_DYNAMIC_VISCOSITY][iN4]) + fmax(0.0,  potent[MUT][iN4] / sigmafN4)); else GN4 = ((prop_b[MU_DYNAMIC_VISCOSITY][iN4 - maxelm]) + fmax(0.0,  potent[MUT][iN4] / sigmaf));
	}
	if (iT4 > -1) {
		if (!bT4) GT4 = ((prop[MU_DYNAMIC_VISCOSITY][iT4]) + fmax(0.0,  potent[MUT][iT4] / sigmafT4)); else GT4 = ((prop_b[MU_DYNAMIC_VISCOSITY][iT4 - maxelm]) + fmax(0.0,  potent[MUT][iT4] / sigmaf));
	}
	if (iW4 > -1) {
		if (!bW4) GW4 = ((prop[MU_DYNAMIC_VISCOSITY][iW4]) + fmax(0.0,  potent[MUT][iW4] / sigmafW4)); else GW4 = ((prop_b[MU_DYNAMIC_VISCOSITY][iW4 - maxelm]) + fmax(0.0,  potent[MUT][iW4] / sigmaf));
	}
	if (iS4 > -1) {
		if (!bS4) GS4 = ((prop[MU_DYNAMIC_VISCOSITY][iS4]) + fmax(0.0,  potent[MUT][iS4] / sigmafS4)); else GS4 = ((prop_b[MU_DYNAMIC_VISCOSITY][iS4 - maxelm]) + fmax(0.0,  potent[MUT][iS4] / sigmaf));
	}
	if (iB4 > -1) {
		if (!bB4) GB4 = ((prop[MU_DYNAMIC_VISCOSITY][iB4]) + fmax(0.0,  potent[MUT][iB4] / sigmafB4)); else GB4 = ((prop_b[MU_DYNAMIC_VISCOSITY][iB4 - maxelm]) + fmax(0.0,  potent[MUT][iB4] / sigmaf));
	}

	doublereal Ge = GP, Gw = GP, Gn = GP, Gs = GP, Gt = GP, Gb = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE > -1) {
		Ge = GE * GP / (feplus * GE + (1 - feplus) * GP); // ���������.
	}
	if (iW > -1) {
		Gw = GW * GP / (fwplus * GW + (1 - fwplus) * GP);
	}
	if (iN > -1) {
		Gn = GN * GP / (fnplus * GN + (1 - fnplus) * GP);
	}
	if (iS > -1) {
		Gs = GS * GP / (fsplus * GS + (1 - fsplus) * GP);
	}
	if (iT > -1) {
		Gt = GT * GP / (ftplus * GT + (1 - ftplus) * GP);
	}
	if (iB > -1) {
		Gb = GB * GP / (fbplus * GB + (1 - fbplus) * GP);
	}

	doublereal Ge2 = GP, Gw2 = GP, Gn2 = GP, Gs2 = GP, Gt2 = GP, Gb2 = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE2 > -1) {
		Ge2 = GE2 * GP / (feplus2 * GE2 + (1 - feplus2) * GP); // ���������.
	}
	if (iW2 > -1) {
		Gw2 = GW2 * GP / (fwplus2 * GW2 + (1 - fwplus2) * GP);
	}
	if (iN2 > -1) {
		Gn2 = GN2 * GP / (fnplus2 * GN2 + (1 - fnplus2) * GP);
	}
	if (iS2 > -1) {
		Gs2 = GS2 * GP / (fsplus2 * GS2 + (1 - fsplus2) * GP);
	}
	if (iT2 > -1) {
		Gt2 = GT2 * GP / (ftplus2 * GT2 + (1 - ftplus2) * GP);
	}
	if (iB2 > -1) {
		Gb2 = GB2 * GP / (fbplus2 * GB2 + (1 - fbplus2) * GP);
	}


	doublereal Ge3 = GP, Gw3 = GP, Gn3 = GP, Gs3 = GP, Gt3 = GP, Gb3 = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE3 > -1) {
		Ge3 = GE3 * GP / (feplus3 * GE3 + (1 - feplus3) * GP); // ���������.
	}
	if (iW3 > -1) {
		Gw3 = GW3 * GP / (fwplus3 * GW3 + (1 - fwplus3) * GP);
	}
	if (iN3 > -1) {
		Gn3 = GN3 * GP / (fnplus3 * GN3 + (1 - fnplus3) * GP);
	}
	if (iS3 > -1) {
		Gs3 = GS3 * GP / (fsplus3 * GS3 + (1 - fsplus3) * GP);
	}
	if (iT3 > -1) {
		Gt3 = GT3 * GP / (ftplus3 * GT3 + (1 - ftplus3) * GP);
	}
	if (iB3 > -1) {
		Gb3 = GB3 * GP / (fbplus3 * GB3 + (1 - fbplus3) * GP);
	}

	doublereal Ge4 = GP, Gw4 = GP, Gn4 = GP, Gs4 = GP, Gt4 = GP, Gb4 = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE4 > -1) {
		Ge4 = GE4 * GP / (feplus4 * GE4 + (1 - feplus4) * GP); // ���������.
	}
	if (iW4 > -1) {
		Gw4 = GW4 * GP / (fwplus4 * GW4 + (1 - fwplus4) * GP);
	}
	if (iN4 > -1) {
		Gn4 = GN4 * GP / (fnplus4 * GN4 + (1 - fnplus4) * GP);
	}
	if (iS4 > -1) {
		Gs4 = GS4 * GP / (fsplus4 * GS4 + (1 - fsplus4) * GP);
	}
	if (iT4 > -1) {
		Gt4 = GT4 * GP / (ftplus4 * GT4 + (1 - ftplus4) * GP);
	}
	if (iB4 > -1) {
		Gb4 = GB4 * GP / (fbplus4 * GB4 + (1 - fbplus4) * GP);
	}




	const doublereal ZeroDiffusion = 0.0;// 1.0e-30;
	// ������������ ������������ ������:
	doublereal De = ZeroDiffusion, Dw = ZeroDiffusion, Dn = ZeroDiffusion, Ds = ZeroDiffusion, Dt = ZeroDiffusion, Db = ZeroDiffusion; // �������������
	if (iE > -1) {
		if (bE) {
			// ��������� ����.
			De = Ge * border_neighbor[iE - maxelm].dS / dxe;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
				De = Ge * dy * dz / dxe;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				De = Ge * dy_loc * dz_loc / dxe;
			}
		}

	}
	if (iW > -1) {
		if (bW) {
			// ��������� ����
			Dw = Gw * border_neighbor[iW - maxelm].dS / dxw;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				Dw = Gw * dy * dz / dxw;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw = Gw * dy_loc * dz_loc / dxw;
			}
		}

	}
	if (iN > -1) {
		if (bN) {
			// ��������� ����.
			Dn = Gn * border_neighbor[iN - maxelm].dS / dyn;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
				Dn = Gn * dx * dz / dyn;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn = Gn * dx_loc * dz_loc / dyn;
			}
		}
	}
	if (iS > -1) {
		if (bS) {
			// ��������� ����
			Ds = Gs * border_neighbor[iS - maxelm].dS / dys;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
				Ds = Gs * dx * dz / dys;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds = Gs * dx_loc * dz_loc / dys;
			}
		}
	}
	if (iT > -1) {
		if (bT) {
			// ��������� ����.
			Dt = Gt * border_neighbor[iT - maxelm].dS / dzt;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
				Dt = Gt * dx * dy / dzt;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt = Gt * dx_loc * dy_loc / dzt;
			}
		}


	}
	if (iB > -1) {
		if (bB) {
			// ��������� ����
			Db = Gb * border_neighbor[iB - maxelm].dS / dzb;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
				Db = Gb * dx * dy / dzb;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;


				Db = Gb * dx_loc * dy_loc / dzb;
			}
		}


	}

	// ������������ ������������ ������:
	doublereal De2 = ZeroDiffusion, Dw2 = ZeroDiffusion, Dn2 = ZeroDiffusion, Ds2 = ZeroDiffusion, Dt2 = ZeroDiffusion, Db2 = ZeroDiffusion; // �������������
	if (iE2 > -1) {
		if (bE2) {
			// ��������� ����.
			De2 = Ge2 * border_neighbor[iE2 - maxelm].dS / dxe2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
				De2 = Ge2 * dy * dz / dxe2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;


				De2 = Ge2 * dy_loc * dz_loc / dxe2;
			}
		}

	}
	if (iW2 > -1) {
		if (bW2) {
			// ��������� ����
			Dw2 = Gw2 * border_neighbor[iW2 - maxelm].dS / dxw2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW2]]) {
				Dw2 = Gw2 * dy * dz / dxw2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw2 = Gw2 * dy_loc * dz_loc / dxw2;
			}
		}

	}
	if (iN2 > -1) {
		if (bN2) {
			// ��������� ����.
			Dn2 = Gn2 * border_neighbor[iN2 - maxelm].dS / dyn2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
				Dn2 = Gn2 * dx * dz / dyn2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn2 = Gn2 * dx_loc * dz_loc / dyn2;
			}
		}
	}
	if (iS2 > -1) {
		if (bS2) {
			// ��������� ����
			Ds2 = Gs2 * border_neighbor[iS2 - maxelm].dS / dys2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
				Ds2 = Gs2 * dx * dz / dys2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds2 = Gs2 * dx_loc * dz_loc / dys2;
			}
		}
	}
	if (iT2 > -1) {
		if (bT2) {
			// ��������� ����.
			Dt2 = Gt2 * border_neighbor[iT2 - maxelm].dS / dzt2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
				Dt2 = Gt2 * dx * dy / dzt2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt2 = Gt2 * dx_loc * dy_loc / dzt2;
			}
		}


	}
	if (iB2 > -1) {
		if (bB2) {
			// ��������� ����
			Db2 = Gb2 * border_neighbor[iB2 - maxelm].dS / dzb2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
				Db2 = Gb2 * dx * dy / dzb2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Db2 = Gb2 * dx_loc * dy_loc / dzb2;
			}
		}


	}

	// ������������ ������������ ������:
	doublereal De3 = ZeroDiffusion, Dw3 = ZeroDiffusion, Dn3 = ZeroDiffusion, Ds3 = ZeroDiffusion, Dt3 = ZeroDiffusion, Db3 = ZeroDiffusion; // �������������
	if (iE3 > -1) {
		if (bE3) {
			// ��������� ����.
			De3 = Ge3 * border_neighbor[iE3 - maxelm].dS / dxe3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
				De3 = Ge3 * dy * dz / dxe3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				De3 = Ge3 * dy_loc * dz_loc / dxe3;
			}
		}

	}
	if (iW3 > -1) {
		if (bW3) {
			// ��������� ����
			Dw3 = Gw3 * border_neighbor[iW3 - maxelm].dS / dxw3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
				Dw3 = Gw3 * dy * dz / dxw3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw3 = Gw3 * dy_loc * dz_loc / dxw3;
			}
		}

	}
	if (iN3 > -1) {
		if (bN3) {
			// ��������� ����.
			Dn3 = Gn3 * border_neighbor[iN3 - maxelm].dS / dyn3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
				Dn3 = Gn3 * dx * dz / dyn3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn3 = Gn3 * dx_loc * dz_loc / dyn3;
			}
		}
	}
	if (iS3 > -1) {
		if (bS3) {
			// ��������� ����
			Ds3 = Gs3 * border_neighbor[iS3 - maxelm].dS / dys3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
				Ds3 = Gs3 * dx * dz / dys3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds3 = Gs3 * dx_loc * dz_loc / dys3;
			}
		}
	}
	if (iT3 > -1) {
		if (bT3) {
			// ��������� ����.
			Dt3 = Gt3 * border_neighbor[iT3 - maxelm].dS / dzt3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
				Dt3 = Gt3 * dx * dy / dzt3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt3 = Gt3 * dx_loc * dy_loc / dzt3;
			}
		}


	}
	if (iB3 > -1) {
		if (bB3) {
			// ��������� ����
			Db3 = Gb3 * border_neighbor[iB3 - maxelm].dS / dzb3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
				Db3 = Gb3 * dx * dy / dzb3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Db3 = Gb3 * dx_loc * dy_loc / dzb3;
			}
		}


	}

	// ������������ ������������ ������:
	doublereal De4 = ZeroDiffusion, Dw4 = ZeroDiffusion, Dn4 = ZeroDiffusion, Ds4 = ZeroDiffusion, Dt4 = ZeroDiffusion, Db4 = ZeroDiffusion; // �������������
	if (iE4 > -1) {
		if (bE4) {
			// ��������� ����.
			De4 = Ge4 * border_neighbor[iE4 - maxelm].dS / dxe4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
				De4 = Ge4 * dy * dz / dxe4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;


				De4 = Ge4 * dy_loc * dz_loc / dxe4;
			}
		}

	}
	if (iW4 > -1) {
		if (bW4) {
			// ��������� ����
			Dw4 = Gw4 * border_neighbor[iW4 - maxelm].dS / dxw4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
				Dw4 = Gw4 * dy * dz / dxw4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw4 = Gw4 * dy_loc * dz_loc / dxw4;
			}
		}

	}
	if (iN4 > -1) {
		if (bN4) {
			// ��������� ����.
			Dn4 = Gn4 * border_neighbor[iN4 - maxelm].dS / dyn4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
				Dn4 = Gn4 * dx * dz / dyn4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn4 = Gn4 * dx_loc * dz_loc / dyn4;
			}
		}
	}
	if (iS4 > -1) {
		if (bS4) {
			// ��������� ����
			Ds4 = Gs4 * border_neighbor[iS4 - maxelm].dS / dys4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
				Ds4 = Gs4 * dx * dz / dys4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds4 = Gs4 * dx_loc * dz_loc / dys4;
			}
		}
	}
	if (iT4 > -1) {
		if (bT4) {
			// ��������� ����.
			Dt4 = Gt4 * border_neighbor[iT4 - maxelm].dS / dzt4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
				Dt4 = Gt4 * dx * dy / dzt4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt4 = Gt4 * dx_loc * dy_loc / dzt4;
			}
		}


	}
	if (iB4 > -1) {
		if (bB4) {
			// ��������� ����
			Db4 = Gb4 * border_neighbor[iB4 - maxelm].dS / dzb4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
				Db4 = Gb4 * dx * dy / dzb4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Db4 = Gb4 * dx_loc * dy_loc / dzb4;
			}
		}
	}

	// ������������� ��������� ������� ���������.
	De = fmax(De, 0.0);
	Dw = fmax(Dw, 0.0);
	Dn = fmax(Dn, 0.0);
	Ds = fmax(Ds, 0.0);
	Dt = fmax(Dt, 0.0);
	Db = fmax(Db, 0.0);

	De2 = fmax(De2, 0.0);
	Dw2 = fmax(Dw2, 0.0);
	Dn2 = fmax(Dn2, 0.0);
	Ds2 = fmax(Ds2, 0.0);
	Dt2 = fmax(Dt2, 0.0);
	Db2 = fmax(Db2, 0.0);

	De3 = fmax(De3, 0.0);
	Dw3 = fmax(Dw3, 0.0);
	Dn3 = fmax(Dn3, 0.0);
	Ds3 = fmax(Ds3, 0.0);
	Dt3 = fmax(Dt3, 0.0);
	Db3 = fmax(Db3, 0.0);

	De4 = fmax(De4, 0.0);
	Dw4 = fmax(Dw4, 0.0);
	Dn4 = fmax(Dn4, 0.0);
	Ds4 = fmax(Ds4, 0.0);
	Dt4 = fmax(Dt4, 0.0);
	Db4 = fmax(Db4, 0.0);

	// ����� �����:
	doublereal Pe = 0.0, Pw = 0.0, Pn = 0.0, Ps = 0.0, Pt = 0.0, Pb = 0.0;
	if (iE > -1) {
		Pe = (Fe) / De;
	}
	if (iW > -1) {
		Pw = -(Fw) / Dw;
	}
	if (iN > -1) {
		Pn = (Fn) / Dn;
	}
	if (iS > -1) {
		Ps = -(Fs) / Ds;
	}
	if (iT > -1) {
		Pt = (Ft) / Dt;
	}
	if (iB > -1) {
		Pb = -(Fb) / Db;
	}

	// ����� �����:
	doublereal Pe2 = 0.0, Pw2 = 0.0, Pn2 = 0.0, Ps2 = 0.0, Pt2 = 0.0, Pb2 = 0.0;
	if (iE2 > -1) {
		Pe2 = (Fe2) / De2;
	}
	if (iW2 > -1) {
		Pw2 = -(Fw2) / Dw2;
	}
	if (iN2 > -1) {
		Pn2 = (Fn2) / Dn2;
	}
	if (iS2 > -1) {
		Ps2 = -(Fs2) / Ds2;
	}
	if (iT2 > -1) {
		Pt2 = (Ft2) / Dt2;
	}
	if (iB2 > -1) {
		Pb2 = -(Fb2) / Db2;
	}

	// ����� �����:
	doublereal Pe3 = 0.0, Pw3 = 0.0, Pn3 = 0.0, Ps3 = 0.0, Pt3 = 0.0, Pb3 = 0.0;
	if (iE3 > -1) {
		Pe3 = (Fe3) / De3;
	}
	if (iW3 > -1) {
		Pw3 = -(Fw3) / Dw3;
	}
	if (iN3 > -1) {
		Pn3 = (Fn3) / Dn3;
	}
	if (iS3 > -1) {
		Ps3 = -(Fs3) / Ds3;
	}
	if (iT3 > -1) {
		Pt3 = (Ft3) / Dt3;
	}
	if (iB3 > -1) {
		Pb3 = -(Fb3) / Db3;
	}

	// ����� �����:
	doublereal Pe4 = 0.0, Pw4 = 0.0, Pn4 = 0.0, Ps4 = 0.0, Pt4 = 0.0, Pb4 = 0.0;
	if (iE4 > -1) {
		Pe4 = (Fe4) / De4;
	}
	if (iW4 > -1) {
		Pw4 = -(Fw4) / Dw4;
	}
	if (iN4 > -1) {
		Pn4 = (Fn4) / Dn4;
	}
	if (iS4 > -1) {
		Ps4 = -(Fs4) / Ds4;
	}
	if (iT4 > -1) {
		Pt4 = (Ft4) / Dt4;
	}
	if (iB4 > -1) {
		Pb4 = -(Fb4) / Db4;
	}


	// ������� � ������ ����� ��� ������������� ����� �������� QUICK
	// � ���� ������������� ������ ���������� ���������.
	// addition to the right side QUICK Leonard.
	doublereal attrs = 0.0;

	if (b_on_adaptive_local_refinement_mesh) {

		// ���������������� ���������.
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].an = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].as = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].at = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = 0.0;

		sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = 0.0;

		sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = 0.0;

		sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = 0.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = 0.0;

	}

	if (ishconvection < distsheme) {

		if (1) {
			// ������� ��� ����������� ������ � ��������������� � ���������� 7.05.2017. 
			if (b_on_adaptive_local_refinement_mesh) {
				// ���������� ������������� ����������� �������:
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);

				// ���������� ������������� ����������� �������:
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(-(Fe2), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(Fw2, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(-(Fn2), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(Fs2, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(-(Ft2), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(Fb2, 0);

				// ���������� ������������� ����������� �������:
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(-(Fe3), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(Fw3, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(-(Fn3), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(Fs3, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(-(Ft3), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(Fb3, 0);

				// ���������� ������������� ����������� �������:
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(-(Fe4), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(Fw4, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(-(Fn4), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(Fs4, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(-(Ft4), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(Fb4, 0);

			}
			else {
				// 25 07 2015
				// ���������� ������������� ����������� �������:
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);
				//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
			}
		}
		else
		{
			// �������� �� ������ ����������������������� 25 ���� 2015.
			if (!bE) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
			}
			else {
				integer inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fabs(Fe);
				}
				else {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
				}
			}
			if (!bW) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
			}
			else {
				integer inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fabs(Fw);
				}
				else {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
				}
			}
			if (!bN) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
			}
			else {
				integer inumber = iN - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fabs(Fn);
				}
				else {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
				}
			}
			if (!bS) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
			}
			else {
				integer inumber = iS - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fabs(Fs);
				}
				else {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
				}
			}
			if (!bT) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
			}
			else {
				integer inumber = iT - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fabs(Ft);
				}
				else {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
				}
			}
			if (!bB) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);
			}
			else
			{
				integer inumber = iB - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fabs(Fb);
				}
				else {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);
				}
			}


			if (b_on_adaptive_local_refinement_mesh) {

				if (!bE2) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(-(Fe2), 0);
				}
				else {
					integer inumber = iE2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fabs(Fe2);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(-(Fe2), 0);
					}
				}
				if (!bW2) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(Fw2, 0);
				}
				else {
					integer inumber = iW2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fabs(Fw2);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(Fw2, 0);
					}
				}
				if (!bN2) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(-(Fn2), 0);
				}
				else {
					integer inumber = iN2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fabs(Fn2);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(-(Fn2), 0);
					}
				}
				if (!bS2) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(Fs2, 0);
				}
				else {
					integer inumber = iS2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fabs(Fs2);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(Fs2, 0);
					}
				}
				if (!bT2) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(-(Ft2), 0);
				}
				else {
					integer inumber = iT2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fabs(Ft2);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(-(Ft2), 0);
					}
				}
				if (!bB2) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(Fb2, 0);
				}
				else
				{
					integer inumber = iB2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fabs(Fb2);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(Fb2, 0);
					}
				}

				if (!bE3) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(-(Fe3), 0);
				}
				else {
					integer inumber = iE3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fabs(Fe3);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(-(Fe3), 0);
					}
				}
				if (!bW3) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(Fw3, 0);
				}
				else {
					integer inumber = iW3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fabs(Fw3);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(Fw3, 0);
					}
				}
				if (!bN3) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(-(Fn3), 0);
				}
				else {
					integer inumber = iN3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fabs(Fn3);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(-(Fn3), 0);
					}
				}
				if (!bS3) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(Fs3, 0);
				}
				else {
					integer inumber = iS3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fabs(Fs3);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(Fs3, 0);
					}
				}
				if (!bT3) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(-(Ft3), 0);
				}
				else {
					integer inumber = iT3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fabs(Ft3);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(-(Ft3), 0);
					}
				}
				if (!bB3) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(Fb3, 0);
				}
				else
				{
					integer inumber = iB3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fabs(Fb3);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(Fb3, 0);
					}
				}

				if (!bE4) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(-(Fe4), 0);
				}
				else {
					integer inumber = iE4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fabs(Fe4);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(-(Fe4), 0);
					}
				}
				if (!bW4) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(Fw4, 0);
				}
				else {
					integer inumber = iW4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fabs(Fw4);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(Fw4, 0);
					}
				}
				if (!bN4) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(-(Fn4), 0);
				}
				else {
					integer inumber = iN4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fabs(Fn4);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(-(Fn4), 0);
					}
				}
				if (!bS4) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(Fs4, 0);
				}
				else {
					integer inumber = iS4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fabs(Fs4);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(Fs4, 0);
					}
				}
				if (!bT4) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(-(Ft4), 0);
				}
				else {
					integer inumber = iT4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fabs(Ft4);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(-(Ft4), 0);
					}
				}
				if (!bB4) {
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(Fb4, 0);
				}
				else
				{
					integer inumber = iB4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fabs(Fb4);
					}
					else {
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(Fb4, 0);
					}
				}
			}


		}

		// ������ ��� ����������� ������ � ��������� � ����������. 7.05.2017.
		//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
		// ��� ���������:
		// ����� ������������� !!! (���������� �� ����).
		if (b_on_adaptive_local_refinement_mesh) {
			// ���� 

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(+(Fe), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-(Fw), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(+(Fn), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-(Fs), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(+(Ft), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(-(Fb), 0);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(+(Fe2), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(-(Fw2), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(+(Fn2), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(-(Fs2), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(+(Ft2), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(-(Fb2), 0);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(+(Fe3), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(-(Fw3), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(+(Fn3), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(-(Fs3), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(+(Ft3), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(-(Fb3), 0);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(+(Fe4), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(-(Fw4), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(+(Fn4), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(-(Fs4), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(+(Ft4), 0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(-(Fb4), 0);
			/*
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = sl[GAMMA_LANGTRY_MENTER_SL][iP].ae +  sl[GAMMA_LANGTRY_MENTER_SL][iP].aw + sl[GAMMA_LANGTRY_MENTER_SL][iP].an + sl[GAMMA_LANGTRY_MENTER_SL][iP].as + sl[GAMMA_LANGTRY_MENTER_SL][iP].at + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
			if (b_on_adaptive_local_refinement_mesh) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2;
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3;
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4;
			}
			*/
		}
		else {
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(+(Fe), 0);
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-(Fw), 0);
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(+(Fn), 0);
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-(Fs), 0);
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(+(Ft), 0);
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(-(Fb), 0);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = sl[GAMMA_LANGTRY_MENTER_SL][iP].ae +
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw + sl[GAMMA_LANGTRY_MENTER_SL][iP].an +
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as + sl[GAMMA_LANGTRY_MENTER_SL][iP].at +
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;

		}

		// 13 ������� 2016
		// ��� ��������. ��� ����� �� �������� � ����������. �� ������� ��� ��������� ������������ ������������, ��
		// ������������� ���������� ���� � �������, �� �� ��������������� ��� ����� ���������� ������������ ����������� ������ ������.
		//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);




		if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ap < 1.0e-36) {
			printf("Zero diagonal coefficient in internal volume in my_elmatr_quad_turbulent_Gamma_Langtry_Menter_SST_3D.\n");
#if doubleintprecision == 1
			printf("ap=%e iP=%lld\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ap, iP);
#else
			printf("ap=%e iP=%d\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ap, iP);
#endif
			if (b_on_adaptive_local_refinement_mesh) {
				printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ae, sl[GAMMA_LANGTRY_MENTER_SL][iP].aw, sl[GAMMA_LANGTRY_MENTER_SL][iP].an, sl[GAMMA_LANGTRY_MENTER_SL][iP].as, sl[GAMMA_LANGTRY_MENTER_SL][iP].at, sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);
				printf("ae2=%e aw2=%e an2=%e as2=%e at2=%e ab2=%e\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2, sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2, sl[GAMMA_LANGTRY_MENTER_SL][iP].an2, sl[GAMMA_LANGTRY_MENTER_SL][iP].as2, sl[GAMMA_LANGTRY_MENTER_SL][iP].at2, sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2);
				printf("ae3=%e aw3=%e an3=%e as3=%e at3=%e ab3=%e\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3, sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3, sl[GAMMA_LANGTRY_MENTER_SL][iP].an3, sl[GAMMA_LANGTRY_MENTER_SL][iP].as3, sl[GAMMA_LANGTRY_MENTER_SL][iP].at3, sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3);
				printf("ae4=%e aw4=%e an4=%e as4=%e at4=%e ab4=%e\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4, sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4, sl[GAMMA_LANGTRY_MENTER_SL][iP].an4, sl[GAMMA_LANGTRY_MENTER_SL][iP].as4, sl[GAMMA_LANGTRY_MENTER_SL][iP].at4, sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4);
			}
			else {
				printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ae, sl[GAMMA_LANGTRY_MENTER_SL][iP].aw, sl[GAMMA_LANGTRY_MENTER_SL][iP].an, sl[GAMMA_LANGTRY_MENTER_SL][iP].as, sl[GAMMA_LANGTRY_MENTER_SL][iP].at, sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);
			}
			system("pause");

			// ��� ������ ������ �� ������� �� ��������
			// �.�. �� ������ ����� ����� ������� ������������ �������.
			// ������� ����� �������� �������� ����� ������ ����� ��������� �����������, ��� ����� �����������.
			// �� �������� �� ��-���� ��������� �������� ����� ���������.
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), BULG);//+fmax(-(Fe),0); 
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), BULG);//+fmax(Fw,0); 
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), BULG);//+fmax(-(Fn),0); 
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), BULG);//+fmax(Fs,0); 
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), BULG);//+fmax(-(Ft),0); 
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), BULG);//+fmax(Fb,0); 
			if (b_on_adaptive_local_refinement_mesh) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), BULG);//+fmax(-(Fe2),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), BULG);//+fmax(Fw2,0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), BULG);//+fmax(-(Fn2),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), BULG);//+fmax(Fs2,0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), BULG);//+fmax(-(Ft2),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), BULG);//+fmax(Fb2,0); 

				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), BULG);//+fmax(-(Fe3),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), BULG);//+fmax(Fw3,0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), BULG);//+fmax(-(Fn3),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), BULG);//+fmax(Fs3,0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), BULG);//+fmax(-(Ft3),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), BULG);//+fmax(Fb3,0); 

				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), BULG);//+fmax(-(Fe4),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), BULG);//+fmax(Fw4,0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), BULG);//+fmax(-(Fn4),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), BULG);//+fmax(Fs4,0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), BULG);//+fmax(-(Ft4),0); 
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), BULG);//+fmax(Fb4,0); 
			}
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = sl[GAMMA_LANGTRY_MENTER_SL][iP].ae + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw + sl[GAMMA_LANGTRY_MENTER_SL][iP].an + sl[GAMMA_LANGTRY_MENTER_SL][iP].as + sl[GAMMA_LANGTRY_MENTER_SL][iP].at + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
			if (b_on_adaptive_local_refinement_mesh) {
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2;
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3;
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4;
			}
		}


		// ������ ��� ����������� ������ � ��������� � ����������. 7.05.2017.
		//sumanb=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
		// ��� ���������:
		// ����� ������������� !!! (���������� �� ����).
		/*
		sumanb = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(+(Fe), 0);
		sumanb += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-(Fw), 0);
		sumanb += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(+(Fn), 0);
		sumanb += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-(Fs), 0);
		sumanb += Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(+(Ft), 0);
		sumanb += Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(-(Fb), 0);
		if (b_on_adaptive_local_refinement_mesh) {
			sumanb += De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(+(Fe2), 0);
			sumanb += Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(-(Fw2), 0);
			sumanb += Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(+(Fn2), 0);
			sumanb += Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(-(Fs2), 0);
			sumanb += Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(+(Ft2), 0);
			sumanb += Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(-(Fb2), 0);

			sumanb += De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(+(Fe3), 0);
			sumanb += Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(-(Fw3), 0);
			sumanb += Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(+(Fn3), 0);
			sumanb += Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(-(Fs3), 0);
			sumanb += Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(+(Ft3), 0);
			sumanb += Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(-(Fb3), 0);

			sumanb += De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(+(Fe4), 0);
			sumanb += Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(-(Fw4), 0);
			sumanb += Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(+(Fn4), 0);
			sumanb += Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(-(Fs4), 0);
			sumanb += Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(+(Ft4), 0);
			sumanb += Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(-(Fb4), 0);
		}
		*/
		//13 ������� 2016.
		//sumanb = fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);
		/*sumanb = sl[GAMMA_LANGTRY_MENTER_SL][iP].ae + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw + sl[GAMMA_LANGTRY_MENTER_SL][iP].an + sl[GAMMA_LANGTRY_MENTER_SL][iP].as + sl[GAMMA_LANGTRY_MENTER_SL][iP].at + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
		if (b_on_adaptive_local_refinement_mesh) {
			sumanb += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2;
			sumanb += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3;
			sumanb += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4;
		}*/

	}
	else if (ishconvection < QUICK)
	{

		if (b_on_adaptive_local_refinement_mesh) {

			// ���������� ������������� ����������� �������:
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = -(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = (Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an = -(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as = (Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at = -(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = (Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			// ���������� ������������� ����������� �������:
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = -(Fe2)*fC(Pe2, ishconvection, true, feplus2) + De2 * fD(Pe2, ishconvection, true, feplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = (Fw2)*fC(Pw2, ishconvection, true, fwplus2) + Dw2 * fD(Pw2, ishconvection, true, fwplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = -(Fn2)*fC(Pn2, ishconvection, true, fnplus2) + Dn2 * fD(Pn2, ishconvection, true, fnplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = (Fs2)*fC(Ps2, ishconvection, true, fsplus2) + Ds2 * fD(Ps2, ishconvection, true, fsplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = -(Ft2)*fC(Pt2, ishconvection, true, ftplus2) + Dt2 * fD(Pt2, ishconvection, true, ftplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = (Fb2)*fC(Pb2, ishconvection, true, fbplus2) + Db2 * fD(Pb2, ishconvection, true, fbplus2);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = -(Fe3)*fC(Pe3, ishconvection, true, feplus3) + De3 * fD(Pe3, ishconvection, true, feplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = (Fw3)*fC(Pw3, ishconvection, true, fwplus3) + Dw3 * fD(Pw3, ishconvection, true, fwplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = -(Fn3)*fC(Pn3, ishconvection, true, fnplus3) + Dn3 * fD(Pn3, ishconvection, true, fnplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = (Fs3)*fC(Ps3, ishconvection, true, fsplus3) + Ds3 * fD(Ps3, ishconvection, true, fsplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = -(Ft3)*fC(Pt3, ishconvection, true, ftplus3) + Dt3 * fD(Pt3, ishconvection, true, ftplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = (Fb3)*fC(Pb3, ishconvection, true, fbplus3) + Db3 * fD(Pb3, ishconvection, true, fbplus3);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = -(Fe4)*fC(Pe4, ishconvection, true, feplus4) + De4 * fD(Pe4, ishconvection, true, feplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = (Fw4)*fC(Pw4, ishconvection, true, fwplus4) + Dw4 * fD(Pw4, ishconvection, true, fwplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = -(Fn4)*fC(Pn4, ishconvection, true, fnplus4) + Dn4 * fD(Pn4, ishconvection, true, fnplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = (Fs4)*fC(Ps4, ishconvection, true, fsplus4) + Ds4 * fD(Ps4, ishconvection, true, fsplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = -(Ft4)*fC(Pt4, ishconvection, true, ftplus4) + Dt4 * fD(Pt4, ishconvection, true, ftplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = (Fb4)*fC(Pb4, ishconvection, true, fbplus4) + Db4 * fD(Pb4, ishconvection, true, fbplus4);

			// 08.05.2017.
			// ��� ���������:
			// ����� ������������� !!! (���������� �� ����).

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fe2)*fC(Pe2, ishconvection, true, feplus2) + De2 * fD(Pe2, ishconvection, true, feplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fw2)*fC(Pw2, ishconvection, true, fwplus2) + Dw2 * fD(Pw2, ishconvection, true, fwplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fn2)*fC(Pn2, ishconvection, true, fnplus2) + Dn2 * fD(Pn2, ishconvection, true, fnplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fs2)*fC(Ps2, ishconvection, true, fsplus2) + Ds2 * fD(Ps2, ishconvection, true, fsplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Ft2)*fC(Pt2, ishconvection, true, ftplus2) + Dt2 * fD(Pt2, ishconvection, true, ftplus2);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fb2)*fC(Pb2, ishconvection, true, fbplus2) + Db2 * fD(Pb2, ishconvection, true, fbplus2);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fe3)*fC(Pe3, ishconvection, true, feplus3) + De3 * fD(Pe3, ishconvection, true, feplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fw3)*fC(Pw3, ishconvection, true, fwplus3) + Dw3 * fD(Pw3, ishconvection, true, fwplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fn3)*fC(Pn3, ishconvection, true, fnplus3) + Dn3 * fD(Pn3, ishconvection, true, fnplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fs3)*fC(Ps3, ishconvection, true, fsplus3) + Ds3 * fD(Ps3, ishconvection, true, fsplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Ft3)*fC(Pt3, ishconvection, true, ftplus3) + Dt3 * fD(Pt3, ishconvection, true, ftplus3);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fb3)*fC(Pb3, ishconvection, true, fbplus3) + Db3 * fD(Pb3, ishconvection, true, fbplus3);

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fe4)*fC(Pe4, ishconvection, true, feplus4) + De4 * fD(Pe4, ishconvection, true, feplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fw4)*fC(Pw4, ishconvection, true, fwplus4) + Dw4 * fD(Pw4, ishconvection, true, fwplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fn4)*fC(Pn4, ishconvection, true, fnplus4) + Dn4 * fD(Pn4, ishconvection, true, fnplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fs4)*fC(Ps4, ishconvection, true, fsplus4) + Ds4 * fD(Ps4, ishconvection, true, fsplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Ft4)*fC(Pt4, ishconvection, true, ftplus4) + Dt4 * fD(Pt4, ishconvection, true, ftplus4);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fb4)*fC(Pb4, ishconvection, true, fbplus4) + Db4 * fD(Pb4, ishconvection, true, fbplus4);

			/*
			sumanb = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sumanb += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sumanb += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sumanb += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sumanb += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sumanb += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			sumanb += +(Fe2)*fC(Pe2, ishconvection, true, feplus2) + De2 * fD(Pe2, ishconvection, true, feplus2);
			sumanb += -(Fw2)*fC(Pw2, ishconvection, true, fwplus2) + Dw2 * fD(Pw2, ishconvection, true, fwplus2);
			sumanb += +(Fn2)*fC(Pn2, ishconvection, true, fnplus2) + Dn2 * fD(Pn2, ishconvection, true, fnplus2);
			sumanb += -(Fs2)*fC(Ps2, ishconvection, true, fsplus2) + Ds2 * fD(Ps2, ishconvection, true, fsplus2);
			sumanb += +(Ft2)*fC(Pt2, ishconvection, true, ftplus2) + Dt2 * fD(Pt2, ishconvection, true, ftplus2);
			sumanb += -(Fb2)*fC(Pb2, ishconvection, true, fbplus2) + Db2 * fD(Pb2, ishconvection, true, fbplus2);

			sumanb += +(Fe3)*fC(Pe3, ishconvection, true, feplus3) + De3 * fD(Pe3, ishconvection, true, feplus3);
			sumanb += -(Fw3)*fC(Pw3, ishconvection, true, fwplus3) + Dw3 * fD(Pw3, ishconvection, true, fwplus3);
			sumanb += +(Fn3)*fC(Pn3, ishconvection, true, fnplus3) + Dn3 * fD(Pn3, ishconvection, true, fnplus3);
			sumanb += -(Fs3)*fC(Ps3, ishconvection, true, fsplus3) + Ds3 * fD(Ps3, ishconvection, true, fsplus3);
			sumanb += +(Ft3)*fC(Pt3, ishconvection, true, ftplus3) + Dt3 * fD(Pt3, ishconvection, true, ftplus3);
			sumanb += -(Fb3)*fC(Pb3, ishconvection, true, fbplus3) + Db3 * fD(Pb3, ishconvection, true, fbplus3);

			sumanb += +(Fe4)*fC(Pe4, ishconvection, true, feplus4) + De4 * fD(Pe4, ishconvection, true, feplus4);
			sumanb += -(Fw4)*fC(Pw4, ishconvection, true, fwplus4) + Dw4 * fD(Pw4, ishconvection, true, fwplus4);
			sumanb += +(Fn4)*fC(Pn4, ishconvection, true, fnplus4) + Dn4 * fD(Pn4, ishconvection, true, fnplus4);
			sumanb += -(Fs4)*fC(Ps4, ishconvection, true, fsplus4) + Ds4 * fD(Ps4, ishconvection, true, fsplus4);
			sumanb += +(Ft4)*fC(Pt4, ishconvection, true, ftplus4) + Dt4 * fD(Pt4, ishconvection, true, ftplus4);
			sumanb += -(Fb4)*fC(Pb4, ishconvection, true, fbplus4) + Db4 * fD(Pb4, ishconvection, true, fbplus4);
			*/

		}
		else
		{
			// ���������� ������������� ����������� �������:
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = -(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = (Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an = -(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as = (Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at = -(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = (Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);
			//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;

			// ������ ��� ����������� ������ � ��������� � ����������. 7.05.2017.
			//sumanb=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;
			//13 ������� 2016.
			//sumanb = fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);

			// 08.05.2017.
			// ��� ���������:
			// ����� ������������� !!! (���������� �� ����).

			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			/*
			sumanb = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sumanb += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sumanb += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sumanb += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sumanb += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sumanb += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);
			*/
		}
	}
	else if (ishconvection >= QUICK)
	{
		// � 3D ������������ ������ ����� ������������ �� ��� ���������� �����.
		// ����� �������� ����� ������ ������� � ����������� � ������� ��������� ���������� ���������.

		// X - direction
		doublereal positionxP, positionxE, positionxW, positionxEE, positionxWW, positionxe, positionxw;
		doublereal SpeedP = 0.0, SpeedE = 0.0, SpeedW = 0.0, SpeedEE = 0.0, SpeedWW = 0.0, SpeedN = 0.0, SpeedS = 0.0;
		doublereal SpeedNN = 0.0, SpeedSS = 0.0, SpeedT = 0.0, SpeedB = 0.0, SpeedTT = 0.0, SpeedBB = 0.0;
		doublereal Speede = 0.0, Speedw = 0.0, Speedn = 0.0, Speeds = 0.0, Speedt = 0.0, Speedb = 0.0;
		// Y - direction
		doublereal positionyP, positionyN, positionyS, positionyNN, positionySS, positionyn, positionys;
		// Z - direction
		doublereal positionzP, positionzT, positionzB, positionzTT, positionzBB, positionzt, positionzb;

		doublereal SpeedE2 = 0.0, SpeedW2 = 0.0, SpeedEE2 = 0.0, SpeedWW2 = 0.0, SpeedN2 = 0.0, SpeedS2 = 0.0;
		doublereal SpeedNN2 = 0.0, SpeedSS2 = 0.0, SpeedT2 = 0.0, SpeedB2 = 0.0, SpeedTT2 = 0.0, SpeedBB2 = 0.0;
		doublereal Speede2 = 0.0, Speedw2 = 0.0, Speedn2 = 0.0, Speeds2 = 0.0, Speedt2 = 0.0, Speedb2 = 0.0;

		doublereal  SpeedE3 = 0.0, SpeedW3 = 0.0, SpeedEE3 = 0.0, SpeedWW3 = 0.0, SpeedN3 = 0.0, SpeedS3 = 0.0;
		doublereal SpeedNN3 = 0.0, SpeedSS3 = 0.0, SpeedT3 = 0.0, SpeedB3 = 0.0, SpeedTT3 = 0.0, SpeedBB3 = 0.0;
		doublereal Speede3 = 0.0, Speedw3 = 0.0, Speedn3 = 0.0, Speeds3 = 0.0, Speedt3 = 0.0, Speedb3 = 0.0;

		doublereal SpeedE4 = 0.0, SpeedW4 = 0.0, SpeedEE4 = 0.0, SpeedWW4 = 0.0, SpeedN4 = 0.0, SpeedS4 = 0.0;
		doublereal SpeedNN4 = 0.0, SpeedSS4 = 0.0, SpeedT4 = 0.0, SpeedB4 = 0.0, SpeedTT4 = 0.0, SpeedBB4 = 0.0;
		doublereal Speede4 = 0.0, Speedw4 = 0.0, Speedn4 = 0.0, Speeds4 = 0.0, Speedt4 = 0.0, Speedb4 = 0.0;

		// X - direction
		doublereal  positionxE2, positionxW2, positionxEE2, positionxWW2, positionxe2, positionxw2;
		// Y - direction
		doublereal  positionyN2, positionyS2, positionyNN2, positionySS2, positionyn2, positionys2;
		// Z - direction
		doublereal  positionzT2, positionzB2, positionzTT2, positionzBB2, positionzt2, positionzb2;

		// X - direction
		doublereal  positionxE3, positionxW3, positionxEE3, positionxWW3, positionxe3, positionxw3;
		// Y - direction
		doublereal  positionyN3, positionyS3, positionyNN3, positionySS3, positionyn3, positionys3;
		// Z - direction
		doublereal  positionzT3, positionzB3, positionzTT3, positionzBB3, positionzt3, positionzb3;

		// X - direction
		doublereal  positionxE4, positionxW4, positionxEE4, positionxWW4, positionxe4, positionxw4;
		// Y - direction
		doublereal  positionyN4, positionyS4, positionyNN4, positionySS4, positionyn4, positionys4;
		// Z - direction
		doublereal  positionzT4, positionzB4, positionzTT4, positionzBB4, positionzt4, positionzb4;

		TOCHKA pointP;
		//center_cord3D(iP, nvtx, pa, pointP, 100);
		pointP = center_coord_loc[iP];

		positionxP = pointP.x; positionyP = pointP.y; positionzP = pointP.z;
		SpeedP = potent[GAMMA_LANGTRY_MENTER][iP];
		// X - direction
		if (!bE) {
			SpeedE = potent[GAMMA_LANGTRY_MENTER][iE];
			//center_cord3D(iE,nvtx,pa,pointP,E_SIDE);
			pointP = center_coord_loc[iE];

			positionxE = pointP.x;
			positionxe = positionxP + 0.5 * dx;

			integer iEE = neighbors_for_the_internal_node[E_SIDE][0][iE];
			if (iEE < 0) {
				iEE = neighbors_for_the_internal_node[E_SIDE][1][iE];
			}
			if (iEE < 0) {
				iEE = neighbors_for_the_internal_node[E_SIDE][2][iE];
			}
			if (iEE < 0) {
				iEE = neighbors_for_the_internal_node[E_SIDE][3][iE];
			}

			if ((iEE >= 0) && (iEE < maxelm)) {
				// ���������� ����
				SpeedEE = potent[GAMMA_LANGTRY_MENTER][iEE];
				//center_cord3D(iEE,nvtx,pa,pointP,EE_SIDE);
				pointP = center_coord_loc[iEE];
				positionxEE = pointP.x;
			}
			else
			{
				// ��������� ����
				if ((iEE >= maxelm) && (iEE < maxelm + maxbound)) {
					SpeedEE = potent[GAMMA_LANGTRY_MENTER][iEE];
				}
				else {
					SpeedEE = SpeedE;
				}
				//volume3D(iE, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iE];

				positionxEE = positionxE + 0.5 * pointP.x;
			}
		}
		else {
			// ��� ��������� ����
			SpeedE = potent[GAMMA_LANGTRY_MENTER][iE];
			SpeedEE = potent[GAMMA_LANGTRY_MENTER][iE];
			positionxe = positionxP + 0.5 * dx;
			positionxE = positionxP + 0.5 * dx;
			positionxEE = positionxP + dx; // ����� ���� �� ���������� !
		}

		if (!bW) {
			//center_cord3D(iW,nvtx,pa,pointP,W_SIDE);
			pointP = center_coord_loc[iW];
			positionxW = pointP.x;
			positionxw = positionxP - 0.5 * dx;
			SpeedW = potent[GAMMA_LANGTRY_MENTER][iW];

			integer iWW = neighbors_for_the_internal_node[W_SIDE][0][iW];
			if (iWW < 0) {
				iWW = neighbors_for_the_internal_node[W_SIDE][1][iW];
			}
			if (iWW < 0) {
				iWW = neighbors_for_the_internal_node[W_SIDE][2][iW];
			}
			if (iWW < 0) {
				iWW = neighbors_for_the_internal_node[W_SIDE][3][iW];
			}

			if ((iWW >= 0) && (iWW < maxelm)) {
				// ���������� ����
				SpeedWW = potent[GAMMA_LANGTRY_MENTER][iWW];
				//center_cord3D(iWW,nvtx,pa,pointP,WW_SIDE);
				pointP = center_coord_loc[iWW];
				positionxWW = pointP.x;
			}
			else
			{
				// ��������� ����
				if ((iWW >= maxelm) && (iWW < maxelm + maxbound)) {
					SpeedWW = potent[GAMMA_LANGTRY_MENTER][iWW];
				}
				else {
					SpeedWW = SpeedW;
				}
				//volume3D(iW, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iW];

				positionxWW = positionxW - 0.5 * pointP.x;
			}
		}
		else {
			// ��� ��������� ����
			SpeedW = potent[GAMMA_LANGTRY_MENTER][iW]; // Attantion !! Debug
			SpeedWW = potent[GAMMA_LANGTRY_MENTER][iW];
			//printf("SpeedW==%e\n",SpeedW); getchar();
			positionxw = positionxP - 0.5 * dx;
			positionxW = positionxP - 0.5 * dx;
			positionxWW = positionxP - dx; // ����� ���� �� ���������� !
		}

		// Y - direction
		if (!bN) {
			SpeedN = potent[GAMMA_LANGTRY_MENTER][iN];
			//center_cord3D(iN,nvtx,pa,pointP,N_SIDE);
			pointP = center_coord_loc[iN];
			positionyN = pointP.y;
			positionyn = positionxP + 0.5 * dy;

			integer iNN = neighbors_for_the_internal_node[N_SIDE][0][iN];
			if (iNN < 0) {
				iNN = neighbors_for_the_internal_node[N_SIDE][1][iN];
			}
			if (iNN < 0) {
				iNN = neighbors_for_the_internal_node[N_SIDE][2][iN];
			}
			if (iNN < 0) {
				iNN = neighbors_for_the_internal_node[N_SIDE][3][iN];
			}

			if ((iNN >= 0) && (iNN < maxelm)) {
				// ���������� ����
				SpeedNN = potent[GAMMA_LANGTRY_MENTER][iNN];
				//center_cord3D(iNN,nvtx,pa,pointP,NN_SIDE);
				pointP = center_coord_loc[iNN];
				positionyNN = pointP.y;
			}
			else
			{
				// ��������� ����
				if ((iNN >= maxelm) && (iNN < maxelm + maxbound)) {
					SpeedNN = potent[GAMMA_LANGTRY_MENTER][iNN];
				}
				else {
					SpeedNN = SpeedN;
				}
				//volume3D(iN, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iN];

				positionyNN = positionyN + 0.5 * pointP.y;
			}
		}
		else {
			// ��� ��������� ����
			SpeedN = potent[GAMMA_LANGTRY_MENTER][iN];
			SpeedNN = potent[GAMMA_LANGTRY_MENTER][iN];
			positionyn = positionyP + 0.5 * dy;
			positionyN = positionyP + 0.5 * dy;
			positionyNN = positionyP + dy; // ����� ���� �� ���������� !
		}

		if (!bS) {
			SpeedS = potent[GAMMA_LANGTRY_MENTER][iS];
			//center_cord3D(iS,nvtx,pa,pointP,S_SIDE);
			pointP = center_coord_loc[iS];
			positionyS = pointP.y;
			positionys = positionyP - 0.5 * dy;

			integer iSS = neighbors_for_the_internal_node[S_SIDE][0][iS];
			if (iSS < 0) {
				iSS = neighbors_for_the_internal_node[S_SIDE][1][iS];
			}
			if (iSS < 0) {
				iSS = neighbors_for_the_internal_node[S_SIDE][2][iS];
			}
			if (iSS < 0) {
				iSS = neighbors_for_the_internal_node[S_SIDE][3][iS];
			}

			if ((iSS >= 0) && (iSS < maxelm)) {
				// ���������� ����
				SpeedSS = potent[GAMMA_LANGTRY_MENTER][iSS];
				//center_cord3D(iSS,nvtx,pa,pointP,SS_SIDE);
				pointP = center_coord_loc[iSS];
				positionySS = pointP.y;
			}
			else
			{
				// ��������� ����
				if ((iSS >= maxelm) && (iSS < maxelm + maxbound)) {
					SpeedSS = potent[GAMMA_LANGTRY_MENTER][iSS];
				}
				else {
					SpeedSS = SpeedS;
				}
				//volume3D(iS, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iS];

				positionySS = positionyS - 0.5 * pointP.y;
			}
		}
		else {
			// ��� ��������� ����
			SpeedS = potent[GAMMA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
			SpeedSS = potent[GAMMA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
			positionys = positionyP - 0.5 * dy;
			positionyS = positionyP - 0.5 * dy;
			positionySS = positionyP - dy; // ����� ���� �� ���������� !
		}

		// Z - direction
		if (!bT) {
			SpeedT = potent[GAMMA_LANGTRY_MENTER][iT];
			//center_cord3D(iT,nvtx,pa,pointP,T_SIDE);
			pointP = center_coord_loc[iT];
			positionzT = pointP.z;
			positionzt = positionzP + 0.5 * dz;

			integer iTT = neighbors_for_the_internal_node[T_SIDE][0][iT];
			if (iTT < 0) {
				iTT = neighbors_for_the_internal_node[T_SIDE][1][iT];
			}
			if (iTT < 0) {
				iTT = neighbors_for_the_internal_node[T_SIDE][2][iT];
			}
			if (iTT < 0) {
				iTT = neighbors_for_the_internal_node[T_SIDE][3][iT];
			}

			if ((iTT >= 0) && (iTT < maxelm)) {
				// ���������� ����
				SpeedTT = potent[GAMMA_LANGTRY_MENTER][iTT];
				//center_cord3D(iTT,nvtx,pa,pointP,TT_SIDE);
				pointP = center_coord_loc[iTT];
				positionzTT = pointP.z;
			}
			else
			{
				// ��������� ����
				if ((iTT >= maxelm) && (iTT < maxelm + maxbound)) {

					SpeedTT = potent[GAMMA_LANGTRY_MENTER][iTT];
				}
				else {
					SpeedTT = SpeedT;
				}
				//volume3D(iT, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iT];

				positionzTT = positionzT + 0.5 * pointP.z;
			}
		}
		else {
			// ��� ��������� ����
			SpeedT = potent[GAMMA_LANGTRY_MENTER][iT];
			SpeedTT = potent[GAMMA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
			positionzt = positionzP + 0.5 * dz;
			positionzT = positionzP + 0.5 * dz;
			positionzTT = positionzP + dz; // ����� ���� �� ���������� !
		}

		if (!bB) {
			SpeedB = potent[GAMMA_LANGTRY_MENTER][iB];
			//center_cord3D(iB,nvtx,pa,pointP,B_SIDE);
			pointP = center_coord_loc[iB];
			positionzB = pointP.z;
			positionzb = positionzP - 0.5 * dz;

			integer iBB = neighbors_for_the_internal_node[B_SIDE][0][iB];
			if (iBB < 0) {
				iBB = neighbors_for_the_internal_node[B_SIDE][1][iB];
			}
			if (iBB < 0) {
				iBB = neighbors_for_the_internal_node[B_SIDE][2][iB];
			}
			if (iBB < 0) {
				iBB = neighbors_for_the_internal_node[B_SIDE][3][iB];
			}

			if ((iBB >= 0) && (iBB < maxelm)) {
				// ���������� ����
				SpeedBB = potent[GAMMA_LANGTRY_MENTER][iBB];
				//center_cord3D(iBB,nvtx,pa,pointP,BB_SIDE);
				pointP = center_coord_loc[iBB];
				positionzBB = pointP.z;
			}
			else
			{
				// ��������� ����
				if ((iBB >= maxelm) && (iBB < maxelm + maxbound)) {
					SpeedBB = potent[GAMMA_LANGTRY_MENTER][iBB];
				}
				else {
					SpeedBB = SpeedB;
				}
				//volume3D(iB, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iB];

				positionzBB = positionzB - 0.5 * pointP.z;
			}
		}
		else {
			// ��� ��������� ����
			SpeedB = potent[GAMMA_LANGTRY_MENTER][iB];
			SpeedBB = potent[GAMMA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
			positionzb = positionzP - 0.5 * dz;
			positionzB = positionzP - 0.5 * dz;
			positionzBB = positionzP - dz; // ����� ���� �� ���������� !
		}

		if (b_on_adaptive_local_refinement_mesh)
		{

			// X - direction
			if ((!bE2) && (iE2 > -1)) {
				SpeedE2 = potent[GAMMA_LANGTRY_MENTER][iE2];
				//center_cord3D(iE,nvtx,pa,pointP,E_SIDE);
				pointP = center_coord_loc[iE2];

				positionxE2 = pointP.x;
				positionxe2 = positionxP + 0.5 * dx;

				integer iEE2 = neighbors_for_the_internal_node[E_SIDE][0][iE2];
				if (iEE2 < 0) {
					iEE2 = neighbors_for_the_internal_node[E_SIDE][1][iE2];
				}
				if (iEE2 < 0) {
					iEE2 = neighbors_for_the_internal_node[E_SIDE][2][iE2];
				}
				if (iEE2 < 0) {
					iEE2 = neighbors_for_the_internal_node[E_SIDE][3][iE2];
				}
				if ((iEE2 >= 0) && (iEE2 < maxelm)) {
					// ���������� ����
					SpeedEE2 = potent[GAMMA_LANGTRY_MENTER][iEE2];
					//center_cord3D(iEE2,nvtx,pa,pointP,EE_SIDE);
					pointP = center_coord_loc[iEE2];
					positionxEE2 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iEE2 >= maxelm) && (iEE2 < maxelm + maxbound)) {
						SpeedEE2 = potent[GAMMA_LANGTRY_MENTER][iEE2];
					}
					else {
						SpeedEE2 = SpeedE2;
						//std::cout << "iEE2 =" << iEE2 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iE, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iE2];

					positionxEE2 = positionxE2 + 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iE2 > -1) {
					SpeedE2 = potent[GAMMA_LANGTRY_MENTER][iE2];
					SpeedEE2 = potent[GAMMA_LANGTRY_MENTER][iE2];
				}
				else {
					SpeedE2 = potent[GAMMA_LANGTRY_MENTER][iE];
					SpeedEE2 = potent[GAMMA_LANGTRY_MENTER][iE];
				}
				positionxe2 = positionxP + 0.5 * dx;
				positionxE2 = positionxP + 0.5 * dx;
				positionxEE2 = positionxP + dx; // ����� ���� �� ���������� !
			}

			if ((!bW2) && ((iW2 > -1))) {
				//center_cord3D(iW,nvtx,pa,pointP,W_SIDE);
				pointP = center_coord_loc[iW2];
				positionxW2 = pointP.x;
				positionxw2 = positionxP - 0.5 * dx;
				SpeedW2 = potent[GAMMA_LANGTRY_MENTER][iW2];

				integer iWW2 = neighbors_for_the_internal_node[W_SIDE][0][iW2];
				if (iWW2 < 0) {
					iWW2 = neighbors_for_the_internal_node[W_SIDE][1][iW2];
				}
				if (iWW2 < 0) {
					iWW2 = neighbors_for_the_internal_node[W_SIDE][2][iW2];
				}
				if (iWW2 < 0) {
					iWW2 = neighbors_for_the_internal_node[W_SIDE][3][iW2];
				}

				if ((iWW2 >= 0) && (iWW2 < maxelm)) {
					// ���������� ����
					SpeedWW2 = potent[GAMMA_LANGTRY_MENTER][iWW2];
					//center_cord3D(iWW2,nvtx,pa,pointP,WW_SIDE);
					pointP = center_coord_loc[iWW2];
					positionxWW2 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iWW2 >= maxelm) && (iWW2 < maxelm + maxbound)) {
						SpeedWW2 = potent[GAMMA_LANGTRY_MENTER][iWW2];
					}
					else {

						SpeedWW2 = SpeedW2;
						//std::cout << "iWW2 =" << iWW2 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iW2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iW2];

					positionxWW2 = positionxW2 - 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iW2 > -1) {
					SpeedW2 = potent[GAMMA_LANGTRY_MENTER][iW2];
					SpeedWW2 = potent[GAMMA_LANGTRY_MENTER][iW2];
				}
				else {
					SpeedW2 = potent[GAMMA_LANGTRY_MENTER][iW];
					SpeedWW2 = potent[GAMMA_LANGTRY_MENTER][iW];
				}
				//printf("SpeedW2==%e\n",SpeedW2); getchar();
				positionxw2 = positionxP - 0.5 * dx;
				positionxW2 = positionxP - 0.5 * dx;
				positionxWW2 = positionxP - dx; // ����� ���� �� ���������� !
			}

			// Y - direction
			if ((!bN2) && (iN2 > -1)) {
				SpeedN2 = potent[GAMMA_LANGTRY_MENTER][iN2];
				//center_cord3D(iN2,nvtx,pa,pointP,N_SIDE);
				pointP = center_coord_loc[iN2];
				positionyN2 = pointP.y;
				positionyn2 = positionxP + 0.5 * dy;

				integer iNN2 = neighbors_for_the_internal_node[N_SIDE][0][iN2];
				if (iNN2 < 0) {
					iNN2 = neighbors_for_the_internal_node[N_SIDE][1][iN2];
				}
				if (iNN2 < 0) {
					iNN2 = neighbors_for_the_internal_node[N_SIDE][2][iN2];
				}
				if (iNN2 < 0) {
					iNN2 = neighbors_for_the_internal_node[N_SIDE][3][iN2];
				}


				if ((iNN2 >= 0) && (iNN2 < maxelm)) {
					// ���������� ����
					SpeedNN2 = potent[GAMMA_LANGTRY_MENTER][iNN2];
					//center_cord3D(iNN2,nvtx,pa,pointP,NN_SIDE);
					pointP = center_coord_loc[iNN2];
					positionyNN2 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iNN2 >= maxelm) && (iNN2 < maxelm + maxbound)) {
						SpeedNN2 = potent[GAMMA_LANGTRY_MENTER][iNN2];
					}
					else {
						SpeedNN2 = SpeedN2;
					}
					//volume3D(iN2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iN2];

					positionyNN2 = positionyN2 + 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iN2 > -1) {
					SpeedN2 = potent[GAMMA_LANGTRY_MENTER][iN2];
					SpeedNN2 = potent[GAMMA_LANGTRY_MENTER][iN2];
				}
				else {
					SpeedN2 = potent[GAMMA_LANGTRY_MENTER][iN];
					SpeedNN2 = potent[GAMMA_LANGTRY_MENTER][iN];
				}
				positionyn2 = positionyP + 0.5 * dy;
				positionyN2 = positionyP + 0.5 * dy;
				positionyNN2 = positionyP + dy; // ����� ���� �� ���������� !
			}

			if ((!bS2) && (iS2 > -1)) {
				SpeedS2 = potent[GAMMA_LANGTRY_MENTER][iS2];
				//center_cord3D(iS2,nvtx,pa,pointP,S_SIDE);
				pointP = center_coord_loc[iS2];
				positionyS2 = pointP.y;
				positionys2 = positionyP - 0.5 * dy;

				integer iSS2 = neighbors_for_the_internal_node[S_SIDE][0][iS2];
				if (iSS2 < 0) {
					iSS2 = neighbors_for_the_internal_node[S_SIDE][1][iS2];
				}
				if (iSS2 < 0) {
					iSS2 = neighbors_for_the_internal_node[S_SIDE][2][iS2];
				}
				if (iSS2 < 0) {
					iSS2 = neighbors_for_the_internal_node[S_SIDE][3][iS2];
				}

				if ((iSS2 >= 0) && (iSS2 < maxelm)) {
					// ���������� ����
					SpeedSS2 = potent[GAMMA_LANGTRY_MENTER][iSS2];
					//center_cord3D(iSS2,nvtx,pa,pointP,SS_SIDE);
					pointP = center_coord_loc[iSS2];
					positionySS2 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iSS2 >= maxelm) && (iSS2 < maxelm + maxbound)) {
						SpeedSS2 = potent[GAMMA_LANGTRY_MENTER][iSS2];
					}
					else {
						SpeedSS2 = SpeedS2;
					}
					//volume3D(iS2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iS2];

					positionySS2 = positionyS2 - 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iS2 > -1) {
					SpeedS2 = potent[GAMMA_LANGTRY_MENTER][iS2]; // ATTANTION !!!!
					SpeedSS2 = potent[GAMMA_LANGTRY_MENTER][iS2]; // ������� �������� ������ ������� ����.
				}
				else {
					SpeedS2 = potent[GAMMA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
					SpeedSS2 = potent[GAMMA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
				}
				positionys2 = positionyP - 0.5 * dy;
				positionyS2 = positionyP - 0.5 * dy;
				positionySS2 = positionyP - dy; // ����� ���� �� ���������� !
			}

			// Z - direction
			if ((!bT2) && (iT2 > -1)) {
				SpeedT2 = potent[GAMMA_LANGTRY_MENTER][iT2];
				//center_cord3D(iT2,nvtx,pa,pointP,T_SIDE);
				pointP = center_coord_loc[iT2];
				positionzT2 = pointP.z;
				positionzt2 = positionzP + 0.5 * dz;

				integer iTT2 = neighbors_for_the_internal_node[T_SIDE][0][iT2];
				if (iTT2 < 0) {
					iTT2 = neighbors_for_the_internal_node[T_SIDE][1][iT2];
				}
				if (iTT2 < 0) {
					iTT2 = neighbors_for_the_internal_node[T_SIDE][2][iT2];
				}
				if (iTT2 < 0) {
					iTT2 = neighbors_for_the_internal_node[T_SIDE][3][iT2];
				}


				if ((iTT2 >= 0) && (iTT2 < maxelm)) {
					// ���������� ����
					SpeedTT2 = potent[GAMMA_LANGTRY_MENTER][iTT2];
					//center_cord3D(iTT2,nvtx,pa,pointP,TT_SIDE);
					pointP = center_coord_loc[iTT2];
					positionzTT2 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iTT2 >= maxelm) && (iTT2 < maxelm + maxbound)) {
						SpeedTT2 = potent[GAMMA_LANGTRY_MENTER][iTT2];
					}
					else {
						SpeedTT2 = SpeedT2;
					}
					//volume3D(iT2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iT2];

					positionzTT2 = positionzT2 + 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iT2 > -1) {
					SpeedT2 = potent[GAMMA_LANGTRY_MENTER][iT2];
					SpeedTT2 = potent[GAMMA_LANGTRY_MENTER][iT2]; // �������� ������ ������� ����
				}
				else {
					SpeedT2 = potent[GAMMA_LANGTRY_MENTER][iT];
					SpeedTT2 = potent[GAMMA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
				}
				positionzt2 = positionzP + 0.5 * dz;
				positionzT2 = positionzP + 0.5 * dz;
				positionzTT2 = positionzP + dz; // ����� ���� �� ���������� !
			}

			if ((!bB2) && (iB2 > -1)) {
				SpeedB2 = potent[GAMMA_LANGTRY_MENTER][iB2];
				//center_cord3D(iB2,nvtx,pa,pointP,B_SIDE);
				pointP = center_coord_loc[iB2];
				positionzB2 = pointP.z;
				positionzb2 = positionzP - 0.5 * dz;

				integer iBB2 = neighbors_for_the_internal_node[B_SIDE][0][iB2];
				if (iBB2 < 0) {
					iBB2 = neighbors_for_the_internal_node[B_SIDE][1][iB2];
				}
				if (iBB2 < 0) {
					iBB2 = neighbors_for_the_internal_node[B_SIDE][2][iB2];
				}
				if (iBB2 < 0) {
					iBB2 = neighbors_for_the_internal_node[B_SIDE][3][iB2];
				}

				if ((iBB2 >= 0) && (iBB2 < maxelm)) {
					// ���������� ����
					SpeedBB2 = potent[GAMMA_LANGTRY_MENTER][iBB2];
					//center_cord3D(iBB2,nvtx,pa,pointP,BB_SIDE);
					pointP = center_coord_loc[iBB2];
					positionzBB2 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iBB2 >= maxelm) && (iBB2 < maxelm + maxbound)) {
						SpeedBB2 = potent[GAMMA_LANGTRY_MENTER][iBB2];
					}
					else {
						SpeedBB2 = SpeedB2;
					}
					//volume3D(iB2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iB2];

					positionzBB2 = positionzB2 - 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iB2 > -1) {
					SpeedB2 = potent[GAMMA_LANGTRY_MENTER][iB2];
					SpeedBB2 = potent[GAMMA_LANGTRY_MENTER][iB2]; // �������� ������ ������� ����
				}
				else {
					SpeedB2 = potent[GAMMA_LANGTRY_MENTER][iB];
					SpeedBB2 = potent[GAMMA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
				}
				positionzb2 = positionzP - 0.5 * dz;
				positionzB2 = positionzP - 0.5 * dz;
				positionzBB2 = positionzP - dz; // ����� ���� �� ���������� !
			}


			// X - direction
			if ((!bE3) && (iE3 > -1)) {
				SpeedE3 = potent[GAMMA_LANGTRY_MENTER][iE3];
				//center_cord3D(iE3,nvtx,pa,pointP,E_SIDE);
				pointP = center_coord_loc[iE3];

				positionxE3 = pointP.x;
				positionxe3 = positionxP + 0.5 * dx;

				integer iEE3 = neighbors_for_the_internal_node[E_SIDE][0][iE3];
				if (iEE3 < 0) {
					iEE3 = neighbors_for_the_internal_node[E_SIDE][1][iE3];
				}
				if (iEE3 < 0) {
					iEE3 = neighbors_for_the_internal_node[E_SIDE][2][iE3];
				}
				if (iEE3 < 0) {
					iEE3 = neighbors_for_the_internal_node[E_SIDE][3][iE3];
				}

				if ((iEE3 >= 0) && (iEE3 < maxelm)) {
					// ���������� ����
					SpeedEE3 = potent[GAMMA_LANGTRY_MENTER][iEE3];
					//center_cord3D(iEE3,nvtx,pa,pointP,EE_SIDE);
					pointP = center_coord_loc[iEE3];
					positionxEE3 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iEE3 >= maxelm) && (iEE3 < maxelm + maxbound)) {
						SpeedEE3 = potent[GAMMA_LANGTRY_MENTER][iEE3];
					}
					else {
						SpeedEE3 = SpeedE3;
					}
					//volume3D(iE3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iE3];

					positionxEE3 = positionxE3 + 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iE3 > -1) {
					SpeedE3 = potent[GAMMA_LANGTRY_MENTER][iE3];
					SpeedEE3 = potent[GAMMA_LANGTRY_MENTER][iE3];
				}
				else {
					SpeedE3 = potent[GAMMA_LANGTRY_MENTER][iE];
					SpeedEE3 = potent[GAMMA_LANGTRY_MENTER][iE];
				}
				positionxe3 = positionxP + 0.5 * dx;
				positionxE3 = positionxP + 0.5 * dx;
				positionxEE3 = positionxP + dx; // ����� ���� �� ���������� !
			}

			if ((!bW3) && (iW3 > -1)) {
				//center_cord3D(iW3,nvtx,pa,pointP,W_SIDE);
				pointP = center_coord_loc[iW3];
				positionxW3 = pointP.x;
				positionxw3 = positionxP - 0.5 * dx;
				SpeedW3 = potent[GAMMA_LANGTRY_MENTER][iW3];

				integer iWW3 = neighbors_for_the_internal_node[W_SIDE][0][iW3];
				if (iWW3 < 0) {
					iWW3 = neighbors_for_the_internal_node[W_SIDE][1][iW3];
				}
				if (iWW3 < 0) {
					iWW3 = neighbors_for_the_internal_node[W_SIDE][2][iW3];
				}
				if (iWW3 < 0) {
					iWW3 = neighbors_for_the_internal_node[W_SIDE][3][iW3];
				}

				if ((iWW3 >= 0) && (iWW3 < maxelm)) {
					// ���������� ����
					SpeedWW3 = potent[GAMMA_LANGTRY_MENTER][iWW3];
					//center_cord3D(iWW3,nvtx,pa,pointP,WW_SIDE);
					pointP = center_coord_loc[iWW3];
					positionxWW3 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iWW3 >= maxelm) && (iWW3 < maxelm + maxbound)) {
						SpeedWW3 = potent[GAMMA_LANGTRY_MENTER][iWW3];
					}
					else {
						SpeedWW3 = SpeedW3;
						//std::cout << "iWW3 =" << iWW3 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iW3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iW3];

					positionxWW3 = positionxW3 - 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iW3 > -1) {
					SpeedW3 = potent[GAMMA_LANGTRY_MENTER][iW3]; // Attantion !! Debug
					SpeedWW3 = potent[GAMMA_LANGTRY_MENTER][iW3];
				}
				else {
					SpeedW3 = potent[GAMMA_LANGTRY_MENTER][iW]; // Attantion !! Debug
					SpeedWW3 = potent[GAMMA_LANGTRY_MENTER][iW];
				}
				//printf("SpeedW3==%e\n",SpeedW3); getchar();
				positionxw3 = positionxP - 0.5 * dx;
				positionxW3 = positionxP - 0.5 * dx;
				positionxWW3 = positionxP - dx; // ����� ���� �� ���������� !
			}

			// Y - direction
			if ((!bN3) && (iN3 > -1)) {
				SpeedN3 = potent[GAMMA_LANGTRY_MENTER][iN3];
				//center_cord3D(iN3,nvtx,pa,pointP,N_SIDE);
				pointP = center_coord_loc[iN3];
				positionyN3 = pointP.y;
				positionyn3 = positionxP + 0.5 * dy;

				integer iNN3 = neighbors_for_the_internal_node[N_SIDE][0][iN3];
				if (iNN3 < 0) {
					iNN3 = neighbors_for_the_internal_node[N_SIDE][1][iN3];
				}
				if (iNN3 < 0) {
					iNN3 = neighbors_for_the_internal_node[N_SIDE][2][iN3];
				}
				if (iNN3 < 0) {
					iNN3 = neighbors_for_the_internal_node[N_SIDE][3][iN3];
				}

				if ((iNN3 >= 0) && (iNN3 < maxelm)) {
					// ���������� ����
					SpeedNN3 = potent[GAMMA_LANGTRY_MENTER][iNN3];
					//center_cord3D(iNN3,nvtx,pa,pointP,NN_SIDE);
					pointP = center_coord_loc[iNN3];
					positionyNN3 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iNN3 >= maxelm) && (iNN3 < maxelm + maxbound)) {
						SpeedNN3 = potent[GAMMA_LANGTRY_MENTER][iNN3];
					}
					else {
						SpeedNN3 = SpeedN3;
					}
					//volume3D(iN3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iN3];

					positionyNN3 = positionyN3 + 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iN3 > -1) {
					SpeedN3 = potent[GAMMA_LANGTRY_MENTER][iN3];
					SpeedNN3 = potent[GAMMA_LANGTRY_MENTER][iN3];
				}
				else {
					SpeedN3 = potent[GAMMA_LANGTRY_MENTER][iN];
					SpeedNN3 = potent[GAMMA_LANGTRY_MENTER][iN];
				}
				positionyn3 = positionyP + 0.5 * dy;
				positionyN3 = positionyP + 0.5 * dy;
				positionyNN3 = positionyP + dy; // ����� ���� �� ���������� !
			}

			if ((!bS3) && (iS3 > -1)) {
				SpeedS3 = potent[GAMMA_LANGTRY_MENTER][iS3];
				//center_cord3D(iS3,nvtx,pa,pointP,S_SIDE);
				pointP = center_coord_loc[iS3];
				positionyS3 = pointP.y;
				positionys3 = positionyP - 0.5 * dy;

				integer iSS3 = neighbors_for_the_internal_node[S_SIDE][0][iS3];
				if (iSS3 < 0) {
					iSS3 = neighbors_for_the_internal_node[S_SIDE][1][iS3];
				}
				if (iSS3 < 0) {
					iSS3 = neighbors_for_the_internal_node[S_SIDE][2][iS3];
				}
				if (iSS3 < 0) {
					iSS3 = neighbors_for_the_internal_node[S_SIDE][3][iS3];
				}

				if ((iSS3 >= 0) && (iSS3 < maxelm)) {
					// ���������� ����
					SpeedSS3 = potent[GAMMA_LANGTRY_MENTER][iSS3];
					//center_cord3D(iSS3,nvtx,pa,pointP,SS_SIDE);
					pointP = center_coord_loc[iSS3];
					positionySS3 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iSS3 >= maxelm) && (iSS3 < maxelm + maxbound)) {
						SpeedSS3 = potent[GAMMA_LANGTRY_MENTER][iSS3];
					}
					else {
						SpeedSS3 = SpeedS3;
					}
					//volume3D(iS3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iS3];

					positionySS3 = positionyS3 - 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iS3 > -1) {
					SpeedS3 = potent[GAMMA_LANGTRY_MENTER][iS3]; // ATTANTION !!!!
					SpeedSS3 = potent[GAMMA_LANGTRY_MENTER][iS3]; // ������� �������� ������ ������� ����.
				}
				else {
					SpeedS3 = potent[GAMMA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
					SpeedSS3 = potent[GAMMA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
				}
				positionys3 = positionyP - 0.5 * dy;
				positionyS3 = positionyP - 0.5 * dy;
				positionySS3 = positionyP - dy; // ����� ���� �� ���������� !
			}

			// Z - direction
			if ((!bT3) && (iT3 > -1)) {
				SpeedT3 = potent[GAMMA_LANGTRY_MENTER][iT3];
				//center_cord3D(iT3,nvtx,pa,pointP,T_SIDE);
				pointP = center_coord_loc[iT3];
				positionzT3 = pointP.z;
				positionzt3 = positionzP + 0.5 * dz;

				integer iTT3 = neighbors_for_the_internal_node[T_SIDE][0][iT3];
				if (iTT3 < 0) {
					iTT3 = neighbors_for_the_internal_node[T_SIDE][1][iT3];
				}
				if (iTT3 < 0) {
					iTT3 = neighbors_for_the_internal_node[T_SIDE][2][iT3];
				}
				if (iTT3 < 0) {
					iTT3 = neighbors_for_the_internal_node[T_SIDE][3][iT3];
				}

				if ((iTT3 >= 0) && (iTT3 < maxelm)) {
					// ���������� ����
					SpeedTT3 = potent[GAMMA_LANGTRY_MENTER][iTT3];
					//center_cord3D(iTT3,nvtx,pa,pointP,TT_SIDE);
					pointP = center_coord_loc[iTT3];
					positionzTT3 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iTT3 >= maxelm) && (iTT3 < maxelm + maxbound)) {
						SpeedTT3 = potent[GAMMA_LANGTRY_MENTER][iTT3];
					}
					else {
						SpeedTT3 = SpeedT3;
					}
					//volume3D(iT3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iT3];

					positionzTT3 = positionzT3 + 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iT3 > -1) {
					SpeedT3 = potent[GAMMA_LANGTRY_MENTER][iT3];
					SpeedTT3 = potent[GAMMA_LANGTRY_MENTER][iT3]; // �������� ������ ������� ����
				}
				else {
					SpeedT3 = potent[GAMMA_LANGTRY_MENTER][iT];
					SpeedTT3 = potent[GAMMA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
				}
				positionzt3 = positionzP + 0.5 * dz;
				positionzT3 = positionzP + 0.5 * dz;
				positionzTT3 = positionzP + dz; // ����� ���� �� ���������� !
			}

			if ((!bB3) && (iB3 > -1)) {
				SpeedB3 = potent[GAMMA_LANGTRY_MENTER][iB3];
				//center_cord3D(iB3,nvtx,pa,pointP,B_SIDE);
				pointP = center_coord_loc[iB3];
				positionzB3 = pointP.z;
				positionzb3 = positionzP - 0.5 * dz;

				integer iBB3 = neighbors_for_the_internal_node[B_SIDE][0][iB3];
				if (iBB3 < 0) {
					iBB3 = neighbors_for_the_internal_node[B_SIDE][1][iB3];
				}
				if (iBB3 < 0) {
					iBB3 = neighbors_for_the_internal_node[B_SIDE][2][iB3];
				}
				if (iBB3 < 0) {
					iBB3 = neighbors_for_the_internal_node[B_SIDE][3][iB3];
				}

				if ((iBB3 >= 0) && (iBB3 < maxelm)) {
					// ���������� ����
					SpeedBB3 = potent[GAMMA_LANGTRY_MENTER][iBB3];
					//center_cord3D(iBB3,nvtx,pa,pointP,BB_SIDE);
					pointP = center_coord_loc[iBB3];
					positionzBB3 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iBB3 >= maxelm) && (iBB3 < maxelm + maxbound)) {
						SpeedBB3 = potent[GAMMA_LANGTRY_MENTER][iBB3];
					}
					else {
						SpeedBB3 = SpeedB3;
					}
					//volume3D(iB3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iB3];

					positionzBB3 = positionzB3 - 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iB3 > -1) {
					SpeedB3 = potent[GAMMA_LANGTRY_MENTER][iB3];
					SpeedBB3 = potent[GAMMA_LANGTRY_MENTER][iB3]; // �������� ������ ������� ����
				}
				else {
					SpeedB3 = potent[GAMMA_LANGTRY_MENTER][iB];
					SpeedBB3 = potent[GAMMA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
				}
				positionzb3 = positionzP - 0.5 * dz;
				positionzB3 = positionzP - 0.5 * dz;
				positionzBB3 = positionzP - dz; // ����� ���� �� ���������� !
			}

			// X - direction
			if ((!bE4) && (iE4 > -1)) {
				SpeedE4 = potent[GAMMA_LANGTRY_MENTER][iE4];
				//center_cord3D(iE4,nvtx,pa,pointP,E_SIDE);
				pointP = center_coord_loc[iE4];

				positionxE4 = pointP.x;
				positionxe4 = positionxP + 0.5 * dx;

				integer iEE4 = neighbors_for_the_internal_node[E_SIDE][0][iE4];
				if (iEE4 < 0) {
					iEE4 = neighbors_for_the_internal_node[E_SIDE][1][iE4];
				}
				if (iEE4 < 0) {
					iEE4 = neighbors_for_the_internal_node[E_SIDE][2][iE4];
				}
				if (iEE4 < 0) {
					iEE4 = neighbors_for_the_internal_node[E_SIDE][3][iE4];
				}

				if ((iEE4 >= 0) && (iEE4 < maxelm)) {
					// ���������� ����
					SpeedEE4 = potent[GAMMA_LANGTRY_MENTER][iEE4];
					//center_cord3D(iEE4,nvtx,pa,pointP,EE_SIDE);
					pointP = center_coord_loc[iEE4];
					positionxEE4 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iEE4 >= maxelm) && (iEE4 < maxelm + maxbound)) {
						SpeedEE4 = potent[GAMMA_LANGTRY_MENTER][iEE4];
					}
					else {
						SpeedEE4 = SpeedE4;
					}
					//volume3D(iE4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iE4];

					positionxEE4 = positionxE4 + 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iE4 > -1) {
					SpeedE4 = potent[GAMMA_LANGTRY_MENTER][iE4];
					SpeedEE4 = potent[GAMMA_LANGTRY_MENTER][iE4];
				}
				else {
					SpeedE4 = potent[GAMMA_LANGTRY_MENTER][iE];
					SpeedEE4 = potent[GAMMA_LANGTRY_MENTER][iE];
				}
				positionxe4 = positionxP + 0.5 * dx;
				positionxE4 = positionxP + 0.5 * dx;
				positionxEE4 = positionxP + dx; // ����� ���� �� ���������� !
			}

			if ((!bW4) && (iW4 > -1)) {
				//center_cord3D(iW4,nvtx,pa,pointP,W_SIDE);
				pointP = center_coord_loc[iW4];
				positionxW4 = pointP.x;
				positionxw4 = positionxP - 0.5 * dx;
				SpeedW4 = potent[GAMMA_LANGTRY_MENTER][iW4];

				integer iWW4 = neighbors_for_the_internal_node[W_SIDE][0][iW4];
				if (iWW4 < 0) {
					iWW4 = neighbors_for_the_internal_node[W_SIDE][1][iW4];
				}
				if (iWW4 < 0) {
					iWW4 = neighbors_for_the_internal_node[W_SIDE][2][iW4];
				}
				if (iWW4 < 0) {
					iWW4 = neighbors_for_the_internal_node[W_SIDE][3][iW4];
				}

				if ((iWW4 >= 0) && (iWW4 < maxelm)) {
					// ���������� ����
					SpeedWW4 = potent[GAMMA_LANGTRY_MENTER][iWW4];
					//center_cord3D(iWW4,nvtx,pa,pointP,WW_SIDE);
					pointP = center_coord_loc[iWW4];
					positionxWW4 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iWW4 >= maxelm) && (iWW4 < maxelm + maxbound)) {
						SpeedWW4 = potent[GAMMA_LANGTRY_MENTER][iWW4];
					}
					else {
						SpeedWW4 = SpeedW4;
						//std::cout << "iWW4 =" << iWW4 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iW4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iW4];

					positionxWW4 = positionxW4 - 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iW4 > -1) {
					SpeedW4 = potent[GAMMA_LANGTRY_MENTER][iW4]; // Attantion !! Debug
					SpeedWW4 = potent[GAMMA_LANGTRY_MENTER][iW4];
				}
				else {
					SpeedW4 = potent[GAMMA_LANGTRY_MENTER][iW]; // Attantion !! Debug
					SpeedWW4 = potent[GAMMA_LANGTRY_MENTER][iW];
				}
				//printf("SpeedW4==%e\n",SpeedW4); getchar();
				positionxw4 = positionxP - 0.5 * dx;
				positionxW4 = positionxP - 0.5 * dx;
				positionxWW4 = positionxP - dx; // ����� ���� �� ���������� !
			}

			// Y - direction
			if ((!bN4) && (iN4 > -1)) {
				SpeedN4 = potent[GAMMA_LANGTRY_MENTER][iN4];
				//center_cord3D(iN4,nvtx,pa,pointP,N_SIDE);
				pointP = center_coord_loc[iN4];
				positionyN4 = pointP.y;
				positionyn4 = positionxP + 0.5 * dy;

				integer iNN4 = neighbors_for_the_internal_node[N_SIDE][0][iN4];
				if (iNN4 < 0) {
					iNN4 = neighbors_for_the_internal_node[N_SIDE][1][iN4];
				}
				if (iNN4 < 0) {
					iNN4 = neighbors_for_the_internal_node[N_SIDE][2][iN4];
				}
				if (iNN4 < 0) {
					iNN4 = neighbors_for_the_internal_node[N_SIDE][3][iN4];
				}

				if ((iNN4 >= 0) && (iNN4 < maxelm)) {
					// ���������� ����
					SpeedNN4 = potent[GAMMA_LANGTRY_MENTER][iNN4];
					//center_cord3D(iNN4,nvtx,pa,pointP,NN_SIDE);
					pointP = center_coord_loc[iNN4];
					positionyNN4 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iNN4 >= maxelm) && (iNN4 < maxelm + maxbound)) {
						SpeedNN4 = potent[GAMMA_LANGTRY_MENTER][iNN4];
					}
					else {
						SpeedNN4 = SpeedN4;
					}
					//volume3D(iN4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iN4];

					positionyNN4 = positionyN4 + 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iN4 > -1) {
					SpeedN4 = potent[GAMMA_LANGTRY_MENTER][iN4];
					SpeedNN4 = potent[GAMMA_LANGTRY_MENTER][iN4];
				}
				else {
					SpeedN4 = potent[GAMMA_LANGTRY_MENTER][iN];
					SpeedNN4 = potent[GAMMA_LANGTRY_MENTER][iN];
				}
				positionyn4 = positionyP + 0.5 * dy;
				positionyN4 = positionyP + 0.5 * dy;
				positionyNN4 = positionyP + dy; // ����� ���� �� ���������� !
			}

			if ((!bS4) && (iS4 > -1)) {
				SpeedS4 = potent[GAMMA_LANGTRY_MENTER][iS4];
				//center_cord3D(iS4,nvtx,pa,pointP,S_SIDE);
				pointP = center_coord_loc[iS4];
				positionyS4 = pointP.y;
				positionys4 = positionyP - 0.5 * dy;

				integer iSS4 = neighbors_for_the_internal_node[S_SIDE][0][iS4];
				if (iSS4 < 0) {
					iSS4 = neighbors_for_the_internal_node[S_SIDE][1][iS4];
				}
				if (iSS4 < 0) {
					iSS4 = neighbors_for_the_internal_node[S_SIDE][2][iS4];
				}
				if (iSS4 < 0) {
					iSS4 = neighbors_for_the_internal_node[S_SIDE][3][iS4];
				}

				if ((iSS4 >= 0) && (iSS4 < maxelm)) {
					// ���������� ����
					SpeedSS4 = potent[GAMMA_LANGTRY_MENTER][iSS4];
					//center_cord3D(iSS4,nvtx,pa,pointP,SS_SIDE);
					pointP = center_coord_loc[iSS4];
					positionySS4 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iSS4 >= maxelm) && (iSS4 < maxelm + maxbound)) {
						SpeedSS4 = potent[GAMMA_LANGTRY_MENTER][iSS4];
					}
					else {
						SpeedSS4 = SpeedS4;
					}
					//volume3D(iS4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iS4];

					positionySS4 = positionyS4 - 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iS4 > -1) {
					SpeedS4 = potent[GAMMA_LANGTRY_MENTER][iS4]; // ATTANTION !!!!
					SpeedSS4 = potent[GAMMA_LANGTRY_MENTER][iS4]; // ������� �������� ������ ������� ����.
				}
				else {
					SpeedS4 = potent[GAMMA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
					SpeedSS4 = potent[GAMMA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
				}
				positionys4 = positionyP - 0.5 * dy;
				positionyS4 = positionyP - 0.5 * dy;
				positionySS4 = positionyP - dy; // ����� ���� �� ���������� !
			}

			// Z - direction
			if ((!bT4) && (iT4 > -1)) {
				SpeedT4 = potent[GAMMA_LANGTRY_MENTER][iT4];
				//center_cord3D(iT4,nvtx,pa,pointP,T_SIDE);
				pointP = center_coord_loc[iT4];
				positionzT4 = pointP.z;
				positionzt4 = positionzP + 0.5 * dz;

				integer iTT4 = neighbors_for_the_internal_node[T_SIDE][0][iT4];
				if (iTT4 < 0) {
					iTT4 = neighbors_for_the_internal_node[T_SIDE][1][iT4];
				}
				if (iTT4 < 0) {
					iTT4 = neighbors_for_the_internal_node[T_SIDE][2][iT4];
				}
				if (iTT4 < 0) {
					iTT4 = neighbors_for_the_internal_node[T_SIDE][3][iT4];
				}

				if ((iTT4 >= 0) && (iTT4 < maxelm)) {
					// ���������� ����
					SpeedTT4 = potent[GAMMA_LANGTRY_MENTER][iTT4];
					//center_cord3D(iTT4,nvtx,pa,pointP,TT_SIDE);
					pointP = center_coord_loc[iTT4];
					positionzTT4 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iTT4 >= maxelm) && (iTT4 < maxelm + maxbound)) {
						SpeedTT4 = potent[GAMMA_LANGTRY_MENTER][iTT4];
					}
					else {
						SpeedTT4 = SpeedT4;
					}
					//volume3D(iT4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iT4];

					positionzTT4 = positionzT4 + 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iT4 > -1) {
					SpeedT4 = potent[GAMMA_LANGTRY_MENTER][iT4];
					SpeedTT4 = potent[GAMMA_LANGTRY_MENTER][iT4]; // �������� ������ ������� ����
				}
				else {
					SpeedT4 = potent[GAMMA_LANGTRY_MENTER][iT];
					SpeedTT4 = potent[GAMMA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
				}
				positionzt4 = positionzP + 0.5 * dz;
				positionzT4 = positionzP + 0.5 * dz;
				positionzTT4 = positionzP + dz; // ����� ���� �� ���������� !
			}

			if ((!bB4) && (iB4 > -1)) {
				SpeedB4 = potent[GAMMA_LANGTRY_MENTER][iB4];
				//center_cord3D(iB4,nvtx,pa,pointP,B_SIDE);
				pointP = center_coord_loc[iB4];
				positionzB4 = pointP.z;
				positionzb4 = positionzP - 0.5 * dz;

				integer iBB4 = neighbors_for_the_internal_node[B_SIDE][0][iB4];
				if (iBB4 < 0) {
					iBB4 = neighbors_for_the_internal_node[B_SIDE][1][iB4];
				}
				if (iBB4 < 0) {
					iBB4 = neighbors_for_the_internal_node[B_SIDE][2][iB4];
				}
				if (iBB4 < 0) {
					iBB4 = neighbors_for_the_internal_node[B_SIDE][3][iB4];
				}

				if ((iBB4 >= 0) && (iBB4 < maxelm)) {
					// ���������� ����
					SpeedBB4 = potent[GAMMA_LANGTRY_MENTER][iBB4];
					//center_cord3D(iBB4,nvtx,pa,pointP,BB_SIDE);
					pointP = center_coord_loc[iBB4];
					positionzBB4 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iBB4 >= maxelm) && (iBB4 < maxelm + maxbound)) {
						SpeedBB4 = potent[GAMMA_LANGTRY_MENTER][iBB4];
					}
					else {
						SpeedBB4 = SpeedB4;
					}
					//volume3D(iB4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iB4];

					positionzBB4 = positionzB4 - 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iB4 > -1) {
					SpeedB4 = potent[GAMMA_LANGTRY_MENTER][iB4];
					SpeedBB4 = potent[GAMMA_LANGTRY_MENTER][iB4]; // �������� ������ ������� ����
				}
				else {
					SpeedB4 = potent[GAMMA_LANGTRY_MENTER][iB];
					SpeedBB4 = potent[GAMMA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
				}
				positionzb4 = positionzP - 0.5 * dz;
				positionzB4 = positionzP - 0.5 * dz;
				positionzBB4 = positionzP - dz; // ����� ���� �� ���������� !
			}

		}



		if ((ishconvection >= QUICK) && (ishconvection <= FROMM)) {
			// ������ ����� ������������ �� ��������� �. ��������� PHOENICS.
		   // ������������� ishconvection ������ ��������� �������� �� ����� QUICK (������������) � ������ �� ����� UNEVENQUICK (�� ������� ���������).

			// X - direction
			Speede = cell_face_value_global(ishconvection, (Fe), SpeedW, SpeedP, SpeedE, SpeedEE);
			Speedw = cell_face_value_global(ishconvection, (Fw), SpeedWW, SpeedW, SpeedP, SpeedE);
			// Y - direction
			Speedn = cell_face_value_global(ishconvection, (Fn), SpeedS, SpeedP, SpeedN, SpeedNN);
			Speeds = cell_face_value_global(ishconvection, (Fs), SpeedSS, SpeedS, SpeedP, SpeedN);
			// Z - direction
			Speedt = cell_face_value_global(ishconvection, (Ft), SpeedB, SpeedP, SpeedT, SpeedTT);
			Speedb = cell_face_value_global(ishconvection, (Fb), SpeedBB, SpeedB, SpeedP, SpeedT);

			if (b_on_adaptive_local_refinement_mesh) {

				// X - direction
				Speede2 = cell_face_value_global(ishconvection, (Fe2), SpeedW2, SpeedP, SpeedE2, SpeedEE2);
				Speedw2 = cell_face_value_global(ishconvection, (Fw2), SpeedWW2, SpeedW2, SpeedP, SpeedE2);
				// Y - direction
				Speedn2 = cell_face_value_global(ishconvection, (Fn2), SpeedS2, SpeedP, SpeedN2, SpeedNN2);
				Speeds2 = cell_face_value_global(ishconvection, (Fs2), SpeedSS2, SpeedS2, SpeedP, SpeedN2);
				// Z - direction
				Speedt2 = cell_face_value_global(ishconvection, (Ft2), SpeedB2, SpeedP, SpeedT2, SpeedTT2);
				Speedb2 = cell_face_value_global(ishconvection, (Fb2), SpeedBB2, SpeedB2, SpeedP, SpeedT2);


				// X - direction
				Speede3 = cell_face_value_global(ishconvection, (Fe3), SpeedW3, SpeedP, SpeedE3, SpeedEE3);
				Speedw3 = cell_face_value_global(ishconvection, (Fw3), SpeedWW3, SpeedW3, SpeedP, SpeedE3);
				// Y - direction
				Speedn3 = cell_face_value_global(ishconvection, (Fn3), SpeedS3, SpeedP, SpeedN3, SpeedNN3);
				Speeds3 = cell_face_value_global(ishconvection, (Fs3), SpeedSS3, SpeedS3, SpeedP, SpeedN3);
				// Z - direction
				Speedt3 = cell_face_value_global(ishconvection, (Ft3), SpeedB3, SpeedP, SpeedT3, SpeedTT3);
				Speedb3 = cell_face_value_global(ishconvection, (Fb3), SpeedBB3, SpeedB3, SpeedP, SpeedT3);


				// X - direction
				Speede4 = cell_face_value_global(ishconvection, (Fe4), SpeedW4, SpeedP, SpeedE4, SpeedEE4);
				Speedw4 = cell_face_value_global(ishconvection, (Fw4), SpeedWW4, SpeedW4, SpeedP, SpeedE4);
				// Y - direction
				Speedn4 = cell_face_value_global(ishconvection, (Fn4), SpeedS4, SpeedP, SpeedN4, SpeedNN4);
				Speeds4 = cell_face_value_global(ishconvection, (Fs4), SpeedSS4, SpeedS4, SpeedP, SpeedN4);
				// Z - direction
				Speedt4 = cell_face_value_global(ishconvection, (Ft4), SpeedB4, SpeedP, SpeedT4, SpeedTT4);
				Speedb4 = cell_face_value_global(ishconvection, (Fb4), SpeedBB4, SpeedB4, SpeedP, SpeedT4);

			}

		}

		if (ishconvection >= UNEVENQUICK) {


			/*
			// ������������������ �������� ��������� � ����� ���������� ���������� ����� QUICK �� ������������� �����.
			// ���������� ���� ������������ �� ������: ...
			// � ������ ������ ������ ���������� �� ������������.
			//doublereal gamma1E, gamma2E, gamma1W, gamma2W, delta1E, delta2E, delta1W, delta2W;
			//doublereal gamma1N, gamma2N, gamma1S, gamma2S, delta1N, delta2N, delta1S, delta2S;
			//doublereal gamma1T, gamma2T, gamma1B, gamma2B, delta1T, delta2T, delta1B, delta2B;
			// X - direction
			// gamma
			//gamma1E=((positionxe-positionxE)*(positionxe-positionxP))/((positionxW-positionxE)*(positionxW-positionxP));
			//gamma2E=((positionxe-positionxP)*(positionxe-positionxW))/((positionxE-positionxP)*(positionxE-positionxW));
			//gamma1W=((positionxw-positionxP)*(positionxw-positionxW))/((positionxWW-positionxP)*(positionxWW-positionxW));
			//gamma2W=((positionxw-positionxW)*(positionxw-positionxWW))/((positionxP-positionxW)*(positionxP-positionxWW));
			// delta
			//delta1E=((positionxe-positionxEE)*(positionxe-positionxE))/((positionxP-positionxEE)*(positionxP-positionxE));
			//delta2E=((positionxe-positionxE)*(positionxe-positionxP))/((positionxEE-positionxE)*(positionxEE-positionxP));
			//delta1W=((positionxw-positionxE)*(positionxw-positionxP))/((positionxW-positionxE)*(positionxW-positionxP));
			//delta2W=((positionxw-positionxP)*(positionxw-positionxW))/((positionxE-positionxP)*(positionxE-positionxW));
			// Y - direction
			// gamma
			//gamma1N=((positionyn-positionyN)*(positionyn-positionyP))/((positionyS-positionyN)*(positionyS-positionyP));
			//gamma2N=((positionyn-positionyP)*(positionyn-positionyS))/((positionyN-positionyP)*(positionyN-positionyS));
			//gamma1S=((positionys-positionyP)*(positionys-positionyS))/((positionySS-positionyP)*(positionySS-positionyS));
			//gamma2S=((positionys-positionyS)*(positionys-positionySS))/((positionyP-positionyS)*(positionyP-positionySS));
			// delta
			//delta1N=((positionyn-positionyNN)*(positionyn-positionyN))/((positionyP-positionyNN)*(positionyP-positionyN));
			//delta2N=((positionyn-positionyN)*(positionyn-positionyP))/((positionyNN-positionyN)*(positionyNN-positionyP));
			//delta1S=((positionys-positionyN)*(positionys-positionyP))/((positionyS-positionyN)*(positionyS-positionyP));
			//delta2S=((positionys-positionyP)*(positionys-positionyS))/((positionyN-positionyP)*(positionyN-positionyS));
			// Z - direction
			// gamma
			//gamma1T=((positionzt-positionzT)*(positionzt-positionzP))/((positionzB-positionzT)*(positionzB-positionzP));
			//gamma2T=((positionzt-positionzP)*(positionzt-positionzB))/((positionzT-positionzP)*(positionzT-positionzB));
			//gamma1B=((positionzb-positionzP)*(positionzb-positionzB))/((positionzBB-positionzP)*(positionzBB-positionzB));
			//gamma2B=((positionzb-positionzB)*(positionzb-positionzBB))/((positionzP-positionzB)*(positionzP-positionzBB));
			// delta
			//delta1T=((positionzt-positionzTT)*(positionzt-positionzT))/((positionzP-positionzTT)*(positionzP-positionzT));
			//delta2T=((positionzt-positionzT)*(positionzt-positionzP))/((positionzTT-positionzT)*(positionzTT-positionzP));
			//delta1B=((positionzb-positionzT)*(positionzb-positionzP))/((positionzB-positionzT)*(positionzB-positionzP));
			//delta2B=((positionzb-positionzP)*(positionzb-positionzB))/((positionzT-positionzP)*(positionzT-positionzB));
			*/




			// ���������� ������� �������� �� ����� ��
			// ������������ ����� �������� QUICK.
			/* ������� ������������:
			 *  A	B	C	D	e	+/-
			 *  W	P	E	-	e	+
			 *  -	P	E	EE  e   -
			 *  WW   W	P	-	w	+
			 *  -	W	P	E	w	-
			 *  S	P	N	-	n	+
			 *  -	P	N	NN  n	-
			 *  SS   S	P	-	s	+
			 *  -	S	P	N	s	-
			 *  B	P	T	-	t	+
			 *  -	P	T	TT  t	-
			 *  BB   B	P	-	b	+
			 *  -	B	P	T	b	-
			 */

			if (ishconvection == UNEVENQUICK) {
				// X - direction
				Speede = workQUICK(dx, 2.0 * (positionxE - positionxe), positionxW, positionxP, positionxE, positionxEE, SpeedW, SpeedP, SpeedE, SpeedEE, (Fe));
				Speedw = workQUICK(2.0 * (positionxw - positionxW), dx, positionxWW, positionxW, positionxP, positionxE, SpeedWW, SpeedW, SpeedP, SpeedE, (Fw));
				// Y - direction
				Speedn = workQUICK(dy, 2.0 * (positionyN - positionyn), positionyS, positionyP, positionyN, positionyNN, SpeedS, SpeedP, SpeedN, SpeedNN, (Fn));
				Speeds = workQUICK(2.0 * (positionys - positionyS), dy, positionySS, positionyS, positionyP, positionyN, SpeedSS, SpeedS, SpeedP, SpeedN, (Fs));
				// Z - direction
				Speedt = workQUICK(dz, 2.0 * (positionzT - positionzt), positionzB, positionzP, positionzT, positionzTT, SpeedB, SpeedP, SpeedT, SpeedTT, (Ft));
				Speedb = workQUICK(2.0 * (positionzb - positionzB), dz, positionzBB, positionzB, positionzP, positionzT, SpeedBB, SpeedB, SpeedP, SpeedT, (Fb));


				if (b_on_adaptive_local_refinement_mesh) {

					// X - direction
					Speede2 = workQUICK(dx, 2.0 * (positionxE2 - positionxe2), positionxW2, positionxP, positionxE2, positionxEE2, SpeedW2, SpeedP, SpeedE2, SpeedEE2, (Fe2));
					Speedw2 = workQUICK(2.0 * (positionxw2 - positionxW2), dx, positionxWW2, positionxW2, positionxP, positionxE2, SpeedWW2, SpeedW2, SpeedP, SpeedE2, (Fw2));
					// Y - direction
					Speedn2 = workQUICK(dy, 2.0 * (positionyN2 - positionyn2), positionyS2, positionyP, positionyN2, positionyNN2, SpeedS2, SpeedP, SpeedN2, SpeedNN2, (Fn2));
					Speeds2 = workQUICK(2.0 * (positionys2 - positionyS2), dy, positionySS2, positionyS2, positionyP, positionyN2, SpeedSS2, SpeedS2, SpeedP, SpeedN2, (Fs2));
					// Z - direction
					Speedt2 = workQUICK(dz, 2.0 * (positionzT2 - positionzt2), positionzB2, positionzP, positionzT2, positionzTT2, SpeedB2, SpeedP, SpeedT2, SpeedTT2, (Ft2));
					Speedb2 = workQUICK(2.0 * (positionzb2 - positionzB2), dz, positionzBB2, positionzB2, positionzP, positionzT2, SpeedBB2, SpeedB2, SpeedP, SpeedT2, (Fb2));


					// X - direction
					Speede3 = workQUICK(dx, 2.0 * (positionxE3 - positionxe3), positionxW3, positionxP, positionxE3, positionxEE3, SpeedW3, SpeedP, SpeedE3, SpeedEE3, (Fe3));
					Speedw3 = workQUICK(2.0 * (positionxw3 - positionxW3), dx, positionxWW3, positionxW3, positionxP, positionxE3, SpeedWW3, SpeedW3, SpeedP, SpeedE3, (Fw3));
					// Y - direction
					Speedn3 = workQUICK(dy, 2.0 * (positionyN3 - positionyn3), positionyS3, positionyP, positionyN3, positionyNN3, SpeedS3, SpeedP, SpeedN3, SpeedNN3, (Fn3));
					Speeds3 = workQUICK(2.0 * (positionys3 - positionyS3), dy, positionySS3, positionyS3, positionyP, positionyN3, SpeedSS3, SpeedS3, SpeedP, SpeedN3, (Fs3));
					// Z - direction
					Speedt3 = workQUICK(dz, 2.0 * (positionzT3 - positionzt3), positionzB3, positionzP, positionzT3, positionzTT3, SpeedB3, SpeedP, SpeedT3, SpeedTT3, (Ft3));
					Speedb3 = workQUICK(2.0 * (positionzb3 - positionzB3), dz, positionzBB3, positionzB3, positionzP, positionzT3, SpeedBB3, SpeedB3, SpeedP, SpeedT3, (Fb3));


					// X - direction
					Speede4 = workQUICK(dx, 2.0 * (positionxE4 - positionxe4), positionxW4, positionxP, positionxE4, positionxEE4, SpeedW4, SpeedP, SpeedE4, SpeedEE4, (Fe4));
					Speedw4 = workQUICK(2.0 * (positionxw4 - positionxW4), dx, positionxWW4, positionxW4, positionxP, positionxE4, SpeedWW4, SpeedW4, SpeedP, SpeedE4, (Fw4));
					// Y - direction
					Speedn4 = workQUICK(dy, 2.0 * (positionyN4 - positionyn4), positionyS4, positionyP, positionyN4, positionyNN4, SpeedS4, SpeedP, SpeedN4, SpeedNN4, (Fn4));
					Speeds4 = workQUICK(2.0 * (positionys4 - positionyS4), dy, positionySS4, positionyS4, positionyP, positionyN4, SpeedSS4, SpeedS4, SpeedP, SpeedN4, (Fs4));
					// Z - direction
					Speedt4 = workQUICK(dz, 2.0 * (positionzT4 - positionzt4), positionzB4, positionzP, positionzT4, positionzTT4, SpeedB4, SpeedP, SpeedT4, SpeedTT4, (Ft4));
					Speedb4 = workQUICK(2.0 * (positionzb4 - positionzB4), dz, positionzBB4, positionzB4, positionzP, positionzT4, SpeedBB4, SpeedB4, SpeedP, SpeedT4, (Fb4));

				}

			}

			if ((ishconvection > UNEVENQUICK) && (ishconvection <= UNEVEN_CUBISTA)) {
				// ���� �� ������ ������ ������������� ����������� ������������ ������ ������ ������ �����:
				// 1. UNEVEN_MUSCL, 2. UNEVEN_SOUCUP, 3. UNEVEN_HLPA, 4. UNEVEN_SMART.
				// ������������� ����� ������ ��������������� ��������.

				// X - direction
				Speede = workKN_VOLKOV(positionxW, positionxP, positionxE, positionxEE, SpeedW, SpeedP, SpeedE, SpeedEE, (Fe), ishconvection);
				Speedw = workKN_VOLKOV(positionxWW, positionxW, positionxP, positionxE, SpeedWW, SpeedW, SpeedP, SpeedE, (Fw), ishconvection);
				// Y - direction
				Speedn = workKN_VOLKOV(positionyS, positionyP, positionyN, positionyNN, SpeedS, SpeedP, SpeedN, SpeedNN, (Fn), ishconvection);
				Speeds = workKN_VOLKOV(positionySS, positionyS, positionyP, positionyN, SpeedSS, SpeedS, SpeedP, SpeedN, (Fs), ishconvection);
				// Z - direction
				Speedt = workKN_VOLKOV(positionzB, positionzP, positionzT, positionzTT, SpeedB, SpeedP, SpeedT, SpeedTT, (Ft), ishconvection);
				Speedb = workKN_VOLKOV(positionzBB, positionzB, positionzP, positionzT, SpeedBB, SpeedB, SpeedP, SpeedT, (Fb), ishconvection);

				// debug ������ �������� ������.
				//printf("%f, %f, %f, %f, %f, %f\n",Speede,Speedw,Speedn,Speeds,Speedt,Speedb);
				//getchar();


				if (b_on_adaptive_local_refinement_mesh) {

					// X - direction
					Speede2 = workKN_VOLKOV(positionxW2, positionxP, positionxE2, positionxEE2, SpeedW2, SpeedP, SpeedE2, SpeedEE2, (Fe2), ishconvection);
					Speedw2 = workKN_VOLKOV(positionxWW2, positionxW2, positionxP, positionxE2, SpeedWW2, SpeedW2, SpeedP, SpeedE2, (Fw2), ishconvection);
					// Y - direction
					Speedn2 = workKN_VOLKOV(positionyS2, positionyP, positionyN2, positionyNN2, SpeedS2, SpeedP, SpeedN2, SpeedNN2, (Fn2), ishconvection);
					Speeds2 = workKN_VOLKOV(positionySS2, positionyS2, positionyP, positionyN2, SpeedSS2, SpeedS2, SpeedP, SpeedN2, (Fs2), ishconvection);
					// Z - direction
					Speedt2 = workKN_VOLKOV(positionzB2, positionzP, positionzT2, positionzTT2, SpeedB2, SpeedP, SpeedT2, SpeedTT2, (Ft2), ishconvection);
					Speedb2 = workKN_VOLKOV(positionzBB2, positionzB2, positionzP, positionzT2, SpeedBB2, SpeedB2, SpeedP, SpeedT2, (Fb2), ishconvection);


					// X - direction
					Speede3 = workKN_VOLKOV(positionxW3, positionxP, positionxE3, positionxEE3, SpeedW3, SpeedP, SpeedE3, SpeedEE3, (Fe3), ishconvection);
					Speedw3 = workKN_VOLKOV(positionxWW3, positionxW3, positionxP, positionxE3, SpeedWW3, SpeedW3, SpeedP, SpeedE3, (Fw3), ishconvection);
					// Y - direction
					Speedn3 = workKN_VOLKOV(positionyS3, positionyP, positionyN3, positionyNN3, SpeedS3, SpeedP, SpeedN3, SpeedNN3, (Fn3), ishconvection);
					Speeds3 = workKN_VOLKOV(positionySS3, positionyS3, positionyP, positionyN3, SpeedSS3, SpeedS3, SpeedP, SpeedN3, (Fs3), ishconvection);
					// Z - direction
					Speedt3 = workKN_VOLKOV(positionzB3, positionzP, positionzT3, positionzTT3, SpeedB3, SpeedP, SpeedT3, SpeedTT3, (Ft3), ishconvection);
					Speedb3 = workKN_VOLKOV(positionzBB3, positionzB3, positionzP, positionzT3, SpeedBB3, SpeedB3, SpeedP, SpeedT3, (Fb3), ishconvection);


					// X - direction
					Speede4 = workKN_VOLKOV(positionxW4, positionxP, positionxE4, positionxEE4, SpeedW4, SpeedP, SpeedE4, SpeedEE4, (Fe4), ishconvection);
					Speedw4 = workKN_VOLKOV(positionxWW4, positionxW4, positionxP, positionxE4, SpeedWW4, SpeedW4, SpeedP, SpeedE4, (Fw4), ishconvection);
					// Y - direction
					Speedn4 = workKN_VOLKOV(positionyS4, positionyP, positionyN4, positionyNN4, SpeedS4, SpeedP, SpeedN4, SpeedNN4, (Fn4), ishconvection);
					Speeds4 = workKN_VOLKOV(positionySS4, positionyS4, positionyP, positionyN4, SpeedSS4, SpeedS4, SpeedP, SpeedN4, (Fs4), ishconvection);
					// Z - direction
					Speedt4 = workKN_VOLKOV(positionzB4, positionzP, positionzT4, positionzTT4, SpeedB4, SpeedP, SpeedT4, SpeedTT4, (Ft4), ishconvection);
					Speedb4 = workKN_VOLKOV(positionzBB4, positionzB4, positionzP, positionzT4, SpeedBB4, SpeedB4, SpeedP, SpeedT4, (Fb4), ishconvection);

				}

			}

		} // endif (ishconvection >= UNEVENQUICK)



		// ������: SIMPLE method for the solution of incompressible flows on non-staggered grids
		// I. Sezai - Eastern Mediterranean University, Mechanical Engineering Department, Mersin 10-Turkey Revised in January, 2011.

		// ���������� ������������� ����������� �������:
		// ����������� ����� ���������� ���������:
		// ������ ����������� ������ ��������������� �����, 
		// � ���������� ����� ������������ � ������ ����� 
		// �������� ������� ���������.
		if (1) {

			/*
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae=De*fD(Pe, EXP2, true, feplus) + fmax(-(Fe),0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw=Dw*fD(Pw, EXP2, true, fwplus) + fmax((Fw),0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an=Dn*fD(Pn, EXP2, true, fnplus) + fmax(-(Fn),0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as=Ds*fD(Ps, EXP2, true, fsplus) + fmax((Fs),0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at=Dt*fD(Pt, EXP2, true, ftplus) + fmax(-(Ft),0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab=Db*fD(Pb, EXP2, true, fbplus) + fmax((Fb),0.0);
			*/

			// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
			// ����� ������ UDS.
			// ��� ����������� � ����������.
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De + fmax(-(Fe), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw + fmax((Fw), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn + fmax(-(Fn), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds + fmax((Fs), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt + fmax(-(Ft), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db + fmax((Fb), 0.0);


			// 08.05.2017.
			// ��� ���������:
			// ����� ������������� !!! (���������� �� ����).
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = De + fmax(+(Fe), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw + fmax(-(Fw), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn + fmax(+(Fn), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds + fmax(-(Fs), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt + fmax(+(Ft), 0.0);
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db + fmax(-(Fb), 0.0);

			/*
			sumanb = De + fmax(+(Fe), 0.0);
			sumanb += Dw + fmax(-(Fw), 0.0);
			sumanb += Dn + fmax(+(Fn), 0.0);
			sumanb += Ds + fmax(-(Fs), 0.0);
			sumanb += Dt + fmax(+(Ft), 0.0);
			sumanb += Db + fmax(-(Fb), 0.0);
			*/


			if (b_on_adaptive_local_refinement_mesh) {

				// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
			 // ����� ������ UDS.
			 // ��� ����������� � ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fmax(-(Fe2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fmax((Fw2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fmax(-(Fn2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fmax((Fs2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fmax(-(Ft2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fmax((Fb2), 0.0);


				// 08.05.2017.
				// ��� ���������:
				// ����� ������������� !!! (���������� �� ����).
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += De2 + fmax(+(Fe2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw2 + fmax(-(Fw2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn2 + fmax(+(Fn2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds2 + fmax(-(Fs2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt2 + fmax(+(Ft2), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db2 + fmax(-(Fb2), 0.0);

				/*
				sumanb += De2 + fmax(+(Fe2), 0.0);
				sumanb += Dw2 + fmax(-(Fw2), 0.0);
				sumanb += Dn2 + fmax(+(Fn2), 0.0);
				sumanb += Ds2 + fmax(-(Fs2), 0.0);
				sumanb += Dt2 + fmax(+(Ft2), 0.0);
				sumanb += Db2 + fmax(-(Fb2), 0.0);
				*/


				// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
		   // ����� ������ UDS.
		   // ��� ����������� � ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fmax(-(Fe3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fmax((Fw3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fmax(-(Fn3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fmax((Fs3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fmax(-(Ft3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fmax((Fb3), 0.0);


				// 08.05.2017.
				// ��� ���������:
				// ����� ������������� !!! (���������� �� ����).
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += De3 + fmax(+(Fe3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw3 + fmax(-(Fw3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn3 + fmax(+(Fn3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds3 + fmax(-(Fs3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt3 + fmax(+(Ft3), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db3 + fmax(-(Fb3), 0.0);

				/*
				sumanb += De3 + fmax(+(Fe3), 0.0);
				sumanb += Dw3 + fmax(-(Fw3), 0.0);
				sumanb += Dn3 + fmax(+(Fn3), 0.0);
				sumanb += Ds3 + fmax(-(Fs3), 0.0);
				sumanb += Dt3 + fmax(+(Ft3), 0.0);
				sumanb += Db3 + fmax(-(Fb3), 0.0);
				*/


				// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
		   // ����� ������ UDS.
		   // ��� ����������� � ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fmax(-(Fe4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fmax((Fw4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fmax(-(Fn4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fmax((Fs4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fmax(-(Ft4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fmax((Fb4), 0.0);


				// 08.05.2017.
				// ��� ���������:
				// ����� ������������� !!! (���������� �� ����).
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += De4 + fmax(+(Fe4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dw4 + fmax(-(Fw4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dn4 + fmax(+(Fn4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Ds4 + fmax(-(Fs4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Dt4 + fmax(+(Ft4), 0.0);
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += Db4 + fmax(-(Fb4), 0.0);

				/*
				sumanb += De4 + fmax(+(Fe4), 0.0);
				sumanb += Dw4 + fmax(-(Fw4), 0.0);
				sumanb += Dn4 + fmax(+(Fn4), 0.0);
				sumanb += Ds4 + fmax(-(Fs4), 0.0);
				sumanb += Dt4 + fmax(+(Ft4), 0.0);
				sumanb += Db4 + fmax(-(Fb4), 0.0);
				*/

			}

		}
		else {

			// 8.05.2017
			// �� � ���� ������ �� �������� ��� ����� ����.


			// ��� ������ ������ �� ��������, ������� ���� � ���������� � ���� ������� ����������,
			// �� ��� ����� ���������� ������������ ��� �� ������������� ����������� ������.
		   // 30 07 2015
		   // 
		   // ������ ������ ������� ����� ���������� �� UDS.
			if (!bE) {
				// ������ ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De + fmax(-(Fe), 0.0);
			}
			else {
				integer inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De + fabs(Fe);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae = De + fabs(Fe);
				}
			}

			if (!bW) {
				// ������ ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw + fmax((Fw), 0.0);
			}
			else {
				integer inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw + fabs(Fw);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw = Dw + fabs(Fw);
				}
			}

			if (!bN) {
				// ������ ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn + fmax(-(Fn), 0.0);
			}
			else {
				integer inumber = iN - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn + fabs(Fn);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an = Dn + fabs(Fn);
				}
			}

			if (!bS) {
				// ������ ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds + fmax((Fs), 0.0);
			}
			else {
				integer inumber = iS - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds + fabs(Fs);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as = Ds + fabs(Fs);
				}
			}

			if (!bT) {
				// ������ ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt + fmax(-(Ft), 0.0);
			}
			else {
				integer inumber = iT - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt + fabs(Ft);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at = Dt + fabs(Ft);
				}
			}

			if (!bB) {
				// ������ ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db + fmax((Fb), 0.0);
			}
			else {
				integer inumber = iB - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db + fabs(Fb);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab = Db + fabs(Fb);
				}
			}


			if (b_on_adaptive_local_refinement_mesh) {

				// ������ ������ ������� ����� ���������� �� UDS.
				if (!bE2) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fmax(-(Fe2), 0.0);
				}
				else {
					integer inumber = iE2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fabs(Fe2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fabs(Fe2);
					}
				}

				if (!bW2) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fmax((Fw2), 0.0);
				}
				else {
					integer inumber = iW2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fabs(Fw2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fabs(Fw2);
					}
				}

				if (!bN2) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fmax(-(Fn2), 0.0);
				}
				else {
					integer inumber = iN2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fabs(Fn2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fabs(Fn2);
					}
				}

				if (!bS2) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fmax((Fs2), 0.0);
				}
				else {
					integer inumber = iS2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fabs(Fs2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fabs(Fs2);
					}
				}

				if (!bT2) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fmax(-(Ft2), 0.0);
				}
				else {
					integer inumber = iT2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fabs(Ft2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fabs(Ft2);
					}
				}

				if (!bB2) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fmax((Fb2), 0.0);
				}
				else {
					integer inumber = iB2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fabs(Fb2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fabs(Fb2);
					}
				}


				// ������ ������ ������� ����� ���������� �� UDS.
				if (!bE3) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fmax(-(Fe3), 0.0);
				}
				else {
					integer inumber = iE3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fabs(Fe3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fabs(Fe3);
					}
				}

				if (!bW3) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fmax((Fw3), 0.0);
				}
				else {
					integer inumber = iW3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fabs(Fw3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fabs(Fw3);
					}
				}

				if (!bN3) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fmax(-(Fn3), 0.0);
				}
				else {
					integer inumber = iN3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fabs(Fn3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fabs(Fn3);
					}
				}

				if (!bS3) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fmax((Fs3), 0.0);
				}
				else {
					integer inumber = iS3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fabs(Fs3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fabs(Fs3);
					}
				}

				if (!bT3) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fmax(-(Ft3), 0.0);
				}
				else {
					integer inumber = iT3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fabs(Ft3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fabs(Ft3);
					}
				}

				if (!bB3) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fmax((Fb3), 0.0);
				}
				else {
					integer inumber = iB3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fabs(Fb3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fabs(Fb3);
					}
				}


				// ������ ������ ������� ����� ���������� �� UDS.
				if (!bE4) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fmax(-(Fe4), 0.0);
				}
				else {
					integer inumber = iE4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fabs(Fe4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fabs(Fe4);
					}
				}

				if (!bW4) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fmax((Fw4), 0.0);
				}
				else {
					integer inumber = iW4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fabs(Fw4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fabs(Fw4);
					}
				}

				if (!bN4) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fmax(-(Fn4), 0.0);
				}
				else {
					integer inumber = iN4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fabs(Fn4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fabs(Fn4);
					}
				}

				if (!bS4) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fmax((Fs4), 0.0);
				}
				else {
					integer inumber = iS4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fabs(Fs4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fabs(Fs4);
					}
				}

				if (!bT4) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fmax(-(Ft4), 0.0);
				}
				else {
					integer inumber = iT4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fabs(Ft4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fabs(Ft4);
					}
				}

				if (!bB4) {
					// ������ ����������.
					sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fmax((Fb4), 0.0);
				}
				else {
					integer inumber = iB4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fabs(Fb4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fabs(Fb4);
					}
				}



			}


			// 7.05.2017 ������� ��� ����������� ������ � ��������������� � ����������.
			sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = sl[GAMMA_LANGTRY_MENTER_SL][iP].ae + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw + sl[GAMMA_LANGTRY_MENTER_SL][iP].an + sl[GAMMA_LANGTRY_MENTER_SL][iP].as + sl[GAMMA_LANGTRY_MENTER_SL][iP].at + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;

			//sumanb = sl[GAMMA_LANGTRY_MENTER_SL][iP].ae + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw + sl[GAMMA_LANGTRY_MENTER_SL][iP].an + sl[GAMMA_LANGTRY_MENTER_SL][iP].as + sl[GAMMA_LANGTRY_MENTER_SL][iP].at + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;


			if (b_on_adaptive_local_refinement_mesh) {

				// 7.05.2017 ������� ��� ����������� ������ � ��������������� � ����������.
				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2;

				//sumanb += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2;

				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3;

				//sumanb += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3;

				sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4;

				//sumanb += sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 + sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4;


			}

		}


		// 7.05.2017 ������� ��� ����������� ������ � ��������������� � ����������.
		//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;

		//sumanb=sl[GAMMA_LANGTRY_MENTER_SL][iP].ae+sl[GAMMA_LANGTRY_MENTER_SL][iP].aw+sl[GAMMA_LANGTRY_MENTER_SL][iP].an+sl[GAMMA_LANGTRY_MENTER_SL][iP].as+sl[GAMMA_LANGTRY_MENTER_SL][iP].at+sl[GAMMA_LANGTRY_MENTER_SL][iP].ab;

		//13 ������� 2016.
		//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);
		//sumanb = fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ab);


		if (1) {
			// ����� � ������ ����� (����� ���������� ���������):
			// X - direction
			attrs += -fmax((Fe), 0) * (Speede - SpeedP) + fmax(-(Fe), 0) * (Speede - SpeedE);
			attrs += -fmax(-(Fw), 0) * (Speedw - SpeedP) + fmax((Fw), 0) * (Speedw - SpeedW);
			// Y - direction
			attrs += -fmax((Fn), 0) * (Speedn - SpeedP) + fmax(-(Fn), 0) * (Speedn - SpeedN);
			attrs += -fmax(-(Fs), 0) * (Speeds - SpeedP) + fmax((Fs), 0) * (Speeds - SpeedS);
			// Z - direction
			attrs += -fmax((Ft), 0) * (Speedt - SpeedP) + fmax(-(Ft), 0) * (Speedt - SpeedT);
			attrs += -fmax(-(Fb), 0) * (Speedb - SpeedP) + fmax((Fb), 0) * (Speedb - SpeedB);


			if (b_on_adaptive_local_refinement_mesh) {

				// ����� � ������ ����� (����� ���������� ���������):
				// X - direction
				attrs += -fmax((Fe2), 0) * (Speede2 - SpeedP) + fmax(-(Fe2), 0) * (Speede2 - SpeedE2);
				attrs += -fmax(-(Fw2), 0) * (Speedw2 - SpeedP) + fmax((Fw2), 0) * (Speedw2 - SpeedW2);
				// Y - direction
				attrs += -fmax((Fn2), 0) * (Speedn2 - SpeedP) + fmax(-(Fn2), 0) * (Speedn2 - SpeedN2);
				attrs += -fmax(-(Fs2), 0) * (Speeds2 - SpeedP) + fmax((Fs2), 0) * (Speeds2 - SpeedS2);
				// Z - direction
				attrs += -fmax((Ft2), 0) * (Speedt2 - SpeedP) + fmax(-(Ft2), 0) * (Speedt2 - SpeedT2);
				attrs += -fmax(-(Fb2), 0) * (Speedb2 - SpeedP) + fmax((Fb2), 0) * (Speedb2 - SpeedB2);


				// ����� � ������ ����� (����� ���������� ���������):
				// X - direction
				attrs += -fmax((Fe3), 0) * (Speede3 - SpeedP) + fmax(-(Fe3), 0) * (Speede3 - SpeedE3);
				attrs += -fmax(-(Fw3), 0) * (Speedw3 - SpeedP) + fmax((Fw3), 0) * (Speedw3 - SpeedW3);
				// Y - direction
				attrs += -fmax((Fn3), 0) * (Speedn3 - SpeedP) + fmax(-(Fn3), 0) * (Speedn3 - SpeedN3);
				attrs += -fmax(-(Fs3), 0) * (Speeds3 - SpeedP) + fmax((Fs3), 0) * (Speeds3 - SpeedS3);
				// Z - direction
				attrs += -fmax((Ft3), 0) * (Speedt3 - SpeedP) + fmax(-(Ft3), 0) * (Speedt3 - SpeedT3);
				attrs += -fmax(-(Fb3), 0) * (Speedb3 - SpeedP) + fmax((Fb3), 0) * (Speedb3 - SpeedB3);


				// ����� � ������ ����� (����� ���������� ���������):
				// X - direction
				attrs += -fmax((Fe4), 0) * (Speede4 - SpeedP) + fmax(-(Fe4), 0) * (Speede4 - SpeedE4);
				attrs += -fmax(-(Fw4), 0) * (Speedw4 - SpeedP) + fmax((Fw4), 0) * (Speedw4 - SpeedW4);
				// Y - direction
				attrs += -fmax((Fn4), 0) * (Speedn4 - SpeedP) + fmax(-(Fn4), 0) * (Speedn4 - SpeedN4);
				attrs += -fmax(-(Fs4), 0) * (Speeds4 - SpeedP) + fmax((Fs4), 0) * (Speeds4 - SpeedS4);
				// Z - direction
				attrs += -fmax((Ft4), 0) * (Speedt4 - SpeedP) + fmax(-(Ft4), 0) * (Speedt4 - SpeedT4);
				attrs += -fmax(-(Fb4), 0) * (Speedb4 - SpeedP) + fmax((Fb4), 0) * (Speedb4 - SpeedB4);

			}

		}
		else {
			// �������.
		   // 30 07 2015
		   // 
		   // ������ ������ ������� ����� ���������� �� UDS.
			if (!bE) {
				attrs += -fmax((Fe), 0) * (Speede - SpeedP) + fmax(-(Fe), 0) * (Speede - SpeedE);
			}
			if (!bW) {
				attrs += -fmax(-(Fw), 0) * (Speedw - SpeedP) + fmax((Fw), 0) * (Speedw - SpeedW);
			}
			if (!bN) {
				attrs += -fmax((Fn), 0) * (Speedn - SpeedP) + fmax(-(Fn), 0) * (Speedn - SpeedN);
			}
			if (!bS) {
				attrs += -fmax(-(Fs), 0) * (Speeds - SpeedP) + fmax((Fs), 0) * (Speeds - SpeedS);
			}
			if (!bT) {
				attrs += -fmax((Ft), 0) * (Speedt - SpeedP) + fmax(-(Ft), 0) * (Speedt - SpeedT);
			}
			if (!bB) {
				attrs += -fmax(-(Fb), 0) * (Speedb - SpeedP) + fmax((Fb), 0) * (Speedb - SpeedB);
			}


			if (b_on_adaptive_local_refinement_mesh) {

				if (!bE2) {
					attrs += -fmax((Fe2), 0) * (Speede2 - SpeedP) + fmax(-(Fe2), 0) * (Speede2 - SpeedE2);
				}
				if (!bW2) {
					attrs += -fmax(-(Fw2), 0) * (Speedw2 - SpeedP) + fmax((Fw2), 0) * (Speedw2 - SpeedW2);
				}
				if (!bN2) {
					attrs += -fmax((Fn2), 0) * (Speedn2 - SpeedP) + fmax(-(Fn2), 0) * (Speedn2 - SpeedN2);
				}
				if (!bS2) {
					attrs += -fmax(-(Fs2), 0) * (Speeds2 - SpeedP) + fmax((Fs2), 0) * (Speeds2 - SpeedS2);
				}
				if (!bT2) {
					attrs += -fmax((Ft2), 0) * (Speedt2 - SpeedP) + fmax(-(Ft2), 0) * (Speedt2 - SpeedT2);
				}
				if (!bB2) {
					attrs += -fmax(-(Fb2), 0) * (Speedb2 - SpeedP) + fmax((Fb2), 0) * (Speedb2 - SpeedB2);
				}


				if (!bE3) {
					attrs += -fmax((Fe3), 0) * (Speede3 - SpeedP) + fmax(-(Fe3), 0) * (Speede3 - SpeedE3);
				}
				if (!bW3) {
					attrs += -fmax(-(Fw3), 0) * (Speedw3 - SpeedP) + fmax((Fw3), 0) * (Speedw3 - SpeedW3);
				}
				if (!bN3) {
					attrs += -fmax((Fn3), 0) * (Speedn3 - SpeedP) + fmax(-(Fn3), 0) * (Speedn3 - SpeedN3);
				}
				if (!bS3) {
					attrs += -fmax(-(Fs3), 0) * (Speeds3 - SpeedP) + fmax((Fs3), 0) * (Speeds3 - SpeedS3);
				}
				if (!bT3) {
					attrs += -fmax((Ft3), 0) * (Speedt3 - SpeedP) + fmax(-(Ft3), 0) * (Speedt3 - SpeedT3);
				}
				if (!bB3) {
					attrs += -fmax(-(Fb3), 0) * (Speedb3 - SpeedP) + fmax((Fb3), 0) * (Speedb3 - SpeedB3);
				}


				if (!bE4) {
					attrs += -fmax((Fe4), 0) * (Speede4 - SpeedP) + fmax(-(Fe4), 0) * (Speede4 - SpeedE4);
				}
				if (!bW4) {
					attrs += -fmax(-(Fw4), 0) * (Speedw4 - SpeedP) + fmax((Fw4), 0) * (Speedw4 - SpeedW4);
				}
				if (!bN4) {
					attrs += -fmax((Fn4), 0) * (Speedn4 - SpeedP) + fmax(-(Fn4), 0) * (Speedn4 - SpeedN4);
				}
				if (!bS4) {
					attrs += -fmax(-(Fs4), 0) * (Speeds4 - SpeedP) + fmax((Fs4), 0) * (Speeds4 - SpeedS4);
				}
				if (!bT4) {
					attrs += -fmax((Ft4), 0) * (Speedt4 - SpeedP) + fmax(-(Ft4), 0) * (Speedt4 - SpeedT4);
				}
				if (!bB4) {
					attrs += -fmax(-(Fb4), 0) * (Speedb4 - SpeedP) + fmax((Fb4), 0) * (Speedb4 - SpeedB4);
				}

			}


		}

		//attrs=0.0; // ����� ����� ������� ����������� ����������� (�������� ����� ��������).

	}

	// ��������� ����������� ��������� ���������������.
	// ������ �������. �������������� ���������� �� ������������.
	// 

	// ������������ ����
	doublereal dSc = 0.0, dSp = 0.0;

	// ��� ������������� �������� Fe, Fw, Fn, Fs, Ft, Fb - ���������� �������� ������ ����� ���������� �������������� ���-���.
		// 02.05.2017
		// ��� ������� �.�. �������� � ������������� ������������ �������������.
		// ������ ���� �������: deltaF=(Fe-Fw+Fn-Fs+Ft-Fb);
		// ����������� ����� ����������� � ���������� �� ANSYS Icepak.
		// ��� ����������� �� ���� �������� ����� �� ���������� ���� - ��� ���������� ������ ��������.
		//doublereal deltaF=(Fe-Fw+Fn-Fs+Ft-Fb);
		// ��� ������ ���������� ��������� ������������� ��� ��������� ����� ����.
		// � ������ ���� ����������� ��������� ��� �������� �������� �������� � ������������ ������ (��)
		// ��� ���������� ��������� ��������� ������������ ������������, �� ������ ���������� ������ ���������� ����� 
		// ����� ����� ���������� ��������� ���������.
		// 8.05.2017.
	doublereal deltaF = fabs(Fe - Fw + Fn - Fs + Ft - Fb);
	if (b_on_adaptive_local_refinement_mesh) {
		deltaF = fabs(Fe - Fw + Fn - Fs + Ft - Fb + Fe2 - Fw2 + Fn2 - Fs2 + Ft2 - Fb2 + Fe3 - Fw3 + Fn3 - Fs3 + Ft3 - Fb3 + Fe4 - Fw4 + Fn4 - Fs4 + Ft4 - Fb4);
	}
	if (deltaF != deltaF) {
		printf("Fe=%e Fw=%e Fn=%e Fs=%e Ft=%e Fb=%e\n", Fe, Fw, Fn, Fs, Ft, Fb);
		printf("Fe2=%e Fw2=%e Fn2=%e Fs2=%e Ft2=%e Fb2=%e\n", Fe2, Fw2, Fn2, Fs2, Ft2, Fb2);
		printf("Fe3=%e Fw3=%e Fn3=%e Fs3=%e Ft3=%e Fb3=%e\n", Fe3, Fw3, Fn3, Fs3, Ft3, Fb3);
		printf("Fe4=%e Fw4=%e Fn4=%e Fs4=%e Ft4=%e Fb4=%e\n", Fe4, Fw4, Fn4, Fs4, Ft4, Fb4);
		printf("ERROR deltaF=%e\n", deltaF);
	}

	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ap != sl[GAMMA_LANGTRY_MENTER_SL][iP].ap) {
		printf("ap!=ap assemble bug. Apriory deltaF. iP=%lld ap=%e\n", iP, sl[GAMMA_LANGTRY_MENTER_SL][iP].ap);
		system("pause");
	}
	//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += deltaF;//-->//sl[GAMMA_LANGTRY_MENTER_SL][iP].ap+=apzero1+deltaF;//+deltaF; // ������������ ������� ������� deltaF ������ ��������������.  ���������� ��������� 
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ap != sl[GAMMA_LANGTRY_MENTER_SL][iP].ap) {
		printf("ap!=ap assemble bug. Apost deltaF. iP=%lld ap=%e\n", iP, sl[GAMMA_LANGTRY_MENTER_SL][iP].ap);
		system("pause");
	}


	sl[GAMMA_LANGTRY_MENTER_SL][iP].b = attrs; // ����� ���������� ��������� ��� ����� ������� ����������� �����������.
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].b != sl[GAMMA_LANGTRY_MENTER_SL][iP].b) {
		printf("exptsr+attrs error NAN or INF in control volume %lld GAMMA_LANGTRY_MENTER\n", iP);
		system("pause");
	}


	// ���������

	doublereal Re_theta_C = Re_theta_c(potent[RE_THETA_LANGTRY_MENTER][iP]);
	doublereal ReV = rP * SInvariantStrainRateTensor[iP] * distance_to_wall[iP] * distance_to_wall[iP] / prop[MU_DYNAMIC_VISCOSITY][iP];
	doublereal Fonset1 = ReV / (2.193 * Re_theta_C);
	doublereal Fonset2 = fmin(fmax(Fonset1, Fonset1* Fonset1* Fonset1* Fonset1),2.0);
	doublereal R_T = rP * fmax(K_limiter_min, potent[TURBULENT_KINETIK_ENERGY][iP]) / 
		(prop[MU_DYNAMIC_VISCOSITY][iP]* fmax(Omega_limiter_min, potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]));
	doublereal Fonset3 = fmax(1.0-(R_T/2.5)* (R_T / 2.5)* (R_T / 2.5), 0.0);
	doublereal Fonset = fmax(Fonset2 - Fonset3, 0.0);
	doublereal Fturb = exp(-(R_T/4.0)* (R_T / 4.0)* (R_T / 4.0)* (R_T / 4.0));

	doublereal Re_omega = rP* fmax(Omega_limiter_min, potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP])* distance_to_wall[iP] * distance_to_wall[iP]/ prop[MU_DYNAMIC_VISCOSITY][iP];
	doublereal Fsublayer = exp(-(Re_omega/200.0)* (Re_omega / 200.0));
	doublereal Flength = F_length_1(potent[RE_THETA_LANGTRY_MENTER][iP]) * (1.0 - Fsublayer) + 40.0 * Fsublayer;

	doublereal Pgamma = Flength* eqin.fluidinfo[0].Ca1* rP* SInvariantStrainRateTensor[iP]*pow(fmin(1.0,fmax(0.0, potent[GAMMA_LANGTRY_MENTER][iP]))*Fonset,0.5)*
		(1.0 - eqin.fluidinfo[0].Cepsilon1* fmin(1.0, fmax(0.0, potent[GAMMA_LANGTRY_MENTER][iP])));

	// ����������
	doublereal Egamma = eqin.fluidinfo[0].Ca2 * rP * potent[CURL][iP] * fmin(1.0, fmax(0.0, potent[GAMMA_LANGTRY_MENTER][iP])) *
		Fturb * (eqin.fluidinfo[0].Cepsilon2* fmin(1.0, fmax(0.0, potent[GAMMA_LANGTRY_MENTER][iP])) - 1.0);


	

	if (brthdsd_ON) {
		dSc = Pgamma; // ���������.
		// ����������
		dSc -= Egamma;
	}
	else {
		dSc = 0.0;
	}
	sl[GAMMA_LANGTRY_MENTER_SL][iP].b += dSc * dx * dy * dz; // ��������� ����� ����������.


	if (btimedep) {
		// ���� ��������������� �����.
		doublereal apzero = rP * dx * dy * dz / tauparam;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ap += apzero;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].b += apzero * gamma_Langtry_Menter_old_time_step[iP];
	}


	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].b != sl[GAMMA_LANGTRY_MENTER_SL][iP].b) {
		printf("dSc*dx*dy*dz error NAN or INF in control volume %lld\n", iP);
		printf("GAMMA_LANGTRY_MENTER\n");
		system("pause");
	}

	if (fabs(sl[GAMMA_LANGTRY_MENTER_SL][iP].ap) < 1.0e-30) {
		sl[GAMMA_LANGTRY_MENTER_SL][iP].ap = 1.0;
		sl[GAMMA_LANGTRY_MENTER_SL][iP].b = potent[GAMMA_LANGTRY_MENTER][iP];

		if (sl[GAMMA_LANGTRY_MENTER_SL][iP].b != sl[GAMMA_LANGTRY_MENTER_SL][iP].b) {
			printf("Zero ap in GAMMA_LANGTRY_MENTER component.\n");
		}
	}



	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ap != sl[GAMMA_LANGTRY_MENTER_SL][iP].ap) {
		printf("ap!=ap assemble bug. iP=%lld ap=%e\n", iP, sl[GAMMA_LANGTRY_MENTER_SL][iP].ap);
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ae != sl[GAMMA_LANGTRY_MENTER_SL][iP].ae) {
		printf("ae!=ae assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].aw != sl[GAMMA_LANGTRY_MENTER_SL][iP].aw) {
		printf("aw!=aw assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].an != sl[GAMMA_LANGTRY_MENTER_SL][iP].an) {
		printf("an!=an assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].as != sl[GAMMA_LANGTRY_MENTER_SL][iP].as) {
		printf("as!=as assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].at != sl[GAMMA_LANGTRY_MENTER_SL][iP].at) {
		printf("at!=at assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ab != sl[GAMMA_LANGTRY_MENTER_SL][iP].ab) {
		printf("ab!=ab assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2 != sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2) {
		printf("ae2!=ae2 assemble bug %e %e\n", sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2, sl[GAMMA_LANGTRY_MENTER_SL][iP].ae2);
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2 != sl[GAMMA_LANGTRY_MENTER_SL][iP].aw2) {
		printf("aw2!=aw2 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].an2 != sl[GAMMA_LANGTRY_MENTER_SL][iP].an2) {
		printf("an2!=an2 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].as2 != sl[GAMMA_LANGTRY_MENTER_SL][iP].as2) {
		printf("as2!=as2 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].at2 != sl[GAMMA_LANGTRY_MENTER_SL][iP].at2) {
		printf("at2!=at2 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2 != sl[GAMMA_LANGTRY_MENTER_SL][iP].ab2) {
		printf("ab2!=ab2 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3 != sl[GAMMA_LANGTRY_MENTER_SL][iP].ae3) {
		printf("ae3!=ae3 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3 != sl[GAMMA_LANGTRY_MENTER_SL][iP].aw3) {
		printf("aw3!=aw3 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].an3 != sl[GAMMA_LANGTRY_MENTER_SL][iP].an3) {
		printf("an3!=an3 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].as3 != sl[GAMMA_LANGTRY_MENTER_SL][iP].as3) {
		printf("as3!=as3 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].at3 != sl[GAMMA_LANGTRY_MENTER_SL][iP].at3) {
		printf("at3!=at3 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3 != sl[GAMMA_LANGTRY_MENTER_SL][iP].ab3) {
		printf("ab3!=ab3 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4 != sl[GAMMA_LANGTRY_MENTER_SL][iP].ae4) {
		printf("ae4!=ae4 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4 != sl[GAMMA_LANGTRY_MENTER_SL][iP].aw4) {
		printf("aw4!=aw4 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].an4 != sl[GAMMA_LANGTRY_MENTER_SL][iP].an4) {
		printf("an4!=an4 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].as4 != sl[GAMMA_LANGTRY_MENTER_SL][iP].as4) {
		printf("as4!=as4 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].at4 != sl[GAMMA_LANGTRY_MENTER_SL][iP].at4) {
		printf("at4!=at4 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4 != sl[GAMMA_LANGTRY_MENTER_SL][iP].ab4) {
		printf("ab4!=ab4 assemble bug\n");
		printf("GAMMA_LANGTRY_MENTER \n");
		system("pause");
	}

} // my_elmatr_quad_turbulent_Gamma_Langtry_Menter_SST_3D


// �������� ���� ��������� ������� ���� ��� ����������� ��������� 
// ��������� - ��������, ��� ������������ ����������� ������������ ������.
// ��� ������������� ��������� ������������ Hex �����.
// ��� ������ ����������� ������ ��� �������� Re_theta � ������ ������� ��������.
void my_elmatr_quad_turbulent_Re_Theta_Langtry_MenterSST_3D(
	integer iP,
	BOUND* border_neighbor,
	integer lw,
	integer ls,
	equation3D**& sl,
	equation3D_bon**& slb,
	//doublereal** diag_coef,
	//integer iVar,
	bool btimedep,
	doublereal tauparam,
	int* ptr,
	int** nvtx,
	doublereal** potent,
	//doublereal* potent_temper,
	TOCHKA* pa,
	float** prop,
	float** prop_b,
	integer maxelm,
	int*** neighbors_for_the_internal_node,
	//doublereal* alpha,
	//doublereal dgx,
	//doublereal dgy,
	//doublereal dgz,
	//doublereal dbeta, 
	integer ishconvection,
	doublereal*& Re_Theta_Menter_Lantgry_old_time_step,
	//bool bBussineskApproach,
	//doublereal temp_ref,
	//bool bfirst_start,
	//doublereal RCh, 
	//integer iflowregime,
	//doublereal* speedoldtimestep, 
	//doublereal* toldtimestep,
	BLOCK* b,
	integer lb,
	TPROP* matlist,
	doublereal**& mf,
	//bool bVERYStable,
	//doublereal &sumanb,
	integer* ilevel_alice,
	doublereal*& distance_to_wall,
	doublereal*& SInvariantStrainRateTensor,
	bool brthdsd_ON,
	TOCHKA*& center_coord_loc,
	TOCHKA*& volume_loc,
	integer maxbound
) {

	// iP - ����� ����������� ������������ ������
	// iP ���������� �� 0 �� maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // ������ �������� ����������� �������
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iE = iE; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iN = iN; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iT = iT;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iS = iS; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iW = iW; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iB = iB;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iP = iP;

	// 26.09.2016 ������� ��� ���� �����.
	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // ������ �������� ����������� �������
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // ������ �������� ����������� �������
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // ������ �������� ����������� �������

	if (b_on_adaptive_local_refinement_mesh) {
		// NON_EXISTENT_NODE ���� �� ������������ � [0..maxelm+maxbound-1] ���� ������������.
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}


	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iE2 = iE2; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iN2 = iN2; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iT2 = iT2;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iS2 = iS2; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iW2 = iW2; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iB2 = iB2;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iE3 = iE3; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iN3 = iN3; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iT3 = iT3;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iS3 = iS3; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iW3 = iW3; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iB3 = iB3;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iE4 = iE4; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iN4 = iN4; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iT4 = iT4;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].iS4 = iS4; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iW4 = iW4; sl[RE_THETA_LANGTRY_MENTER_SL][iP].iB4 = iB4;



	// ���������������� ���������.
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = 0.0;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = 0.0;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = 0.0;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = 0.0;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = 0.0;


	// ������� ����������� �����.
	// �� ������� ������ ����� ���������� � ����� �������� ������ ���.
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].bE2 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bW2 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bS2 = false;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].bN2 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bB2 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bT2 = false;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].bE3 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bW3 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bS3 = false;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].bN3 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bB3 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bT3 = false;

	sl[RE_THETA_LANGTRY_MENTER_SL][iP].bE4 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bW4 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bS4 = false;
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].bN4 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bB4 = false; sl[RE_THETA_LANGTRY_MENTER_SL][iP].bT4 = false;

	if (CHECK_NODE_FOR_EXISTENCE(iE2)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bE2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iW2)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bW2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iN2)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bN2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iS2)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bS2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iT2)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bT2 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iB2)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bB2 = true;

	if (CHECK_NODE_FOR_EXISTENCE(iE3)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bE3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iW3)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bW3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iN3)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bN3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iS3)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bS3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iT3)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bT3 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iB3)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bB3 = true;

	if (CHECK_NODE_FOR_EXISTENCE(iE4)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bE4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iW4)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bW4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iN4)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bN4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iS4)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bS4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iT4)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bT4 = true;
	if (CHECK_NODE_FOR_EXISTENCE(iB4)) sl[RE_THETA_LANGTRY_MENTER_SL][iP].bB4 = true;

	// ���������� ��.	

	// ���� � ����� �� ������ ����� ������� ��������� �������
	// �� ��������������� ���������� ����� true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if ((iE >= maxelm) && (iE < maxelm + maxbound)) bE = true;
	if ((iN >= maxelm) && (iN < maxelm + maxbound)) bN = true;
	if ((iT >= maxelm) && (iT < maxelm + maxbound)) bT = true;
	if ((iW >= maxelm) && (iW < maxelm + maxbound)) bW = true;
	if ((iS >= maxelm) && (iS < maxelm + maxbound)) bS = true;
	if ((iB >= maxelm) && (iB < maxelm + maxbound)) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if ((iE2 >= maxelm) && (iE2 < maxelm + maxbound)) bE2 = true;
	if ((iN2 >= maxelm) && (iN2 < maxelm + maxbound)) bN2 = true;
	if ((iT2 >= maxelm) && (iT2 < maxelm + maxbound)) bT2 = true;
	if ((iW2 >= maxelm) && (iW2 < maxelm + maxbound)) bW2 = true;
	if ((iS2 >= maxelm) && (iS2 < maxelm + maxbound)) bS2 = true;
	if ((iB2 >= maxelm) && (iB2 < maxelm + maxbound)) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if ((iE3 >= maxelm) && (iE3 < maxelm + maxbound)) bE3 = true;
	if ((iN3 >= maxelm) && (iN3 < maxelm + maxbound)) bN3 = true;
	if ((iT3 >= maxelm) && (iT3 < maxelm + maxbound)) bT3 = true;
	if ((iW3 >= maxelm) && (iW3 < maxelm + maxbound)) bW3 = true;
	if ((iS3 >= maxelm) && (iS3 < maxelm + maxbound)) bS3 = true;
	if ((iB3 >= maxelm) && (iB3 < maxelm + maxbound)) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if ((iE4 >= maxelm) && (iE4 < maxelm + maxbound)) bE4 = true;
	if ((iN4 >= maxelm) && (iN4 < maxelm + maxbound)) bN4 = true;
	if ((iT4 >= maxelm) && (iT4 < maxelm + maxbound)) bT4 = true;
	if ((iW4 >= maxelm) && (iW4 < maxelm + maxbound)) bW4 = true;
	if ((iS4 >= maxelm) && (iS4 < maxelm + maxbound)) bS4 = true;
	if ((iB4 >= maxelm) && (iB4 < maxelm + maxbound)) bB4 = true;

	// ���������� �������� �������� ������������ ������:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
	//volume3D(iP, nvtx, pa, dx, dy, dz);

	TOCHKA pvol = volume_loc[iP];
	dx = pvol.x;
	dy = pvol.y;
	dz = pvol.z;

	//printf("%.2f %.2f\n",dx,dy); // debug GOOD
	//getchar();

	doublereal dxe = 0.5 * dx, dxw = 0.5 * dx, dyn = 0.5 * dy, dys = 0.5 * dy, dzt = 0.5 * dz, dzb = 0.5 * dz;
	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE)) {
		if (!bE) dxe = 0.5 * (pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW)) {
		if (!bW) dxw = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5 * (pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN)) {
		if (!bN) dyn = 0.5 * (pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS)) {
		if (!bS) dys = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5 * (pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT)) {
		if (!bT) dzt = 0.5 * (pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB)) {
		if (!bB) dzb = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5 * (pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5 * dx, dxw2 = 0.5 * dx, dyn2 = 0.5 * dy, dys2 = 0.5 * dy, dzt2 = 0.5 * dz, dzb2 = 0.5 * dz;
	doublereal dxe3 = 0.5 * dx, dxw3 = 0.5 * dx, dyn3 = 0.5 * dy, dys3 = 0.5 * dy, dzt3 = 0.5 * dz, dzb3 = 0.5 * dz;
	doublereal dxe4 = 0.5 * dx, dxw4 = 0.5 * dx, dyn4 = 0.5 * dy, dys4 = 0.5 * dy, dzt4 = 0.5 * dz, dzb4 = 0.5 * dz;

	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE2)) {
		if (!bE2) dxe2 = 0.5 * (pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW2)) {
		if (!bW2) dxw2 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5 * (pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN2)) {
		if (!bN2) dyn2 = 0.5 * (pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS2)) {
		if (!bS2) dys2 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5 * (pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT2)) {
		if (!bT2) dzt2 = 0.5 * (pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB2)) {
		if (!bB2) dzb2 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5 * (pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE3)) {
		if (!bE3) dxe3 = 0.5 * (pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW3)) {
		if (!bW3) dxw3 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5 * (pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN3)) {
		if (!bN3) dyn3 = 0.5 * (pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS3)) {
		if (!bS3) dys3 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5 * (pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT3)) {
		if (!bT3) dzt3 = 0.5 * (pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB3)) {
		if (!bB3) dzb3 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5 * (pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// �.�. �������� ��������� ������ ����, �� ����� ��� ������������
	// x - direction
	if (CHECK_NODE_FOR_EXISTENCE(iE4)) {
		if (!bE4) dxe4 = 0.5 * (pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iW4)) {
		if (!bW4) dxw4 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5 * (pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (CHECK_NODE_FOR_EXISTENCE(iN4)) {
		if (!bN4) dyn4 = 0.5 * (pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iS4)) {
		if (!bS4) dys4 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5 * (pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (CHECK_NODE_FOR_EXISTENCE(iT4)) {
		if (!bT4) dzt4 = 0.5 * (pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (CHECK_NODE_FOR_EXISTENCE(iB4)) {
		if (!bB4) dzb4 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5 * (pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}


	// ���� ��������������� ��������� �����:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5 * dx / dxe;
	fwplus = 0.5 * dx / dxw;
	// y-direction
	fnplus = 0.5 * dy / dyn;
	fsplus = 0.5 * dy / dys;
	// z-direction
	ftplus = 0.5 * dz / dzt;
	fbplus = 0.5 * dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5 * dx / dxe2;
	fwplus2 = 0.5 * dx / dxw2;
	// y-direction
	fnplus2 = 0.5 * dy / dyn2;
	fsplus2 = 0.5 * dy / dys2;
	// z-direction
	ftplus2 = 0.5 * dz / dzt2;
	fbplus2 = 0.5 * dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5 * dx / dxe3;
	fwplus3 = 0.5 * dx / dxw3;
	// y-direction
	fnplus3 = 0.5 * dy / dyn3;
	fsplus3 = 0.5 * dy / dys3;
	// z-direction
	ftplus3 = 0.5 * dz / dzt3;
	fbplus3 = 0.5 * dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5 * dx / dxe4;
	fwplus4 = 0.5 * dx / dxw4;
	// y-direction
	fnplus4 = 0.5 * dy / dyn4;
	fsplus4 = 0.5 * dy / dys4;
	// z-direction
	ftplus4 = 0.5 * dz / dzt4;
	fbplus4 = 0.5 * dz / dzb4;

	//printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	//getchar();

	// ��������� �� ����� �� ���������������� ������� �������������	
	doublereal rP, rE = 0.0, rN = 0.0, rT = 0.0, rW = 0.0, rS = 0.0, rB = 0.0;
	rP = prop[RHO][iP];
	if (iE > -1) {
		if (!bE) rE = prop[RHO][iE]; else rE = prop_b[RHO][iE - maxelm];
	}
	if (iN > -1) {
		if (!bN) rN = prop[RHO][iN]; else rN = prop_b[RHO][iN - maxelm];
	}
	if (iT > -1) {
		if (!bT) rT = prop[RHO][iT]; else rT = prop_b[RHO][iT - maxelm];
	}
	if (iW > -1) {
		if (!bW) rW = prop[RHO][iW]; else rW = prop_b[RHO][iW - maxelm];
	}
	if (iS > -1) {
		if (!bS) rS = prop[RHO][iS]; else rS = prop_b[RHO][iS - maxelm];
	}
	if (iB > -1) {
		if (!bB) rB = prop[RHO][iB]; else rB = prop_b[RHO][iB - maxelm];
	}

	doublereal  rE2 = 0.0, rN2 = 0.0, rT2 = 0.0, rW2 = 0.0, rS2 = 0.0, rB2 = 0.0;
	if (iE2 > -1) {
		if (!bE2) rE2 = prop[RHO][iE2]; else rE2 = prop_b[RHO][iE2 - maxelm];
	}
	if (iN2 > -1) {
		if (!bN2) rN2 = prop[RHO][iN2]; else rN2 = prop_b[RHO][iN2 - maxelm];
	}
	if (iT2 > -1) {
		if (!bT2) rT2 = prop[RHO][iT2]; else rT2 = prop_b[RHO][iT2 - maxelm];
	}
	if (iW2 > -1) {
		if (!bW2) rW2 = prop[RHO][iW2]; else rW2 = prop_b[RHO][iW2 - maxelm];
	}
	if (iS2 > -1) {
		if (!bS2) rS2 = prop[RHO][iS2]; else rS2 = prop_b[RHO][iS2 - maxelm];
	}
	if (iB2 > -1) {
		if (!bB2) rB2 = prop[RHO][iB2]; else rB2 = prop_b[RHO][iB2 - maxelm];
	}

	doublereal  rE3 = 0.0, rN3 = 0.0, rT3 = 0.0, rW3 = 0.0, rS3 = 0.0, rB3 = 0.0;
	if (iE3 > -1) {
		if (!bE3) rE3 = prop[RHO][iE3]; else rE3 = prop_b[RHO][iE3 - maxelm];
	}
	if (iN3 > -1) {
		if (!bN3) rN3 = prop[RHO][iN3]; else rN3 = prop_b[RHO][iN3 - maxelm];
	}
	if (iT3 > -1) {
		if (!bT3) rT3 = prop[RHO][iT3]; else rT3 = prop_b[RHO][iT3 - maxelm];
	}
	if (iW3 > -1) {
		if (!bW3) rW3 = prop[RHO][iW3]; else rW3 = prop_b[RHO][iW3 - maxelm];
	}
	if (iS3 > -1) {
		if (!bS3) rS3 = prop[RHO][iS3]; else rS3 = prop_b[RHO][iS3 - maxelm];
	}
	if (iB3 > -1) {
		if (!bB3) rB3 = prop[RHO][iB3]; else rB3 = prop_b[RHO][iB3 - maxelm];
	}

	doublereal  rE4 = 0.0, rN4 = 0.0, rT4 = 0.0, rW4 = 0.0, rS4 = 0.0, rB4 = 0.0;
	if (iE4 > -1) {
		if (!bE4) rE4 = prop[RHO][iE4]; else rE4 = prop_b[RHO][iE4 - maxelm];
	}
	if (iN4 > -1) {
		if (!bN4) rN4 = prop[RHO][iN4]; else rN4 = prop_b[RHO][iN4 - maxelm];
	}
	if (iT4 > -1) {
		if (!bT4) rT4 = prop[RHO][iT4]; else rT4 = prop_b[RHO][iT4 - maxelm];
	}
	if (iW4 > -1) {
		if (!bW4) rW4 = prop[RHO][iW4]; else rW4 = prop_b[RHO][iW4 - maxelm];
	}
	if (iS4 > -1) {
		if (!bS4) rS4 = prop[RHO][iS4]; else rS4 = prop_b[RHO][iS4 - maxelm];
	}
	if (iB4 > -1) {
		if (!bB4) rB4 = prop[RHO][iB4]; else rB4 = prop_b[RHO][iB4 - maxelm];
	}

	/*
	doublereal rhoe, rhow, rhon, rhos, rhot, rhob;
	// �������� ���������  �� ����� ��:
	rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP);  // ���������.
	rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP);
	rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP);
	rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP);
	rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP);
	rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP);


	doublereal rhoe2, rhow2, rhon2, rhos2, rhot2, rhob2;
	doublereal rhoe3, rhow3, rhon3, rhos3, rhot3, rhob3;
	doublereal rhoe4, rhow4, rhon4, rhos4, rhot4, rhob4;

	rhoe2 = rE2 * rP / (feplus2*rE2 + (1.0 - feplus2)*rP);
	rhow2 = rW2 * rP / (fwplus2*rW2 + (1.0 - fwplus2)*rP);
	rhon2 = rN2 * rP / (fnplus2*rN2 + (1.0 - fnplus2)*rP);
	rhos2 = rS2 * rP / (fsplus2*rS2 + (1.0 - fsplus2)*rP);
	rhot2 = rT2 * rP / (ftplus2*rT2 + (1.0 - ftplus2)*rP);
	rhob2 = rB2 * rP / (fbplus2*rB2 + (1.0 - fbplus2)*rP);

	rhoe3 = rE3 * rP / (feplus3*rE3 + (1.0 - feplus3)*rP);
	rhow3 = rW3 * rP / (fwplus3*rW3 + (1.0 - fwplus3)*rP);
	rhon3 = rN3 * rP / (fnplus3*rN3 + (1.0 - fnplus3)*rP);
	rhos3 = rS3 * rP / (fsplus3*rS3 + (1.0 - fsplus3)*rP);
	rhot3 = rT3 * rP / (ftplus3*rT3 + (1.0 - ftplus3)*rP);
	rhob3 = rB3 * rP / (fbplus3*rB3 + (1.0 - fbplus3)*rP);

	rhoe4 = rE4 * rP / (feplus4*rE4 + (1.0 - feplus4)*rP);
	rhow4 = rW4 * rP / (fwplus4*rW4 + (1.0 - fwplus4)*rP);
	rhon4 = rN4 * rP / (fnplus4*rN4 + (1.0 - fnplus4)*rP);
	rhos4 = rS4 * rP / (fsplus4*rS4 + (1.0 - fsplus4)*rP);
	rhot4 = rT4 * rP / (ftplus4*rT4 + (1.0 - ftplus4)*rP);
	rhob4 = rB4 * rP / (fbplus4*rB4 + (1.0 - fbplus4)*rP);
	*/
	doublereal rhoe = 0.0, rhow = 0.0, rhon = 0.0, rhos = 0.0, rhot = 0.0, rhob = 0.0;
	// ������������ ��������� ������� ���, ����� ����������� 
	// ���������� �����������.
	if (iE > -1) {
		if (!bE) rhoe = rE * rP / (feplus * rE + (1.0 - feplus) * rP); else rhoe = rE; // ��������� !
	}
	if (iW > -1) {
		if (!bW) rhow = rW * rP / (fwplus * rW + (1.0 - fwplus) * rP); else rhow = rW;
	}
	if (iN > -1) {
		if (!bN) rhon = rN * rP / (fnplus * rN + (1.0 - fnplus) * rP); else rhon = rN;
	}
	if (iS > -1) {
		if (!bS) rhos = rS * rP / (fsplus * rS + (1.0 - fsplus) * rP); else rhos = rS;
	}
	if (iT > -1) {
		if (!bT) rhot = rT * rP / (ftplus * rT + (1.0 - ftplus) * rP); else rhot = rT;
	}
	if (iB > -1) {
		if (!bB) rhob = rB * rP / (fbplus * rB + (1.0 - fbplus) * rP); else rhob = rB;
	}

	doublereal rhoe2 = 0.0, rhow2 = 0.0, rhon2 = 0.0, rhos2 = 0.0, rhot2 = 0.0, rhob2 = 0.0;
	doublereal rhoe3 = 0.0, rhow3 = 0.0, rhon3 = 0.0, rhos3 = 0.0, rhot3 = 0.0, rhob3 = 0.0;
	doublereal rhoe4 = 0.0, rhow4 = 0.0, rhon4 = 0.0, rhos4 = 0.0, rhot4 = 0.0, rhob4 = 0.0;

	if (iE2 > -1) {
		if (!bE2)  rhoe2 = rE2 * rP / (feplus2 * rE2 + (1.0 - feplus2) * rP); else rhoe2 = rE2; // ��������� !
	}
	if (iW2 > -1) {
		if (!bW2)  rhow2 = rW2 * rP / (fwplus2 * rW2 + (1.0 - fwplus2) * rP); else rhow2 = rW2;
	}
	if (iN2 > -1) {
		if (!bN2) rhon2 = rN2 * rP / (fnplus2 * rN2 + (1.0 - fnplus2) * rP); else rhon2 = rN2;
	}
	if (iS2 > -1) {
		if (!bS2)  rhos2 = rS2 * rP / (fsplus2 * rS2 + (1.0 - fsplus2) * rP); else rhos2 = rS2;
	}
	if (iT2 > -1) {
		if (!bT2)  rhot2 = rT2 * rP / (ftplus2 * rT2 + (1.0 - ftplus2) * rP); else rhot2 = rT2;
	}
	if (iB2 > -1) {
		if (!bB2) rhob2 = rB2 * rP / (fbplus2 * rB2 + (1.0 - fbplus2) * rP); else rhob2 = rB2;
	}

	if (iE3 > -1) {
		if (!bE3) rhoe3 = rE3 * rP / (feplus3 * rE3 + (1.0 - feplus3) * rP); else rhoe3 = rE3;
	}
	if (iW3 > -1) {
		if (!bW3) rhow3 = rW3 * rP / (fwplus3 * rW3 + (1.0 - fwplus3) * rP); else rhow3 = rW3;
	}
	if (iN3 > -1) {
		if (!bN3) rhon3 = rN3 * rP / (fnplus3 * rN3 + (1.0 - fnplus3) * rP); else rhon3 = rN3;
	}
	if (iS3 > -1) {
		if (!bS3) rhos3 = rS3 * rP / (fsplus3 * rS3 + (1.0 - fsplus3) * rP); else rhos3 = rS3;
	}
	if (iT3 > -1) {
		if (!bT3) rhot3 = rT3 * rP / (ftplus3 * rT3 + (1.0 - ftplus3) * rP); else rhot3 = rT3;
	}
	if (iB3 > -1) {
		if (!bB3) rhob3 = rB3 * rP / (fbplus3 * rB3 + (1.0 - fbplus3) * rP); else rhob3 = rB3;
	}

	if (iE4 > -1) {
		if (!bE4) rhoe4 = rE4 * rP / (feplus4 * rE4 + (1.0 - feplus4) * rP); else rhoe4 = rE4;
	}
	if (iW4 > -1) {
		if (!bW4) rhow4 = rW4 * rP / (fwplus4 * rW4 + (1.0 - fwplus4) * rP); else rhow4 = rW4;
	}
	if (iN4 > -1) {
		if (!bN4) rhon4 = rN4 * rP / (fnplus4 * rN4 + (1.0 - fnplus4) * rP); else rhon4 = rN4;
	}
	if (iS4 > -1) {
		if (!bS4) rhos4 = rS4 * rP / (fsplus4 * rS4 + (1.0 - fsplus4) * rP); else rhos4 = rS4;
	}
	if (iT4 > -1) {
		if (!bT4) rhot4 = rT4 * rP / (ftplus4 * rT4 + (1.0 - ftplus4) * rP); else rhot4 = rT4;
	}
	if (iB4 > -1) {
		if (!bB4) rhob4 = rB4 * rP / (fbplus4 * rB4 + (1.0 - fbplus4) * rP); else rhob4 = rB4;
	}



	/*
	   ����������� ����������:
	   �� �������� ������������ ���-��� �������� �� ����� ��
	   ������ ���� ��������� � � ���������� ��� ��������� ��������.
	   ��� ������������� ������ ���� ��������� ��� ���������� ������� �� ����� ������������ ������
	   � ���������� ���������� ��������, ���������������� � ������������ ������������� �.�. � ����������
	   ����������� ��� ������ ����� ������������� ��������� �������������. ���� ��� ���������� ������� ��
	   ��������� ������������ ���-��� �� ���� ��������� ���������� � �� ��������� (���� ��� � ������������
	   ������ � ��������� ���������� �������� � ����������������. ��� �������� �������� ��������� �����������
	   ������� ����� ��������� ��������� ����������) �� ������ ����� �� ����� ������������� ��������� �������������
	   � ������������� ������� ����� ��������.
	   ������������ ���������� ������������ �������� �� ��� ��� ������������, �
	   �� ����������� ������ ���. ��� ����������� ����� ��������� ������������� �������� ���� ��� �� ������
	   ����������������� �������� � ��������. ��� ���������� ��������� ������������ ������������ �
	   ���������� ��� ��������� ��������. ��� ������� � ������� �������� ���������
	   SIMPLE (�� ������ ����������������� ���������� ������� ������������ ������� � ������� ��������, ��
	   ���� � ��� ��� ����� �� �� ��������� �������� ���������� ����������� ����� ������ ����� (������� ���� ��������� � ������)
	   � ������� ������� ��� ������������ ������������ �������� � ���������� ��������).
	   ��������� ������������� ��������...
	   �������� ������ �������������� ������ ��������, �.�. �� ��� ������������ ������������
	   ��� ���� ����� ��� �� ���������. ������� ������������ �������� ������������ ���-��� ������
	   �� ������ �������� ��������� SIMPLE. �� ������ �������� ������������� ������� ������������
	   �������� ����� ���������� ������� �������� ������������� ��������.
	*/


	// ������������ ����� ����� ����� ��.
	// � ���������� �������� � ������ ���������������� ��������������� 
	// ���������� ���������������� �������� ���-���.
	doublereal Fe = 0.0, Fw = 0.0, Fn = 0.0, Fs = 0.0, Ft = 0.0, Fb = 0.0;


	// ��� ���� �����.
	doublereal Fe1 = 0.0, Fe2 = 0.0, Fe3 = 0.0, Fe4 = 0.0;
	doublereal Fw1 = 0.0, Fw2 = 0.0, Fw3 = 0.0, Fw4 = 0.0;
	doublereal Fn1 = 0.0, Fn2 = 0.0, Fn3 = 0.0, Fn4 = 0.0;
	doublereal Fs1 = 0.0, Fs2 = 0.0, Fs3 = 0.0, Fs4 = 0.0;
	doublereal Ft1 = 0.0, Ft2 = 0.0, Ft3 = 0.0, Ft4 = 0.0;
	doublereal Fb1 = 0.0, Fb2 = 0.0, Fb3 = 0.0, Fb4 = 0.0;

	// �������� ����� ������������ � �� �����������, 
	// ��� ������ ����������� ��������� �� �������� ������ ������� ����.
	if (!b_on_adaptive_local_refinement_mesh) {
		if (Fe != Fe) {
			printf("Fe=%e\n", Fe);
		}
		Fe = mf[iP][E_SIDE];
		if (Fe != Fe) {
			printf("Fe=%e\n", Fe);
			system("pause");
		}
		Fn = mf[iP][N_SIDE];
		Ft = mf[iP][T_SIDE];
		Fw = mf[iP][W_SIDE];
		Fs = mf[iP][S_SIDE];
		Fb = mf[iP][B_SIDE];
	}
	else {
		// ����� �� ����. 24.11.2018

		if (iE > -1) {
			if (bE) {
				// ��������� ����.
				Fe = mf[iP][E_SIDE] * (border_neighbor[iE - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
					Fe = mf[iP][E_SIDE];
				}
				else {

					Fe = mf[iE][W_SIDE];

				}
			}
		}

		if (iW > -1) {
			if (bW) {
				// ��������� ����.
				Fw = mf[iP][W_SIDE] * (border_neighbor[iW - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
					Fw = mf[iP][W_SIDE];
				}
				else {

					Fw = mf[iW][E_SIDE];

				}
			}
		}

		if (iN > -1) {
			if (bN) {
				// ��������� ����.
				Fn = mf[iP][N_SIDE] * (border_neighbor[iN - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
					Fn = mf[iP][N_SIDE];
				}
				else {

					Fn = mf[iN][S_SIDE];

				}
			}
		}

		if (iS > -1) {
			if (bS) {
				// ��������� ����.
				Fs = mf[iP][S_SIDE] * (border_neighbor[iS - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
					Fs = mf[iP][S_SIDE];
				}
				else {

					Fs = mf[iS][N_SIDE];

				}
			}
		}

		if (iT > -1) {
			if (bT) {
				// ��������� ����.
				Ft = mf[iP][T_SIDE] * (border_neighbor[iT - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
					Ft = mf[iP][T_SIDE];
				}
				else {

					Ft = mf[iT][B_SIDE];

				}
			}
		}

		if (iB > -1) {
			if (bB) {
				// ��������� ����.
				Fb = mf[iP][B_SIDE] * (border_neighbor[iB - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
					Fb = mf[iP][B_SIDE];
				}
				else {

					Fb = mf[iB][T_SIDE];

				}
			}
		}

		if (iE2 > -1) {
			if (bE2) {
				// ��������� ����.
				Fe2 = mf[iP][E_SIDE] * (border_neighbor[iE2 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
					Fe2 = mf[iP][E_SIDE];
				}
				else {

					Fe2 = mf[iE2][W_SIDE];

				}
			}
		}

		if (iW2 > -1) {
			if (bW2) {
				// ��������� ����.
				Fw2 = mf[iP][W_SIDE] * (border_neighbor[iW2 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW2]]) {
					Fw = mf[iP][W_SIDE];
				}
				else {

					Fw = mf[iW2][E_SIDE];

				}
			}
		}

		if (iN2 > -1) {
			if (bN2) {
				// ��������� ����.
				Fn2 = mf[iP][N_SIDE] * (border_neighbor[iN2 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
					Fn2 = mf[iP][N_SIDE];
				}
				else {

					Fn2 = mf[iN2][S_SIDE];

				}
			}
		}

		if (iS2 > -1) {
			if (bS2) {
				// ��������� ����.
				Fs2 = mf[iP][S_SIDE] * (border_neighbor[iS2 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
					Fs2 = mf[iP][S_SIDE];
				}
				else {

					Fs2 = mf[iS2][N_SIDE];

				}
			}
		}

		if (iT2 > -1) {
			if (bT2) {
				// ��������� ����.
				Ft2 = mf[iP][T_SIDE] * (border_neighbor[iT2 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
					Ft2 = mf[iP][T_SIDE];
				}
				else {

					Ft2 = mf[iT2][B_SIDE];

				}
			}
		}

		if (iB2 > -1) {
			if (bB2) {
				// ��������� ����.
				Fb2 = mf[iP][B_SIDE] * (border_neighbor[iB2 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
					Fb2 = mf[iP][B_SIDE];
				}
				else {

					Fb2 = mf[iB2][T_SIDE];

				}
			}
		}


		if (iE3 > -1) {
			if (bE3) {
				// ��������� ����.
				Fe3 = mf[iP][E_SIDE] * (border_neighbor[iE3 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
					Fe3 = mf[iP][E_SIDE];
				}
				else {

					Fe3 = mf[iE3][W_SIDE];

				}
			}
		}

		if (iW3 > -1) {
			if (bW3) {
				// ��������� ����.
				Fw3 = mf[iP][W_SIDE] * (border_neighbor[iW3 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
					Fw3 = mf[iP][W_SIDE];
				}
				else {

					Fw3 = mf[iW3][E_SIDE];

				}
			}
		}

		if (iN3 > -1) {
			if (bN3) {
				// ��������� ����.
				Fn3 = mf[iP][N_SIDE] * (border_neighbor[iN3 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
					Fn3 = mf[iP][N_SIDE];
				}
				else {

					Fn3 = mf[iN3][S_SIDE];

				}
			}
		}

		if (iS3 > -1) {
			if (bS3) {
				// ��������� ����.
				Fs3 = mf[iP][S_SIDE] * (border_neighbor[iS3 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
					Fs3 = mf[iP][S_SIDE];
				}
				else {

					Fs3 = mf[iS3][N_SIDE];

				}
			}
		}

		if (iT3 > -1) {
			if (bT3) {
				// ��������� ����.
				Ft3 = mf[iP][T_SIDE] * (border_neighbor[iT3 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
					Ft3 = mf[iP][T_SIDE];
				}
				else {

					Ft3 = mf[iT3][B_SIDE];

				}
			}
		}

		if (iB3 > -1) {
			if (bB3) {
				// ��������� ����.
				Fb3 = mf[iP][B_SIDE] * (border_neighbor[iB3 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
					Fb3 = mf[iP][B_SIDE];
				}
				else {

					Fb3 = mf[iB3][T_SIDE];

				}
			}
		}

		if (iE4 > -1) {
			if (bE4) {
				// ��������� ����.
				Fe4 = mf[iP][E_SIDE] * (border_neighbor[iE4 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
					Fe4 = mf[iP][E_SIDE];
				}
				else {

					Fe4 = mf[iE4][W_SIDE];

				}
			}
		}

		if (iW4 > -1) {
			if (bW4) {
				// ��������� ����.
				Fw4 = mf[iP][W_SIDE] * (border_neighbor[iW4 - maxelm].dS / (dy * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
					Fw4 = mf[iP][W_SIDE];
				}
				else {

					Fw4 = mf[iW4][E_SIDE];

				}
			}
		}

		if (iN4 > -1) {
			if (bN4) {
				// ��������� ����.
				Fn4 = mf[iP][N_SIDE] * (border_neighbor[iN4 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
					Fn4 = mf[iP][N_SIDE];
				}
				else {

					Fn4 = mf[iN4][S_SIDE];

				}
			}
		}

		if (iS4 > -1) {
			if (bS4) {
				// ��������� ����.
				Fs4 = mf[iP][S_SIDE] * (border_neighbor[iS4 - maxelm].dS / (dx * dz));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
					Fs4 = mf[iP][S_SIDE];
				}
				else {

					Fs4 = mf[iS4][N_SIDE];

				}
			}
		}

		if (iT4 > -1) {
			if (bT4) {
				// ��������� ����.
				Ft4 = mf[iP][T_SIDE] * (border_neighbor[iT4 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
					Ft4 = mf[iP][T_SIDE];
				}
				else {

					Ft4 = mf[iT4][B_SIDE];

				}
			}
		}

		if (iB4 > -1) {
			if (bB4) {
				// ��������� ����.
				Fb4 = mf[iP][B_SIDE] * (border_neighbor[iB4 - maxelm].dS / (dx * dy));
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
					Fb4 = mf[iP][B_SIDE];
				}
				else {

					Fb4 = mf[iB4][T_SIDE];

				}
			}
		}


	}



	doublereal eps = 1e-37; // ��� ��������� ������������� ����.

	//eqin.fluidinfo[0].sigma_nu
	// ������������ ��������:
	doublereal  GP, GE, GW, GN, GS, GT, GB;
	doublereal  GE2, GW2, GN2, GS2, GT2, GB2;
	doublereal  GE3, GW3, GN3, GS3, GT3, GB3;
	doublereal  GE4, GW4, GN4, GS4, GT4, GB4;

	
	doublereal sigma_theta_T = eqin.fluidinfo[0].sigma_theta_T;

	doublereal sigma_theta_TE = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TW = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TN = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TS = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TT = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TB = eqin.fluidinfo[0].sigma_theta_T;

	doublereal sigma_theta_TE2 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TW2 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TN2 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TS2 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TT2 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TB2 = eqin.fluidinfo[0].sigma_theta_T;

	doublereal sigma_theta_TE3 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TW3 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TN3 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TS3 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TT3 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TB3 = eqin.fluidinfo[0].sigma_theta_T;

	doublereal sigma_theta_TE4 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TW4 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TN4 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TS4 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TT4 = eqin.fluidinfo[0].sigma_theta_T;
	doublereal sigma_theta_TB4 = eqin.fluidinfo[0].sigma_theta_T;

	

	// ���������� ������������ ��������:
	GP = ((sigma_theta_T*prop[MU_DYNAMIC_VISCOSITY][iP]) + fmax(0.0, sigma_theta_T * potent[MUT][iP])); // � ������ ����������� ��.
	if (iE > -1) {
		if (!bE) GE = ((sigma_theta_TE*prop[MU_DYNAMIC_VISCOSITY][iE]) + fmax(0.0, sigma_theta_TE * potent[MUT][iE])); else GE = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iE - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iE]));
	}
	if (iN > -1) {
		if (!bN) GN = ((sigma_theta_TN*prop[MU_DYNAMIC_VISCOSITY][iN]) + fmax(0.0, sigma_theta_TN * potent[MUT][iN])); else GN = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iN - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iN]));
	}
	if (iT > -1) {
		if (!bT) GT = ((sigma_theta_TT*prop[MU_DYNAMIC_VISCOSITY][iT]) + fmax(0.0, sigma_theta_TT * potent[MUT][iT])); else GT = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iT - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iT]));
	}
	if (iW > -1) {
		if (!bW) GW = ((sigma_theta_TW*prop[MU_DYNAMIC_VISCOSITY][iW]) + fmax(0.0, sigma_theta_TW * potent[MUT][iW])); else GW = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iW - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iW]));
	}
	if (iS > -1) {
		if (!bS) GS = ((sigma_theta_TS*prop[MU_DYNAMIC_VISCOSITY][iS]) + fmax(0.0, sigma_theta_TS * potent[MUT][iS])); else GS = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iS - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iS]));
	}
	if (iB > -1) {
		if (!bB) GB = ((sigma_theta_TB*prop[MU_DYNAMIC_VISCOSITY][iB]) + fmax(0.0, sigma_theta_TB * potent[MUT][iB])); else GB = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iB - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iB]));
	}

	if (iE2 > -1) {
		if (!bE2) GE2 = ((sigma_theta_TE2 * prop[MU_DYNAMIC_VISCOSITY][iE2]) + fmax(0.0, sigma_theta_TE2 * potent[MUT][iE2])); else GE2 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iE2 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iE2]));
	}
	if (iN2 > -1) {
		if (!bN2) GN2 = ((sigma_theta_TN2 * prop[MU_DYNAMIC_VISCOSITY][iN2]) + fmax(0.0, sigma_theta_TN2 * potent[MUT][iN2])); else GN2 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iN2 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iN2]));
	}
	if (iT2 > -1) {
		if (!bT2) GT2 = ((sigma_theta_TT2 * prop[MU_DYNAMIC_VISCOSITY][iT2]) + fmax(0.0, sigma_theta_TT2 * potent[MUT][iT2])); else GT2 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iT2 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iT2]));
	}
	if (iW2 > -1) {
		if (!bW2) GW2 = ((sigma_theta_TW2 * prop[MU_DYNAMIC_VISCOSITY][iW2]) + fmax(0.0, sigma_theta_TW2 * potent[MUT][iW2])); else GW2 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iW2 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iW2]));
	}
	if (iS2 > -1) {
		if (!bS2) GS2 = ((sigma_theta_TS2 * prop[MU_DYNAMIC_VISCOSITY][iS2]) + fmax(0.0, sigma_theta_TS2 * potent[MUT][iS2])); else GS2 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iS2 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iS2]));
	}
	if (iB2 > -1) {
		if (!bB2) GB2 = ((sigma_theta_TB2 * prop[MU_DYNAMIC_VISCOSITY][iB2]) + fmax(0.0, sigma_theta_TB2 * potent[MUT][iB2])); else GB2 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iB2 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iB2]));
	}

	if (iE3 > -1) {
		if (!bE3) GE3 = ((sigma_theta_TE3 * prop[MU_DYNAMIC_VISCOSITY][iE3]) + fmax(0.0, sigma_theta_TE3 * potent[MUT][iE3])); else GE3 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iE3 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iE3]));
	}
	if (iN3 > -1) {
		if (!bN3) GN3 = ((sigma_theta_TN3 * prop[MU_DYNAMIC_VISCOSITY][iN3]) + fmax(0.0, sigma_theta_TN3 * potent[MUT][iN3])); else GN3 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iN3 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iN3]));
	}
	if (iT3 > -1) {
		if (!bT3) GT3 = ((sigma_theta_TT3 * prop[MU_DYNAMIC_VISCOSITY][iT3]) + fmax(0.0, sigma_theta_TT3 * potent[MUT][iT3])); else GT3 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iT3 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iT3]));
	}
	if (iW3 > -1) {
		if (!bW3) GW3 = ((sigma_theta_TW3 * prop[MU_DYNAMIC_VISCOSITY][iW3]) + fmax(0.0, sigma_theta_TW3 * potent[MUT][iW3])); else GW3 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iW3 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iW3]));
	}
	if (iS3 > -1) {
		if (!bS3) GS3 = ((sigma_theta_TS3 * prop[MU_DYNAMIC_VISCOSITY][iS3]) + fmax(0.0, sigma_theta_TS3 * potent[MUT][iS3])); else GS3 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iS3 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iS3]));
	}
	if (iB3 > -1) {
		if (!bB3) GB3 = ((sigma_theta_TB3 * prop[MU_DYNAMIC_VISCOSITY][iB3]) + fmax(0.0, sigma_theta_TB3 * potent[MUT][iB3])); else GB3 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iB3 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iB3]));
	}

	if (iE4 > -1) {
		if (!bE4) GE4 = ((sigma_theta_TE4 * prop[MU_DYNAMIC_VISCOSITY][iE4]) + fmax(0.0, sigma_theta_TE4 * potent[MUT][iE4])); else GE4 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iE4 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iE4]));
	}
	if (iN4 > -1) {
		if (!bN4) GN4 = ((sigma_theta_TN4 * prop[MU_DYNAMIC_VISCOSITY][iN4]) + fmax(0.0, sigma_theta_TN4 * potent[MUT][iN4])); else GN4 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iN4 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iN4]));
	}
	if (iT4 > -1) {
		if (!bT4) GT4 = ((sigma_theta_TT4 * prop[MU_DYNAMIC_VISCOSITY][iT4]) + fmax(0.0, sigma_theta_TT4 * potent[MUT][iT4])); else GT4 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iT4 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iT4]));
	}
	if (iW4 > -1) {
		if (!bW4) GW4 = ((sigma_theta_TW4 * prop[MU_DYNAMIC_VISCOSITY][iW4]) + fmax(0.0, sigma_theta_TW4 * potent[MUT][iW4])); else GW4 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iW4 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iW4]));
	}
	if (iS4 > -1) {
		if (!bS4) GS4 = ((sigma_theta_TS4 * prop[MU_DYNAMIC_VISCOSITY][iS4]) + fmax(0.0, sigma_theta_TS4 * potent[MUT][iS4])); else GS4 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iS4 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iS4]));
	}
	if (iB4 > -1) {
		if (!bB4) GB4 = ((sigma_theta_TB4 * prop[MU_DYNAMIC_VISCOSITY][iB4]) + fmax(0.0, sigma_theta_TB4 * potent[MUT][iB4])); else GB4 = ((sigma_theta_T * prop_b[MU_DYNAMIC_VISCOSITY][iB4 - maxelm]) + fmax(0.0, sigma_theta_T * potent[MUT][iB4]));
	}

	doublereal Ge = GP, Gw = GP, Gn = GP, Gs = GP, Gt = GP, Gb = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE > -1) {
		Ge = GE * GP / (feplus * GE + (1 - feplus) * GP); // ���������.
	}
	if (iW > -1) {
		Gw = GW * GP / (fwplus * GW + (1 - fwplus) * GP);
	}
	if (iN > -1) {
		Gn = GN * GP / (fnplus * GN + (1 - fnplus) * GP);
	}
	if (iS > -1) {
		Gs = GS * GP / (fsplus * GS + (1 - fsplus) * GP);
	}
	if (iT > -1) {
		Gt = GT * GP / (ftplus * GT + (1 - ftplus) * GP);
	}
	if (iB > -1) {
		Gb = GB * GP / (fbplus * GB + (1 - fbplus) * GP);
	}

	doublereal Ge2 = GP, Gw2 = GP, Gn2 = GP, Gs2 = GP, Gt2 = GP, Gb2 = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE2 > -1) {
		Ge2 = GE2 * GP / (feplus2 * GE2 + (1 - feplus2) * GP); // ���������.
	}
	if (iW2 > -1) {
		Gw2 = GW2 * GP / (fwplus2 * GW2 + (1 - fwplus2) * GP);
	}
	if (iN2 > -1) {
		Gn2 = GN2 * GP / (fnplus2 * GN2 + (1 - fnplus2) * GP);
	}
	if (iS2 > -1) {
		Gs2 = GS2 * GP / (fsplus2 * GS2 + (1 - fsplus2) * GP);
	}
	if (iT2 > -1) {
		Gt2 = GT2 * GP / (ftplus2 * GT2 + (1 - ftplus2) * GP);
	}
	if (iB2 > -1) {
		Gb2 = GB2 * GP / (fbplus2 * GB2 + (1 - fbplus2) * GP);
	}


	doublereal Ge3 = GP, Gw3 = GP, Gn3 = GP, Gs3 = GP, Gt3 = GP, Gb3 = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE3 > -1) {
		Ge3 = GE3 * GP / (feplus3 * GE3 + (1 - feplus3) * GP); // ���������.
	}
	if (iW3 > -1) {
		Gw3 = GW3 * GP / (fwplus3 * GW3 + (1 - fwplus3) * GP);
	}
	if (iN3 > -1) {
		Gn3 = GN3 * GP / (fnplus3 * GN3 + (1 - fnplus3) * GP);
	}
	if (iS3 > -1) {
		Gs3 = GS3 * GP / (fsplus3 * GS3 + (1 - fsplus3) * GP);
	}
	if (iT3 > -1) {
		Gt3 = GT3 * GP / (ftplus3 * GT3 + (1 - ftplus3) * GP);
	}
	if (iB3 > -1) {
		Gb3 = GB3 * GP / (fbplus3 * GB3 + (1 - fbplus3) * GP);
	}

	doublereal Ge4 = GP, Gw4 = GP, Gn4 = GP, Gs4 = GP, Gt4 = GP, Gb4 = GP;
	// �������� ������������ �������� �� ����� ��.
	if (iE4 > -1) {
		Ge4 = GE4 * GP / (feplus4 * GE4 + (1 - feplus4) * GP); // ���������.
	}
	if (iW4 > -1) {
		Gw4 = GW4 * GP / (fwplus4 * GW4 + (1 - fwplus4) * GP);
	}
	if (iN4 > -1) {
		Gn4 = GN4 * GP / (fnplus4 * GN4 + (1 - fnplus4) * GP);
	}
	if (iS4 > -1) {
		Gs4 = GS4 * GP / (fsplus4 * GS4 + (1 - fsplus4) * GP);
	}
	if (iT4 > -1) {
		Gt4 = GT4 * GP / (ftplus4 * GT4 + (1 - ftplus4) * GP);
	}
	if (iB4 > -1) {
		Gb4 = GB4 * GP / (fbplus4 * GB4 + (1 - fbplus4) * GP);
	}




	const doublereal ZeroDiffusion = 0.0;// 1.0e-30;
	// ������������ ������������ ������:
	doublereal De = ZeroDiffusion, Dw = ZeroDiffusion, Dn = ZeroDiffusion, Ds = ZeroDiffusion, Dt = ZeroDiffusion, Db = ZeroDiffusion; // �������������
	if (iE > -1) {
		if (bE) {
			// ��������� ����.
			De = Ge * border_neighbor[iE - maxelm].dS / dxe;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
				De = Ge * dy * dz / dxe;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				De = Ge * dy_loc * dz_loc / dxe;
			}
		}

	}
	if (iW > -1) {
		if (bW) {
			// ��������� ����
			Dw = Gw * border_neighbor[iW - maxelm].dS / dxw;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				Dw = Gw * dy * dz / dxw;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw = Gw * dy_loc * dz_loc / dxw;
			}
		}

	}
	if (iN > -1) {
		if (bN) {
			// ��������� ����.
			Dn = Gn * border_neighbor[iN - maxelm].dS / dyn;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
				Dn = Gn * dx * dz / dyn;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn = Gn * dx_loc * dz_loc / dyn;
			}
		}
	}
	if (iS > -1) {
		if (bS) {
			// ��������� ����
			Ds = Gs * border_neighbor[iS - maxelm].dS / dys;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
				Ds = Gs * dx * dz / dys;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds = Gs * dx_loc * dz_loc / dys;
			}
		}
	}
	if (iT > -1) {
		if (bT) {
			// ��������� ����.
			Dt = Gt * border_neighbor[iT - maxelm].dS / dzt;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
				Dt = Gt * dx * dy / dzt;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt = Gt * dx_loc * dy_loc / dzt;
			}
		}


	}
	if (iB > -1) {
		if (bB) {
			// ��������� ����
			Db = Gb * border_neighbor[iB - maxelm].dS / dzb;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
				Db = Gb * dx * dy / dzb;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;


				Db = Gb * dx_loc * dy_loc / dzb;
			}
		}


	}

	// ������������ ������������ ������:
	doublereal De2 = ZeroDiffusion, Dw2 = ZeroDiffusion, Dn2 = ZeroDiffusion, Ds2 = ZeroDiffusion, Dt2 = ZeroDiffusion, Db2 = ZeroDiffusion; // �������������
	if (iE2 > -1) {
		if (bE2) {
			// ��������� ����.
			De2 = Ge2 * border_neighbor[iE2 - maxelm].dS / dxe2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
				De2 = Ge2 * dy * dz / dxe2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;


				De2 = Ge2 * dy_loc * dz_loc / dxe2;
			}
		}

	}
	if (iW2 > -1) {
		if (bW2) {
			// ��������� ����
			Dw2 = Gw2 * border_neighbor[iW2 - maxelm].dS / dxw2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW2]]) {
				Dw2 = Gw2 * dy * dz / dxw2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw2 = Gw2 * dy_loc * dz_loc / dxw2;
			}
		}

	}
	if (iN2 > -1) {
		if (bN2) {
			// ��������� ����.
			Dn2 = Gn2 * border_neighbor[iN2 - maxelm].dS / dyn2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
				Dn2 = Gn2 * dx * dz / dyn2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn2 = Gn2 * dx_loc * dz_loc / dyn2;
			}
		}
	}
	if (iS2 > -1) {
		if (bS2) {
			// ��������� ����
			Ds2 = Gs2 * border_neighbor[iS2 - maxelm].dS / dys2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
				Ds2 = Gs2 * dx * dz / dys2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds2 = Gs2 * dx_loc * dz_loc / dys2;
			}
		}
	}
	if (iT2 > -1) {
		if (bT2) {
			// ��������� ����.
			Dt2 = Gt2 * border_neighbor[iT2 - maxelm].dS / dzt2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
				Dt2 = Gt2 * dx * dy / dzt2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt2 = Gt2 * dx_loc * dy_loc / dzt2;
			}
		}


	}
	if (iB2 > -1) {
		if (bB2) {
			// ��������� ����
			Db2 = Gb2 * border_neighbor[iB2 - maxelm].dS / dzb2;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
				Db2 = Gb2 * dx * dy / dzb2;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB2];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Db2 = Gb2 * dx_loc * dy_loc / dzb2;
			}
		}


	}

	// ������������ ������������ ������:
	doublereal De3 = ZeroDiffusion, Dw3 = ZeroDiffusion, Dn3 = ZeroDiffusion, Ds3 = ZeroDiffusion, Dt3 = ZeroDiffusion, Db3 = ZeroDiffusion; // �������������
	if (iE3 > -1) {
		if (bE3) {
			// ��������� ����.
			De3 = Ge3 * border_neighbor[iE3 - maxelm].dS / dxe3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
				De3 = Ge3 * dy * dz / dxe3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				De3 = Ge3 * dy_loc * dz_loc / dxe3;
			}
		}

	}
	if (iW3 > -1) {
		if (bW3) {
			// ��������� ����
			Dw3 = Gw3 * border_neighbor[iW3 - maxelm].dS / dxw3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
				Dw3 = Gw3 * dy * dz / dxw3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw3 = Gw3 * dy_loc * dz_loc / dxw3;
			}
		}

	}
	if (iN3 > -1) {
		if (bN3) {
			// ��������� ����.
			Dn3 = Gn3 * border_neighbor[iN3 - maxelm].dS / dyn3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
				Dn3 = Gn3 * dx * dz / dyn3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn3 = Gn3 * dx_loc * dz_loc / dyn3;
			}
		}
	}
	if (iS3 > -1) {
		if (bS3) {
			// ��������� ����
			Ds3 = Gs3 * border_neighbor[iS3 - maxelm].dS / dys3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
				Ds3 = Gs3 * dx * dz / dys3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds3 = Gs3 * dx_loc * dz_loc / dys3;
			}
		}
	}
	if (iT3 > -1) {
		if (bT3) {
			// ��������� ����.
			Dt3 = Gt3 * border_neighbor[iT3 - maxelm].dS / dzt3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
				Dt3 = Gt3 * dx * dy / dzt3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt3 = Gt3 * dx_loc * dy_loc / dzt3;
			}
		}


	}
	if (iB3 > -1) {
		if (bB3) {
			// ��������� ����
			Db3 = Gb3 * border_neighbor[iB3 - maxelm].dS / dzb3;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
				Db3 = Gb3 * dx * dy / dzb3;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB3];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Db3 = Gb3 * dx_loc * dy_loc / dzb3;
			}
		}


	}

	// ������������ ������������ ������:
	doublereal De4 = ZeroDiffusion, Dw4 = ZeroDiffusion, Dn4 = ZeroDiffusion, Ds4 = ZeroDiffusion, Dt4 = ZeroDiffusion, Db4 = ZeroDiffusion; // �������������
	if (iE4 > -1) {
		if (bE4) {
			// ��������� ����.
			De4 = Ge4 * border_neighbor[iE4 - maxelm].dS / dxe4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
				De4 = Ge4 * dy * dz / dxe4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iE4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;


				De4 = Ge4 * dy_loc * dz_loc / dxe4;
			}
		}

	}
	if (iW4 > -1) {
		if (bW4) {
			// ��������� ����
			Dw4 = Gw4 * border_neighbor[iW4 - maxelm].dS / dxw4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
				Dw4 = Gw4 * dy * dz / dxw4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iW4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dw4 = Gw4 * dy_loc * dz_loc / dxw4;
			}
		}

	}
	if (iN4 > -1) {
		if (bN4) {
			// ��������� ����.
			Dn4 = Gn4 * border_neighbor[iN4 - maxelm].dS / dyn4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
				Dn4 = Gn4 * dx * dz / dyn4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iN4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dn4 = Gn4 * dx_loc * dz_loc / dyn4;
			}
		}
	}
	if (iS4 > -1) {
		if (bS4) {
			// ��������� ����
			Ds4 = Gs4 * border_neighbor[iS4 - maxelm].dS / dys4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
				Ds4 = Gs4 * dx * dz / dys4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iS4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Ds4 = Gs4 * dx_loc * dz_loc / dys4;
			}
		}
	}
	if (iT4 > -1) {
		if (bT4) {
			// ��������� ����.
			Dt4 = Gt4 * border_neighbor[iT4 - maxelm].dS / dzt4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
				Dt4 = Gt4 * dx * dy / dzt4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iT4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Dt4 = Gt4 * dx_loc * dy_loc / dzt4;
			}
		}


	}
	if (iB4 > -1) {
		if (bB4) {
			// ��������� ����
			Db4 = Gb4 * border_neighbor[iB4 - maxelm].dS / dzb4;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
				Db4 = Gb4 * dx * dy / dzb4;
			}
			else {
				// ���������� �������� ��������� ������������ ������:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// ����� �������� ������������ ������
				//volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				pvol = volume_loc[iB4];
				dx_loc = pvol.x;
				dy_loc = pvol.y;
				dz_loc = pvol.z;

				Db4 = Gb4 * dx_loc * dy_loc / dzb4;
			}
		}
	}

	// ������������� ��������� ������� ���������.
	De = fmax(De, 0.0);
	Dw = fmax(Dw, 0.0);
	Dn = fmax(Dn, 0.0);
	Ds = fmax(Ds, 0.0);
	Dt = fmax(Dt, 0.0);
	Db = fmax(Db, 0.0);

	De2 = fmax(De2, 0.0);
	Dw2 = fmax(Dw2, 0.0);
	Dn2 = fmax(Dn2, 0.0);
	Ds2 = fmax(Ds2, 0.0);
	Dt2 = fmax(Dt2, 0.0);
	Db2 = fmax(Db2, 0.0);

	De3 = fmax(De3, 0.0);
	Dw3 = fmax(Dw3, 0.0);
	Dn3 = fmax(Dn3, 0.0);
	Ds3 = fmax(Ds3, 0.0);
	Dt3 = fmax(Dt3, 0.0);
	Db3 = fmax(Db3, 0.0);

	De4 = fmax(De4, 0.0);
	Dw4 = fmax(Dw4, 0.0);
	Dn4 = fmax(Dn4, 0.0);
	Ds4 = fmax(Ds4, 0.0);
	Dt4 = fmax(Dt4, 0.0);
	Db4 = fmax(Db4, 0.0);

	// ����� �����:
	doublereal Pe = 0.0, Pw = 0.0, Pn = 0.0, Ps = 0.0, Pt = 0.0, Pb = 0.0;
	if (iE > -1) {
		Pe = (Fe) / De;
	}
	if (iW > -1) {
		Pw = -(Fw) / Dw;
	}
	if (iN > -1) {
		Pn = (Fn) / Dn;
	}
	if (iS > -1) {
		Ps = -(Fs) / Ds;
	}
	if (iT > -1) {
		Pt = (Ft) / Dt;
	}
	if (iB > -1) {
		Pb = -(Fb) / Db;
	}

	// ����� �����:
	doublereal Pe2 = 0.0, Pw2 = 0.0, Pn2 = 0.0, Ps2 = 0.0, Pt2 = 0.0, Pb2 = 0.0;
	if (iE2 > -1) {
		Pe2 = (Fe2) / De2;
	}
	if (iW2 > -1) {
		Pw2 = -(Fw2) / Dw2;
	}
	if (iN2 > -1) {
		Pn2 = (Fn2) / Dn2;
	}
	if (iS2 > -1) {
		Ps2 = -(Fs2) / Ds2;
	}
	if (iT2 > -1) {
		Pt2 = (Ft2) / Dt2;
	}
	if (iB2 > -1) {
		Pb2 = -(Fb2) / Db2;
	}

	// ����� �����:
	doublereal Pe3 = 0.0, Pw3 = 0.0, Pn3 = 0.0, Ps3 = 0.0, Pt3 = 0.0, Pb3 = 0.0;
	if (iE3 > -1) {
		Pe3 = (Fe3) / De3;
	}
	if (iW3 > -1) {
		Pw3 = -(Fw3) / Dw3;
	}
	if (iN3 > -1) {
		Pn3 = (Fn3) / Dn3;
	}
	if (iS3 > -1) {
		Ps3 = -(Fs3) / Ds3;
	}
	if (iT3 > -1) {
		Pt3 = (Ft3) / Dt3;
	}
	if (iB3 > -1) {
		Pb3 = -(Fb3) / Db3;
	}

	// ����� �����:
	doublereal Pe4 = 0.0, Pw4 = 0.0, Pn4 = 0.0, Ps4 = 0.0, Pt4 = 0.0, Pb4 = 0.0;
	if (iE4 > -1) {
		Pe4 = (Fe4) / De4;
	}
	if (iW4 > -1) {
		Pw4 = -(Fw4) / Dw4;
	}
	if (iN4 > -1) {
		Pn4 = (Fn4) / Dn4;
	}
	if (iS4 > -1) {
		Ps4 = -(Fs4) / Ds4;
	}
	if (iT4 > -1) {
		Pt4 = (Ft4) / Dt4;
	}
	if (iB4 > -1) {
		Pb4 = -(Fb4) / Db4;
	}


	// ������� � ������ ����� ��� ������������� ����� �������� QUICK
	// � ���� ������������� ������ ���������� ���������.
	// addition to the right side QUICK Leonard.
	doublereal attrs = 0.0;

	if (b_on_adaptive_local_refinement_mesh) {

		// ���������������� ���������.
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = 0.0;

		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = 0.0;

		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = 0.0;

		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = 0.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = 0.0;

	}

	if (ishconvection < distsheme) {

		if (1) {
			// ������� ��� ����������� ������ � ��������������� � ���������� 7.05.2017. 
			if (b_on_adaptive_local_refinement_mesh) {
				// ���������� ������������� ����������� �������:
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);

				// ���������� ������������� ����������� �������:
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(-(Fe2), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(Fw2, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(-(Fn2), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(Fs2, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(-(Ft2), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(Fb2, 0);

				// ���������� ������������� ����������� �������:
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(-(Fe3), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(Fw3, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(-(Fn3), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(Fs3, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(-(Ft3), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(Fb3, 0);

				// ���������� ������������� ����������� �������:
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(-(Fe4), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(Fw4, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(-(Fn4), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(Fs4, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(-(Ft4), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(Fb4, 0);

			}
			else {
				// 25 07 2015
				// ���������� ������������� ����������� �������:
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);
				//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
			}
		}
		else
		{
			// �������� �� ������ ����������������������� 25 ���� 2015.
			if (!bE) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
			}
			else {
				integer inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fabs(Fe);
				}
				else {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-(Fe), 0);
				}
			}
			if (!bW) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
			}
			else {
				integer inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fabs(Fw);
				}
				else {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
				}
			}
			if (!bN) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
			}
			else {
				integer inumber = iN - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fabs(Fn);
				}
				else {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-(Fn), 0);
				}
			}
			if (!bS) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
			}
			else {
				integer inumber = iS - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fabs(Fs);
				}
				else {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
				}
			}
			if (!bT) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
			}
			else {
				integer inumber = iT - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fabs(Ft);
				}
				else {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(-(Ft), 0);
				}
			}
			if (!bB) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);
			}
			else
			{
				integer inumber = iB - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fabs(Fb);
				}
				else {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(Fb, 0);
				}
			}


			if (b_on_adaptive_local_refinement_mesh) {

				if (!bE2) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(-(Fe2), 0);
				}
				else {
					integer inumber = iE2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fabs(Fe2);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(-(Fe2), 0);
					}
				}
				if (!bW2) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(Fw2, 0);
				}
				else {
					integer inumber = iW2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fabs(Fw2);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(Fw2, 0);
					}
				}
				if (!bN2) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(-(Fn2), 0);
				}
				else {
					integer inumber = iN2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fabs(Fn2);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(-(Fn2), 0);
					}
				}
				if (!bS2) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(Fs2, 0);
				}
				else {
					integer inumber = iS2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fabs(Fs2);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(Fs2, 0);
					}
				}
				if (!bT2) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(-(Ft2), 0);
				}
				else {
					integer inumber = iT2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fabs(Ft2);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(-(Ft2), 0);
					}
				}
				if (!bB2) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(Fb2, 0);
				}
				else
				{
					integer inumber = iB2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fabs(Fb2);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(Fb2, 0);
					}
				}

				if (!bE3) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(-(Fe3), 0);
				}
				else {
					integer inumber = iE3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fabs(Fe3);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(-(Fe3), 0);
					}
				}
				if (!bW3) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(Fw3, 0);
				}
				else {
					integer inumber = iW3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fabs(Fw3);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(Fw3, 0);
					}
				}
				if (!bN3) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(-(Fn3), 0);
				}
				else {
					integer inumber = iN3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fabs(Fn3);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(-(Fn3), 0);
					}
				}
				if (!bS3) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(Fs3, 0);
				}
				else {
					integer inumber = iS3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fabs(Fs3);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(Fs3, 0);
					}
				}
				if (!bT3) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(-(Ft3), 0);
				}
				else {
					integer inumber = iT3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fabs(Ft3);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(-(Ft3), 0);
					}
				}
				if (!bB3) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(Fb3, 0);
				}
				else
				{
					integer inumber = iB3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fabs(Fb3);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(Fb3, 0);
					}
				}

				if (!bE4) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(-(Fe4), 0);
				}
				else {
					integer inumber = iE4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fabs(Fe4);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(-(Fe4), 0);
					}
				}
				if (!bW4) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(Fw4, 0);
				}
				else {
					integer inumber = iW4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fabs(Fw4);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(Fw4, 0);
					}
				}
				if (!bN4) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(-(Fn4), 0);
				}
				else {
					integer inumber = iN4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fabs(Fn4);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(-(Fn4), 0);
					}
				}
				if (!bS4) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(Fs4, 0);
				}
				else {
					integer inumber = iS4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fabs(Fs4);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(Fs4, 0);
					}
				}
				if (!bT4) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(-(Ft4), 0);
				}
				else {
					integer inumber = iT4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fabs(Ft4);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(-(Ft4), 0);
					}
				}
				if (!bB4) {
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(Fb4, 0);
				}
				else
				{
					integer inumber = iB4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fabs(Fb4);
					}
					else {
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(Fb4, 0);
					}
				}
			}


		}

		// ������ ��� ����������� ������ � ��������� � ����������. 7.05.2017.
		//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
		// ��� ���������:
		// ����� ������������� !!! (���������� �� ����).
		if (b_on_adaptive_local_refinement_mesh) {
			// ���� 

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(+(Fe), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-(Fw), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(+(Fn), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-(Fs), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(+(Ft), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(-(Fb), 0);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(+(Fe2), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(-(Fw2), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(+(Fn2), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(-(Fs2), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(+(Ft2), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(-(Fb2), 0);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(+(Fe3), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(-(Fw3), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(+(Fn3), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(-(Fs3), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(+(Ft3), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(-(Fb3), 0);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(+(Fe4), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(-(Fw4), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(+(Fn4), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(-(Fs4), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(+(Ft4), 0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(-(Fb4), 0);
			/*
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae +  sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
			if (b_on_adaptive_local_refinement_mesh) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2;
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3;
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4;
			}
			*/
		}
		else {
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(+(Fe), 0);
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-(Fw), 0);
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(+(Fn), 0);
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-(Fs), 0);
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(+(Ft), 0);
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(-(Fb), 0);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae +
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an +
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at +
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;

		}

		// 13 ������� 2016
		// ��� ��������. ��� ����� �� �������� � ����������. �� ������� ��� ��������� ������������ ������������, ��
		// ������������� ���������� ���� � �������, �� �� ��������������� ��� ����� ���������� ������������ ����������� ������ ������.
		//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);




		if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap < 1.0e-36) {
			printf("Zero diagonal coefficient in internal volume in my_elmatr_quad_turbulent_Re_Theta_Langtry_MenterSST_3D.\n");
#if doubleintprecision == 1
			printf("ap=%e iP=%lld\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap, iP);
#else
			printf("ap=%e iP=%d\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap, iP);
#endif
			if (b_on_adaptive_local_refinement_mesh) {
				printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae, sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw, sl[RE_THETA_LANGTRY_MENTER_SL][iP].an, sl[RE_THETA_LANGTRY_MENTER_SL][iP].as, sl[RE_THETA_LANGTRY_MENTER_SL][iP].at, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);
				printf("ae2=%e aw2=%e an2=%e as2=%e at2=%e ab2=%e\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2, sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2, sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2, sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2, sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2);
				printf("ae3=%e aw3=%e an3=%e as3=%e at3=%e ab3=%e\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3, sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3, sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3, sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3, sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3);
				printf("ae4=%e aw4=%e an4=%e as4=%e at4=%e ab4=%e\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4, sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4, sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4, sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4, sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4);
			}
			else {
				printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae, sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw, sl[RE_THETA_LANGTRY_MENTER_SL][iP].an, sl[RE_THETA_LANGTRY_MENTER_SL][iP].as, sl[RE_THETA_LANGTRY_MENTER_SL][iP].at, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);
			}
			system("pause");

			// ��� ������ ������ �� ������� �� ��������
			// �.�. �� ������ ����� ����� ������� ������������ �������.
			// ������� ����� �������� �������� ����� ������ ����� ��������� �����������, ��� ����� �����������.
			// �� �������� �� ��-���� ��������� �������� ����� ���������.
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De * ApproxConvective(fabs(Pe), BULG);//+fmax(-(Fe),0); 
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw * ApproxConvective(fabs(Pw), BULG);//+fmax(Fw,0); 
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn * ApproxConvective(fabs(Pn), BULG);//+fmax(-(Fn),0); 
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds * ApproxConvective(fabs(Ps), BULG);//+fmax(Fs,0); 
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt * ApproxConvective(fabs(Pt), BULG);//+fmax(-(Ft),0); 
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db * ApproxConvective(fabs(Pb), BULG);//+fmax(Fb,0); 
			if (b_on_adaptive_local_refinement_mesh) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 * ApproxConvective(fabs(Pe2), BULG);//+fmax(-(Fe2),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 * ApproxConvective(fabs(Pw2), BULG);//+fmax(Fw2,0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 * ApproxConvective(fabs(Pn2), BULG);//+fmax(-(Fn2),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 * ApproxConvective(fabs(Ps2), BULG);//+fmax(Fs2,0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 * ApproxConvective(fabs(Pt2), BULG);//+fmax(-(Ft2),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 * ApproxConvective(fabs(Pb2), BULG);//+fmax(Fb2,0); 

				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 * ApproxConvective(fabs(Pe3), BULG);//+fmax(-(Fe3),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 * ApproxConvective(fabs(Pw3), BULG);//+fmax(Fw3,0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 * ApproxConvective(fabs(Pn3), BULG);//+fmax(-(Fn3),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 * ApproxConvective(fabs(Ps3), BULG);//+fmax(Fs3,0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 * ApproxConvective(fabs(Pt3), BULG);//+fmax(-(Ft3),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 * ApproxConvective(fabs(Pb3), BULG);//+fmax(Fb3,0); 

				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 * ApproxConvective(fabs(Pe4), BULG);//+fmax(-(Fe4),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 * ApproxConvective(fabs(Pw4), BULG);//+fmax(Fw4,0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 * ApproxConvective(fabs(Pn4), BULG);//+fmax(-(Fn4),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 * ApproxConvective(fabs(Ps4), BULG);//+fmax(Fs4,0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 * ApproxConvective(fabs(Pt4), BULG);//+fmax(-(Ft4),0); 
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 * ApproxConvective(fabs(Pb4), BULG);//+fmax(Fb4,0); 
			}
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
			if (b_on_adaptive_local_refinement_mesh) {
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2;
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3;
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4;
			}
		}


		// ������ ��� ����������� ������ � ��������� � ����������. 7.05.2017.
		//sumanb=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
		// ��� ���������:
		// ����� ������������� !!! (���������� �� ����).
		/*
		sumanb = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(+(Fe), 0);
		sumanb += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-(Fw), 0);
		sumanb += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(+(Fn), 0);
		sumanb += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-(Fs), 0);
		sumanb += Dt * ApproxConvective(fabs(Pt), ishconvection) + fmax(+(Ft), 0);
		sumanb += Db * ApproxConvective(fabs(Pb), ishconvection) + fmax(-(Fb), 0);
		if (b_on_adaptive_local_refinement_mesh) {
			sumanb += De2 * ApproxConvective(fabs(Pe2), ishconvection) + fmax(+(Fe2), 0);
			sumanb += Dw2 * ApproxConvective(fabs(Pw2), ishconvection) + fmax(-(Fw2), 0);
			sumanb += Dn2 * ApproxConvective(fabs(Pn2), ishconvection) + fmax(+(Fn2), 0);
			sumanb += Ds2 * ApproxConvective(fabs(Ps2), ishconvection) + fmax(-(Fs2), 0);
			sumanb += Dt2 * ApproxConvective(fabs(Pt2), ishconvection) + fmax(+(Ft2), 0);
			sumanb += Db2 * ApproxConvective(fabs(Pb2), ishconvection) + fmax(-(Fb2), 0);

			sumanb += De3 * ApproxConvective(fabs(Pe3), ishconvection) + fmax(+(Fe3), 0);
			sumanb += Dw3 * ApproxConvective(fabs(Pw3), ishconvection) + fmax(-(Fw3), 0);
			sumanb += Dn3 * ApproxConvective(fabs(Pn3), ishconvection) + fmax(+(Fn3), 0);
			sumanb += Ds3 * ApproxConvective(fabs(Ps3), ishconvection) + fmax(-(Fs3), 0);
			sumanb += Dt3 * ApproxConvective(fabs(Pt3), ishconvection) + fmax(+(Ft3), 0);
			sumanb += Db3 * ApproxConvective(fabs(Pb3), ishconvection) + fmax(-(Fb3), 0);

			sumanb += De4 * ApproxConvective(fabs(Pe4), ishconvection) + fmax(+(Fe4), 0);
			sumanb += Dw4 * ApproxConvective(fabs(Pw4), ishconvection) + fmax(-(Fw4), 0);
			sumanb += Dn4 * ApproxConvective(fabs(Pn4), ishconvection) + fmax(+(Fn4), 0);
			sumanb += Ds4 * ApproxConvective(fabs(Ps4), ishconvection) + fmax(-(Fs4), 0);
			sumanb += Dt4 * ApproxConvective(fabs(Pt4), ishconvection) + fmax(+(Ft4), 0);
			sumanb += Db4 * ApproxConvective(fabs(Pb4), ishconvection) + fmax(-(Fb4), 0);
		}
		*/
		//13 ������� 2016.
		//sumanb = fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);
		/*sumanb = sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
		if (b_on_adaptive_local_refinement_mesh) {
			sumanb += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2;
			sumanb += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3;
			sumanb += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4;
		}*/

	}
	else if (ishconvection < QUICK)
	{

		if (b_on_adaptive_local_refinement_mesh) {

			// ���������� ������������� ����������� �������:
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = -(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = (Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = -(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = (Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = -(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = (Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			// ���������� ������������� ����������� �������:
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = -(Fe2)*fC(Pe2, ishconvection, true, feplus2) + De2 * fD(Pe2, ishconvection, true, feplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = (Fw2)*fC(Pw2, ishconvection, true, fwplus2) + Dw2 * fD(Pw2, ishconvection, true, fwplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = -(Fn2)*fC(Pn2, ishconvection, true, fnplus2) + Dn2 * fD(Pn2, ishconvection, true, fnplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = (Fs2)*fC(Ps2, ishconvection, true, fsplus2) + Ds2 * fD(Ps2, ishconvection, true, fsplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = -(Ft2)*fC(Pt2, ishconvection, true, ftplus2) + Dt2 * fD(Pt2, ishconvection, true, ftplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = (Fb2)*fC(Pb2, ishconvection, true, fbplus2) + Db2 * fD(Pb2, ishconvection, true, fbplus2);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = -(Fe3)*fC(Pe3, ishconvection, true, feplus3) + De3 * fD(Pe3, ishconvection, true, feplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = (Fw3)*fC(Pw3, ishconvection, true, fwplus3) + Dw3 * fD(Pw3, ishconvection, true, fwplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = -(Fn3)*fC(Pn3, ishconvection, true, fnplus3) + Dn3 * fD(Pn3, ishconvection, true, fnplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = (Fs3)*fC(Ps3, ishconvection, true, fsplus3) + Ds3 * fD(Ps3, ishconvection, true, fsplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = -(Ft3)*fC(Pt3, ishconvection, true, ftplus3) + Dt3 * fD(Pt3, ishconvection, true, ftplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = (Fb3)*fC(Pb3, ishconvection, true, fbplus3) + Db3 * fD(Pb3, ishconvection, true, fbplus3);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = -(Fe4)*fC(Pe4, ishconvection, true, feplus4) + De4 * fD(Pe4, ishconvection, true, feplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = (Fw4)*fC(Pw4, ishconvection, true, fwplus4) + Dw4 * fD(Pw4, ishconvection, true, fwplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = -(Fn4)*fC(Pn4, ishconvection, true, fnplus4) + Dn4 * fD(Pn4, ishconvection, true, fnplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = (Fs4)*fC(Ps4, ishconvection, true, fsplus4) + Ds4 * fD(Ps4, ishconvection, true, fsplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = -(Ft4)*fC(Pt4, ishconvection, true, ftplus4) + Dt4 * fD(Pt4, ishconvection, true, ftplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = (Fb4)*fC(Pb4, ishconvection, true, fbplus4) + Db4 * fD(Pb4, ishconvection, true, fbplus4);

			// 08.05.2017.
			// ��� ���������:
			// ����� ������������� !!! (���������� �� ����).

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fe2)*fC(Pe2, ishconvection, true, feplus2) + De2 * fD(Pe2, ishconvection, true, feplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fw2)*fC(Pw2, ishconvection, true, fwplus2) + Dw2 * fD(Pw2, ishconvection, true, fwplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fn2)*fC(Pn2, ishconvection, true, fnplus2) + Dn2 * fD(Pn2, ishconvection, true, fnplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fs2)*fC(Ps2, ishconvection, true, fsplus2) + Ds2 * fD(Ps2, ishconvection, true, fsplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Ft2)*fC(Pt2, ishconvection, true, ftplus2) + Dt2 * fD(Pt2, ishconvection, true, ftplus2);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fb2)*fC(Pb2, ishconvection, true, fbplus2) + Db2 * fD(Pb2, ishconvection, true, fbplus2);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fe3)*fC(Pe3, ishconvection, true, feplus3) + De3 * fD(Pe3, ishconvection, true, feplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fw3)*fC(Pw3, ishconvection, true, fwplus3) + Dw3 * fD(Pw3, ishconvection, true, fwplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fn3)*fC(Pn3, ishconvection, true, fnplus3) + Dn3 * fD(Pn3, ishconvection, true, fnplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fs3)*fC(Ps3, ishconvection, true, fsplus3) + Ds3 * fD(Ps3, ishconvection, true, fsplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Ft3)*fC(Pt3, ishconvection, true, ftplus3) + Dt3 * fD(Pt3, ishconvection, true, ftplus3);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fb3)*fC(Pb3, ishconvection, true, fbplus3) + Db3 * fD(Pb3, ishconvection, true, fbplus3);

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fe4)*fC(Pe4, ishconvection, true, feplus4) + De4 * fD(Pe4, ishconvection, true, feplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fw4)*fC(Pw4, ishconvection, true, fwplus4) + Dw4 * fD(Pw4, ishconvection, true, fwplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fn4)*fC(Pn4, ishconvection, true, fnplus4) + Dn4 * fD(Pn4, ishconvection, true, fnplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fs4)*fC(Ps4, ishconvection, true, fsplus4) + Ds4 * fD(Ps4, ishconvection, true, fsplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Ft4)*fC(Pt4, ishconvection, true, ftplus4) + Dt4 * fD(Pt4, ishconvection, true, ftplus4);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fb4)*fC(Pb4, ishconvection, true, fbplus4) + Db4 * fD(Pb4, ishconvection, true, fbplus4);

			/*
			sumanb = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sumanb += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sumanb += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sumanb += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sumanb += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sumanb += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			sumanb += +(Fe2)*fC(Pe2, ishconvection, true, feplus2) + De2 * fD(Pe2, ishconvection, true, feplus2);
			sumanb += -(Fw2)*fC(Pw2, ishconvection, true, fwplus2) + Dw2 * fD(Pw2, ishconvection, true, fwplus2);
			sumanb += +(Fn2)*fC(Pn2, ishconvection, true, fnplus2) + Dn2 * fD(Pn2, ishconvection, true, fnplus2);
			sumanb += -(Fs2)*fC(Ps2, ishconvection, true, fsplus2) + Ds2 * fD(Ps2, ishconvection, true, fsplus2);
			sumanb += +(Ft2)*fC(Pt2, ishconvection, true, ftplus2) + Dt2 * fD(Pt2, ishconvection, true, ftplus2);
			sumanb += -(Fb2)*fC(Pb2, ishconvection, true, fbplus2) + Db2 * fD(Pb2, ishconvection, true, fbplus2);

			sumanb += +(Fe3)*fC(Pe3, ishconvection, true, feplus3) + De3 * fD(Pe3, ishconvection, true, feplus3);
			sumanb += -(Fw3)*fC(Pw3, ishconvection, true, fwplus3) + Dw3 * fD(Pw3, ishconvection, true, fwplus3);
			sumanb += +(Fn3)*fC(Pn3, ishconvection, true, fnplus3) + Dn3 * fD(Pn3, ishconvection, true, fnplus3);
			sumanb += -(Fs3)*fC(Ps3, ishconvection, true, fsplus3) + Ds3 * fD(Ps3, ishconvection, true, fsplus3);
			sumanb += +(Ft3)*fC(Pt3, ishconvection, true, ftplus3) + Dt3 * fD(Pt3, ishconvection, true, ftplus3);
			sumanb += -(Fb3)*fC(Pb3, ishconvection, true, fbplus3) + Db3 * fD(Pb3, ishconvection, true, fbplus3);

			sumanb += +(Fe4)*fC(Pe4, ishconvection, true, feplus4) + De4 * fD(Pe4, ishconvection, true, feplus4);
			sumanb += -(Fw4)*fC(Pw4, ishconvection, true, fwplus4) + Dw4 * fD(Pw4, ishconvection, true, fwplus4);
			sumanb += +(Fn4)*fC(Pn4, ishconvection, true, fnplus4) + Dn4 * fD(Pn4, ishconvection, true, fnplus4);
			sumanb += -(Fs4)*fC(Ps4, ishconvection, true, fsplus4) + Ds4 * fD(Ps4, ishconvection, true, fsplus4);
			sumanb += +(Ft4)*fC(Pt4, ishconvection, true, ftplus4) + Dt4 * fD(Pt4, ishconvection, true, ftplus4);
			sumanb += -(Fb4)*fC(Pb4, ishconvection, true, fbplus4) + Db4 * fD(Pb4, ishconvection, true, fbplus4);
			*/

		}
		else
		{
			// ���������� ������������� ����������� �������:
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = -(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = (Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = -(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = (Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = -(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = (Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);
			//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;

			// ������ ��� ����������� ������ � ��������� � ����������. 7.05.2017.
			//sumanb=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;
			//13 ������� 2016.
			//sumanb = fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);

			// 08.05.2017.
			// ��� ���������:
			// ����� ������������� !!! (���������� �� ����).

			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);

			/*
			sumanb = +(Fe)*fC(Pe, ishconvection, true, feplus) + De * fD(Pe, ishconvection, true, feplus);
			sumanb += -(Fw)*fC(Pw, ishconvection, true, fwplus) + Dw * fD(Pw, ishconvection, true, fwplus);
			sumanb += +(Fn)*fC(Pn, ishconvection, true, fnplus) + Dn * fD(Pn, ishconvection, true, fnplus);
			sumanb += -(Fs)*fC(Ps, ishconvection, true, fsplus) + Ds * fD(Ps, ishconvection, true, fsplus);
			sumanb += +(Ft)*fC(Pt, ishconvection, true, ftplus) + Dt * fD(Pt, ishconvection, true, ftplus);
			sumanb += -(Fb)*fC(Pb, ishconvection, true, fbplus) + Db * fD(Pb, ishconvection, true, fbplus);
			*/
		}
	}
	else if (ishconvection >= QUICK)
	{
		// � 3D ������������ ������ ����� ������������ �� ��� ���������� �����.
		// ����� �������� ����� ������ ������� � ����������� � ������� ��������� ���������� ���������.

		// X - direction
		doublereal positionxP, positionxE, positionxW, positionxEE, positionxWW, positionxe, positionxw;
		doublereal SpeedP = 0.0, SpeedE = 0.0, SpeedW = 0.0, SpeedEE = 0.0, SpeedWW = 0.0, SpeedN = 0.0, SpeedS = 0.0;
		doublereal SpeedNN = 0.0, SpeedSS = 0.0, SpeedT = 0.0, SpeedB = 0.0, SpeedTT = 0.0, SpeedBB = 0.0;
		doublereal Speede = 0.0, Speedw = 0.0, Speedn = 0.0, Speeds = 0.0, Speedt = 0.0, Speedb = 0.0;
		// Y - direction
		doublereal positionyP, positionyN, positionyS, positionyNN, positionySS, positionyn, positionys;
		// Z - direction
		doublereal positionzP, positionzT, positionzB, positionzTT, positionzBB, positionzt, positionzb;

		doublereal SpeedE2 = 0.0, SpeedW2 = 0.0, SpeedEE2 = 0.0, SpeedWW2 = 0.0, SpeedN2 = 0.0, SpeedS2 = 0.0;
		doublereal SpeedNN2 = 0.0, SpeedSS2 = 0.0, SpeedT2 = 0.0, SpeedB2 = 0.0, SpeedTT2 = 0.0, SpeedBB2 = 0.0;
		doublereal Speede2 = 0.0, Speedw2 = 0.0, Speedn2 = 0.0, Speeds2 = 0.0, Speedt2 = 0.0, Speedb2 = 0.0;

		doublereal  SpeedE3 = 0.0, SpeedW3 = 0.0, SpeedEE3 = 0.0, SpeedWW3 = 0.0, SpeedN3 = 0.0, SpeedS3 = 0.0;
		doublereal SpeedNN3 = 0.0, SpeedSS3 = 0.0, SpeedT3 = 0.0, SpeedB3 = 0.0, SpeedTT3 = 0.0, SpeedBB3 = 0.0;
		doublereal Speede3 = 0.0, Speedw3 = 0.0, Speedn3 = 0.0, Speeds3 = 0.0, Speedt3 = 0.0, Speedb3 = 0.0;

		doublereal SpeedE4 = 0.0, SpeedW4 = 0.0, SpeedEE4 = 0.0, SpeedWW4 = 0.0, SpeedN4 = 0.0, SpeedS4 = 0.0;
		doublereal SpeedNN4 = 0.0, SpeedSS4 = 0.0, SpeedT4 = 0.0, SpeedB4 = 0.0, SpeedTT4 = 0.0, SpeedBB4 = 0.0;
		doublereal Speede4 = 0.0, Speedw4 = 0.0, Speedn4 = 0.0, Speeds4 = 0.0, Speedt4 = 0.0, Speedb4 = 0.0;

		// X - direction
		doublereal  positionxE2, positionxW2, positionxEE2, positionxWW2, positionxe2, positionxw2;
		// Y - direction
		doublereal  positionyN2, positionyS2, positionyNN2, positionySS2, positionyn2, positionys2;
		// Z - direction
		doublereal  positionzT2, positionzB2, positionzTT2, positionzBB2, positionzt2, positionzb2;

		// X - direction
		doublereal  positionxE3, positionxW3, positionxEE3, positionxWW3, positionxe3, positionxw3;
		// Y - direction
		doublereal  positionyN3, positionyS3, positionyNN3, positionySS3, positionyn3, positionys3;
		// Z - direction
		doublereal  positionzT3, positionzB3, positionzTT3, positionzBB3, positionzt3, positionzb3;

		// X - direction
		doublereal  positionxE4, positionxW4, positionxEE4, positionxWW4, positionxe4, positionxw4;
		// Y - direction
		doublereal  positionyN4, positionyS4, positionyNN4, positionySS4, positionyn4, positionys4;
		// Z - direction
		doublereal  positionzT4, positionzB4, positionzTT4, positionzBB4, positionzt4, positionzb4;

		TOCHKA pointP;
		//center_cord3D(iP, nvtx, pa, pointP, 100);
		pointP = center_coord_loc[iP];

		positionxP = pointP.x; positionyP = pointP.y; positionzP = pointP.z;
		SpeedP = potent[RE_THETA_LANGTRY_MENTER][iP];
		// X - direction
		if (!bE) {
			SpeedE = potent[RE_THETA_LANGTRY_MENTER][iE];
			//center_cord3D(iE,nvtx,pa,pointP,E_SIDE);
			pointP = center_coord_loc[iE];

			positionxE = pointP.x;
			positionxe = positionxP + 0.5 * dx;

			integer iEE = neighbors_for_the_internal_node[E_SIDE][0][iE];
			if (iEE < 0) {
				iEE = neighbors_for_the_internal_node[E_SIDE][1][iE];
			}
			if (iEE < 0) {
				iEE = neighbors_for_the_internal_node[E_SIDE][2][iE];
			}
			if (iEE < 0) {
				iEE = neighbors_for_the_internal_node[E_SIDE][3][iE];
			}

			if ((iEE >= 0) && (iEE < maxelm)) {
				// ���������� ����
				SpeedEE = potent[RE_THETA_LANGTRY_MENTER][iEE];
				//center_cord3D(iEE,nvtx,pa,pointP,EE_SIDE);
				pointP = center_coord_loc[iEE];
				positionxEE = pointP.x;
			}
			else
			{
				// ��������� ����
				if ((iEE >= maxelm) && (iEE < maxelm + maxbound)) {
					SpeedEE = potent[RE_THETA_LANGTRY_MENTER][iEE];
				}
				else {
					SpeedEE = SpeedE;
				}
				//volume3D(iE, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iE];

				positionxEE = positionxE + 0.5 * pointP.x;
			}
		}
		else {
			// ��� ��������� ����
			SpeedE = potent[RE_THETA_LANGTRY_MENTER][iE];
			SpeedEE = potent[RE_THETA_LANGTRY_MENTER][iE];
			positionxe = positionxP + 0.5 * dx;
			positionxE = positionxP + 0.5 * dx;
			positionxEE = positionxP + dx; // ����� ���� �� ���������� !
		}

		if (!bW) {
			//center_cord3D(iW,nvtx,pa,pointP,W_SIDE);
			pointP = center_coord_loc[iW];
			positionxW = pointP.x;
			positionxw = positionxP - 0.5 * dx;
			SpeedW = potent[RE_THETA_LANGTRY_MENTER][iW];

			integer iWW = neighbors_for_the_internal_node[W_SIDE][0][iW];
			if (iWW < 0) {
				iWW = neighbors_for_the_internal_node[W_SIDE][1][iW];
			}
			if (iWW < 0) {
				iWW = neighbors_for_the_internal_node[W_SIDE][2][iW];
			}
			if (iWW < 0) {
				iWW = neighbors_for_the_internal_node[W_SIDE][3][iW];
			}

			if ((iWW >= 0) && (iWW < maxelm)) {
				// ���������� ����
				SpeedWW = potent[RE_THETA_LANGTRY_MENTER][iWW];
				//center_cord3D(iWW,nvtx,pa,pointP,WW_SIDE);
				pointP = center_coord_loc[iWW];
				positionxWW = pointP.x;
			}
			else
			{
				// ��������� ����
				if ((iWW >= maxelm) && (iWW < maxelm + maxbound)) {
					SpeedWW = potent[RE_THETA_LANGTRY_MENTER][iWW];
				}
				else {
					SpeedWW = SpeedW;
				}
				//volume3D(iW, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iW];

				positionxWW = positionxW - 0.5 * pointP.x;
			}
		}
		else {
			// ��� ��������� ����
			SpeedW = potent[RE_THETA_LANGTRY_MENTER][iW]; // Attantion !! Debug
			SpeedWW = potent[RE_THETA_LANGTRY_MENTER][iW];
			//printf("SpeedW==%e\n",SpeedW); getchar();
			positionxw = positionxP - 0.5 * dx;
			positionxW = positionxP - 0.5 * dx;
			positionxWW = positionxP - dx; // ����� ���� �� ���������� !
		}

		// Y - direction
		if (!bN) {
			SpeedN = potent[RE_THETA_LANGTRY_MENTER][iN];
			//center_cord3D(iN,nvtx,pa,pointP,N_SIDE);
			pointP = center_coord_loc[iN];
			positionyN = pointP.y;
			positionyn = positionxP + 0.5 * dy;

			integer iNN = neighbors_for_the_internal_node[N_SIDE][0][iN];
			if (iNN < 0) {
				iNN = neighbors_for_the_internal_node[N_SIDE][1][iN];
			}
			if (iNN < 0) {
				iNN = neighbors_for_the_internal_node[N_SIDE][2][iN];
			}
			if (iNN < 0) {
				iNN = neighbors_for_the_internal_node[N_SIDE][3][iN];
			}

			if ((iNN >= 0) && (iNN < maxelm)) {
				// ���������� ����
				SpeedNN = potent[RE_THETA_LANGTRY_MENTER][iNN];
				//center_cord3D(iNN,nvtx,pa,pointP,NN_SIDE);
				pointP = center_coord_loc[iNN];
				positionyNN = pointP.y;
			}
			else
			{
				// ��������� ����
				if ((iNN >= maxelm) && (iNN < maxelm + maxbound)) {
					SpeedNN = potent[RE_THETA_LANGTRY_MENTER][iNN];
				}
				else {
					SpeedNN = SpeedN;
				}
				//volume3D(iN, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iN];

				positionyNN = positionyN + 0.5 * pointP.y;
			}
		}
		else {
			// ��� ��������� ����
			SpeedN = potent[RE_THETA_LANGTRY_MENTER][iN];
			SpeedNN = potent[RE_THETA_LANGTRY_MENTER][iN];
			positionyn = positionyP + 0.5 * dy;
			positionyN = positionyP + 0.5 * dy;
			positionyNN = positionyP + dy; // ����� ���� �� ���������� !
		}

		if (!bS) {
			SpeedS = potent[RE_THETA_LANGTRY_MENTER][iS];
			//center_cord3D(iS,nvtx,pa,pointP,S_SIDE);
			pointP = center_coord_loc[iS];
			positionyS = pointP.y;
			positionys = positionyP - 0.5 * dy;

			integer iSS = neighbors_for_the_internal_node[S_SIDE][0][iS];
			if (iSS < 0) {
				iSS = neighbors_for_the_internal_node[S_SIDE][1][iS];
			}
			if (iSS < 0) {
				iSS = neighbors_for_the_internal_node[S_SIDE][2][iS];
			}
			if (iSS < 0) {
				iSS = neighbors_for_the_internal_node[S_SIDE][3][iS];
			}

			if ((iSS >= 0) && (iSS < maxelm)) {
				// ���������� ����
				SpeedSS = potent[RE_THETA_LANGTRY_MENTER][iSS];
				//center_cord3D(iSS,nvtx,pa,pointP,SS_SIDE);
				pointP = center_coord_loc[iSS];
				positionySS = pointP.y;
			}
			else
			{
				// ��������� ����
				if ((iSS >= maxelm) && (iSS < maxelm + maxbound)) {
					SpeedSS = potent[RE_THETA_LANGTRY_MENTER][iSS];
				}
				else {
					SpeedSS = SpeedS;
				}
				//volume3D(iS, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iS];

				positionySS = positionyS - 0.5 * pointP.y;
			}
		}
		else {
			// ��� ��������� ����
			SpeedS = potent[RE_THETA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
			SpeedSS = potent[RE_THETA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
			positionys = positionyP - 0.5 * dy;
			positionyS = positionyP - 0.5 * dy;
			positionySS = positionyP - dy; // ����� ���� �� ���������� !
		}

		// Z - direction
		if (!bT) {
			SpeedT = potent[RE_THETA_LANGTRY_MENTER][iT];
			//center_cord3D(iT,nvtx,pa,pointP,T_SIDE);
			pointP = center_coord_loc[iT];
			positionzT = pointP.z;
			positionzt = positionzP + 0.5 * dz;

			integer iTT = neighbors_for_the_internal_node[T_SIDE][0][iT];
			if (iTT < 0) {
				iTT = neighbors_for_the_internal_node[T_SIDE][1][iT];
			}
			if (iTT < 0) {
				iTT = neighbors_for_the_internal_node[T_SIDE][2][iT];
			}
			if (iTT < 0) {
				iTT = neighbors_for_the_internal_node[T_SIDE][3][iT];
			}

			if ((iTT >= 0) && (iTT < maxelm)) {
				// ���������� ����
				SpeedTT = potent[RE_THETA_LANGTRY_MENTER][iTT];
				//center_cord3D(iTT,nvtx,pa,pointP,TT_SIDE);
				pointP = center_coord_loc[iTT];
				positionzTT = pointP.z;
			}
			else
			{
				// ��������� ����
				if ((iTT >= maxelm) && (iTT < maxelm + maxbound)) {

					SpeedTT = potent[RE_THETA_LANGTRY_MENTER][iTT];
				}
				else {
					SpeedTT = SpeedT;
				}
				//volume3D(iT, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iT];

				positionzTT = positionzT + 0.5 * pointP.z;
			}
		}
		else {
			// ��� ��������� ����
			SpeedT = potent[RE_THETA_LANGTRY_MENTER][iT];
			SpeedTT = potent[RE_THETA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
			positionzt = positionzP + 0.5 * dz;
			positionzT = positionzP + 0.5 * dz;
			positionzTT = positionzP + dz; // ����� ���� �� ���������� !
		}

		if (!bB) {
			SpeedB = potent[RE_THETA_LANGTRY_MENTER][iB];
			//center_cord3D(iB,nvtx,pa,pointP,B_SIDE);
			pointP = center_coord_loc[iB];
			positionzB = pointP.z;
			positionzb = positionzP - 0.5 * dz;

			integer iBB = neighbors_for_the_internal_node[B_SIDE][0][iB];
			if (iBB < 0) {
				iBB = neighbors_for_the_internal_node[B_SIDE][1][iB];
			}
			if (iBB < 0) {
				iBB = neighbors_for_the_internal_node[B_SIDE][2][iB];
			}
			if (iBB < 0) {
				iBB = neighbors_for_the_internal_node[B_SIDE][3][iB];
			}

			if ((iBB >= 0) && (iBB < maxelm)) {
				// ���������� ����
				SpeedBB = potent[RE_THETA_LANGTRY_MENTER][iBB];
				//center_cord3D(iBB,nvtx,pa,pointP,BB_SIDE);
				pointP = center_coord_loc[iBB];
				positionzBB = pointP.z;
			}
			else
			{
				// ��������� ����
				if ((iBB >= maxelm) && (iBB < maxelm + maxbound)) {
					SpeedBB = potent[RE_THETA_LANGTRY_MENTER][iBB];
				}
				else {
					SpeedBB = SpeedB;
				}
				//volume3D(iB, nvtx, pa, pointP.x, pointP.y, pointP.z);
				pointP = volume_loc[iB];

				positionzBB = positionzB - 0.5 * pointP.z;
			}
		}
		else {
			// ��� ��������� ����
			SpeedB = potent[RE_THETA_LANGTRY_MENTER][iB];
			SpeedBB = potent[RE_THETA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
			positionzb = positionzP - 0.5 * dz;
			positionzB = positionzP - 0.5 * dz;
			positionzBB = positionzP - dz; // ����� ���� �� ���������� !
		}

		if (b_on_adaptive_local_refinement_mesh)
		{

			// X - direction
			if ((!bE2) && (iE2 > -1)) {
				SpeedE2 = potent[RE_THETA_LANGTRY_MENTER][iE2];
				//center_cord3D(iE,nvtx,pa,pointP,E_SIDE);
				pointP = center_coord_loc[iE2];

				positionxE2 = pointP.x;
				positionxe2 = positionxP + 0.5 * dx;

				integer iEE2 = neighbors_for_the_internal_node[E_SIDE][0][iE2];
				if (iEE2 < 0) {
					iEE2 = neighbors_for_the_internal_node[E_SIDE][1][iE2];
				}
				if (iEE2 < 0) {
					iEE2 = neighbors_for_the_internal_node[E_SIDE][2][iE2];
				}
				if (iEE2 < 0) {
					iEE2 = neighbors_for_the_internal_node[E_SIDE][3][iE2];
				}
				if ((iEE2 >= 0) && (iEE2 < maxelm)) {
					// ���������� ����
					SpeedEE2 = potent[RE_THETA_LANGTRY_MENTER][iEE2];
					//center_cord3D(iEE2,nvtx,pa,pointP,EE_SIDE);
					pointP = center_coord_loc[iEE2];
					positionxEE2 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iEE2 >= maxelm) && (iEE2 < maxelm + maxbound)) {
						SpeedEE2 = potent[RE_THETA_LANGTRY_MENTER][iEE2];
					}
					else {
						SpeedEE2 = SpeedE2;
						//std::cout << "iEE2 =" << iEE2 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iE, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iE2];

					positionxEE2 = positionxE2 + 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iE2 > -1) {
					SpeedE2 = potent[RE_THETA_LANGTRY_MENTER][iE2];
					SpeedEE2 = potent[RE_THETA_LANGTRY_MENTER][iE2];
				}
				else {
					SpeedE2 = potent[RE_THETA_LANGTRY_MENTER][iE];
					SpeedEE2 = potent[RE_THETA_LANGTRY_MENTER][iE];
				}
				positionxe2 = positionxP + 0.5 * dx;
				positionxE2 = positionxP + 0.5 * dx;
				positionxEE2 = positionxP + dx; // ����� ���� �� ���������� !
			}

			if ((!bW2) && ((iW2 > -1))) {
				//center_cord3D(iW,nvtx,pa,pointP,W_SIDE);
				pointP = center_coord_loc[iW2];
				positionxW2 = pointP.x;
				positionxw2 = positionxP - 0.5 * dx;
				SpeedW2 = potent[RE_THETA_LANGTRY_MENTER][iW2];

				integer iWW2 = neighbors_for_the_internal_node[W_SIDE][0][iW2];
				if (iWW2 < 0) {
					iWW2 = neighbors_for_the_internal_node[W_SIDE][1][iW2];
				}
				if (iWW2 < 0) {
					iWW2 = neighbors_for_the_internal_node[W_SIDE][2][iW2];
				}
				if (iWW2 < 0) {
					iWW2 = neighbors_for_the_internal_node[W_SIDE][3][iW2];
				}

				if ((iWW2 >= 0) && (iWW2 < maxelm)) {
					// ���������� ����
					SpeedWW2 = potent[RE_THETA_LANGTRY_MENTER][iWW2];
					//center_cord3D(iWW2,nvtx,pa,pointP,WW_SIDE);
					pointP = center_coord_loc[iWW2];
					positionxWW2 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iWW2 >= maxelm) && (iWW2 < maxelm + maxbound)) {
						SpeedWW2 = potent[RE_THETA_LANGTRY_MENTER][iWW2];
					}
					else {

						SpeedWW2 = SpeedW2;
						//std::cout << "iWW2 =" << iWW2 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iW2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iW2];

					positionxWW2 = positionxW2 - 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iW2 > -1) {
					SpeedW2 = potent[RE_THETA_LANGTRY_MENTER][iW2];
					SpeedWW2 = potent[RE_THETA_LANGTRY_MENTER][iW2];
				}
				else {
					SpeedW2 = potent[RE_THETA_LANGTRY_MENTER][iW];
					SpeedWW2 = potent[RE_THETA_LANGTRY_MENTER][iW];
				}
				//printf("SpeedW2==%e\n",SpeedW2); getchar();
				positionxw2 = positionxP - 0.5 * dx;
				positionxW2 = positionxP - 0.5 * dx;
				positionxWW2 = positionxP - dx; // ����� ���� �� ���������� !
			}

			// Y - direction
			if ((!bN2) && (iN2 > -1)) {
				SpeedN2 = potent[RE_THETA_LANGTRY_MENTER][iN2];
				//center_cord3D(iN2,nvtx,pa,pointP,N_SIDE);
				pointP = center_coord_loc[iN2];
				positionyN2 = pointP.y;
				positionyn2 = positionxP + 0.5 * dy;

				integer iNN2 = neighbors_for_the_internal_node[N_SIDE][0][iN2];
				if (iNN2 < 0) {
					iNN2 = neighbors_for_the_internal_node[N_SIDE][1][iN2];
				}
				if (iNN2 < 0) {
					iNN2 = neighbors_for_the_internal_node[N_SIDE][2][iN2];
				}
				if (iNN2 < 0) {
					iNN2 = neighbors_for_the_internal_node[N_SIDE][3][iN2];
				}


				if ((iNN2 >= 0) && (iNN2 < maxelm)) {
					// ���������� ����
					SpeedNN2 = potent[RE_THETA_LANGTRY_MENTER][iNN2];
					//center_cord3D(iNN2,nvtx,pa,pointP,NN_SIDE);
					pointP = center_coord_loc[iNN2];
					positionyNN2 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iNN2 >= maxelm) && (iNN2 < maxelm + maxbound)) {
						SpeedNN2 = potent[RE_THETA_LANGTRY_MENTER][iNN2];
					}
					else {
						SpeedNN2 = SpeedN2;
					}
					//volume3D(iN2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iN2];

					positionyNN2 = positionyN2 + 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iN2 > -1) {
					SpeedN2 = potent[RE_THETA_LANGTRY_MENTER][iN2];
					SpeedNN2 = potent[RE_THETA_LANGTRY_MENTER][iN2];
				}
				else {
					SpeedN2 = potent[RE_THETA_LANGTRY_MENTER][iN];
					SpeedNN2 = potent[RE_THETA_LANGTRY_MENTER][iN];
				}
				positionyn2 = positionyP + 0.5 * dy;
				positionyN2 = positionyP + 0.5 * dy;
				positionyNN2 = positionyP + dy; // ����� ���� �� ���������� !
			}

			if ((!bS2) && (iS2 > -1)) {
				SpeedS2 = potent[RE_THETA_LANGTRY_MENTER][iS2];
				//center_cord3D(iS2,nvtx,pa,pointP,S_SIDE);
				pointP = center_coord_loc[iS2];
				positionyS2 = pointP.y;
				positionys2 = positionyP - 0.5 * dy;

				integer iSS2 = neighbors_for_the_internal_node[S_SIDE][0][iS2];
				if (iSS2 < 0) {
					iSS2 = neighbors_for_the_internal_node[S_SIDE][1][iS2];
				}
				if (iSS2 < 0) {
					iSS2 = neighbors_for_the_internal_node[S_SIDE][2][iS2];
				}
				if (iSS2 < 0) {
					iSS2 = neighbors_for_the_internal_node[S_SIDE][3][iS2];
				}

				if ((iSS2 >= 0) && (iSS2 < maxelm)) {
					// ���������� ����
					SpeedSS2 = potent[RE_THETA_LANGTRY_MENTER][iSS2];
					//center_cord3D(iSS2,nvtx,pa,pointP,SS_SIDE);
					pointP = center_coord_loc[iSS2];
					positionySS2 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iSS2 >= maxelm) && (iSS2 < maxelm + maxbound)) {
						SpeedSS2 = potent[RE_THETA_LANGTRY_MENTER][iSS2];
					}
					else {
						SpeedSS2 = SpeedS2;
					}
					//volume3D(iS2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iS2];

					positionySS2 = positionyS2 - 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iS2 > -1) {
					SpeedS2 = potent[RE_THETA_LANGTRY_MENTER][iS2]; // ATTANTION !!!!
					SpeedSS2 = potent[RE_THETA_LANGTRY_MENTER][iS2]; // ������� �������� ������ ������� ����.
				}
				else {
					SpeedS2 = potent[RE_THETA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
					SpeedSS2 = potent[RE_THETA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
				}
				positionys2 = positionyP - 0.5 * dy;
				positionyS2 = positionyP - 0.5 * dy;
				positionySS2 = positionyP - dy; // ����� ���� �� ���������� !
			}

			// Z - direction
			if ((!bT2) && (iT2 > -1)) {
				SpeedT2 = potent[RE_THETA_LANGTRY_MENTER][iT2];
				//center_cord3D(iT2,nvtx,pa,pointP,T_SIDE);
				pointP = center_coord_loc[iT2];
				positionzT2 = pointP.z;
				positionzt2 = positionzP + 0.5 * dz;

				integer iTT2 = neighbors_for_the_internal_node[T_SIDE][0][iT2];
				if (iTT2 < 0) {
					iTT2 = neighbors_for_the_internal_node[T_SIDE][1][iT2];
				}
				if (iTT2 < 0) {
					iTT2 = neighbors_for_the_internal_node[T_SIDE][2][iT2];
				}
				if (iTT2 < 0) {
					iTT2 = neighbors_for_the_internal_node[T_SIDE][3][iT2];
				}


				if ((iTT2 >= 0) && (iTT2 < maxelm)) {
					// ���������� ����
					SpeedTT2 = potent[RE_THETA_LANGTRY_MENTER][iTT2];
					//center_cord3D(iTT2,nvtx,pa,pointP,TT_SIDE);
					pointP = center_coord_loc[iTT2];
					positionzTT2 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iTT2 >= maxelm) && (iTT2 < maxelm + maxbound)) {
						SpeedTT2 = potent[RE_THETA_LANGTRY_MENTER][iTT2];
					}
					else {
						SpeedTT2 = SpeedT2;
					}
					//volume3D(iT2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iT2];

					positionzTT2 = positionzT2 + 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iT2 > -1) {
					SpeedT2 = potent[RE_THETA_LANGTRY_MENTER][iT2];
					SpeedTT2 = potent[RE_THETA_LANGTRY_MENTER][iT2]; // �������� ������ ������� ����
				}
				else {
					SpeedT2 = potent[RE_THETA_LANGTRY_MENTER][iT];
					SpeedTT2 = potent[RE_THETA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
				}
				positionzt2 = positionzP + 0.5 * dz;
				positionzT2 = positionzP + 0.5 * dz;
				positionzTT2 = positionzP + dz; // ����� ���� �� ���������� !
			}

			if ((!bB2) && (iB2 > -1)) {
				SpeedB2 = potent[RE_THETA_LANGTRY_MENTER][iB2];
				//center_cord3D(iB2,nvtx,pa,pointP,B_SIDE);
				pointP = center_coord_loc[iB2];
				positionzB2 = pointP.z;
				positionzb2 = positionzP - 0.5 * dz;

				integer iBB2 = neighbors_for_the_internal_node[B_SIDE][0][iB2];
				if (iBB2 < 0) {
					iBB2 = neighbors_for_the_internal_node[B_SIDE][1][iB2];
				}
				if (iBB2 < 0) {
					iBB2 = neighbors_for_the_internal_node[B_SIDE][2][iB2];
				}
				if (iBB2 < 0) {
					iBB2 = neighbors_for_the_internal_node[B_SIDE][3][iB2];
				}

				if ((iBB2 >= 0) && (iBB2 < maxelm)) {
					// ���������� ����
					SpeedBB2 = potent[RE_THETA_LANGTRY_MENTER][iBB2];
					//center_cord3D(iBB2,nvtx,pa,pointP,BB_SIDE);
					pointP = center_coord_loc[iBB2];
					positionzBB2 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iBB2 >= maxelm) && (iBB2 < maxelm + maxbound)) {
						SpeedBB2 = potent[RE_THETA_LANGTRY_MENTER][iBB2];
					}
					else {
						SpeedBB2 = SpeedB2;
					}
					//volume3D(iB2, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iB2];

					positionzBB2 = positionzB2 - 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iB2 > -1) {
					SpeedB2 = potent[RE_THETA_LANGTRY_MENTER][iB2];
					SpeedBB2 = potent[RE_THETA_LANGTRY_MENTER][iB2]; // �������� ������ ������� ����
				}
				else {
					SpeedB2 = potent[RE_THETA_LANGTRY_MENTER][iB];
					SpeedBB2 = potent[RE_THETA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
				}
				positionzb2 = positionzP - 0.5 * dz;
				positionzB2 = positionzP - 0.5 * dz;
				positionzBB2 = positionzP - dz; // ����� ���� �� ���������� !
			}


			// X - direction
			if ((!bE3) && (iE3 > -1)) {
				SpeedE3 = potent[RE_THETA_LANGTRY_MENTER][iE3];
				//center_cord3D(iE3,nvtx,pa,pointP,E_SIDE);
				pointP = center_coord_loc[iE3];

				positionxE3 = pointP.x;
				positionxe3 = positionxP + 0.5 * dx;

				integer iEE3 = neighbors_for_the_internal_node[E_SIDE][0][iE3];
				if (iEE3 < 0) {
					iEE3 = neighbors_for_the_internal_node[E_SIDE][1][iE3];
				}
				if (iEE3 < 0) {
					iEE3 = neighbors_for_the_internal_node[E_SIDE][2][iE3];
				}
				if (iEE3 < 0) {
					iEE3 = neighbors_for_the_internal_node[E_SIDE][3][iE3];
				}

				if ((iEE3 >= 0) && (iEE3 < maxelm)) {
					// ���������� ����
					SpeedEE3 = potent[RE_THETA_LANGTRY_MENTER][iEE3];
					//center_cord3D(iEE3,nvtx,pa,pointP,EE_SIDE);
					pointP = center_coord_loc[iEE3];
					positionxEE3 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iEE3 >= maxelm) && (iEE3 < maxelm + maxbound)) {
						SpeedEE3 = potent[RE_THETA_LANGTRY_MENTER][iEE3];
					}
					else {
						SpeedEE3 = SpeedE3;
					}
					//volume3D(iE3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iE3];

					positionxEE3 = positionxE3 + 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iE3 > -1) {
					SpeedE3 = potent[RE_THETA_LANGTRY_MENTER][iE3];
					SpeedEE3 = potent[RE_THETA_LANGTRY_MENTER][iE3];
				}
				else {
					SpeedE3 = potent[RE_THETA_LANGTRY_MENTER][iE];
					SpeedEE3 = potent[RE_THETA_LANGTRY_MENTER][iE];
				}
				positionxe3 = positionxP + 0.5 * dx;
				positionxE3 = positionxP + 0.5 * dx;
				positionxEE3 = positionxP + dx; // ����� ���� �� ���������� !
			}

			if ((!bW3) && (iW3 > -1)) {
				//center_cord3D(iW3,nvtx,pa,pointP,W_SIDE);
				pointP = center_coord_loc[iW3];
				positionxW3 = pointP.x;
				positionxw3 = positionxP - 0.5 * dx;
				SpeedW3 = potent[RE_THETA_LANGTRY_MENTER][iW3];

				integer iWW3 = neighbors_for_the_internal_node[W_SIDE][0][iW3];
				if (iWW3 < 0) {
					iWW3 = neighbors_for_the_internal_node[W_SIDE][1][iW3];
				}
				if (iWW3 < 0) {
					iWW3 = neighbors_for_the_internal_node[W_SIDE][2][iW3];
				}
				if (iWW3 < 0) {
					iWW3 = neighbors_for_the_internal_node[W_SIDE][3][iW3];
				}

				if ((iWW3 >= 0) && (iWW3 < maxelm)) {
					// ���������� ����
					SpeedWW3 = potent[RE_THETA_LANGTRY_MENTER][iWW3];
					//center_cord3D(iWW3,nvtx,pa,pointP,WW_SIDE);
					pointP = center_coord_loc[iWW3];
					positionxWW3 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iWW3 >= maxelm) && (iWW3 < maxelm + maxbound)) {
						SpeedWW3 = potent[RE_THETA_LANGTRY_MENTER][iWW3];
					}
					else {
						SpeedWW3 = SpeedW3;
						//std::cout << "iWW3 =" << iWW3 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iW3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iW3];

					positionxWW3 = positionxW3 - 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iW3 > -1) {
					SpeedW3 = potent[RE_THETA_LANGTRY_MENTER][iW3]; // Attantion !! Debug
					SpeedWW3 = potent[RE_THETA_LANGTRY_MENTER][iW3];
				}
				else {
					SpeedW3 = potent[RE_THETA_LANGTRY_MENTER][iW]; // Attantion !! Debug
					SpeedWW3 = potent[RE_THETA_LANGTRY_MENTER][iW];
				}
				//printf("SpeedW3==%e\n",SpeedW3); getchar();
				positionxw3 = positionxP - 0.5 * dx;
				positionxW3 = positionxP - 0.5 * dx;
				positionxWW3 = positionxP - dx; // ����� ���� �� ���������� !
			}

			// Y - direction
			if ((!bN3) && (iN3 > -1)) {
				SpeedN3 = potent[RE_THETA_LANGTRY_MENTER][iN3];
				//center_cord3D(iN3,nvtx,pa,pointP,N_SIDE);
				pointP = center_coord_loc[iN3];
				positionyN3 = pointP.y;
				positionyn3 = positionxP + 0.5 * dy;

				integer iNN3 = neighbors_for_the_internal_node[N_SIDE][0][iN3];
				if (iNN3 < 0) {
					iNN3 = neighbors_for_the_internal_node[N_SIDE][1][iN3];
				}
				if (iNN3 < 0) {
					iNN3 = neighbors_for_the_internal_node[N_SIDE][2][iN3];
				}
				if (iNN3 < 0) {
					iNN3 = neighbors_for_the_internal_node[N_SIDE][3][iN3];
				}

				if ((iNN3 >= 0) && (iNN3 < maxelm)) {
					// ���������� ����
					SpeedNN3 = potent[RE_THETA_LANGTRY_MENTER][iNN3];
					//center_cord3D(iNN3,nvtx,pa,pointP,NN_SIDE);
					pointP = center_coord_loc[iNN3];
					positionyNN3 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iNN3 >= maxelm) && (iNN3 < maxelm + maxbound)) {
						SpeedNN3 = potent[RE_THETA_LANGTRY_MENTER][iNN3];
					}
					else {
						SpeedNN3 = SpeedN3;
					}
					//volume3D(iN3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iN3];

					positionyNN3 = positionyN3 + 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iN3 > -1) {
					SpeedN3 = potent[RE_THETA_LANGTRY_MENTER][iN3];
					SpeedNN3 = potent[RE_THETA_LANGTRY_MENTER][iN3];
				}
				else {
					SpeedN3 = potent[RE_THETA_LANGTRY_MENTER][iN];
					SpeedNN3 = potent[RE_THETA_LANGTRY_MENTER][iN];
				}
				positionyn3 = positionyP + 0.5 * dy;
				positionyN3 = positionyP + 0.5 * dy;
				positionyNN3 = positionyP + dy; // ����� ���� �� ���������� !
			}

			if ((!bS3) && (iS3 > -1)) {
				SpeedS3 = potent[RE_THETA_LANGTRY_MENTER][iS3];
				//center_cord3D(iS3,nvtx,pa,pointP,S_SIDE);
				pointP = center_coord_loc[iS3];
				positionyS3 = pointP.y;
				positionys3 = positionyP - 0.5 * dy;

				integer iSS3 = neighbors_for_the_internal_node[S_SIDE][0][iS3];
				if (iSS3 < 0) {
					iSS3 = neighbors_for_the_internal_node[S_SIDE][1][iS3];
				}
				if (iSS3 < 0) {
					iSS3 = neighbors_for_the_internal_node[S_SIDE][2][iS3];
				}
				if (iSS3 < 0) {
					iSS3 = neighbors_for_the_internal_node[S_SIDE][3][iS3];
				}

				if ((iSS3 >= 0) && (iSS3 < maxelm)) {
					// ���������� ����
					SpeedSS3 = potent[RE_THETA_LANGTRY_MENTER][iSS3];
					//center_cord3D(iSS3,nvtx,pa,pointP,SS_SIDE);
					pointP = center_coord_loc[iSS3];
					positionySS3 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iSS3 >= maxelm) && (iSS3 < maxelm + maxbound)) {
						SpeedSS3 = potent[RE_THETA_LANGTRY_MENTER][iSS3];
					}
					else {
						SpeedSS3 = SpeedS3;
					}
					//volume3D(iS3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iS3];

					positionySS3 = positionyS3 - 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iS3 > -1) {
					SpeedS3 = potent[RE_THETA_LANGTRY_MENTER][iS3]; // ATTANTION !!!!
					SpeedSS3 = potent[RE_THETA_LANGTRY_MENTER][iS3]; // ������� �������� ������ ������� ����.
				}
				else {
					SpeedS3 = potent[RE_THETA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
					SpeedSS3 = potent[RE_THETA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
				}
				positionys3 = positionyP - 0.5 * dy;
				positionyS3 = positionyP - 0.5 * dy;
				positionySS3 = positionyP - dy; // ����� ���� �� ���������� !
			}

			// Z - direction
			if ((!bT3) && (iT3 > -1)) {
				SpeedT3 = potent[RE_THETA_LANGTRY_MENTER][iT3];
				//center_cord3D(iT3,nvtx,pa,pointP,T_SIDE);
				pointP = center_coord_loc[iT3];
				positionzT3 = pointP.z;
				positionzt3 = positionzP + 0.5 * dz;

				integer iTT3 = neighbors_for_the_internal_node[T_SIDE][0][iT3];
				if (iTT3 < 0) {
					iTT3 = neighbors_for_the_internal_node[T_SIDE][1][iT3];
				}
				if (iTT3 < 0) {
					iTT3 = neighbors_for_the_internal_node[T_SIDE][2][iT3];
				}
				if (iTT3 < 0) {
					iTT3 = neighbors_for_the_internal_node[T_SIDE][3][iT3];
				}

				if ((iTT3 >= 0) && (iTT3 < maxelm)) {
					// ���������� ����
					SpeedTT3 = potent[RE_THETA_LANGTRY_MENTER][iTT3];
					//center_cord3D(iTT3,nvtx,pa,pointP,TT_SIDE);
					pointP = center_coord_loc[iTT3];
					positionzTT3 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iTT3 >= maxelm) && (iTT3 < maxelm + maxbound)) {
						SpeedTT3 = potent[RE_THETA_LANGTRY_MENTER][iTT3];
					}
					else {
						SpeedTT3 = SpeedT3;
					}
					//volume3D(iT3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iT3];

					positionzTT3 = positionzT3 + 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iT3 > -1) {
					SpeedT3 = potent[RE_THETA_LANGTRY_MENTER][iT3];
					SpeedTT3 = potent[RE_THETA_LANGTRY_MENTER][iT3]; // �������� ������ ������� ����
				}
				else {
					SpeedT3 = potent[RE_THETA_LANGTRY_MENTER][iT];
					SpeedTT3 = potent[RE_THETA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
				}
				positionzt3 = positionzP + 0.5 * dz;
				positionzT3 = positionzP + 0.5 * dz;
				positionzTT3 = positionzP + dz; // ����� ���� �� ���������� !
			}

			if ((!bB3) && (iB3 > -1)) {
				SpeedB3 = potent[RE_THETA_LANGTRY_MENTER][iB3];
				//center_cord3D(iB3,nvtx,pa,pointP,B_SIDE);
				pointP = center_coord_loc[iB3];
				positionzB3 = pointP.z;
				positionzb3 = positionzP - 0.5 * dz;

				integer iBB3 = neighbors_for_the_internal_node[B_SIDE][0][iB3];
				if (iBB3 < 0) {
					iBB3 = neighbors_for_the_internal_node[B_SIDE][1][iB3];
				}
				if (iBB3 < 0) {
					iBB3 = neighbors_for_the_internal_node[B_SIDE][2][iB3];
				}
				if (iBB3 < 0) {
					iBB3 = neighbors_for_the_internal_node[B_SIDE][3][iB3];
				}

				if ((iBB3 >= 0) && (iBB3 < maxelm)) {
					// ���������� ����
					SpeedBB3 = potent[RE_THETA_LANGTRY_MENTER][iBB3];
					//center_cord3D(iBB3,nvtx,pa,pointP,BB_SIDE);
					pointP = center_coord_loc[iBB3];
					positionzBB3 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iBB3 >= maxelm) && (iBB3 < maxelm + maxbound)) {
						SpeedBB3 = potent[RE_THETA_LANGTRY_MENTER][iBB3];
					}
					else {
						SpeedBB3 = SpeedB3;
					}
					//volume3D(iB3, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iB3];

					positionzBB3 = positionzB3 - 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iB3 > -1) {
					SpeedB3 = potent[RE_THETA_LANGTRY_MENTER][iB3];
					SpeedBB3 = potent[RE_THETA_LANGTRY_MENTER][iB3]; // �������� ������ ������� ����
				}
				else {
					SpeedB3 = potent[RE_THETA_LANGTRY_MENTER][iB];
					SpeedBB3 = potent[RE_THETA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
				}
				positionzb3 = positionzP - 0.5 * dz;
				positionzB3 = positionzP - 0.5 * dz;
				positionzBB3 = positionzP - dz; // ����� ���� �� ���������� !
			}

			// X - direction
			if ((!bE4) && (iE4 > -1)) {
				SpeedE4 = potent[RE_THETA_LANGTRY_MENTER][iE4];
				//center_cord3D(iE4,nvtx,pa,pointP,E_SIDE);
				pointP = center_coord_loc[iE4];

				positionxE4 = pointP.x;
				positionxe4 = positionxP + 0.5 * dx;

				integer iEE4 = neighbors_for_the_internal_node[E_SIDE][0][iE4];
				if (iEE4 < 0) {
					iEE4 = neighbors_for_the_internal_node[E_SIDE][1][iE4];
				}
				if (iEE4 < 0) {
					iEE4 = neighbors_for_the_internal_node[E_SIDE][2][iE4];
				}
				if (iEE4 < 0) {
					iEE4 = neighbors_for_the_internal_node[E_SIDE][3][iE4];
				}

				if ((iEE4 >= 0) && (iEE4 < maxelm)) {
					// ���������� ����
					SpeedEE4 = potent[RE_THETA_LANGTRY_MENTER][iEE4];
					//center_cord3D(iEE4,nvtx,pa,pointP,EE_SIDE);
					pointP = center_coord_loc[iEE4];
					positionxEE4 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iEE4 >= maxelm) && (iEE4 < maxelm + maxbound)) {
						SpeedEE4 = potent[RE_THETA_LANGTRY_MENTER][iEE4];
					}
					else {
						SpeedEE4 = SpeedE4;
					}
					//volume3D(iE4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iE4];

					positionxEE4 = positionxE4 + 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iE4 > -1) {
					SpeedE4 = potent[RE_THETA_LANGTRY_MENTER][iE4];
					SpeedEE4 = potent[RE_THETA_LANGTRY_MENTER][iE4];
				}
				else {
					SpeedE4 = potent[RE_THETA_LANGTRY_MENTER][iE];
					SpeedEE4 = potent[RE_THETA_LANGTRY_MENTER][iE];
				}
				positionxe4 = positionxP + 0.5 * dx;
				positionxE4 = positionxP + 0.5 * dx;
				positionxEE4 = positionxP + dx; // ����� ���� �� ���������� !
			}

			if ((!bW4) && (iW4 > -1)) {
				//center_cord3D(iW4,nvtx,pa,pointP,W_SIDE);
				pointP = center_coord_loc[iW4];
				positionxW4 = pointP.x;
				positionxw4 = positionxP - 0.5 * dx;
				SpeedW4 = potent[RE_THETA_LANGTRY_MENTER][iW4];

				integer iWW4 = neighbors_for_the_internal_node[W_SIDE][0][iW4];
				if (iWW4 < 0) {
					iWW4 = neighbors_for_the_internal_node[W_SIDE][1][iW4];
				}
				if (iWW4 < 0) {
					iWW4 = neighbors_for_the_internal_node[W_SIDE][2][iW4];
				}
				if (iWW4 < 0) {
					iWW4 = neighbors_for_the_internal_node[W_SIDE][3][iW4];
				}

				if ((iWW4 >= 0) && (iWW4 < maxelm)) {
					// ���������� ����
					SpeedWW4 = potent[RE_THETA_LANGTRY_MENTER][iWW4];
					//center_cord3D(iWW4,nvtx,pa,pointP,WW_SIDE);
					pointP = center_coord_loc[iWW4];
					positionxWW4 = pointP.x;
				}
				else
				{
					// ��������� ����
					if ((iWW4 >= maxelm) && (iWW4 < maxelm + maxbound)) {
						SpeedWW4 = potent[RE_THETA_LANGTRY_MENTER][iWW4];
					}
					else {
						SpeedWW4 = SpeedW4;
						//std::cout << "iWW4 =" << iWW4 << " maxelm=" << maxelm << " maxbound=" << maxbound << std::endl;
						//getchar();
					}
					//volume3D(iW4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iW4];

					positionxWW4 = positionxW4 - 0.5 * pointP.x;
				}
			}
			else {
				// ��� ��������� ����
				if (iW4 > -1) {
					SpeedW4 = potent[RE_THETA_LANGTRY_MENTER][iW4]; // Attantion !! Debug
					SpeedWW4 = potent[RE_THETA_LANGTRY_MENTER][iW4];
				}
				else {
					SpeedW4 = potent[RE_THETA_LANGTRY_MENTER][iW]; // Attantion !! Debug
					SpeedWW4 = potent[RE_THETA_LANGTRY_MENTER][iW];
				}
				//printf("SpeedW4==%e\n",SpeedW4); getchar();
				positionxw4 = positionxP - 0.5 * dx;
				positionxW4 = positionxP - 0.5 * dx;
				positionxWW4 = positionxP - dx; // ����� ���� �� ���������� !
			}

			// Y - direction
			if ((!bN4) && (iN4 > -1)) {
				SpeedN4 = potent[RE_THETA_LANGTRY_MENTER][iN4];
				//center_cord3D(iN4,nvtx,pa,pointP,N_SIDE);
				pointP = center_coord_loc[iN4];
				positionyN4 = pointP.y;
				positionyn4 = positionxP + 0.5 * dy;

				integer iNN4 = neighbors_for_the_internal_node[N_SIDE][0][iN4];
				if (iNN4 < 0) {
					iNN4 = neighbors_for_the_internal_node[N_SIDE][1][iN4];
				}
				if (iNN4 < 0) {
					iNN4 = neighbors_for_the_internal_node[N_SIDE][2][iN4];
				}
				if (iNN4 < 0) {
					iNN4 = neighbors_for_the_internal_node[N_SIDE][3][iN4];
				}

				if ((iNN4 >= 0) && (iNN4 < maxelm)) {
					// ���������� ����
					SpeedNN4 = potent[RE_THETA_LANGTRY_MENTER][iNN4];
					//center_cord3D(iNN4,nvtx,pa,pointP,NN_SIDE);
					pointP = center_coord_loc[iNN4];
					positionyNN4 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iNN4 >= maxelm) && (iNN4 < maxelm + maxbound)) {
						SpeedNN4 = potent[RE_THETA_LANGTRY_MENTER][iNN4];
					}
					else {
						SpeedNN4 = SpeedN4;
					}
					//volume3D(iN4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iN4];

					positionyNN4 = positionyN4 + 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iN4 > -1) {
					SpeedN4 = potent[RE_THETA_LANGTRY_MENTER][iN4];
					SpeedNN4 = potent[RE_THETA_LANGTRY_MENTER][iN4];
				}
				else {
					SpeedN4 = potent[RE_THETA_LANGTRY_MENTER][iN];
					SpeedNN4 = potent[RE_THETA_LANGTRY_MENTER][iN];
				}
				positionyn4 = positionyP + 0.5 * dy;
				positionyN4 = positionyP + 0.5 * dy;
				positionyNN4 = positionyP + dy; // ����� ���� �� ���������� !
			}

			if ((!bS4) && (iS4 > -1)) {
				SpeedS4 = potent[RE_THETA_LANGTRY_MENTER][iS4];
				//center_cord3D(iS4,nvtx,pa,pointP,S_SIDE);
				pointP = center_coord_loc[iS4];
				positionyS4 = pointP.y;
				positionys4 = positionyP - 0.5 * dy;

				integer iSS4 = neighbors_for_the_internal_node[S_SIDE][0][iS4];
				if (iSS4 < 0) {
					iSS4 = neighbors_for_the_internal_node[S_SIDE][1][iS4];
				}
				if (iSS4 < 0) {
					iSS4 = neighbors_for_the_internal_node[S_SIDE][2][iS4];
				}
				if (iSS4 < 0) {
					iSS4 = neighbors_for_the_internal_node[S_SIDE][3][iS4];
				}

				if ((iSS4 >= 0) && (iSS4 < maxelm)) {
					// ���������� ����
					SpeedSS4 = potent[RE_THETA_LANGTRY_MENTER][iSS4];
					//center_cord3D(iSS4,nvtx,pa,pointP,SS_SIDE);
					pointP = center_coord_loc[iSS4];
					positionySS4 = pointP.y;
				}
				else
				{
					// ��������� ����
					if ((iSS4 >= maxelm) && (iSS4 < maxelm + maxbound)) {
						SpeedSS4 = potent[RE_THETA_LANGTRY_MENTER][iSS4];
					}
					else {
						SpeedSS4 = SpeedS4;
					}
					//volume3D(iS4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iS4];

					positionySS4 = positionyS4 - 0.5 * pointP.y;
				}
			}
			else {
				// ��� ��������� ����
				if (iS4 > -1) {
					SpeedS4 = potent[RE_THETA_LANGTRY_MENTER][iS4]; // ATTANTION !!!!
					SpeedSS4 = potent[RE_THETA_LANGTRY_MENTER][iS4]; // ������� �������� ������ ������� ����.
				}
				else {
					SpeedS4 = potent[RE_THETA_LANGTRY_MENTER][iS]; // ATTANTION !!!!
					SpeedSS4 = potent[RE_THETA_LANGTRY_MENTER][iS]; // ������� �������� ������ ������� ����.
				}
				positionys4 = positionyP - 0.5 * dy;
				positionyS4 = positionyP - 0.5 * dy;
				positionySS4 = positionyP - dy; // ����� ���� �� ���������� !
			}

			// Z - direction
			if ((!bT4) && (iT4 > -1)) {
				SpeedT4 = potent[RE_THETA_LANGTRY_MENTER][iT4];
				//center_cord3D(iT4,nvtx,pa,pointP,T_SIDE);
				pointP = center_coord_loc[iT4];
				positionzT4 = pointP.z;
				positionzt4 = positionzP + 0.5 * dz;

				integer iTT4 = neighbors_for_the_internal_node[T_SIDE][0][iT4];
				if (iTT4 < 0) {
					iTT4 = neighbors_for_the_internal_node[T_SIDE][1][iT4];
				}
				if (iTT4 < 0) {
					iTT4 = neighbors_for_the_internal_node[T_SIDE][2][iT4];
				}
				if (iTT4 < 0) {
					iTT4 = neighbors_for_the_internal_node[T_SIDE][3][iT4];
				}

				if ((iTT4 >= 0) && (iTT4 < maxelm)) {
					// ���������� ����
					SpeedTT4 = potent[RE_THETA_LANGTRY_MENTER][iTT4];
					//center_cord3D(iTT4,nvtx,pa,pointP,TT_SIDE);
					pointP = center_coord_loc[iTT4];
					positionzTT4 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iTT4 >= maxelm) && (iTT4 < maxelm + maxbound)) {
						SpeedTT4 = potent[RE_THETA_LANGTRY_MENTER][iTT4];
					}
					else {
						SpeedTT4 = SpeedT4;
					}
					//volume3D(iT4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iT4];

					positionzTT4 = positionzT4 + 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iT4 > -1) {
					SpeedT4 = potent[RE_THETA_LANGTRY_MENTER][iT4];
					SpeedTT4 = potent[RE_THETA_LANGTRY_MENTER][iT4]; // �������� ������ ������� ����
				}
				else {
					SpeedT4 = potent[RE_THETA_LANGTRY_MENTER][iT];
					SpeedTT4 = potent[RE_THETA_LANGTRY_MENTER][iT]; // �������� ������ ������� ����
				}
				positionzt4 = positionzP + 0.5 * dz;
				positionzT4 = positionzP + 0.5 * dz;
				positionzTT4 = positionzP + dz; // ����� ���� �� ���������� !
			}

			if ((!bB4) && (iB4 > -1)) {
				SpeedB4 = potent[RE_THETA_LANGTRY_MENTER][iB4];
				//center_cord3D(iB4,nvtx,pa,pointP,B_SIDE);
				pointP = center_coord_loc[iB4];
				positionzB4 = pointP.z;
				positionzb4 = positionzP - 0.5 * dz;

				integer iBB4 = neighbors_for_the_internal_node[B_SIDE][0][iB4];
				if (iBB4 < 0) {
					iBB4 = neighbors_for_the_internal_node[B_SIDE][1][iB4];
				}
				if (iBB4 < 0) {
					iBB4 = neighbors_for_the_internal_node[B_SIDE][2][iB4];
				}
				if (iBB4 < 0) {
					iBB4 = neighbors_for_the_internal_node[B_SIDE][3][iB4];
				}

				if ((iBB4 >= 0) && (iBB4 < maxelm)) {
					// ���������� ����
					SpeedBB4 = potent[RE_THETA_LANGTRY_MENTER][iBB4];
					//center_cord3D(iBB4,nvtx,pa,pointP,BB_SIDE);
					pointP = center_coord_loc[iBB4];
					positionzBB4 = pointP.z;
				}
				else
				{
					// ��������� ����
					if ((iBB4 >= maxelm) && (iBB4 < maxelm + maxbound)) {
						SpeedBB4 = potent[RE_THETA_LANGTRY_MENTER][iBB4];
					}
					else {
						SpeedBB4 = SpeedB4;
					}
					//volume3D(iB4, nvtx, pa, pointP.x, pointP.y, pointP.z);
					pointP = volume_loc[iB4];

					positionzBB4 = positionzB4 - 0.5 * pointP.z;
				}
			}
			else {
				// ��� ��������� ����
				if (iB4 > -1) {
					SpeedB4 = potent[RE_THETA_LANGTRY_MENTER][iB4];
					SpeedBB4 = potent[RE_THETA_LANGTRY_MENTER][iB4]; // �������� ������ ������� ����
				}
				else {
					SpeedB4 = potent[RE_THETA_LANGTRY_MENTER][iB];
					SpeedBB4 = potent[RE_THETA_LANGTRY_MENTER][iB]; // �������� ������ ������� ����
				}
				positionzb4 = positionzP - 0.5 * dz;
				positionzB4 = positionzP - 0.5 * dz;
				positionzBB4 = positionzP - dz; // ����� ���� �� ���������� !
			}

		}



		if ((ishconvection >= QUICK) && (ishconvection <= FROMM)) {
			// ������ ����� ������������ �� ��������� �. ��������� PHOENICS.
		   // ������������� ishconvection ������ ��������� �������� �� ����� QUICK (������������) � ������ �� ����� UNEVENQUICK (�� ������� ���������).

			// X - direction
			Speede = cell_face_value_global(ishconvection, (Fe), SpeedW, SpeedP, SpeedE, SpeedEE);
			Speedw = cell_face_value_global(ishconvection, (Fw), SpeedWW, SpeedW, SpeedP, SpeedE);
			// Y - direction
			Speedn = cell_face_value_global(ishconvection, (Fn), SpeedS, SpeedP, SpeedN, SpeedNN);
			Speeds = cell_face_value_global(ishconvection, (Fs), SpeedSS, SpeedS, SpeedP, SpeedN);
			// Z - direction
			Speedt = cell_face_value_global(ishconvection, (Ft), SpeedB, SpeedP, SpeedT, SpeedTT);
			Speedb = cell_face_value_global(ishconvection, (Fb), SpeedBB, SpeedB, SpeedP, SpeedT);

			if (b_on_adaptive_local_refinement_mesh) {

				// X - direction
				Speede2 = cell_face_value_global(ishconvection, (Fe2), SpeedW2, SpeedP, SpeedE2, SpeedEE2);
				Speedw2 = cell_face_value_global(ishconvection, (Fw2), SpeedWW2, SpeedW2, SpeedP, SpeedE2);
				// Y - direction
				Speedn2 = cell_face_value_global(ishconvection, (Fn2), SpeedS2, SpeedP, SpeedN2, SpeedNN2);
				Speeds2 = cell_face_value_global(ishconvection, (Fs2), SpeedSS2, SpeedS2, SpeedP, SpeedN2);
				// Z - direction
				Speedt2 = cell_face_value_global(ishconvection, (Ft2), SpeedB2, SpeedP, SpeedT2, SpeedTT2);
				Speedb2 = cell_face_value_global(ishconvection, (Fb2), SpeedBB2, SpeedB2, SpeedP, SpeedT2);


				// X - direction
				Speede3 = cell_face_value_global(ishconvection, (Fe3), SpeedW3, SpeedP, SpeedE3, SpeedEE3);
				Speedw3 = cell_face_value_global(ishconvection, (Fw3), SpeedWW3, SpeedW3, SpeedP, SpeedE3);
				// Y - direction
				Speedn3 = cell_face_value_global(ishconvection, (Fn3), SpeedS3, SpeedP, SpeedN3, SpeedNN3);
				Speeds3 = cell_face_value_global(ishconvection, (Fs3), SpeedSS3, SpeedS3, SpeedP, SpeedN3);
				// Z - direction
				Speedt3 = cell_face_value_global(ishconvection, (Ft3), SpeedB3, SpeedP, SpeedT3, SpeedTT3);
				Speedb3 = cell_face_value_global(ishconvection, (Fb3), SpeedBB3, SpeedB3, SpeedP, SpeedT3);


				// X - direction
				Speede4 = cell_face_value_global(ishconvection, (Fe4), SpeedW4, SpeedP, SpeedE4, SpeedEE4);
				Speedw4 = cell_face_value_global(ishconvection, (Fw4), SpeedWW4, SpeedW4, SpeedP, SpeedE4);
				// Y - direction
				Speedn4 = cell_face_value_global(ishconvection, (Fn4), SpeedS4, SpeedP, SpeedN4, SpeedNN4);
				Speeds4 = cell_face_value_global(ishconvection, (Fs4), SpeedSS4, SpeedS4, SpeedP, SpeedN4);
				// Z - direction
				Speedt4 = cell_face_value_global(ishconvection, (Ft4), SpeedB4, SpeedP, SpeedT4, SpeedTT4);
				Speedb4 = cell_face_value_global(ishconvection, (Fb4), SpeedBB4, SpeedB4, SpeedP, SpeedT4);

			}

		}

		if (ishconvection >= UNEVENQUICK) {


			/*
			// ������������������ �������� ��������� � ����� ���������� ���������� ����� QUICK �� ������������� �����.
			// ���������� ���� ������������ �� ������: ...
			// � ������ ������ ������ ���������� �� ������������.
			//doublereal gamma1E, gamma2E, gamma1W, gamma2W, delta1E, delta2E, delta1W, delta2W;
			//doublereal gamma1N, gamma2N, gamma1S, gamma2S, delta1N, delta2N, delta1S, delta2S;
			//doublereal gamma1T, gamma2T, gamma1B, gamma2B, delta1T, delta2T, delta1B, delta2B;
			// X - direction
			// gamma
			//gamma1E=((positionxe-positionxE)*(positionxe-positionxP))/((positionxW-positionxE)*(positionxW-positionxP));
			//gamma2E=((positionxe-positionxP)*(positionxe-positionxW))/((positionxE-positionxP)*(positionxE-positionxW));
			//gamma1W=((positionxw-positionxP)*(positionxw-positionxW))/((positionxWW-positionxP)*(positionxWW-positionxW));
			//gamma2W=((positionxw-positionxW)*(positionxw-positionxWW))/((positionxP-positionxW)*(positionxP-positionxWW));
			// delta
			//delta1E=((positionxe-positionxEE)*(positionxe-positionxE))/((positionxP-positionxEE)*(positionxP-positionxE));
			//delta2E=((positionxe-positionxE)*(positionxe-positionxP))/((positionxEE-positionxE)*(positionxEE-positionxP));
			//delta1W=((positionxw-positionxE)*(positionxw-positionxP))/((positionxW-positionxE)*(positionxW-positionxP));
			//delta2W=((positionxw-positionxP)*(positionxw-positionxW))/((positionxE-positionxP)*(positionxE-positionxW));
			// Y - direction
			// gamma
			//gamma1N=((positionyn-positionyN)*(positionyn-positionyP))/((positionyS-positionyN)*(positionyS-positionyP));
			//gamma2N=((positionyn-positionyP)*(positionyn-positionyS))/((positionyN-positionyP)*(positionyN-positionyS));
			//gamma1S=((positionys-positionyP)*(positionys-positionyS))/((positionySS-positionyP)*(positionySS-positionyS));
			//gamma2S=((positionys-positionyS)*(positionys-positionySS))/((positionyP-positionyS)*(positionyP-positionySS));
			// delta
			//delta1N=((positionyn-positionyNN)*(positionyn-positionyN))/((positionyP-positionyNN)*(positionyP-positionyN));
			//delta2N=((positionyn-positionyN)*(positionyn-positionyP))/((positionyNN-positionyN)*(positionyNN-positionyP));
			//delta1S=((positionys-positionyN)*(positionys-positionyP))/((positionyS-positionyN)*(positionyS-positionyP));
			//delta2S=((positionys-positionyP)*(positionys-positionyS))/((positionyN-positionyP)*(positionyN-positionyS));
			// Z - direction
			// gamma
			//gamma1T=((positionzt-positionzT)*(positionzt-positionzP))/((positionzB-positionzT)*(positionzB-positionzP));
			//gamma2T=((positionzt-positionzP)*(positionzt-positionzB))/((positionzT-positionzP)*(positionzT-positionzB));
			//gamma1B=((positionzb-positionzP)*(positionzb-positionzB))/((positionzBB-positionzP)*(positionzBB-positionzB));
			//gamma2B=((positionzb-positionzB)*(positionzb-positionzBB))/((positionzP-positionzB)*(positionzP-positionzBB));
			// delta
			//delta1T=((positionzt-positionzTT)*(positionzt-positionzT))/((positionzP-positionzTT)*(positionzP-positionzT));
			//delta2T=((positionzt-positionzT)*(positionzt-positionzP))/((positionzTT-positionzT)*(positionzTT-positionzP));
			//delta1B=((positionzb-positionzT)*(positionzb-positionzP))/((positionzB-positionzT)*(positionzB-positionzP));
			//delta2B=((positionzb-positionzP)*(positionzb-positionzB))/((positionzT-positionzP)*(positionzT-positionzB));
			*/




			// ���������� ������� �������� �� ����� ��
			// ������������ ����� �������� QUICK.
			/* ������� ������������:
			 *  A	B	C	D	e	+/-
			 *  W	P	E	-	e	+
			 *  -	P	E	EE  e   -
			 *  WW   W	P	-	w	+
			 *  -	W	P	E	w	-
			 *  S	P	N	-	n	+
			 *  -	P	N	NN  n	-
			 *  SS   S	P	-	s	+
			 *  -	S	P	N	s	-
			 *  B	P	T	-	t	+
			 *  -	P	T	TT  t	-
			 *  BB   B	P	-	b	+
			 *  -	B	P	T	b	-
			 */

			if (ishconvection == UNEVENQUICK) {
				// X - direction
				Speede = workQUICK(dx, 2.0 * (positionxE - positionxe), positionxW, positionxP, positionxE, positionxEE, SpeedW, SpeedP, SpeedE, SpeedEE, (Fe));
				Speedw = workQUICK(2.0 * (positionxw - positionxW), dx, positionxWW, positionxW, positionxP, positionxE, SpeedWW, SpeedW, SpeedP, SpeedE, (Fw));
				// Y - direction
				Speedn = workQUICK(dy, 2.0 * (positionyN - positionyn), positionyS, positionyP, positionyN, positionyNN, SpeedS, SpeedP, SpeedN, SpeedNN, (Fn));
				Speeds = workQUICK(2.0 * (positionys - positionyS), dy, positionySS, positionyS, positionyP, positionyN, SpeedSS, SpeedS, SpeedP, SpeedN, (Fs));
				// Z - direction
				Speedt = workQUICK(dz, 2.0 * (positionzT - positionzt), positionzB, positionzP, positionzT, positionzTT, SpeedB, SpeedP, SpeedT, SpeedTT, (Ft));
				Speedb = workQUICK(2.0 * (positionzb - positionzB), dz, positionzBB, positionzB, positionzP, positionzT, SpeedBB, SpeedB, SpeedP, SpeedT, (Fb));


				if (b_on_adaptive_local_refinement_mesh) {

					// X - direction
					Speede2 = workQUICK(dx, 2.0 * (positionxE2 - positionxe2), positionxW2, positionxP, positionxE2, positionxEE2, SpeedW2, SpeedP, SpeedE2, SpeedEE2, (Fe2));
					Speedw2 = workQUICK(2.0 * (positionxw2 - positionxW2), dx, positionxWW2, positionxW2, positionxP, positionxE2, SpeedWW2, SpeedW2, SpeedP, SpeedE2, (Fw2));
					// Y - direction
					Speedn2 = workQUICK(dy, 2.0 * (positionyN2 - positionyn2), positionyS2, positionyP, positionyN2, positionyNN2, SpeedS2, SpeedP, SpeedN2, SpeedNN2, (Fn2));
					Speeds2 = workQUICK(2.0 * (positionys2 - positionyS2), dy, positionySS2, positionyS2, positionyP, positionyN2, SpeedSS2, SpeedS2, SpeedP, SpeedN2, (Fs2));
					// Z - direction
					Speedt2 = workQUICK(dz, 2.0 * (positionzT2 - positionzt2), positionzB2, positionzP, positionzT2, positionzTT2, SpeedB2, SpeedP, SpeedT2, SpeedTT2, (Ft2));
					Speedb2 = workQUICK(2.0 * (positionzb2 - positionzB2), dz, positionzBB2, positionzB2, positionzP, positionzT2, SpeedBB2, SpeedB2, SpeedP, SpeedT2, (Fb2));


					// X - direction
					Speede3 = workQUICK(dx, 2.0 * (positionxE3 - positionxe3), positionxW3, positionxP, positionxE3, positionxEE3, SpeedW3, SpeedP, SpeedE3, SpeedEE3, (Fe3));
					Speedw3 = workQUICK(2.0 * (positionxw3 - positionxW3), dx, positionxWW3, positionxW3, positionxP, positionxE3, SpeedWW3, SpeedW3, SpeedP, SpeedE3, (Fw3));
					// Y - direction
					Speedn3 = workQUICK(dy, 2.0 * (positionyN3 - positionyn3), positionyS3, positionyP, positionyN3, positionyNN3, SpeedS3, SpeedP, SpeedN3, SpeedNN3, (Fn3));
					Speeds3 = workQUICK(2.0 * (positionys3 - positionyS3), dy, positionySS3, positionyS3, positionyP, positionyN3, SpeedSS3, SpeedS3, SpeedP, SpeedN3, (Fs3));
					// Z - direction
					Speedt3 = workQUICK(dz, 2.0 * (positionzT3 - positionzt3), positionzB3, positionzP, positionzT3, positionzTT3, SpeedB3, SpeedP, SpeedT3, SpeedTT3, (Ft3));
					Speedb3 = workQUICK(2.0 * (positionzb3 - positionzB3), dz, positionzBB3, positionzB3, positionzP, positionzT3, SpeedBB3, SpeedB3, SpeedP, SpeedT3, (Fb3));


					// X - direction
					Speede4 = workQUICK(dx, 2.0 * (positionxE4 - positionxe4), positionxW4, positionxP, positionxE4, positionxEE4, SpeedW4, SpeedP, SpeedE4, SpeedEE4, (Fe4));
					Speedw4 = workQUICK(2.0 * (positionxw4 - positionxW4), dx, positionxWW4, positionxW4, positionxP, positionxE4, SpeedWW4, SpeedW4, SpeedP, SpeedE4, (Fw4));
					// Y - direction
					Speedn4 = workQUICK(dy, 2.0 * (positionyN4 - positionyn4), positionyS4, positionyP, positionyN4, positionyNN4, SpeedS4, SpeedP, SpeedN4, SpeedNN4, (Fn4));
					Speeds4 = workQUICK(2.0 * (positionys4 - positionyS4), dy, positionySS4, positionyS4, positionyP, positionyN4, SpeedSS4, SpeedS4, SpeedP, SpeedN4, (Fs4));
					// Z - direction
					Speedt4 = workQUICK(dz, 2.0 * (positionzT4 - positionzt4), positionzB4, positionzP, positionzT4, positionzTT4, SpeedB4, SpeedP, SpeedT4, SpeedTT4, (Ft4));
					Speedb4 = workQUICK(2.0 * (positionzb4 - positionzB4), dz, positionzBB4, positionzB4, positionzP, positionzT4, SpeedBB4, SpeedB4, SpeedP, SpeedT4, (Fb4));

				}

			}

			if ((ishconvection > UNEVENQUICK) && (ishconvection <= UNEVEN_CUBISTA)) {
				// ���� �� ������ ������ ������������� ����������� ������������ ������ ������ ������ �����:
				// 1. UNEVEN_MUSCL, 2. UNEVEN_SOUCUP, 3. UNEVEN_HLPA, 4. UNEVEN_SMART.
				// ������������� ����� ������ ��������������� ��������.

				// X - direction
				Speede = workKN_VOLKOV(positionxW, positionxP, positionxE, positionxEE, SpeedW, SpeedP, SpeedE, SpeedEE, (Fe), ishconvection);
				Speedw = workKN_VOLKOV(positionxWW, positionxW, positionxP, positionxE, SpeedWW, SpeedW, SpeedP, SpeedE, (Fw), ishconvection);
				// Y - direction
				Speedn = workKN_VOLKOV(positionyS, positionyP, positionyN, positionyNN, SpeedS, SpeedP, SpeedN, SpeedNN, (Fn), ishconvection);
				Speeds = workKN_VOLKOV(positionySS, positionyS, positionyP, positionyN, SpeedSS, SpeedS, SpeedP, SpeedN, (Fs), ishconvection);
				// Z - direction
				Speedt = workKN_VOLKOV(positionzB, positionzP, positionzT, positionzTT, SpeedB, SpeedP, SpeedT, SpeedTT, (Ft), ishconvection);
				Speedb = workKN_VOLKOV(positionzBB, positionzB, positionzP, positionzT, SpeedBB, SpeedB, SpeedP, SpeedT, (Fb), ishconvection);

				// debug ������ �������� ������.
				//printf("%f, %f, %f, %f, %f, %f\n",Speede,Speedw,Speedn,Speeds,Speedt,Speedb);
				//getchar();


				if (b_on_adaptive_local_refinement_mesh) {

					// X - direction
					Speede2 = workKN_VOLKOV(positionxW2, positionxP, positionxE2, positionxEE2, SpeedW2, SpeedP, SpeedE2, SpeedEE2, (Fe2), ishconvection);
					Speedw2 = workKN_VOLKOV(positionxWW2, positionxW2, positionxP, positionxE2, SpeedWW2, SpeedW2, SpeedP, SpeedE2, (Fw2), ishconvection);
					// Y - direction
					Speedn2 = workKN_VOLKOV(positionyS2, positionyP, positionyN2, positionyNN2, SpeedS2, SpeedP, SpeedN2, SpeedNN2, (Fn2), ishconvection);
					Speeds2 = workKN_VOLKOV(positionySS2, positionyS2, positionyP, positionyN2, SpeedSS2, SpeedS2, SpeedP, SpeedN2, (Fs2), ishconvection);
					// Z - direction
					Speedt2 = workKN_VOLKOV(positionzB2, positionzP, positionzT2, positionzTT2, SpeedB2, SpeedP, SpeedT2, SpeedTT2, (Ft2), ishconvection);
					Speedb2 = workKN_VOLKOV(positionzBB2, positionzB2, positionzP, positionzT2, SpeedBB2, SpeedB2, SpeedP, SpeedT2, (Fb2), ishconvection);


					// X - direction
					Speede3 = workKN_VOLKOV(positionxW3, positionxP, positionxE3, positionxEE3, SpeedW3, SpeedP, SpeedE3, SpeedEE3, (Fe3), ishconvection);
					Speedw3 = workKN_VOLKOV(positionxWW3, positionxW3, positionxP, positionxE3, SpeedWW3, SpeedW3, SpeedP, SpeedE3, (Fw3), ishconvection);
					// Y - direction
					Speedn3 = workKN_VOLKOV(positionyS3, positionyP, positionyN3, positionyNN3, SpeedS3, SpeedP, SpeedN3, SpeedNN3, (Fn3), ishconvection);
					Speeds3 = workKN_VOLKOV(positionySS3, positionyS3, positionyP, positionyN3, SpeedSS3, SpeedS3, SpeedP, SpeedN3, (Fs3), ishconvection);
					// Z - direction
					Speedt3 = workKN_VOLKOV(positionzB3, positionzP, positionzT3, positionzTT3, SpeedB3, SpeedP, SpeedT3, SpeedTT3, (Ft3), ishconvection);
					Speedb3 = workKN_VOLKOV(positionzBB3, positionzB3, positionzP, positionzT3, SpeedBB3, SpeedB3, SpeedP, SpeedT3, (Fb3), ishconvection);


					// X - direction
					Speede4 = workKN_VOLKOV(positionxW4, positionxP, positionxE4, positionxEE4, SpeedW4, SpeedP, SpeedE4, SpeedEE4, (Fe4), ishconvection);
					Speedw4 = workKN_VOLKOV(positionxWW4, positionxW4, positionxP, positionxE4, SpeedWW4, SpeedW4, SpeedP, SpeedE4, (Fw4), ishconvection);
					// Y - direction
					Speedn4 = workKN_VOLKOV(positionyS4, positionyP, positionyN4, positionyNN4, SpeedS4, SpeedP, SpeedN4, SpeedNN4, (Fn4), ishconvection);
					Speeds4 = workKN_VOLKOV(positionySS4, positionyS4, positionyP, positionyN4, SpeedSS4, SpeedS4, SpeedP, SpeedN4, (Fs4), ishconvection);
					// Z - direction
					Speedt4 = workKN_VOLKOV(positionzB4, positionzP, positionzT4, positionzTT4, SpeedB4, SpeedP, SpeedT4, SpeedTT4, (Ft4), ishconvection);
					Speedb4 = workKN_VOLKOV(positionzBB4, positionzB4, positionzP, positionzT4, SpeedBB4, SpeedB4, SpeedP, SpeedT4, (Fb4), ishconvection);

				}

			}

		} // endif (ishconvection >= UNEVENQUICK)



		// ������: SIMPLE method for the solution of incompressible flows on non-staggered grids
		// I. Sezai - Eastern Mediterranean University, Mechanical Engineering Department, Mersin 10-Turkey Revised in January, 2011.

		// ���������� ������������� ����������� �������:
		// ����������� ����� ���������� ���������:
		// ������ ����������� ������ ��������������� �����, 
		// � ���������� ����� ������������ � ������ ����� 
		// �������� ������� ���������.
		if (1) {

			/*
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae=De*fD(Pe, EXP2, true, feplus) + fmax(-(Fe),0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw=Dw*fD(Pw, EXP2, true, fwplus) + fmax((Fw),0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an=Dn*fD(Pn, EXP2, true, fnplus) + fmax(-(Fn),0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as=Ds*fD(Ps, EXP2, true, fsplus) + fmax((Fs),0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at=Dt*fD(Pt, EXP2, true, ftplus) + fmax(-(Ft),0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab=Db*fD(Pb, EXP2, true, fbplus) + fmax((Fb),0.0);
			*/

			// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
			// ����� ������ UDS.
			// ��� ����������� � ����������.
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De + fmax(-(Fe), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw + fmax((Fw), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn + fmax(-(Fn), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds + fmax((Fs), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt + fmax(-(Ft), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db + fmax((Fb), 0.0);


			// 08.05.2017.
			// ��� ���������:
			// ����� ������������� !!! (���������� �� ����).
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = De + fmax(+(Fe), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw + fmax(-(Fw), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn + fmax(+(Fn), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds + fmax(-(Fs), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt + fmax(+(Ft), 0.0);
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db + fmax(-(Fb), 0.0);

			/*
			sumanb = De + fmax(+(Fe), 0.0);
			sumanb += Dw + fmax(-(Fw), 0.0);
			sumanb += Dn + fmax(+(Fn), 0.0);
			sumanb += Ds + fmax(-(Fs), 0.0);
			sumanb += Dt + fmax(+(Ft), 0.0);
			sumanb += Db + fmax(-(Fb), 0.0);
			*/


			if (b_on_adaptive_local_refinement_mesh) {

				// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
			 // ����� ������ UDS.
			 // ��� ����������� � ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fmax(-(Fe2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fmax((Fw2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fmax(-(Fn2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fmax((Fs2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fmax(-(Ft2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fmax((Fb2), 0.0);


				// 08.05.2017.
				// ��� ���������:
				// ����� ������������� !!! (���������� �� ����).
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += De2 + fmax(+(Fe2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw2 + fmax(-(Fw2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn2 + fmax(+(Fn2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds2 + fmax(-(Fs2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt2 + fmax(+(Ft2), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db2 + fmax(-(Fb2), 0.0);

				/*
				sumanb += De2 + fmax(+(Fe2), 0.0);
				sumanb += Dw2 + fmax(-(Fw2), 0.0);
				sumanb += Dn2 + fmax(+(Fn2), 0.0);
				sumanb += Ds2 + fmax(-(Fs2), 0.0);
				sumanb += Dt2 + fmax(+(Ft2), 0.0);
				sumanb += Db2 + fmax(-(Fb2), 0.0);
				*/


				// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
		   // ����� ������ UDS.
		   // ��� ����������� � ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fmax(-(Fe3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fmax((Fw3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fmax(-(Fn3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fmax((Fs3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fmax(-(Ft3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fmax((Fb3), 0.0);


				// 08.05.2017.
				// ��� ���������:
				// ����� ������������� !!! (���������� �� ����).
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += De3 + fmax(+(Fe3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw3 + fmax(-(Fw3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn3 + fmax(+(Fn3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds3 + fmax(-(Fs3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt3 + fmax(+(Ft3), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db3 + fmax(-(Fb3), 0.0);

				/*
				sumanb += De3 + fmax(+(Fe3), 0.0);
				sumanb += Dw3 + fmax(-(Fw3), 0.0);
				sumanb += Dn3 + fmax(+(Fn3), 0.0);
				sumanb += Ds3 + fmax(-(Fs3), 0.0);
				sumanb += Dt3 + fmax(+(Ft3), 0.0);
				sumanb += Db3 + fmax(-(Fb3), 0.0);
				*/


				// ������� ��� ����������� ������ � ������������� � ���������� 7.05.2017.
		   // ����� ������ UDS.
		   // ��� ����������� � ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fmax(-(Fe4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fmax((Fw4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fmax(-(Fn4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fmax((Fs4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fmax(-(Ft4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fmax((Fb4), 0.0);


				// 08.05.2017.
				// ��� ���������:
				// ����� ������������� !!! (���������� �� ����).
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += De4 + fmax(+(Fe4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dw4 + fmax(-(Fw4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dn4 + fmax(+(Fn4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Ds4 + fmax(-(Fs4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Dt4 + fmax(+(Ft4), 0.0);
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += Db4 + fmax(-(Fb4), 0.0);

				/*
				sumanb += De4 + fmax(+(Fe4), 0.0);
				sumanb += Dw4 + fmax(-(Fw4), 0.0);
				sumanb += Dn4 + fmax(+(Fn4), 0.0);
				sumanb += Ds4 + fmax(-(Fs4), 0.0);
				sumanb += Dt4 + fmax(+(Ft4), 0.0);
				sumanb += Db4 + fmax(-(Fb4), 0.0);
				*/

			}

		}
		else {

			// 8.05.2017
			// �� � ���� ������ �� �������� ��� ����� ����.


			// ��� ������ ������ �� ��������, ������� ���� � ���������� � ���� ������� ����������,
			// �� ��� ����� ���������� ������������ ��� �� ������������� ����������� ������.
		   // 30 07 2015
		   // 
		   // ������ ������ ������� ����� ���������� �� UDS.
			if (!bE) {
				// ������ ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De + fmax(-(Fe), 0.0);
			}
			else {
				integer inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De + fabs(Fe);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae = De + fabs(Fe);
				}
			}

			if (!bW) {
				// ������ ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw + fmax((Fw), 0.0);
			}
			else {
				integer inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw + fabs(Fw);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw = Dw + fabs(Fw);
				}
			}

			if (!bN) {
				// ������ ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn + fmax(-(Fn), 0.0);
			}
			else {
				integer inumber = iN - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn + fabs(Fn);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an = Dn + fabs(Fn);
				}
			}

			if (!bS) {
				// ������ ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds + fmax((Fs), 0.0);
			}
			else {
				integer inumber = iS - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds + fabs(Fs);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as = Ds + fabs(Fs);
				}
			}

			if (!bT) {
				// ������ ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt + fmax(-(Ft), 0.0);
			}
			else {
				integer inumber = iT - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt + fabs(Ft);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at = Dt + fabs(Ft);
				}
			}

			if (!bB) {
				// ������ ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db + fmax((Fb), 0.0);
			}
			else {
				integer inumber = iB - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					// ������� �� ���������: ������ ������.
					// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db + fabs(Fb);
				}
				else {
					// �� ���� ��������� ������� ����� ������� ������� �� �������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab = Db + fabs(Fb);
				}
			}


			if (b_on_adaptive_local_refinement_mesh) {

				// ������ ������ ������� ����� ���������� �� UDS.
				if (!bE2) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fmax(-(Fe2), 0.0);
				}
				else {
					integer inumber = iE2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fabs(Fe2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 = De2 + fabs(Fe2);
					}
				}

				if (!bW2) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fmax((Fw2), 0.0);
				}
				else {
					integer inumber = iW2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fabs(Fw2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 = Dw2 + fabs(Fw2);
					}
				}

				if (!bN2) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fmax(-(Fn2), 0.0);
				}
				else {
					integer inumber = iN2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fabs(Fn2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 = Dn2 + fabs(Fn2);
					}
				}

				if (!bS2) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fmax((Fs2), 0.0);
				}
				else {
					integer inumber = iS2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fabs(Fs2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 = Ds2 + fabs(Fs2);
					}
				}

				if (!bT2) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fmax(-(Ft2), 0.0);
				}
				else {
					integer inumber = iT2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fabs(Ft2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 = Dt2 + fabs(Ft2);
					}
				}

				if (!bB2) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fmax((Fb2), 0.0);
				}
				else {
					integer inumber = iB2 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fabs(Fb2);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 = Db2 + fabs(Fb2);
					}
				}


				// ������ ������ ������� ����� ���������� �� UDS.
				if (!bE3) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fmax(-(Fe3), 0.0);
				}
				else {
					integer inumber = iE3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fabs(Fe3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 = De3 + fabs(Fe3);
					}
				}

				if (!bW3) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fmax((Fw3), 0.0);
				}
				else {
					integer inumber = iW3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fabs(Fw3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 = Dw3 + fabs(Fw3);
					}
				}

				if (!bN3) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fmax(-(Fn3), 0.0);
				}
				else {
					integer inumber = iN3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fabs(Fn3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 = Dn3 + fabs(Fn3);
					}
				}

				if (!bS3) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fmax((Fs3), 0.0);
				}
				else {
					integer inumber = iS3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fabs(Fs3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 = Ds3 + fabs(Fs3);
					}
				}

				if (!bT3) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fmax(-(Ft3), 0.0);
				}
				else {
					integer inumber = iT3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fabs(Ft3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 = Dt3 + fabs(Ft3);
					}
				}

				if (!bB3) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fmax((Fb3), 0.0);
				}
				else {
					integer inumber = iB3 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fabs(Fb3);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 = Db3 + fabs(Fb3);
					}
				}


				// ������ ������ ������� ����� ���������� �� UDS.
				if (!bE4) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fmax(-(Fe4), 0.0);
				}
				else {
					integer inumber = iE4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fabs(Fe4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 = De4 + fabs(Fe4);
					}
				}

				if (!bW4) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fmax((Fw4), 0.0);
				}
				else {
					integer inumber = iW4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fabs(Fw4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 = Dw4 + fabs(Fw4);
					}
				}

				if (!bN4) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fmax(-(Fn4), 0.0);
				}
				else {
					integer inumber = iN4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fabs(Fn4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 = Dn4 + fabs(Fn4);
					}
				}

				if (!bS4) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fmax((Fs4), 0.0);
				}
				else {
					integer inumber = iS4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fabs(Fs4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 = Ds4 + fabs(Fs4);
					}
				}

				if (!bT4) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fmax(-(Ft4), 0.0);
				}
				else {
					integer inumber = iT4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fabs(Ft4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 = Dt4 + fabs(Ft4);
					}
				}

				if (!bB4) {
					// ������ ����������.
					sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fmax((Fb4), 0.0);
				}
				else {
					integer inumber = iB4 - maxelm;
					if (border_neighbor[inumber].MCB == (ls + lw)) {
						// ������� �� ���������: ������ ������.
						// �������� ������� ���� �� �������, ��� �� ����� ������� ������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fabs(Fb4);
					}
					else {
						// �� ���� ��������� ������� ����� ������� ������� �� �������.
						sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 = Db4 + fabs(Fb4);
					}
				}



			}


			// 7.05.2017 ������� ��� ����������� ������ � ��������������� � ����������.
			sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;

			//sumanb = sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;


			if (b_on_adaptive_local_refinement_mesh) {

				// 7.05.2017 ������� ��� ����������� ������ � ��������������� � ����������.
				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2;

				//sumanb += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2;

				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3;

				//sumanb += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3;

				sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4;

				//sumanb += sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 + sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4;


			}

		}


		// 7.05.2017 ������� ��� ����������� ������ � ��������������� � ����������.
		//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;

		//sumanb=sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae+sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw+sl[RE_THETA_LANGTRY_MENTER_SL][iP].an+sl[RE_THETA_LANGTRY_MENTER_SL][iP].as+sl[RE_THETA_LANGTRY_MENTER_SL][iP].at+sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab;

		//13 ������� 2016.
		//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);
		//sumanb = fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].an) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].as) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].at) + fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab);


		if (1) {
			// ����� � ������ ����� (����� ���������� ���������):
			// X - direction
			attrs += -fmax((Fe), 0) * (Speede - SpeedP) + fmax(-(Fe), 0) * (Speede - SpeedE);
			attrs += -fmax(-(Fw), 0) * (Speedw - SpeedP) + fmax((Fw), 0) * (Speedw - SpeedW);
			// Y - direction
			attrs += -fmax((Fn), 0) * (Speedn - SpeedP) + fmax(-(Fn), 0) * (Speedn - SpeedN);
			attrs += -fmax(-(Fs), 0) * (Speeds - SpeedP) + fmax((Fs), 0) * (Speeds - SpeedS);
			// Z - direction
			attrs += -fmax((Ft), 0) * (Speedt - SpeedP) + fmax(-(Ft), 0) * (Speedt - SpeedT);
			attrs += -fmax(-(Fb), 0) * (Speedb - SpeedP) + fmax((Fb), 0) * (Speedb - SpeedB);


			if (b_on_adaptive_local_refinement_mesh) {

				// ����� � ������ ����� (����� ���������� ���������):
				// X - direction
				attrs += -fmax((Fe2), 0) * (Speede2 - SpeedP) + fmax(-(Fe2), 0) * (Speede2 - SpeedE2);
				attrs += -fmax(-(Fw2), 0) * (Speedw2 - SpeedP) + fmax((Fw2), 0) * (Speedw2 - SpeedW2);
				// Y - direction
				attrs += -fmax((Fn2), 0) * (Speedn2 - SpeedP) + fmax(-(Fn2), 0) * (Speedn2 - SpeedN2);
				attrs += -fmax(-(Fs2), 0) * (Speeds2 - SpeedP) + fmax((Fs2), 0) * (Speeds2 - SpeedS2);
				// Z - direction
				attrs += -fmax((Ft2), 0) * (Speedt2 - SpeedP) + fmax(-(Ft2), 0) * (Speedt2 - SpeedT2);
				attrs += -fmax(-(Fb2), 0) * (Speedb2 - SpeedP) + fmax((Fb2), 0) * (Speedb2 - SpeedB2);


				// ����� � ������ ����� (����� ���������� ���������):
				// X - direction
				attrs += -fmax((Fe3), 0) * (Speede3 - SpeedP) + fmax(-(Fe3), 0) * (Speede3 - SpeedE3);
				attrs += -fmax(-(Fw3), 0) * (Speedw3 - SpeedP) + fmax((Fw3), 0) * (Speedw3 - SpeedW3);
				// Y - direction
				attrs += -fmax((Fn3), 0) * (Speedn3 - SpeedP) + fmax(-(Fn3), 0) * (Speedn3 - SpeedN3);
				attrs += -fmax(-(Fs3), 0) * (Speeds3 - SpeedP) + fmax((Fs3), 0) * (Speeds3 - SpeedS3);
				// Z - direction
				attrs += -fmax((Ft3), 0) * (Speedt3 - SpeedP) + fmax(-(Ft3), 0) * (Speedt3 - SpeedT3);
				attrs += -fmax(-(Fb3), 0) * (Speedb3 - SpeedP) + fmax((Fb3), 0) * (Speedb3 - SpeedB3);


				// ����� � ������ ����� (����� ���������� ���������):
				// X - direction
				attrs += -fmax((Fe4), 0) * (Speede4 - SpeedP) + fmax(-(Fe4), 0) * (Speede4 - SpeedE4);
				attrs += -fmax(-(Fw4), 0) * (Speedw4 - SpeedP) + fmax((Fw4), 0) * (Speedw4 - SpeedW4);
				// Y - direction
				attrs += -fmax((Fn4), 0) * (Speedn4 - SpeedP) + fmax(-(Fn4), 0) * (Speedn4 - SpeedN4);
				attrs += -fmax(-(Fs4), 0) * (Speeds4 - SpeedP) + fmax((Fs4), 0) * (Speeds4 - SpeedS4);
				// Z - direction
				attrs += -fmax((Ft4), 0) * (Speedt4 - SpeedP) + fmax(-(Ft4), 0) * (Speedt4 - SpeedT4);
				attrs += -fmax(-(Fb4), 0) * (Speedb4 - SpeedP) + fmax((Fb4), 0) * (Speedb4 - SpeedB4);

			}

		}
		else {
			// �������.
		   // 30 07 2015
		   // 
		   // ������ ������ ������� ����� ���������� �� UDS.
			if (!bE) {
				attrs += -fmax((Fe), 0) * (Speede - SpeedP) + fmax(-(Fe), 0) * (Speede - SpeedE);
			}
			if (!bW) {
				attrs += -fmax(-(Fw), 0) * (Speedw - SpeedP) + fmax((Fw), 0) * (Speedw - SpeedW);
			}
			if (!bN) {
				attrs += -fmax((Fn), 0) * (Speedn - SpeedP) + fmax(-(Fn), 0) * (Speedn - SpeedN);
			}
			if (!bS) {
				attrs += -fmax(-(Fs), 0) * (Speeds - SpeedP) + fmax((Fs), 0) * (Speeds - SpeedS);
			}
			if (!bT) {
				attrs += -fmax((Ft), 0) * (Speedt - SpeedP) + fmax(-(Ft), 0) * (Speedt - SpeedT);
			}
			if (!bB) {
				attrs += -fmax(-(Fb), 0) * (Speedb - SpeedP) + fmax((Fb), 0) * (Speedb - SpeedB);
			}


			if (b_on_adaptive_local_refinement_mesh) {

				if (!bE2) {
					attrs += -fmax((Fe2), 0) * (Speede2 - SpeedP) + fmax(-(Fe2), 0) * (Speede2 - SpeedE2);
				}
				if (!bW2) {
					attrs += -fmax(-(Fw2), 0) * (Speedw2 - SpeedP) + fmax((Fw2), 0) * (Speedw2 - SpeedW2);
				}
				if (!bN2) {
					attrs += -fmax((Fn2), 0) * (Speedn2 - SpeedP) + fmax(-(Fn2), 0) * (Speedn2 - SpeedN2);
				}
				if (!bS2) {
					attrs += -fmax(-(Fs2), 0) * (Speeds2 - SpeedP) + fmax((Fs2), 0) * (Speeds2 - SpeedS2);
				}
				if (!bT2) {
					attrs += -fmax((Ft2), 0) * (Speedt2 - SpeedP) + fmax(-(Ft2), 0) * (Speedt2 - SpeedT2);
				}
				if (!bB2) {
					attrs += -fmax(-(Fb2), 0) * (Speedb2 - SpeedP) + fmax((Fb2), 0) * (Speedb2 - SpeedB2);
				}


				if (!bE3) {
					attrs += -fmax((Fe3), 0) * (Speede3 - SpeedP) + fmax(-(Fe3), 0) * (Speede3 - SpeedE3);
				}
				if (!bW3) {
					attrs += -fmax(-(Fw3), 0) * (Speedw3 - SpeedP) + fmax((Fw3), 0) * (Speedw3 - SpeedW3);
				}
				if (!bN3) {
					attrs += -fmax((Fn3), 0) * (Speedn3 - SpeedP) + fmax(-(Fn3), 0) * (Speedn3 - SpeedN3);
				}
				if (!bS3) {
					attrs += -fmax(-(Fs3), 0) * (Speeds3 - SpeedP) + fmax((Fs3), 0) * (Speeds3 - SpeedS3);
				}
				if (!bT3) {
					attrs += -fmax((Ft3), 0) * (Speedt3 - SpeedP) + fmax(-(Ft3), 0) * (Speedt3 - SpeedT3);
				}
				if (!bB3) {
					attrs += -fmax(-(Fb3), 0) * (Speedb3 - SpeedP) + fmax((Fb3), 0) * (Speedb3 - SpeedB3);
				}


				if (!bE4) {
					attrs += -fmax((Fe4), 0) * (Speede4 - SpeedP) + fmax(-(Fe4), 0) * (Speede4 - SpeedE4);
				}
				if (!bW4) {
					attrs += -fmax(-(Fw4), 0) * (Speedw4 - SpeedP) + fmax((Fw4), 0) * (Speedw4 - SpeedW4);
				}
				if (!bN4) {
					attrs += -fmax((Fn4), 0) * (Speedn4 - SpeedP) + fmax(-(Fn4), 0) * (Speedn4 - SpeedN4);
				}
				if (!bS4) {
					attrs += -fmax(-(Fs4), 0) * (Speeds4 - SpeedP) + fmax((Fs4), 0) * (Speeds4 - SpeedS4);
				}
				if (!bT4) {
					attrs += -fmax((Ft4), 0) * (Speedt4 - SpeedP) + fmax(-(Ft4), 0) * (Speedt4 - SpeedT4);
				}
				if (!bB4) {
					attrs += -fmax(-(Fb4), 0) * (Speedb4 - SpeedP) + fmax((Fb4), 0) * (Speedb4 - SpeedB4);
				}

			}


		}

		//attrs=0.0; // ����� ����� ������� ����������� ����������� (�������� ����� ��������).

	}

	// ��������� ����������� ��������� ���������������.
	// ������ �������. �������������� ���������� �� ������������.
	// 

	// ������������ ����
	doublereal dSc = 0.0, dSp = 0.0;

	// ��� ������������� �������� Fe, Fw, Fn, Fs, Ft, Fb - ���������� �������� ������ ����� ���������� �������������� ���-���.
		// 02.05.2017
		// ��� ������� �.�. �������� � ������������� ������������ �������������.
		// ������ ���� �������: deltaF=(Fe-Fw+Fn-Fs+Ft-Fb);
		// ����������� ����� ����������� � ���������� �� ANSYS Icepak.
		// ��� ����������� �� ���� �������� ����� �� ���������� ���� - ��� ���������� ������ ��������.
		//doublereal deltaF=(Fe-Fw+Fn-Fs+Ft-Fb);
		// ��� ������ ���������� ��������� ������������� ��� ��������� ����� ����.
		// � ������ ���� ����������� ��������� ��� �������� �������� �������� � ������������ ������ (��)
		// ��� ���������� ��������� ��������� ������������ ������������, �� ������ ���������� ������ ���������� ����� 
		// ����� ����� ���������� ��������� ���������.
		// 8.05.2017.
	doublereal deltaF = fabs(Fe - Fw + Fn - Fs + Ft - Fb);
	if (b_on_adaptive_local_refinement_mesh) {
		deltaF = fabs(Fe - Fw + Fn - Fs + Ft - Fb + Fe2 - Fw2 + Fn2 - Fs2 + Ft2 - Fb2 + Fe3 - Fw3 + Fn3 - Fs3 + Ft3 - Fb3 + Fe4 - Fw4 + Fn4 - Fs4 + Ft4 - Fb4);
	}
	if (deltaF != deltaF) {
		printf("Fe=%e Fw=%e Fn=%e Fs=%e Ft=%e Fb=%e\n", Fe, Fw, Fn, Fs, Ft, Fb);
		printf("Fe2=%e Fw2=%e Fn2=%e Fs2=%e Ft2=%e Fb2=%e\n", Fe2, Fw2, Fn2, Fs2, Ft2, Fb2);
		printf("Fe3=%e Fw3=%e Fn3=%e Fs3=%e Ft3=%e Fb3=%e\n", Fe3, Fw3, Fn3, Fs3, Ft3, Fb3);
		printf("Fe4=%e Fw4=%e Fn4=%e Fs4=%e Ft4=%e Fb4=%e\n", Fe4, Fw4, Fn4, Fs4, Ft4, Fb4);
		printf("ERROR deltaF=%e\n", deltaF);
	}

	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap) {
		printf("ap!=ap assemble bug. Apriory deltaF. iP=%lld ap=%e\n", iP, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap);
		system("pause");
	}
	//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += deltaF;//-->//sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap+=apzero1+deltaF;//+deltaF; // ������������ ������� ������� deltaF ������ ��������������.  ���������� ��������� 
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap) {
		printf("ap!=ap assemble bug. Apost deltaF. iP=%lld ap=%e\n", iP, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap);
		system("pause");
	}


	sl[RE_THETA_LANGTRY_MENTER_SL][iP].b = attrs; // ����� ���������� ��������� ��� ����� ������� ����������� �����������.
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].b != sl[RE_THETA_LANGTRY_MENTER_SL][iP].b) {
		printf("exptsr+attrs error NAN or INF in control volume %lld RE_THETA_LANGTRY_MENTER\n", iP);
		system("pause");
	}


	// ���������

	doublereal P_Theta_T=0.0;

	

	doublereal Usquare = fmax(1.0e-6,potent[VXCOR][iP] * potent[VXCOR][iP] + potent[VYCOR][iP] * potent[VYCOR][iP] + potent[VZCOR][iP] * potent[VZCOR][iP]);
	doublereal T_val = (500.0 * prop[MU_DYNAMIC_VISCOSITY][iP]) / (rP*Usquare);

	doublereal delta = fmax(1.0e-9, 375.0 * potent[CURL][iP] * prop[MU_DYNAMIC_VISCOSITY][iP] * potent[RE_THETA_LANGTRY_MENTER][iP] * distance_to_wall[iP] / (rP* Usquare));

	doublereal Re_omega = rP * fmax(Omega_limiter_min, potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]) * distance_to_wall[iP] * distance_to_wall[iP] / prop[MU_DYNAMIC_VISCOSITY][iP];
	doublereal Fwake = exp(-(Re_omega /1.0e5)*(Re_omega /1.0e5));

	doublereal F_Theta_T = fmin(fmax(Fwake*exp(-(distance_to_wall[iP]/ delta)* (distance_to_wall[iP] / delta)* (distance_to_wall[iP] / delta)* (distance_to_wall[iP] / delta)),1.0-((eqin.fluidinfo[0].Cepsilon2* fmin(1.0, fmax(0.0, potent[GAMMA_LANGTRY_MENTER][iP])) -1.0)/(eqin.fluidinfo[0].Cepsilon2  -1.0))*((eqin.fluidinfo[0].Cepsilon2 * fmin(1.0, fmax(0.0, potent[GAMMA_LANGTRY_MENTER][iP])) - 1.0) / (eqin.fluidinfo[0].Cepsilon2 - 1.0))),1.0);

	doublereal Tu = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][iP] / (3.0 * Usquare));
	if (Tu < 0.027) {
		Tu = 0.027;
	}

	doublereal dU_na_ds = (1.0 / Usquare) * (potent[VELOCITY_X_COMPONENT][iP]*potent[VELOCITY_X_COMPONENT][iP]*potent[GRADXVX][iP]+
		potent[VELOCITY_X_COMPONENT][iP] * potent[VELOCITY_Y_COMPONENT][iP] * potent[GRADYVX][iP] + 
		potent[VELOCITY_X_COMPONENT][iP] * potent[VELOCITY_Z_COMPONENT][iP] * potent[GRADZVX][iP] + 
		potent[VELOCITY_Y_COMPONENT][iP] * potent[VELOCITY_X_COMPONENT][iP] * potent[GRADXVY][iP] +
		potent[VELOCITY_Y_COMPONENT][iP] * potent[VELOCITY_Y_COMPONENT][iP] * potent[GRADYVY][iP] +
		potent[VELOCITY_Y_COMPONENT][iP] * potent[VELOCITY_Z_COMPONENT][iP] * potent[GRADZVY][iP] + 
		potent[VELOCITY_Z_COMPONENT][iP] * potent[VELOCITY_X_COMPONENT][iP] * potent[GRADXVZ][iP] +
		potent[VELOCITY_Z_COMPONENT][iP] * potent[VELOCITY_Y_COMPONENT][iP] * potent[GRADYVZ][iP] +
		potent[VELOCITY_Z_COMPONENT][iP] * potent[VELOCITY_Z_COMPONENT][iP] * potent[GRADZVZ][iP]);

	doublereal thetaT = 0.1; // ������� ������������� ������� ��� ����������.
	doublereal lambda_theta, Re_theta_T_eq;

	// 100 ��������.
	for (integer iter_loc = 0; iter_loc < 100; iter_loc++) 
	{
		lambda_theta = (rP * thetaT * thetaT / prop[MU_DYNAMIC_VISCOSITY][iP]) * dU_na_ds;
		if (lambda_theta > 0.1) {
			lambda_theta = 0.1;
		}
		if (lambda_theta < -0.1) {
			lambda_theta = -0.1;
		}

		Re_theta_T_eq = Re_theta_eq(lambda_theta, Tu);
		if (Re_theta_T_eq < 20.0) {
			Re_theta_T_eq = 20.0;
		}

		// ������ ����������: 0,2.
		thetaT = 0.2 * ((Re_theta_T_eq * prop[MU_DYNAMIC_VISCOSITY][iP]) / (rP * sqrt(Usquare))) + 0.8 * thetaT;
	}


	P_Theta_T = eqin.fluidinfo[0].CthetaT * (rP / (T_val)) * (Re_theta_T_eq - potent[RE_THETA_LANGTRY_MENTER][iP]) * (1.0 - F_Theta_T);
	
	if (brthdsd_ON) {
		dSc = P_Theta_T;
		// ���������� �����������.		
	}
	else {
		dSc = 0.0;
	}
	sl[RE_THETA_LANGTRY_MENTER_SL][iP].b += dSc * dx * dy * dz; // ��������� ����� ����������.


	if (btimedep) {
		// ���� ��������������� �����.
		doublereal apzero = rP * dx * dy * dz / tauparam;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap += apzero;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].b += apzero * Re_Theta_Menter_Lantgry_old_time_step[iP];
	}


	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].b != sl[RE_THETA_LANGTRY_MENTER_SL][iP].b) {
		printf("dSc*dx*dy*dz error NAN or INF in control volume %lld\n", iP);
		printf("RE_THETA_LANGTRY_MENTER\n");
		system("pause");
	}

	if (fabs(sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap) < 1.0e-30) {
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap = 1.0;
		sl[RE_THETA_LANGTRY_MENTER_SL][iP].b = potent[RE_THETA_LANGTRY_MENTER][iP];

		if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].b != sl[RE_THETA_LANGTRY_MENTER_SL][iP].b) {
			printf("Zero ap in RE_THETA_LANGTRY_MENTER component.\n");
		}
	}



	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap) {
		printf("ap!=ap assemble bug. iP=%lld ap=%e\n", iP, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ap);
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae) {
		printf("ae!=ae assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw != sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw) {
		printf("aw!=aw assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].an != sl[RE_THETA_LANGTRY_MENTER_SL][iP].an) {
		printf("an!=an assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].as != sl[RE_THETA_LANGTRY_MENTER_SL][iP].as) {
		printf("as!=as assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].at != sl[RE_THETA_LANGTRY_MENTER_SL][iP].at) {
		printf("at!=at assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab) {
		printf("ab!=ab assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2) {
		printf("ae2!=ae2 assemble bug %e %e\n", sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2, sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae2);
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw2) {
		printf("aw2!=aw2 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].an2) {
		printf("an2!=an2 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].as2) {
		printf("as2!=as2 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].at2) {
		printf("at2!=at2 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab2) {
		printf("ab2!=ab2 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae3) {
		printf("ae3!=ae3 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw3) {
		printf("aw3!=aw3 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].an3) {
		printf("an3!=an3 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].as3) {
		printf("as3!=as3 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].at3) {
		printf("at3!=at3 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab3) {
		printf("ab3!=ab3 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ae4) {
		printf("ae4!=ae4 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].aw4) {
		printf("aw4!=aw4 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].an4) {
		printf("an4!=an4 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].as4) {
		printf("as4!=as4 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].at4) {
		printf("at4!=at4 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}
	if (sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4 != sl[RE_THETA_LANGTRY_MENTER_SL][iP].ab4) {
		printf("ab4!=ab4 assemble bug\n");
		printf("RE_THETA_LANGTRY_MENTER \n");
		system("pause");
	}

} // my_elmatr_quad_turbulent_Re_Theta_Langtry_MenterSST_3D



// ���� ��������� ������� ��� ������� �������������� gamma � ������ ������� �������.
void my_elmatr_quad_Gamma_Langtry_Mentor_3D_bound(integer inumber, integer maxelm,
	bool bDirichlet, BOUND* border_neighbor, integer ls, integer lw,
	WALL* w,
	//integer iVar,
	equation3D_bon*& slb,
	TOCHKA* pa, int** nvtx,
	float** prop_b, float** prop,
	doublereal** potent
	//, integer iflowregime
) {



	// bDirichlet   �������������� ������ ������ ��������� ������� �������.
	// bDirichlet == false �������������� ������ ������ ���������� ������� �������.

	// inumber - ����� ���������� ��.
	// inumber ���������� �� 0..maxbound-1

   /*
   for (integer i=0; i<lw; i++) {
	   if (w[i].bsymmetry) {
		 printf("symmetry \n");
	   }
	   if (w[i].bpressure) {
		 printf("bpressure \n");
	   }
   }
   getchar(); // debug;
   */

   // ��������.
   /* 1. ������� �������.
		 ��� ����������� ������ � ������� ����, �.�. �� ������� ������� ������.
		 �� ������� ��������� � �� �������� ������� bpressure.
		 � �� opening � ������ ��������� ������.
	  2. ������� �������. ��� ������ ������, ������� ��������� ��� �������� �������.
   */

   

	// ������� ������� ��������� ������� �������
   //if (bDirichlet && (border_neighbor[inumber].MCB<(ls + lw)) && (border_neighbor[inumber].MCB >= ls) && (!w[border_neighbor[inumber].MCB - ls].bopening) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry) && (!w[border_neighbor[inumber].MCB - ls].bpressure)) {
	if (bDirichlet && (border_neighbor[inumber].MCB < (ls + lw)) && (border_neighbor[inumber].MCB >= ls) &&
		((!w[border_neighbor[inumber].MCB - ls].bpressure) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry))) {

		//system("PAUSE");

		// ��������� ������� �������
		// ������ �������� �� �������
		// ��� �� ������� ��������� � �� �������� �������.

		slb[inumber].aw = 1.0;
		slb[inumber].ai = 0.0;

		//border_neighbor[inumber].Norm - ���������� �������.
		if (w[border_neighbor[inumber].MCB - ls].bopening) {
			if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(potent[VXCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && (fabs(potent[VYCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && (fabs(potent[VZCOR][maxelm + inumber]) > 1.0e-20)))
			{

				
					// �� ������ ������������ �������� (��������� ������� �� ����� ������� �� ������ �������).
					// 10.10.2019

					if ((border_neighbor[inumber].Norm == E_SIDE) && (potent[VXCOR][maxelm + inumber] > 0.0)) {
						// ������� ������� ������
						slb[inumber].b = 1.0;
					}
					if ((border_neighbor[inumber].Norm == W_SIDE) && (potent[VXCOR][maxelm + inumber] < 0.0)) {
						// ������� ������� ������
						slb[inumber].b = 1.0;
					}
					if ((border_neighbor[inumber].Norm == N_SIDE) && (potent[VYCOR][maxelm + inumber] > 0.0)) {
						// ������� ������� ������
						slb[inumber].b = 1.0;
					}
					if ((border_neighbor[inumber].Norm == S_SIDE) && (potent[VYCOR][maxelm + inumber] < 0.0)) {
						// ������� ������� ������
						slb[inumber].b = 1.0;
					}
					if ((border_neighbor[inumber].Norm == T_SIDE) && (potent[VZCOR][maxelm + inumber] > 0.0)) {
						// ������� ������� ������
						slb[inumber].b = 1.0;
					}
					if ((border_neighbor[inumber].Norm == B_SIDE) && (potent[VZCOR][maxelm + inumber] < 0.0)) {
						// ������� ������� ������
						slb[inumber].b = 1.0;
					}
				

			}
		}
		else if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(w[border_neighbor[inumber].MCB - ls].Vx) > 1.0e-20)) ||
			((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vy) > 1.0e-20) ||
			((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vz) > 1.0e-20))
		{

			
				
				// �������������� ������� (gamma=1.0), ������� ��������� ������������.

				if ((border_neighbor[inumber].Norm == E_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx > 0.0)) {
					// ������� ������� ������
					slb[inumber].b = 1.0;
				}
				if ((border_neighbor[inumber].Norm == W_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx < 0.0)) {
					// ������� ������� ������
					slb[inumber].b = 1.0;
				}
				if ((border_neighbor[inumber].Norm == N_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy > 0.0)) {
					// ������� ������� ������
					slb[inumber].b = 1.0;
				}
				if ((border_neighbor[inumber].Norm == S_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy < 0.0)) {
					// ������� ������� ������
					slb[inumber].b = 1.0;
				}
				if ((border_neighbor[inumber].Norm == T_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz > 0.0)) {
					// ������� ������� ������
					slb[inumber].b = 1.0;
				}
				if ((border_neighbor[inumber].Norm == B_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz < 0.0)) {
					// ������� ������� ������
					slb[inumber].b = 1.0;
				}
			

		}
		else {
			// ����������� ������.
			// ������������ ������� ������������ ��������� �� ������� ����������� ������ ����� ����.
			slb[inumber].b = 0.0;
		}

		slb[inumber].iI = NON_EXISTENT_NODE; // �� ������������ � �������
		slb[inumber].iW = border_neighbor[inumber].iB;
#if doubleintprecision == 1
		//printf("%lld, soseddb=%lld\n",inumber, border_neighbor[inumber].iB); getchar(); // debug
#else
		//printf("%d, soseddb=%d\n",inumber, border_neighbor[inumber].iB); getchar(); // debug
#endif


		// ��� ������� �������:
		// ������ ������������ ������� 
		// �� ����� ����.
		slb[inumber].iW1 = NON_EXISTENT_NODE;
		slb[inumber].iW2 = NON_EXISTENT_NODE;
		slb[inumber].iW3 = NON_EXISTENT_NODE;
		slb[inumber].iW4 = NON_EXISTENT_NODE;
	}
	else if (!bDirichlet && ((border_neighbor[inumber].MCB == (ls + lw)) || (border_neighbor[inumber].MCB < ls))) {
		// �������� ���� �������� ������ �������.
		// ���� ������ ������. ������ ������ ����������� �� ������� (border_neighbor[inumber].MCB==(ls+lw)).

		// ��������� ������� ������� ��� �������������� gamma.
		// ������ ������� ���������� �� ������ ������.
		// if (!bDirichlet) �����������
		Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
		
	}
	else if ((border_neighbor[inumber].MCB < (ls + lw)) && (border_neighbor[inumber].MCB >= ls) && ((w[border_neighbor[inumber].MCB - ls].bpressure) || (w[border_neighbor[inumber].MCB - ls].bsymmetry)/*|| ((w[border_neighbor[inumber].MCB - ls].bopening))*/)) {

		// �������� ������� ������.


		if (!bDirichlet) {
			// ��� ����� ���������� �������� ����� ������� �������:
			// ���������. �.�. ��������� ��� ���������� ��������
			// ��������� � ���������� ���������������� � ���������
			// �� ������� ������� � ������������� � ���������, ��
			// ��������� ������� ��� ��������� ����������������
			// ����������� �� ��������� ��� ���������� ��������, � ��������� 
			// �� ������������� � ���������.

			Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);

		}
	}
	else if (!bDirichlet && (border_neighbor[inumber].MCB < (ls + lw)) && (border_neighbor[inumber].MCB >= ls) && ((!w[border_neighbor[inumber].MCB - ls].bpressure) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry))) {
		// ������� �������� �������. ���������� ������� �������.

		if (w[border_neighbor[inumber].MCB - ls].bopening) {
			if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(potent[VXCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && (fabs(potent[VYCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && (fabs(potent[VZCOR][maxelm + inumber]) > 1.0e-20)))
			{

				if ((border_neighbor[inumber].Norm == E_SIDE) && (potent[VXCOR][maxelm + inumber] < 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == W_SIDE) && (potent[VXCOR][maxelm + inumber] > 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == N_SIDE) && (potent[VYCOR][maxelm + inumber] < 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == S_SIDE) && (potent[VYCOR][maxelm + inumber] > 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == T_SIDE) && (potent[VZCOR][maxelm + inumber] < 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == B_SIDE) && (potent[VZCOR][maxelm + inumber] > 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}

			}
		}
		else if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(w[border_neighbor[inumber].MCB - ls].Vx) > 1.0e-20)) ||
			((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vy) > 1.0e-20) ||
			((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vz) > 1.0e-20))
		{

			if ((border_neighbor[inumber].Norm == E_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx < 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == W_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx > 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == N_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy < 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == S_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy > 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == T_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz < 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == B_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz > 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}

		}
	}

} // my_elmatr_quad_Gamma_Langtry_Mentor_3D_bound


// ���� ��������� ������� ��� ������� Re_Theta_T � ������ ������� �������.
void my_elmatr_quad_Re_Theta_T_Langtry_Mentor_3D_bound(integer inumber, integer maxelm,
	bool bDirichlet, BOUND* border_neighbor, integer ls, integer lw,
	WALL* w,
	//integer iVar,
	equation3D_bon*& slb,
	TOCHKA* pa, int** nvtx,
	float** prop_b, float** prop,
	doublereal** potent
	//, integer iflowregime
) {



	// bDirichlet   �������������� ������ ������ ��������� ������� �������.
	// bDirichlet == false �������������� ������ ������ ���������� ������� �������.

	// inumber - ����� ���������� ��.
	// inumber ���������� �� 0..maxbound-1

   /*
   for (integer i=0; i<lw; i++) {
	   if (w[i].bsymmetry) {
		 printf("symmetry \n");
	   }
	   if (w[i].bpressure) {
		 printf("bpressure \n");
	   }
   }
   getchar(); // debug;
   */

   // ��������.
   /* 1. ������� �������.
         ������� ��������� ����������� �� ������� �������.
		 ��� ����������� ������ � ������� ����, �.�. �� ������� ������� ������.
		 �� ������� ��������� � �� �������� ������� bpressure.
		 � �� opening � ������ ��������� ������.
	  2. ������� �������. ��� ������ ������, ������� ��������� ��� �������� �������.
   */



   // ������� ������� ��������� ������� �������
  //if (bDirichlet && (border_neighbor[inumber].MCB<(ls + lw)) && (border_neighbor[inumber].MCB >= ls) && (!w[border_neighbor[inumber].MCB - ls].bopening) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry) && (!w[border_neighbor[inumber].MCB - ls].bpressure)) {
	if (bDirichlet && (border_neighbor[inumber].MCB < (ls + lw)) && (border_neighbor[inumber].MCB >= ls) &&
		((!w[border_neighbor[inumber].MCB - ls].bpressure) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry))) {

		//system("PAUSE");

		// ��������� ������� �������
		// ������ �������� �� �������
		// ��� �� ������� ��������� � �� �������� �������.

		slb[inumber].aw = 1.0;
		slb[inumber].ai = 0.0;

		//border_neighbor[inumber].Norm - ���������� �������.
		if (w[border_neighbor[inumber].MCB - ls].bopening) {
			if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(potent[VXCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && (fabs(potent[VYCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && (fabs(potent[VZCOR][maxelm + inumber]) > 1.0e-20)))
			{


				// �� ������ ������������ �������� (��������� ������� �� ����� ������� �� ������ �������).
				// 10.10.2019

				if ((border_neighbor[inumber].Norm == E_SIDE) && (potent[VXCOR][maxelm + inumber] > 0.0)) {
					// ������� ������� ������
					doublereal Tu_infinity = 100.0*sqrt(2.0*potent[TURBULENT_KINETIK_ENERGY][maxelm+inumber]/3.0)/fabs(potent[VXCOR][maxelm + inumber]);
					slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
				}
				if ((border_neighbor[inumber].Norm == W_SIDE) && (potent[VXCOR][maxelm + inumber] < 0.0)) {
					// ������� ������� ������
					doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(potent[VXCOR][maxelm + inumber]);
					slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
				}
				if ((border_neighbor[inumber].Norm == N_SIDE) && (potent[VYCOR][maxelm + inumber] > 0.0)) {
					// ������� ������� ������
					doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(potent[VYCOR][maxelm + inumber]);
					slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
				}
				if ((border_neighbor[inumber].Norm == S_SIDE) && (potent[VYCOR][maxelm + inumber] < 0.0)) {
					// ������� ������� ������
					doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(potent[VYCOR][maxelm + inumber]);
					slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
				}
				if ((border_neighbor[inumber].Norm == T_SIDE) && (potent[VZCOR][maxelm + inumber] > 0.0)) {
					// ������� ������� ������
					doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(potent[VZCOR][maxelm + inumber]);
					slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
				}
				if ((border_neighbor[inumber].Norm == B_SIDE) && (potent[VZCOR][maxelm + inumber] < 0.0)) {
					// ������� ������� ������
					doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(potent[VZCOR][maxelm + inumber]);
					slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
				}


			}
		}
		else if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(w[border_neighbor[inumber].MCB - ls].Vx) > 1.0e-20)) ||
			((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vy) > 1.0e-20) ||
			((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vz) > 1.0e-20))
		{



			// �������������� ������� (gamma=1.0), ������� ��������� ������������.

			if ((border_neighbor[inumber].Norm == E_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx > 0.0)) {
				// ������� ������� ������
				doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(w[border_neighbor[inumber].MCB - ls].Vx);
				slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
			}
			if ((border_neighbor[inumber].Norm == W_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx < 0.0)) {
				// ������� ������� ������
				doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(w[border_neighbor[inumber].MCB - ls].Vx);
				slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
			}
			if ((border_neighbor[inumber].Norm == N_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy > 0.0)) {
				// ������� ������� ������
				doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(w[border_neighbor[inumber].MCB - ls].Vy);
				slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
			}
			if ((border_neighbor[inumber].Norm == S_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy < 0.0)) {
				// ������� ������� ������
				doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(w[border_neighbor[inumber].MCB - ls].Vy);
				slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
			}
			if ((border_neighbor[inumber].Norm == T_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz > 0.0)) {
				// ������� ������� ������
				doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(w[border_neighbor[inumber].MCB - ls].Vz);
				slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
			}
			if ((border_neighbor[inumber].Norm == B_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz < 0.0)) {
				// ������� ������� ������
				doublereal Tu_infinity = 100.0 * sqrt(2.0 * potent[TURBULENT_KINETIK_ENERGY][maxelm + inumber] / 3.0) / fabs(w[border_neighbor[inumber].MCB - ls].Vz);
				slb[inumber].b = Re_theta_T_farfield(Tu_infinity);
			}


		}
		else {
			// ����������� ������.
			// ������������ ������� ������������ ��������� �� ������� ����������� ������ ����� ����.
			slb[inumber].b = 0.0;
		}

		slb[inumber].iI = NON_EXISTENT_NODE; // �� ������������ � �������
		slb[inumber].iW = border_neighbor[inumber].iB;
#if doubleintprecision == 1
		//printf("%lld, soseddb=%lld\n",inumber, border_neighbor[inumber].iB); getchar(); // debug
#else
		//printf("%d, soseddb=%d\n",inumber, border_neighbor[inumber].iB); getchar(); // debug
#endif


		// ��� ������� �������:
		// ������ ������������ ������� 
		// �� ����� ����.
		slb[inumber].iW1 = NON_EXISTENT_NODE;
		slb[inumber].iW2 = NON_EXISTENT_NODE;
		slb[inumber].iW3 = NON_EXISTENT_NODE;
		slb[inumber].iW4 = NON_EXISTENT_NODE;
	}
	else if (!bDirichlet && ((border_neighbor[inumber].MCB == (ls + lw)) || (border_neighbor[inumber].MCB < ls))) {
		// �������� ���� �������� ������ �������.
		// ���� ������ ������. ������ ������ ����������� �� ������� (border_neighbor[inumber].MCB==(ls+lw)).

		// ��������� ������� ������� ��� �������������� gamma.
		// ������ ������� ���������� �� ������ ������.
		// if (!bDirichlet) // �����������
		Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
		
	}
	else if ((border_neighbor[inumber].MCB < (ls + lw)) && (border_neighbor[inumber].MCB >= ls) && ((w[border_neighbor[inumber].MCB - ls].bpressure) || (w[border_neighbor[inumber].MCB - ls].bsymmetry)/*|| ((w[border_neighbor[inumber].MCB - ls].bopening))*/)) {

		// �������� ������� ������.


		if (!bDirichlet) {
			// ��� ����� ���������� �������� ����� ������� �������:
			// ���������. �.�. ��������� ��� ���������� ��������
			// ��������� � ���������� ���������������� � ���������
			// �� ������� ������� � ������������� � ���������, ��
			// ��������� ������� ��� ��������� ����������������
			// ����������� �� ��������� ��� ���������� ��������, � ��������� 
			// �� ������������� � ���������.

			Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);

		}
	}
	else if (!bDirichlet && (border_neighbor[inumber].MCB < (ls + lw)) && (border_neighbor[inumber].MCB >= ls) && ((!w[border_neighbor[inumber].MCB - ls].bpressure) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry))) {
		// ������� �������� �������. ���������� ������� �������.

		if (w[border_neighbor[inumber].MCB - ls].bopening) {
			if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(potent[VXCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && (fabs(potent[VYCOR][maxelm + inumber]) > 1.0e-20)) ||
				((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && (fabs(potent[VZCOR][maxelm + inumber]) > 1.0e-20)))
			{

				if ((border_neighbor[inumber].Norm == E_SIDE) && (potent[VXCOR][maxelm + inumber] < 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == W_SIDE) && (potent[VXCOR][maxelm + inumber] > 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == N_SIDE) && (potent[VYCOR][maxelm + inumber] < 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == S_SIDE) && (potent[VYCOR][maxelm + inumber] > 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == T_SIDE) && (potent[VZCOR][maxelm + inumber] < 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}
				if ((border_neighbor[inumber].Norm == B_SIDE) && (potent[VZCOR][maxelm + inumber] > 0.0)) {
					// �������� ������� ������
					Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
				}

			}
		}
		else if (((border_neighbor[inumber].Norm == E_SIDE || border_neighbor[inumber].Norm == W_SIDE) && (fabs(w[border_neighbor[inumber].MCB - ls].Vx) > 1.0e-20)) ||
			((border_neighbor[inumber].Norm == N_SIDE || border_neighbor[inumber].Norm == S_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vy) > 1.0e-20) ||
			((border_neighbor[inumber].Norm == T_SIDE || border_neighbor[inumber].Norm == B_SIDE) && fabs(w[border_neighbor[inumber].MCB - ls].Vz) > 1.0e-20))
		{

			if ((border_neighbor[inumber].Norm == E_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx < 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == W_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vx > 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == N_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy < 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == S_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vy > 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == T_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz < 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}
			if ((border_neighbor[inumber].Norm == B_SIDE) && (w[border_neighbor[inumber].MCB - ls].Vz > 0.0)) {
				// �������� ������� ������
				Neiman_Zero_in_Wall_STUB(inumber, slb, border_neighbor);
			}

		}
	}

} // my_elmatr_quad_Re_Theta_T_Langtry_Mentor_3D_bound


#endif
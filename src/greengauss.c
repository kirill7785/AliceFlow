// Файл greengauss.c вычисляет градиенты величин 
// в центрах контрольных объёмов.
// begin 20 октября 2011 г. 
// ссылки: Гаврилов Андрей опыт разработки cfd.
// Основано на теореме Грина-Гаусса.
// begin 15 мая 2012 года. Реализован второй порядок точности.

#pragma once
#ifndef GREEN_GAUSS_C
#define GREEN_GAUSS_C 1

// 14.04.2019 Работает на АЛИС сетке.
// Вычисление градиентов скоростей в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
// Поскольку интерполяция линейная то точность данной формулы O(h). 
// По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
// что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
// совпадает со взвешенным методом наименьших квадратов.
void green_gaussO1(integer iP, doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
								 doublereal* &mf, float* &prop, float* &prop_b,
	BOUND* &border_neighbor, integer *ilevel_alice, TOCHKA* &volume_loc) {


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу 1983г.
	bool bRCh=false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {

		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}
	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;
	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;
	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 >= maxelm) bE2 = true;
		if (iN2 >= maxelm) bN2 = true;
		if (iT2 >= maxelm) bT2 = true;
		if (iW2 >= maxelm) bW2 = true;
		if (iS2 >= maxelm) bS2 = true;
		if (iB2 >= maxelm) bB2 = true;

		

		if (iE3 >= maxelm) bE3 = true;
		if (iN3 >= maxelm) bN3 = true;
		if (iT3 >= maxelm) bT3 = true;
		if (iW3 >= maxelm) bW3 = true;
		if (iS3 >= maxelm) bS3 = true;
		if (iB3 >= maxelm) bB3 = true;

		

		if (iE4 >= maxelm) bE4 = true;
		if (iN4 >= maxelm) bN4 = true;
		if (iT4 >= maxelm) bT4 = true;
		if (iW4 >= maxelm) bW4 = true;
		if (iS4 >= maxelm) bS4 = true;
		if (iB4 >= maxelm) bB4 = true;
	}

	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	//volume3D(iP, nvtx, pa, dx, dy, dz);
	//dx = fabs(dx);
	//dy = fabs(dy);
	//dz = fabs(dz);

	TOCHKA point_loc = volume_loc[iP];
	dx = point_loc.x;
	dy = point_loc.y;
	dz = point_loc.z;


	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	if (b_on_adaptive_local_refinement_mesh) {

		// т.к. известна нумерация вершин куба, то здесь она используется
		// x - direction
		if (iE2 > -1) {
			if (!bE2) dxe2 = 0.5 * (pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
			if (!bE2) dxe2 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		}
		if (iW2 > -1) {
			if (!bW2) dxw2 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
			if (!bW2) dxw2 -= 0.5 * (pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
		}
		// y - direction
		if (iN2 > -1) {
			if (!bN2) dyn2 = 0.5 * (pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
			if (!bN2) dyn2 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		}
		if (iS2 > -1) {
			if (!bS2) dys2 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
			if (!bS2) dys2 -= 0.5 * (pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
		}
		// z - direction
		if (iT2 > -1) {
			if (!bT2) dzt2 = 0.5 * (pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
			if (!bT2) dzt2 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		}
		if (iB2 > -1) {
			if (!bB2) dzb2 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
			if (!bB2) dzb2 -= 0.5 * (pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
		}

		// т.к. известна нумерация вершин куба, то здесь она используется
		// x - direction
		if (iE3 > -1) {
			if (!bE3) dxe3 = 0.5 * (pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
			if (!bE3) dxe3 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		}
		if (iW3 > -1) {
			if (!bW3) dxw3 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
			if (!bW3) dxw3 -= 0.5 * (pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
		}
		// y - direction
		if (iN3 > -1) {
			if (!bN3) dyn3 = 0.5 * (pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
			if (!bN3) dyn3 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		}
		if (iS3 > -1) {
			if (!bS3) dys3 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
			if (!bS3) dys3 -= 0.5 * (pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
		}
		// z - direction
		if (iT3 > -1) {
			if (!bT3) dzt3 = 0.5 * (pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
			if (!bT3) dzt3 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		}
		if (iB3 > -1) {
			if (!bB3) dzb3 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
			if (!bB3) dzb3 -= 0.5 * (pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
		}

		// т.к. известна нумерация вершин куба, то здесь она используется
		// x - direction
		if (iE4 > -1) {
			if (!bE4) dxe4 = 0.5 * (pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
			if (!bE4) dxe4 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		}
		if (iW4 > -1) {
			if (!bW4) dxw4 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
			if (!bW4) dxw4 -= 0.5 * (pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
		}
		// y - direction
		if (iN4 > -1) {
			if (!bN4) dyn4 = 0.5 * (pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
			if (!bN4) dyn4 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		}
		if (iS4 > -1) {
			if (!bS4) dys4 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
			if (!bS4) dys4 -= 0.5 * (pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
		}
		// z - direction
		if (iT4 > -1) {
			if (!bT4) dzt4 = 0.5 * (pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
			if (!bT4) dzt4 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		}
		if (iB4 > -1) {
			if (!bB4) dzb4 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
			if (!bB4) dzb4 -= 0.5 * (pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
		}

	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus=0.5*dx/dxe;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB];
				doublereal dx_loc = point_loc.x;
				doublereal dy_loc = point_loc.y;
				//dz_loc = point_loc.z;

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.

	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 > -1) {

			dSqe2 = dy * dz;

			if (bE2) {
				// граничный узел.
				dSqe2 = border_neighbor[iE2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
					dSqe2 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE2];
					//dx_loc = point_loc.x;
					doublereal  dy_loc = point_loc.y;
					doublereal  dz_loc = point_loc.z;

					dSqe2 = dy_loc * dz_loc;
				}
			}


		}


		if (iW2 > -1) {
			dSqw2 = dy * dz;

			if (bW2) {
				// граничный узел.
				dSqw2 = border_neighbor[iW2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
					dSqw2 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW2];
					//dx_loc = point_loc.x;
					doublereal dy_loc = point_loc.y;
					doublereal dz_loc = point_loc.z;

					dSqw2 = dy_loc * dz_loc;
				}
			}


		}


		if (iN2 > -1) {

			dSqn2 = dx * dz;

			if (bN2) {
				// граничный узел.
				dSqn2 = border_neighbor[iN2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
					dSqn2 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN2];
					doublereal dx_loc = point_loc.x;
					//dy_loc = point_loc.y;
					doublereal dz_loc = point_loc.z;


					dSqn2 = dx_loc * dz_loc;
				}
			}


		}


		if (iS2 > -1) {

			dSqs2 = dx * dz;

			if (bS2) {
				// граничный узел.
				dSqs2 = border_neighbor[iS2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
					dSqs2 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS2];
					doublereal dx_loc = point_loc.x;
					//dy_loc = point_loc.y;
					doublereal dz_loc = point_loc.z;

					dSqs2 = dx_loc * dz_loc;
				}
			}


		}


		if (iT2 > -1) {

			dSqt2 = dx * dy;

			if (bT2) {
				// граничный узел.
				dSqt2 = border_neighbor[iT2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
					dSqt2 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT2];
					doublereal dx_loc = point_loc.x;
					doublereal dy_loc = point_loc.y;
					//dz_loc = point_loc.z;

					dSqt2 = dx_loc * dy_loc;
				}
			}


		}


		if (iB2 > -1) {

			dSqb2 = dx * dy;

			if (bB2) {
				// граничный узел.
				dSqb2 = border_neighbor[iB2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
					dSqb2 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					//doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB2];
					doublereal dx_loc = point_loc.x;
					doublereal dy_loc = point_loc.y;
					//dz_loc = point_loc.z;

					dSqb2 = dx_loc * dy_loc;
				}
			}


		}




		if (iE3 > -1) {

			dSqe3 = dy * dz;

			if (bE3) {
				// граничный узел.
				dSqe3 = border_neighbor[iE3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
					dSqe3 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqe3 = dy_loc * dz_loc;
				}
			}


		}


		if (iW3 > -1) {

			dSqw3 = dy * dz;

			if (bW3) {
				// граничный узел.
				dSqw3 = border_neighbor[iW3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
					dSqw3 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqw3 = dy_loc * dz_loc;
				}
			}


		}


		if (iN3 > -1) {

			dSqn3 = dx * dz;

			if (bN3) {
				// граничный узел.
				dSqn3 = border_neighbor[iN3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
					dSqn3 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqn3 = dx_loc * dz_loc;
				}
			}


		}


		if (iS3 > -1) {

			dSqs3 = dx * dz;

			if (bS3) {
				// граничный узел.
				dSqs3 = border_neighbor[iS3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
					dSqs3 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqs3 = dx_loc * dz_loc;
				}
			}


		}


		if (iT3 > -1) {

			dSqt3 = dx * dy;

			if (bT3) {
				// граничный узел.
				dSqt3 = border_neighbor[iT3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
					dSqt3 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqt3 = dx_loc * dy_loc;
				}
			}


		}


		if (iB3 > -1) {

			dSqb3 = dx * dy;

			if (bB3) {
				// граничный узел.
				dSqb3 = border_neighbor[iB3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
					dSqb3 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqb3 = dx_loc * dy_loc;
				}
			}


		}

		



		if (iE4 > -1) {

			dSqe4 = dy * dz;

			if (bE4) {
				// граничный узел.
				dSqe4 = border_neighbor[iE4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
					dSqe4 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqe4 = dy_loc * dz_loc;
				}
			}


		}


		if (iW4 > -1) {

			dSqw4 = dy * dz;

			if (bW4) {
				// граничный узел.
				dSqw4 = border_neighbor[iW4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
					dSqw4 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqw4 = dy_loc * dz_loc;
				}
			}


		}


		if (iN4 > -1) {

			dSqn4 = dx * dz;

			if (bN4) {
				// граничный узел.
				dSqn4 = border_neighbor[iN4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
					dSqn4 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;


					dSqn4 = dx_loc * dz_loc;
				}
			}


		}


		if (iS4 > -1) {

			dSqs4 = dx * dz;

			if (bS4) {
				// граничный узел.
				dSqs4 = border_neighbor[iS4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
					dSqs4 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqs4 = dx_loc * dz_loc;
				}
			}


		}


		if (iT4 > -1) {

			dSqt4 = dx * dy;

			if (bT4) {
				// граничный узел.
				dSqt4 = border_neighbor[iT4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
					dSqt4 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqt4 = dx_loc * dy_loc;
				}
			}


		}


		if (iB4 > -1) {

			dSqb4 = dx * dy;

			if (bB4) {
				// граничный узел.
				dSqb4 = border_neighbor[iB4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
					dSqb4 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqb4 = dx_loc * dy_loc;
				}
			}


		}
	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}
	

	// плотность на грани КО аппроксимируется средним гармоническим
	doublereal rP=1.0, rE=1.0, rN=1.0, rT=1.0, rW=1.0, rS=1.0, rB=1.0;
	doublereal rE2 = 1.0, rN2 = 1.0, rT2 = 1.0, rW2 = 1.0, rS2 = 1.0, rB2 = 1.0;
	doublereal rE3 = 1.0, rN3 = 1.0, rT3 = 1.0, rW3 = 1.0, rS3 = 1.0, rB3 = 1.0;
	doublereal rE4 = 1.0, rN4 = 1.0, rT4 = 1.0, rW4 = 1.0, rS4 = 1.0, rB4 = 1.0;
	doublereal rhoe = 1.0, rhow = 1.0, rhon = 1.0, rhos = 1.0, rhot = 1.0, rhob = 1.0;
	doublereal rhoe2 = 1.0, rhow2 = 1.0, rhon2 = 1.0, rhos2 = 1.0, rhot2 = 1.0, rhob2 = 1.0;
	doublereal rhoe3 = 1.0, rhow3 = 1.0, rhon3 = 1.0, rhos3 = 1.0, rhot3 = 1.0, rhob3 = 1.0;
	doublereal rhoe4 = 1.0, rhow4 = 1.0, rhon4 = 1.0, rhos4 = 1.0, rhot4 = 1.0, rhob4 = 1.0;


	if (bRCh) {
		rP=prop[iP];
		if (iE > -1) {
			if (!bE) rE = prop[iE]; else rE = prop_b[iE - maxelm];
		}
		if (iN > -1) {
			if (!bN) rN = prop[iN]; else rN = prop_b[iN - maxelm];
		}
		if (iT > -1) {
			if (!bT) rT = prop[iT]; else rT = prop_b[iT - maxelm];
		}
		if (iW > -1) {
			if (!bW) rW = prop[iW]; else rW = prop_b[iW - maxelm];
		}
		if (iS > -1) {
			if (!bS) rS = prop[iS]; else rS = prop_b[iS - maxelm];
		}
		if (iB > -1) {
			if (!bB) rB = prop[iB]; else rB = prop_b[iB - maxelm];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				if (!bE2) rE2 = prop[iE2]; else rE2 = prop_b[iE2 - maxelm];
			}
			if (iN2 > -1) {
				if (!bN2) rN2 = prop[iN2]; else rN2 = prop_b[iN2 - maxelm];
			}
			if (iT2 > -1) {
				if (!bT2) rT2 = prop[iT2]; else rT2 = prop_b[iT2 - maxelm];
			}
			if (iW2 > -1) {
				if (!bW2) rW2 = prop[iW2]; else rW2 = prop_b[iW2 - maxelm];
			}
			if (iS2 > -1) {
				if (!bS2) rS2 = prop[iS2]; else rS2 = prop_b[iS2 - maxelm];
			}
			if (iB2 > -1) {
				if (!bB2) rB2 = prop[iB2]; else rB2 = prop_b[iB2 - maxelm];
			}

			if (iE3 > -1) {
				if (!bE3) rE3 = prop[iE3]; else rE3 = prop_b[iE3 - maxelm];
			}
			if (iN3 > -1) {
				if (!bN3) rN3 = prop[iN3]; else rN3 = prop_b[iN3 - maxelm];
			}
			if (iT3 > -1) {
				if (!bT3) rT3 = prop[iT3]; else rT3 = prop_b[iT3 - maxelm];
			}
			if (iW3 > -1) {
				if (!bW3) rW3 = prop[iW3]; else rW3 = prop_b[iW3 - maxelm];
			}
			if (iS3 > -1) {
				if (!bS3) rS3 = prop[iS3]; else rS3 = prop_b[iS3 - maxelm];
			}
			if (iB3 > -1) {
				if (!bB3) rB3 = prop[iB3]; else rB3 = prop_b[iB3 - maxelm];
			}

			if (iE4 > -1) {
				if (!bE4) rE4 = prop[iE4]; else rE4 = prop_b[iE4 - maxelm];
			}
			if (iN4 > -1) {
				if (!bN4) rN4 = prop[iN4]; else rN4 = prop_b[iN4 - maxelm];
			}
			if (iT4 > -1) {
				if (!bT4) rT4 = prop[iT4]; else rT4 = prop_b[iT4 - maxelm];
			}
			if (iW4 > -1) {
				if (!bW4) rW4 = prop[iW4]; else rW4 = prop_b[iW4 - maxelm];
			}
			if (iS4 > -1) {
				if (!bS4) rS4 = prop[iS4]; else rS4 = prop_b[iS4 - maxelm];
			}
			if (iB4 > -1) {
				if (!bB4) rB4 = prop[iB4]; else rB4 = prop_b[iB4 - maxelm];
			}

		}
	    
		if (iE > -1) {
			rhoe = rE * rP / (feplus*rE + (1.0 - feplus)*rP); // проверено.
		}
		if (iW > -1) {
			rhow = rW * rP / (fwplus*rW + (1.0 - fwplus)*rP);
		}
		if (iN > -1) {
			rhon = rN * rP / (fnplus*rN + (1.0 - fnplus)*rP);
		}
		if (iS > -1) {
			rhos = rS * rP / (fsplus*rS + (1.0 - fsplus)*rP);
		}
		if (iT > -1) {
			rhot = rT * rP / (ftplus*rT + (1.0 - ftplus)*rP);
		}
		if (iB > -1) {
			rhob = rB * rP / (fbplus*rB + (1.0 - fbplus)*rP);
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				rhoe2 = rE2 * rP / (feplus2 * rE2 + (1.0 - feplus2) * rP); // проверено.
			}
			if (iW2 > -1) {
				rhow2 = rW2 * rP / (fwplus2 * rW2 + (1.0 - fwplus2) * rP);
			}
			if (iN2 > -1) {
				rhon2 = rN2 * rP / (fnplus2 * rN2 + (1.0 - fnplus2) * rP);
			}
			if (iS2 > -1) {
				rhos2 = rS2 * rP / (fsplus2 * rS2 + (1.0 - fsplus2) * rP);
			}
			if (iT2 > -1) {
				rhot2 = rT2 * rP / (ftplus2 * rT2 + (1.0 - ftplus2) * rP);
			}
			if (iB2 > -1) {
				rhob2 = rB2 * rP / (fbplus2 * rB2 + (1.0 - fbplus2) * rP);
			}

			if (iE3 > -1) {
				rhoe3 = rE3 * rP / (feplus3 * rE3 + (1.0 - feplus3) * rP); // проверено.
			}
			if (iW3 > -1) {
				rhow3 = rW3 * rP / (fwplus3 * rW3 + (1.0 - fwplus3) * rP);
			}
			if (iN3 > -1) {
				rhon3 = rN3 * rP / (fnplus3 * rN3 + (1.0 - fnplus3) * rP);
			}
			if (iS3 > -1) {
				rhos3 = rS3 * rP / (fsplus3 * rS3 + (1.0 - fsplus3) * rP);
			}
			if (iT3 > -1) {
				rhot3 = rT3 * rP / (ftplus3 * rT3 + (1.0 - ftplus3) * rP);
			}
			if (iB3 > -1) {
				rhob3 = rB3 * rP / (fbplus3 * rB3 + (1.0 - fbplus3) * rP);
			}

			if (iE4 > -1) {
				rhoe4 = rE4 * rP / (feplus4 * rE4 + (1.0 - feplus4) * rP); // проверено.
			}
			if (iW4 > -1) {
				rhow4 = rW4 * rP / (fwplus4 * rW4 + (1.0 - fwplus4) * rP);
			}
			if (iN4 > -1) {
				rhon4 = rN4 * rP / (fnplus4 * rN4 + (1.0 - fnplus4) * rP);
			}
			if (iS4 > -1) {
				rhos4 = rS4 * rP / (fsplus4 * rS4 + (1.0 - fsplus4) * rP);
			}
			if (iT4 > -1) {
				rhot4 = rT4 * rP / (ftplus4 * rT4 + (1.0 - ftplus4) * rP);
			}
			if (iB4 > -1) {
				rhob4 = rB4 * rP / (fbplus4 * rB4 + (1.0 - fbplus4) * rP);
			}

		}
	}

	// линейная интерполяция скорости VX на грань КО.
    doublereal fe=0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// VX
		if (iE > -1) {
			if (!bE) {
				if (bRCh) {
					fe = mf[E_SIDE] / (rhoe*dy*dz); // скорость на грани с учётом поправки Рхи-Чоу 1983г.
				}
				else {
					fe = feplus*potent[VXCOR][iE] + (1.0 - feplus)*potent[VXCOR][iP];
				}
			}
			else fe = potent[VXCOR][iE];
		}
		if (iW > -1) {
			if (!bW) {
				if (bRCh) {
					fw = mf[W_SIDE] / (rhow*dy*dz); // скорость на грани с учётом монотонизирующей поправки Рхи-Чоу 1983г.
				}
				else {
					fw = fwplus*potent[VXCOR][iW] + (1.0 - fwplus)*potent[VXCOR][iP];
				}
			}
			else fw = potent[VXCOR][iW];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				if (!bE2) {
					if (bRCh) {
						fe2 = mf[E_SIDE] / (rhoe2 * dy * dz); // скорость на грани с учётом поправки Рхи-Чоу.
					}
					else {
						fe2 = feplus2 * potent[VXCOR][iE2] + (1.0 - feplus2) * potent[VXCOR][iP];
					}
				}
				else fe2 = potent[VXCOR][iE2];
			}
			if (iW2 > -1) {
				if (!bW2) {
					if (bRCh) {
						fw2 = mf[W_SIDE] / (rhow2 * dy * dz); // скорость на грани с учётом монотонизирующей поправки Рхи-Чоу.
					}
					else {
						fw2 = fwplus2 * potent[VXCOR][iW2] + (1.0 - fwplus2) * potent[VXCOR][iP];
					}
				}
				else fw2 = potent[VXCOR][iW2];
			}

			if (iE3 > -1) {
				if (!bE3) {
					if (bRCh) {
						fe3 = mf[E_SIDE] / (rhoe3 * dy * dz); // скорость на грани с учётом поправки Рхи-Чоу.
					}
					else {
						fe3 = feplus3 * potent[VXCOR][iE3] + (1.0 - feplus3) * potent[VXCOR][iP];
					}
				}
				else fe3 = potent[VXCOR][iE3];
			}
			if (iW3 > -1) {
				if (!bW3) {
					if (bRCh) {
						fw3 = mf[W_SIDE] / (rhow3 * dy * dz); // скорость на грани с учётом монотонизирующей поправки Рхи-Чоу.
					}
					else {
						fw3 = fwplus3 * potent[VXCOR][iW3] + (1.0 - fwplus3) * potent[VXCOR][iP];
					}
				}
				else fw3 = potent[VXCOR][iW3];
			}

			if (iE4 > -1) {
				if (!bE4) {
					if (bRCh) {
						fe4 = mf[E_SIDE] / (rhoe4 * dy * dz); // скорость на грани с учётом поправки Рхи-Чоу.
					}
					else {
						fe4 = feplus4 * potent[VXCOR][iE4] + (1.0 - feplus4) * potent[VXCOR][iP];
					}
				}
				else fe4 = potent[VXCOR][iE4];
			}
			if (iW4 > -1) {
				if (!bW4) {
					if (bRCh) {
						fw4 = mf[W_SIDE] / (rhow4 * dy * dz); // скорость на грани с учётом монотонизирующей поправки Рхи-Чоу.
					}
					else {
						fw4 = fwplus4 * potent[VXCOR][iW4] + (1.0 - fwplus4) * potent[VXCOR][iP];
					}
				}
				else fw4 = potent[VXCOR][iW4];
			}
		}

		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus*potent[VXCOR][iN] + (1.0 - fnplus)*potent[VXCOR][iP]; else fn = potent[VXCOR][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus*potent[VXCOR][iS] + (1.0 - fsplus)*potent[VXCOR][iP]; else fs = potent[VXCOR][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus*potent[VXCOR][iT] + (1.0 - ftplus)*potent[VXCOR][iP]; else ft = potent[VXCOR][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus*potent[VXCOR][iB] + (1.0 - fbplus)*potent[VXCOR][iP]; else fb = potent[VXCOR][iB];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iN2 > -1) {
				if (!bN2) fn2 = fnplus2 * potent[VXCOR][iN2] + (1.0 - fnplus2) * potent[VXCOR][iP]; else fn2 = potent[VXCOR][iN2];
			}
			if (iS2 > -1) {
				if (!bS2) fs2 = fsplus2 * potent[VXCOR][iS2] + (1.0 - fsplus2) * potent[VXCOR][iP]; else fs2 = potent[VXCOR][iS2];
			}
			if (iT2 > -1) {
				if (!bT2) ft2 = ftplus2 * potent[VXCOR][iT2] + (1.0 - ftplus2) * potent[VXCOR][iP]; else ft2 = potent[VXCOR][iT2];
			}
			if (iB2 > -1) {
				if (!bB2) fb2 = fbplus2 * potent[VXCOR][iB2] + (1.0 - fbplus2) * potent[VXCOR][iP]; else fb2 = potent[VXCOR][iB2];
			}

			if (iN3 > -1) {
				if (!bN3) fn3 = fnplus3 * potent[VXCOR][iN3] + (1.0 - fnplus3) * potent[VXCOR][iP]; else fn3 = potent[VXCOR][iN3];
			}
			if (iS3 > -1) {
				if (!bS3) fs3 = fsplus3 * potent[VXCOR][iS3] + (1.0 - fsplus3) * potent[VXCOR][iP]; else fs3 = potent[VXCOR][iS3];
			}
			if (iT3 > -1) {
				if (!bT3) ft3 = ftplus3 * potent[VXCOR][iT3] + (1.0 - ftplus3) * potent[VXCOR][iP]; else ft3 = potent[VXCOR][iT3];
			}
			if (iB3 > -1) {
				if (!bB3) fb3 = fbplus3 * potent[VXCOR][iB3] + (1.0 - fbplus3) * potent[VXCOR][iP]; else fb3 = potent[VXCOR][iB3];
			}

			if (iN4 > -1) {
				if (!bN4) fn4 = fnplus4 * potent[VXCOR][iN4] + (1.0 - fnplus4) * potent[VXCOR][iP]; else fn4 = potent[VXCOR][iN4];
			}
			if (iS4 > -1) {
				if (!bS4) fs4 = fsplus4 * potent[VXCOR][iS4] + (1.0 - fsplus4) * potent[VXCOR][iP]; else fs4 = potent[VXCOR][iS4];
			}
			if (iT4 > -1) {
				if (!bT4) ft4 = ftplus4 * potent[VXCOR][iT4] + (1.0 - ftplus4) * potent[VXCOR][iP]; else ft4 = potent[VXCOR][iT4];
			}
			if (iB4 > -1) {
				if (!bB4) fb4 = fbplus4 * potent[VXCOR][iB4] + (1.0 - fbplus4) * potent[VXCOR][iP]; else fb4 = potent[VXCOR][iB4];
			}
		}
        // градиент VX
	    //potent[GRADXVX][iP]=(fe-fw)/dx;
	    //potent[GRADYVX][iP]=(fn-fs)/dy;
	    //potent[GRADZVX][iP]=(ft-fb)/dz;
		potent[GRADXVX][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYVX][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZVX][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;


		fe = 0.0; fw = 0.0; fn = 0.0; fs = 0.0; ft = 0.0; fb = 0.0;
		fe2 = 0.0; fw2 = 0.0; fn2 = 0.0; fs2 = 0.0; ft2 = 0.0; fb2 = 0.0;
		fe3 = 0.0; fw3 = 0.0; fn3 = 0.0; fs3 = 0.0; ft3 = 0.0; fb3 = 0.0;
		fe4 = 0.0; fw4 = 0.0; fn4 = 0.0; fs4 = 0.0; ft4 = 0.0; fb4 = 0.0;

		// линейная интерполяция скорости VY на грань КО.
		if (iE > -1) {
			if (!bE) fe = feplus*potent[VYCOR][iE] + (1.0 - feplus)*potent[VYCOR][iP]; else fe = potent[VYCOR][iE];
		}
		if (iW > -1) {
			if (!bW) fw = fwplus*potent[VYCOR][iW] + (1.0 - fwplus)*potent[VYCOR][iP]; else fw = potent[VYCOR][iW];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				if (!bE2) fe2 = feplus2 * potent[VYCOR][iE2] + (1.0 - feplus2) * potent[VYCOR][iP]; else fe2 = potent[VYCOR][iE2];
			}
			if (iW2 > -1) {
				if (!bW2) fw2 = fwplus2 * potent[VYCOR][iW2] + (1.0 - fwplus2) * potent[VYCOR][iP]; else fw2 = potent[VYCOR][iW2];
			}
			if (iE3 > -1) {
				if (!bE3) fe3 = feplus3 * potent[VYCOR][iE3] + (1.0 - feplus3) * potent[VYCOR][iP]; else fe3 = potent[VYCOR][iE3];
			}
			if (iW3 > -1) {
				if (!bW3) fw3 = fwplus3 * potent[VYCOR][iW3] + (1.0 - fwplus3) * potent[VYCOR][iP]; else fw3 = potent[VYCOR][iW3];
			}
			if (iE4 > -1) {
				if (!bE4) fe4 = feplus4 * potent[VYCOR][iE4] + (1.0 - feplus4) * potent[VYCOR][iP]; else fe4 = potent[VYCOR][iE4];
			}
			if (iW4 > -1) {
				if (!bW4) fw4 = fwplus4 * potent[VYCOR][iW4] + (1.0 - fwplus4) * potent[VYCOR][iP]; else fw4 = potent[VYCOR][iW4];
			}

		}

		if (iN > -1) {
			if (!bN) {
				if (bRCh) {
					fn = mf[N_SIDE] / (rhon*dx*dz);
				}
				else {
					fn = fnplus*potent[VYCOR][iN] + (1.0 - fnplus)*potent[VYCOR][iP];
				}
			}
			else fn = potent[VYCOR][iN];
		}
		if (iS > -1) {
			if (!bS) {
				if (bRCh) {
					fs = mf[S_SIDE] / (rhos*dx*dz);
				}
				else {
					fs = fsplus*potent[VYCOR][iS] + (1.0 - fsplus)*potent[VYCOR][iP];
				}
			}
			else fs = potent[VYCOR][iS];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iN2 > -1) {
				if (!bN2) {
					if (bRCh) {
						fn2 = mf[N_SIDE] / (rhon2 * dx * dz);
					}
					else {
						fn2 = fnplus2 * potent[VYCOR][iN2] + (1.0 - fnplus2) * potent[VYCOR][iP];
					}
				}
				else fn2 = potent[VYCOR][iN2];
			}
			if (iS2 > -1) {
				if (!bS2) {
					if (bRCh) {
						fs2 = mf[S_SIDE] / (rhos2 * dx * dz);
					}
					else {
						fs2 = fsplus2 * potent[VYCOR][iS2] + (1.0 - fsplus2) * potent[VYCOR][iP];
					}
				}
				else fs2 = potent[VYCOR][iS2];
			}
			if (iN3 > -1) {
				if (!bN3) {
					if (bRCh) {
						fn3 = mf[N_SIDE] / (rhon3 * dx * dz);
					}
					else {
						fn3 = fnplus3 * potent[VYCOR][iN3] + (1.0 - fnplus3) * potent[VYCOR][iP];
					}
				}
				else fn3 = potent[VYCOR][iN3];
			}
			if (iS3 > -1) {
				if (!bS3) {
					if (bRCh) {
						fs3 = mf[S_SIDE] / (rhos3 * dx * dz);
					}
					else {
						fs3 = fsplus3 * potent[VYCOR][iS3] + (1.0 - fsplus3) * potent[VYCOR][iP];
					}
				}
				else fs3 = potent[VYCOR][iS3];
			}
			if (iN4 > -1) {
				if (!bN4) {
					if (bRCh) {
						fn4 = mf[N_SIDE] / (rhon4 * dx * dz);
					}
					else {
						fn4 = fnplus4 * potent[VYCOR][iN4] + (1.0 - fnplus4) * potent[VYCOR][iP];
					}
				}
				else fn4 = potent[VYCOR][iN4];
			}
			if (iS4 > -1) {
				if (!bS4) {
					if (bRCh) {
						fs4 = mf[S_SIDE] / (rhos4 * dx * dz);
					}
					else {
						fs4 = fsplus4 * potent[VYCOR][iS4] + (1.0 - fsplus4) * potent[VYCOR][iP];
					}
				}
				else fs4 = potent[VYCOR][iS4];
			}
		}

		if (iT > -1) {
			if (!bT) ft = ftplus*potent[VYCOR][iT] + (1.0 - ftplus)*potent[VYCOR][iP]; else ft = potent[VYCOR][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus*potent[VYCOR][iB] + (1.0 - fbplus)*potent[VYCOR][iP]; else fb = potent[VYCOR][iB];
		}

		if (b_on_adaptive_local_refinement_mesh) {
			if (iT2 > -1) {
				if (!bT2) ft2 = ftplus2 * potent[VYCOR][iT2] + (1.0 - ftplus2) * potent[VYCOR][iP]; else ft2 = potent[VYCOR][iT2];
			}
			if (iB2 > -1) {
				if (!bB2) fb2 = fbplus2 * potent[VYCOR][iB2] + (1.0 - fbplus2) * potent[VYCOR][iP]; else fb2 = potent[VYCOR][iB2];
			}
			if (iT3 > -1) {
				if (!bT3) ft3 = ftplus3 * potent[VYCOR][iT3] + (1.0 - ftplus3) * potent[VYCOR][iP]; else ft3 = potent[VYCOR][iT3];
			}
			if (iB3 > -1) {
				if (!bB3) fb3 = fbplus3 * potent[VYCOR][iB3] + (1.0 - fbplus3) * potent[VYCOR][iP]; else fb3 = potent[VYCOR][iB3];
			}
			if (iT4 > -1) {
				if (!bT4) ft4 = ftplus4 * potent[VYCOR][iT4] + (1.0 - ftplus4) * potent[VYCOR][iP]; else ft4 = potent[VYCOR][iT4];
			}
			if (iB4 > -1) {
				if (!bB4) fb4 = fbplus4 * potent[VYCOR][iB4] + (1.0 - fbplus4) * potent[VYCOR][iP]; else fb4 = potent[VYCOR][iB4];
			}
		}
		// градиент VY

	    //potent[GRADXVY][iP]=(fe-fw)/dx;
	    //potent[GRADYVY][iP]=(fn-fs)/dy;
	    //potent[GRADZVY][iP]=(ft-fb)/dz;
		potent[GRADXVY][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYVY][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZVY][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;

		fe = 0.0; fw = 0.0; fn = 0.0; fs = 0.0; ft = 0.0; fb = 0.0;
		fe2 = 0.0; fw2 = 0.0; fn2 = 0.0; fs2 = 0.0; ft2 = 0.0; fb2 = 0.0;
		fe3 = 0.0; fw3 = 0.0; fn3 = 0.0; fs3 = 0.0; ft3 = 0.0; fb3 = 0.0;
		fe4 = 0.0; fw4 = 0.0; fn4 = 0.0; fs4 = 0.0; ft4 = 0.0; fb4 = 0.0;

	    // линейная интерполяция скорости VZ на грань КО.
		if (iE > -1) {
			if (!bE) fe = feplus*potent[VZCOR][iE] + (1.0 - feplus)*potent[VZCOR][iP]; else fe = potent[VZCOR][iE];
		}
		if (iW > -1) {
			if (!bW) fw = fwplus*potent[VZCOR][iW] + (1.0 - fwplus)*potent[VZCOR][iP]; else fw = potent[VZCOR][iW];
		}
		if (iN > -1) {
			if (!bN) fn = fnplus*potent[VZCOR][iN] + (1.0 - fnplus)*potent[VZCOR][iP]; else fn = potent[VZCOR][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus*potent[VZCOR][iS] + (1.0 - fsplus)*potent[VZCOR][iP]; else fs = potent[VZCOR][iS];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				if (!bE2) fe2 = feplus2 * potent[VZCOR][iE2] + (1.0 - feplus2) * potent[VZCOR][iP]; else fe2 = potent[VZCOR][iE2];
			}
			if (iW2 > -1) {
				if (!bW2) fw2 = fwplus2 * potent[VZCOR][iW2] + (1.0 - fwplus2) * potent[VZCOR][iP]; else fw2 = potent[VZCOR][iW2];
			}
			if (iN2 > -1) {
				if (!bN2) fn2 = fnplus2 * potent[VZCOR][iN2] + (1.0 - fnplus2) * potent[VZCOR][iP]; else fn2 = potent[VZCOR][iN2];
			}
			if (iS2 > -1) {
				if (!bS2) fs2 = fsplus2 * potent[VZCOR][iS2] + (1.0 - fsplus2) * potent[VZCOR][iP]; else fs2 = potent[VZCOR][iS2];
			}
			if (iE3 > -1) {
				if (!bE3) fe3 = feplus3 * potent[VZCOR][iE3] + (1.0 - feplus3) * potent[VZCOR][iP]; else fe3 = potent[VZCOR][iE3];
			}
			if (iW3 > -1) {
				if (!bW3) fw3 = fwplus3 * potent[VZCOR][iW3] + (1.0 - fwplus3) * potent[VZCOR][iP]; else fw3 = potent[VZCOR][iW3];
			}
			if (iN3 > -1) {
				if (!bN3) fn3 = fnplus3 * potent[VZCOR][iN3] + (1.0 - fnplus3) * potent[VZCOR][iP]; else fn3 = potent[VZCOR][iN3];
			}
			if (iS3 > -1) {
				if (!bS3) fs3 = fsplus3 * potent[VZCOR][iS3] + (1.0 - fsplus3) * potent[VZCOR][iP]; else fs3 = potent[VZCOR][iS3];
			}
			if (iE4 > -1) {
				if (!bE4) fe4 = feplus4 * potent[VZCOR][iE4] + (1.0 - feplus4) * potent[VZCOR][iP]; else fe4 = potent[VZCOR][iE4];
			}
			if (iW4 > -1) {
				if (!bW4) fw4 = fwplus4 * potent[VZCOR][iW4] + (1.0 - fwplus4) * potent[VZCOR][iP]; else fw4 = potent[VZCOR][iW4];
			}
			if (iN4 > -1) {
				if (!bN4) fn4 = fnplus4 * potent[VZCOR][iN4] + (1.0 - fnplus4) * potent[VZCOR][iP]; else fn4 = potent[VZCOR][iN4];
			}
			if (iS4 > -1) {
				if (!bS4) fs4 = fsplus4 * potent[VZCOR][iS4] + (1.0 - fsplus4) * potent[VZCOR][iP]; else fs4 = potent[VZCOR][iS4];
			}
		}

		if (iT > -1) {
			if (!bT) {
				if (bRCh) {
					ft = mf[T_SIDE] / (rhot*dx*dy);
				}
				else {
					ft = ftplus*potent[VZCOR][iT] + (1.0 - ftplus)*potent[VZCOR][iP];
				}
			}
			else ft = potent[VZCOR][iT];
		}
		if (iB > -1) {
			if (!bB) {
				if (bRCh) {
					fb = mf[B_SIDE] / (rhob*dx*dy);
				}
				else {
					fb = fbplus*potent[VZCOR][iB] + (1.0 - fbplus)*potent[VZCOR][iP];
				}
			}
			else fb = potent[VZCOR][iB];
		}

		if (b_on_adaptive_local_refinement_mesh) {

			if (iT2 > -1) {
				if (!bT2) {
					if (bRCh) {
						ft2 = mf[T_SIDE] / (rhot2 * dx * dy);
					}
					else {
						ft2 = ftplus2 * potent[VZCOR][iT2] + (1.0 - ftplus2) * potent[VZCOR][iP];
					}
				}
				else ft2 = potent[VZCOR][iT2];
			}
			if (iB2 > -1) {
				if (!bB2) {
					if (bRCh) {
						fb2 = mf[B_SIDE] / (rhob2 * dx * dy);
					}
					else {
						fb2 = fbplus2 * potent[VZCOR][iB2] + (1.0 - fbplus2) * potent[VZCOR][iP];
					}
				}
				else fb2 = potent[VZCOR][iB2];
			}
			if (iT3 > -1) {
				if (!bT3) {
					if (bRCh) {
						ft3 = mf[T_SIDE] / (rhot3 * dx * dy);
					}
					else {
						ft3 = ftplus3 * potent[VZCOR][iT3] + (1.0 - ftplus3) * potent[VZCOR][iP];
					}
				}
				else ft3 = potent[VZCOR][iT3];
			}
			if (iB3 > -1) {
				if (!bB3) {
					if (bRCh) {
						fb3 = mf[B_SIDE] / (rhob3 * dx * dy);
					}
					else {
						fb3 = fbplus3 * potent[VZCOR][iB3] + (1.0 - fbplus3) * potent[VZCOR][iP];
					}
				}
				else fb3 = potent[VZCOR][iB3];
			}
			if (iT4 > -1) {
				if (!bT4) {
					if (bRCh) {
						ft4 = mf[T_SIDE] / (rhot4 * dx * dy);
					}
					else {
						ft4 = ftplus4 * potent[VZCOR][iT4] + (1.0 - ftplus4) * potent[VZCOR][iP];
					}
				}
				else ft4 = potent[VZCOR][iT4];
			}
			if (iB4 > -1) {
				if (!bB4) {
					if (bRCh) {
						fb4 = mf[B_SIDE] / (rhob4 * dx * dy);
					}
					else {
						fb4 = fbplus4 * potent[VZCOR][iB4] + (1.0 - fbplus4) * potent[VZCOR][iP];
					}
				}
				else fb4 = potent[VZCOR][iB4];
			}
		}

		// градиент VZ
	    //potent[GRADXVZ][iP]=(fe-fw)/dx;
	    //potent[GRADYVZ][iP]=(fn-fs)/dy;
	    //potent[GRADZVZ][iP]=(ft-fb)/dz;
		potent[GRADXVZ][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYVZ][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZVZ][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;
	}
	else {

		if (1) {


			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE])*(potent[VXCOR][iE]) + (potent[VYCOR][iE])*(potent[VYCOR][iE]) + (potent[VZCOR][iE])*(potent[VZCOR][iE]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iE] = 0.0;
						potent[GRADYVX][iE] = 0.0;
						potent[GRADZVX][iE] = 0.0;

						potent[GRADXVY][iE] = 0.0;
						potent[GRADYVY][iE] = 0.0;
						potent[GRADZVY][iE] = 0.0;

						potent[GRADXVZ][iE] = 0.0;
						potent[GRADYVZ][iE] = 0.0;
						potent[GRADZVZ][iE] = 0.0;

						/*
						potent[GRADXVX][iP] = 0.0;
						potent[GRADYVX][iP] = 0.0;
						potent[GRADZVX][iP] = 0.0;

						potent[GRADXVY][iP] = 0.0;
						potent[GRADYVY][iP] = 0.0;
						potent[GRADZVY][iP] = 0.0;

						potent[GRADXVZ][iP] = 0.0;
						potent[GRADYVZ][iP] = 0.0;
						potent[GRADZVZ][iP] = 0.0;
						*/
					}
					else {

						potent[GRADXVX][iE] = potent[GRADXVX][iP];
						potent[GRADYVX][iE] = potent[GRADYVX][iP];
						potent[GRADZVX][iE] = potent[GRADZVX][iP];

						potent[GRADXVY][iE] = potent[GRADXVY][iP];
						potent[GRADYVY][iE] = potent[GRADYVY][iP];
						potent[GRADZVY][iE] = potent[GRADZVY][iP];

						potent[GRADXVZ][iE] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iE] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iE] = potent[GRADZVZ][iP];
					}
				}
			}
			if (iW > -1) {
				if (bW) {

					doublereal dspeed = sqrt((potent[VXCOR][iW]) * (potent[VXCOR][iW]) + (potent[VYCOR][iW]) * (potent[VYCOR][iW]) + (potent[VZCOR][iW]) * (potent[VZCOR][iW]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iW] = 0.0;
						potent[GRADYVX][iW] = 0.0;
						potent[GRADZVX][iW] = 0.0;

						potent[GRADXVY][iW] = 0.0;
						potent[GRADYVY][iW] = 0.0;
						potent[GRADZVY][iW] = 0.0;

						potent[GRADXVZ][iW] = 0.0;
						potent[GRADYVZ][iW] = 0.0;
						potent[GRADZVZ][iW] = 0.0;
					}
					else {

						potent[GRADXVX][iW] = potent[GRADXVX][iP];
						potent[GRADYVX][iW] = potent[GRADYVX][iP];
						potent[GRADZVX][iW] = potent[GRADZVX][iP];

						potent[GRADXVY][iW] = potent[GRADXVY][iP];
						potent[GRADYVY][iW] = potent[GRADYVY][iP];
						potent[GRADZVY][iW] = potent[GRADZVY][iP];

						potent[GRADXVZ][iW] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iW] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iW] = potent[GRADZVZ][iP];
					}
				}
			}
			if (iN > -1) {
				if (bN) {

					doublereal dspeed = sqrt((potent[VXCOR][iN]) * (potent[VXCOR][iN]) + (potent[VYCOR][iN]) * (potent[VYCOR][iN]) + (potent[VZCOR][iN]) * (potent[VZCOR][iN]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iN] = 0.0;
						potent[GRADYVX][iN] = 0.0;
						potent[GRADZVX][iN] = 0.0;

						potent[GRADXVY][iN] = 0.0;
						potent[GRADYVY][iN] = 0.0;
						potent[GRADZVY][iN] = 0.0;

						potent[GRADXVZ][iN] = 0.0;
						potent[GRADYVZ][iN] = 0.0;
						potent[GRADZVZ][iN] = 0.0;
					}
					else {

						potent[GRADXVX][iN] = potent[GRADXVX][iP];
						potent[GRADYVX][iN] = potent[GRADYVX][iP];
						potent[GRADZVX][iN] = potent[GRADZVX][iP];

						potent[GRADXVY][iN] = potent[GRADXVY][iP];
						potent[GRADYVY][iN] = potent[GRADYVY][iP];
						potent[GRADZVY][iN] = potent[GRADZVY][iP];

						potent[GRADXVZ][iN] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iN] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iN] = potent[GRADZVZ][iP];
					}
				}
			}
			if (iS > -1) {
				if (bS) {

					doublereal dspeed = sqrt((potent[VXCOR][iS]) * (potent[VXCOR][iS]) + (potent[VYCOR][iS]) * (potent[VYCOR][iS]) + (potent[VZCOR][iS]) * (potent[VZCOR][iS]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iS] = 0.0;
						potent[GRADYVX][iS] = 0.0;
						potent[GRADZVX][iS] = 0.0;

						potent[GRADXVY][iS] = 0.0;
						potent[GRADYVY][iS] = 0.0;
						potent[GRADZVY][iS] = 0.0;

						potent[GRADXVZ][iS] = 0.0;
						potent[GRADYVZ][iS] = 0.0;
						potent[GRADZVZ][iS] = 0.0;
					}
					else {

						potent[GRADXVX][iS] = potent[GRADXVX][iP];
						potent[GRADYVX][iS] = potent[GRADYVX][iP];
						potent[GRADZVX][iS] = potent[GRADZVX][iP];

						potent[GRADXVY][iS] = potent[GRADXVY][iP];
						potent[GRADYVY][iS] = potent[GRADYVY][iP];
						potent[GRADZVY][iS] = potent[GRADZVY][iP];

						potent[GRADXVZ][iS] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iS] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iS] = potent[GRADZVZ][iP];
					}
				}
			}
			if (iT > -1) {
				if (bT) {

					doublereal dspeed = sqrt((potent[VXCOR][iT]) * (potent[VXCOR][iT]) + (potent[VYCOR][iT]) * (potent[VYCOR][iT]) + (potent[VZCOR][iT]) * (potent[VZCOR][iT]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iT] = 0.0;
						potent[GRADYVX][iT] = 0.0;
						potent[GRADZVX][iT] = 0.0;

						potent[GRADXVY][iT] = 0.0;
						potent[GRADYVY][iT] = 0.0;
						potent[GRADZVY][iT] = 0.0;

						potent[GRADXVZ][iT] = 0.0;
						potent[GRADYVZ][iT] = 0.0;
						potent[GRADZVZ][iT] = 0.0;
					}
					else {

						potent[GRADXVX][iT] = potent[GRADXVX][iP];
						potent[GRADYVX][iT] = potent[GRADYVX][iP];
						potent[GRADZVX][iT] = potent[GRADZVX][iP];

						potent[GRADXVY][iT] = potent[GRADXVY][iP];
						potent[GRADYVY][iT] = potent[GRADYVY][iP];
						potent[GRADZVY][iT] = potent[GRADZVY][iP];

						potent[GRADXVZ][iT] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iT] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iT] = potent[GRADZVZ][iP];
					}
				}
			}
			if (iB > -1) {
				if (bB) {
					doublereal dspeed = sqrt((potent[VXCOR][iB]) * (potent[VXCOR][iB]) + (potent[VYCOR][iB]) * (potent[VYCOR][iB]) + (potent[VZCOR][iB]) * (potent[VZCOR][iB]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iB] = 0.0;
						potent[GRADYVX][iB] = 0.0;
						potent[GRADZVX][iB] = 0.0;

						potent[GRADXVY][iB] = 0.0;
						potent[GRADYVY][iB] = 0.0;
						potent[GRADZVY][iB] = 0.0;

						potent[GRADXVZ][iB] = 0.0;
						potent[GRADYVZ][iB] = 0.0;
						potent[GRADZVZ][iB] = 0.0;
					}
					else {

						potent[GRADXVX][iB] = potent[GRADXVX][iP];
						potent[GRADYVX][iB] = potent[GRADYVX][iP];
						potent[GRADZVX][iB] = potent[GRADZVX][iP];

						potent[GRADXVY][iB] = potent[GRADXVY][iP];
						potent[GRADYVY][iB] = potent[GRADYVY][iP];
						potent[GRADZVY][iB] = potent[GRADZVY][iP];

						potent[GRADXVZ][iB] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iB] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iB] = potent[GRADZVZ][iP];
					}
				}
			}

			if (b_on_adaptive_local_refinement_mesh) {
				if (iE2 > -1) {
					if (bE2) {

						// 10.02.2017
						// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

						doublereal dspeed = sqrt((potent[VXCOR][iE2]) * (potent[VXCOR][iE2]) + (potent[VYCOR][iE2]) * (potent[VYCOR][iE2]) + (potent[VZCOR][iE2]) * (potent[VZCOR][iE2]));

						if (dspeed < 1.0e-10) {
							potent[GRADXVX][iE2] = 0.0;
							potent[GRADYVX][iE2] = 0.0;
							potent[GRADZVX][iE2] = 0.0;

							potent[GRADXVY][iE2] = 0.0;
							potent[GRADYVY][iE2] = 0.0;
							potent[GRADZVY][iE2] = 0.0;

							potent[GRADXVZ][iE2] = 0.0;
							potent[GRADYVZ][iE2] = 0.0;
							potent[GRADZVZ][iE2] = 0.0;


						}
						else {

							potent[GRADXVX][iE2] = potent[GRADXVX][iP];
							potent[GRADYVX][iE2] = potent[GRADYVX][iP];
							potent[GRADZVX][iE2] = potent[GRADZVX][iP];

							potent[GRADXVY][iE2] = potent[GRADXVY][iP];
							potent[GRADYVY][iE2] = potent[GRADYVY][iP];
							potent[GRADZVY][iE2] = potent[GRADZVY][iP];

							potent[GRADXVZ][iE2] = potent[GRADXVZ][iP];
							potent[GRADYVZ][iE2] = potent[GRADYVZ][iP];
							potent[GRADZVZ][iE2] = potent[GRADZVZ][iP];
						}
					}
				}

				if (iE3 > -1) {
					if (bE3) {

						// 10.02.2017
						// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

						doublereal dspeed = sqrt((potent[VXCOR][iE3]) * (potent[VXCOR][iE3]) + (potent[VYCOR][iE3]) * (potent[VYCOR][iE3]) + (potent[VZCOR][iE3]) * (potent[VZCOR][iE3]));

						if (dspeed < 1.0e-10) {
							potent[GRADXVX][iE3] = 0.0;
							potent[GRADYVX][iE3] = 0.0;
							potent[GRADZVX][iE3] = 0.0;

							potent[GRADXVY][iE3] = 0.0;
							potent[GRADYVY][iE3] = 0.0;
							potent[GRADZVY][iE3] = 0.0;

							potent[GRADXVZ][iE3] = 0.0;
							potent[GRADYVZ][iE3] = 0.0;
							potent[GRADZVZ][iE3] = 0.0;


						}
						else {

							potent[GRADXVX][iE3] = potent[GRADXVX][iP];
							potent[GRADYVX][iE3] = potent[GRADYVX][iP];
							potent[GRADZVX][iE3] = potent[GRADZVX][iP];

							potent[GRADXVY][iE3] = potent[GRADXVY][iP];
							potent[GRADYVY][iE3] = potent[GRADYVY][iP];
							potent[GRADZVY][iE3] = potent[GRADZVY][iP];

							potent[GRADXVZ][iE3] = potent[GRADXVZ][iP];
							potent[GRADYVZ][iE3] = potent[GRADYVZ][iP];
							potent[GRADZVZ][iE3] = potent[GRADZVZ][iP];
						}
					}
				}

				if (iE4 > -1) {
					if (bE4) {

						// 10.02.2017
						// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

						doublereal dspeed = sqrt((potent[VXCOR][iE4]) * (potent[VXCOR][iE4]) + (potent[VYCOR][iE4]) * (potent[VYCOR][iE4]) + (potent[VZCOR][iE4]) * (potent[VZCOR][iE4]));

						if (dspeed < 1.0e-10) {
							potent[GRADXVX][iE4] = 0.0;
							potent[GRADYVX][iE4] = 0.0;
							potent[GRADZVX][iE4] = 0.0;

							potent[GRADXVY][iE4] = 0.0;
							potent[GRADYVY][iE4] = 0.0;
							potent[GRADZVY][iE4] = 0.0;

							potent[GRADXVZ][iE4] = 0.0;
							potent[GRADYVZ][iE4] = 0.0;
							potent[GRADZVZ][iE4] = 0.0;


						}
						else {

							potent[GRADXVX][iE4] = potent[GRADXVX][iP];
							potent[GRADYVX][iE4] = potent[GRADYVX][iP];
							potent[GRADZVX][iE4] = potent[GRADZVX][iP];

							potent[GRADXVY][iE4] = potent[GRADXVY][iP];
							potent[GRADYVY][iE4] = potent[GRADYVY][iP];
							potent[GRADZVY][iE4] = potent[GRADZVY][iP];

							potent[GRADXVZ][iE4] = potent[GRADXVZ][iP];
							potent[GRADYVZ][iE4] = potent[GRADYVZ][iP];
							potent[GRADZVZ][iE4] = potent[GRADZVZ][iP];
						}
					}
				}
			

			if (iW2 > -1) {
				if (bW2) {

					doublereal dspeed = sqrt((potent[VXCOR][iW2]) * (potent[VXCOR][iW2]) + (potent[VYCOR][iW2]) * (potent[VYCOR][iW2]) + (potent[VZCOR][iW2]) * (potent[VZCOR][iW2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iW2] = 0.0;
						potent[GRADYVX][iW2] = 0.0;
						potent[GRADZVX][iW2] = 0.0;

						potent[GRADXVY][iW2] = 0.0;
						potent[GRADYVY][iW2] = 0.0;
						potent[GRADZVY][iW2] = 0.0;

						potent[GRADXVZ][iW2] = 0.0;
						potent[GRADYVZ][iW2] = 0.0;
						potent[GRADZVZ][iW2] = 0.0;
					}
					else {

						potent[GRADXVX][iW2] = potent[GRADXVX][iP];
						potent[GRADYVX][iW2] = potent[GRADYVX][iP];
						potent[GRADZVX][iW2] = potent[GRADZVX][iP];

						potent[GRADXVY][iW2] = potent[GRADXVY][iP];
						potent[GRADYVY][iW2] = potent[GRADYVY][iP];
						potent[GRADZVY][iW2] = potent[GRADZVY][iP];

						potent[GRADXVZ][iW2] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iW2] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iW2] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					doublereal dspeed = sqrt((potent[VXCOR][iW3]) * (potent[VXCOR][iW3]) + (potent[VYCOR][iW3]) * (potent[VYCOR][iW3]) + (potent[VZCOR][iW3]) * (potent[VZCOR][iW3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iW3] = 0.0;
						potent[GRADYVX][iW3] = 0.0;
						potent[GRADZVX][iW3] = 0.0;

						potent[GRADXVY][iW3] = 0.0;
						potent[GRADYVY][iW3] = 0.0;
						potent[GRADZVY][iW3] = 0.0;

						potent[GRADXVZ][iW3] = 0.0;
						potent[GRADYVZ][iW3] = 0.0;
						potent[GRADZVZ][iW3] = 0.0;
					}
					else {

						potent[GRADXVX][iW3] = potent[GRADXVX][iP];
						potent[GRADYVX][iW3] = potent[GRADYVX][iP];
						potent[GRADZVX][iW3] = potent[GRADZVX][iP];

						potent[GRADXVY][iW3] = potent[GRADXVY][iP];
						potent[GRADYVY][iW3] = potent[GRADYVY][iP];
						potent[GRADZVY][iW3] = potent[GRADZVY][iP];

						potent[GRADXVZ][iW3] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iW3] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iW3] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					doublereal dspeed = sqrt((potent[VXCOR][iW4]) * (potent[VXCOR][iW4]) + (potent[VYCOR][iW4]) * (potent[VYCOR][iW4]) + (potent[VZCOR][iW4]) * (potent[VZCOR][iW4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iW4] = 0.0;
						potent[GRADYVX][iW4] = 0.0;
						potent[GRADZVX][iW4] = 0.0;

						potent[GRADXVY][iW4] = 0.0;
						potent[GRADYVY][iW4] = 0.0;
						potent[GRADZVY][iW4] = 0.0;

						potent[GRADXVZ][iW4] = 0.0;
						potent[GRADYVZ][iW4] = 0.0;
						potent[GRADZVZ][iW4] = 0.0;
					}
					else {

						potent[GRADXVX][iW4] = potent[GRADXVX][iP];
						potent[GRADYVX][iW4] = potent[GRADYVX][iP];
						potent[GRADZVX][iW4] = potent[GRADZVX][iP];

						potent[GRADXVY][iW4] = potent[GRADXVY][iP];
						potent[GRADYVY][iW4] = potent[GRADYVY][iP];
						potent[GRADZVY][iW4] = potent[GRADZVY][iP];

						potent[GRADXVZ][iW4] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iW4] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iW4] = potent[GRADZVZ][iP];
					}
				}
			}


			if (iN2 > -1) {
				if (bN2) {

					doublereal dspeed = sqrt((potent[VXCOR][iN2]) * (potent[VXCOR][iN2]) + (potent[VYCOR][iN2]) * (potent[VYCOR][iN2]) + (potent[VZCOR][iN2]) * (potent[VZCOR][iN2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iN2] = 0.0;
						potent[GRADYVX][iN2] = 0.0;
						potent[GRADZVX][iN2] = 0.0;

						potent[GRADXVY][iN2] = 0.0;
						potent[GRADYVY][iN2] = 0.0;
						potent[GRADZVY][iN2] = 0.0;

						potent[GRADXVZ][iN2] = 0.0;
						potent[GRADYVZ][iN2] = 0.0;
						potent[GRADZVZ][iN2] = 0.0;
					}
					else {

						potent[GRADXVX][iN2] = potent[GRADXVX][iP];
						potent[GRADYVX][iN2] = potent[GRADYVX][iP];
						potent[GRADZVX][iN2] = potent[GRADZVX][iP];

						potent[GRADXVY][iN2] = potent[GRADXVY][iP];
						potent[GRADYVY][iN2] = potent[GRADYVY][iP];
						potent[GRADZVY][iN2] = potent[GRADZVY][iP];

						potent[GRADXVZ][iN2] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iN2] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iN2] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					doublereal dspeed = sqrt((potent[VXCOR][iN3]) * (potent[VXCOR][iN3]) + (potent[VYCOR][iN3]) * (potent[VYCOR][iN3]) + (potent[VZCOR][iN3]) * (potent[VZCOR][iN3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iN3] = 0.0;
						potent[GRADYVX][iN3] = 0.0;
						potent[GRADZVX][iN3] = 0.0;

						potent[GRADXVY][iN3] = 0.0;
						potent[GRADYVY][iN3] = 0.0;
						potent[GRADZVY][iN3] = 0.0;

						potent[GRADXVZ][iN3] = 0.0;
						potent[GRADYVZ][iN3] = 0.0;
						potent[GRADZVZ][iN3] = 0.0;
					}
					else {

						potent[GRADXVX][iN3] = potent[GRADXVX][iP];
						potent[GRADYVX][iN3] = potent[GRADYVX][iP];
						potent[GRADZVX][iN3] = potent[GRADZVX][iP];

						potent[GRADXVY][iN3] = potent[GRADXVY][iP];
						potent[GRADYVY][iN3] = potent[GRADYVY][iP];
						potent[GRADZVY][iN3] = potent[GRADZVY][iP];

						potent[GRADXVZ][iN3] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iN3] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iN3] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					doublereal dspeed = sqrt((potent[VXCOR][iN4]) * (potent[VXCOR][iN4]) + (potent[VYCOR][iN4]) * (potent[VYCOR][iN4]) + (potent[VZCOR][iN4]) * (potent[VZCOR][iN4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iN4] = 0.0;
						potent[GRADYVX][iN4] = 0.0;
						potent[GRADZVX][iN4] = 0.0;

						potent[GRADXVY][iN4] = 0.0;
						potent[GRADYVY][iN4] = 0.0;
						potent[GRADZVY][iN4] = 0.0;

						potent[GRADXVZ][iN4] = 0.0;
						potent[GRADYVZ][iN4] = 0.0;
						potent[GRADZVZ][iN4] = 0.0;
					}
					else {

						potent[GRADXVX][iN4] = potent[GRADXVX][iP];
						potent[GRADYVX][iN4] = potent[GRADYVX][iP];
						potent[GRADZVX][iN4] = potent[GRADZVX][iP];

						potent[GRADXVY][iN4] = potent[GRADXVY][iP];
						potent[GRADYVY][iN4] = potent[GRADYVY][iP];
						potent[GRADZVY][iN4] = potent[GRADZVY][iP];

						potent[GRADXVZ][iN4] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iN4] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iN4] = potent[GRADZVZ][iP];
					}
				}
			}
		

			if (iS2 > -1) {
				if (bS2) {

					doublereal dspeed = sqrt((potent[VXCOR][iS2]) * (potent[VXCOR][iS2]) + (potent[VYCOR][iS2]) * (potent[VYCOR][iS2]) + (potent[VZCOR][iS2]) * (potent[VZCOR][iS2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iS2] = 0.0;
						potent[GRADYVX][iS2] = 0.0;
						potent[GRADZVX][iS2] = 0.0;

						potent[GRADXVY][iS2] = 0.0;
						potent[GRADYVY][iS2] = 0.0;
						potent[GRADZVY][iS2] = 0.0;

						potent[GRADXVZ][iS2] = 0.0;
						potent[GRADYVZ][iS2] = 0.0;
						potent[GRADZVZ][iS2] = 0.0;
					}
					else {

						potent[GRADXVX][iS2] = potent[GRADXVX][iP];
						potent[GRADYVX][iS2] = potent[GRADYVX][iP];
						potent[GRADZVX][iS2] = potent[GRADZVX][iP];

						potent[GRADXVY][iS2] = potent[GRADXVY][iP];
						potent[GRADYVY][iS2] = potent[GRADYVY][iP];
						potent[GRADZVY][iS2] = potent[GRADZVY][iP];

						potent[GRADXVZ][iS2] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iS2] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iS2] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					doublereal dspeed = sqrt((potent[VXCOR][iS3]) * (potent[VXCOR][iS3]) + (potent[VYCOR][iS3]) * (potent[VYCOR][iS3]) + (potent[VZCOR][iS3]) * (potent[VZCOR][iS3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iS3] = 0.0;
						potent[GRADYVX][iS3] = 0.0;
						potent[GRADZVX][iS3] = 0.0;

						potent[GRADXVY][iS3] = 0.0;
						potent[GRADYVY][iS3] = 0.0;
						potent[GRADZVY][iS3] = 0.0;

						potent[GRADXVZ][iS3] = 0.0;
						potent[GRADYVZ][iS3] = 0.0;
						potent[GRADZVZ][iS3] = 0.0;
					}
					else {

						potent[GRADXVX][iS3] = potent[GRADXVX][iP];
						potent[GRADYVX][iS3] = potent[GRADYVX][iP];
						potent[GRADZVX][iS3] = potent[GRADZVX][iP];

						potent[GRADXVY][iS3] = potent[GRADXVY][iP];
						potent[GRADYVY][iS3] = potent[GRADYVY][iP];
						potent[GRADZVY][iS3] = potent[GRADZVY][iP];

						potent[GRADXVZ][iS3] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iS3] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iS3] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					doublereal dspeed = sqrt((potent[VXCOR][iS4]) * (potent[VXCOR][iS4]) + (potent[VYCOR][iS4]) * (potent[VYCOR][iS4]) + (potent[VZCOR][iS4]) * (potent[VZCOR][iS4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iS4] = 0.0;
						potent[GRADYVX][iS4] = 0.0;
						potent[GRADZVX][iS4] = 0.0;

						potent[GRADXVY][iS4] = 0.0;
						potent[GRADYVY][iS4] = 0.0;
						potent[GRADZVY][iS4] = 0.0;

						potent[GRADXVZ][iS4] = 0.0;
						potent[GRADYVZ][iS4] = 0.0;
						potent[GRADZVZ][iS4] = 0.0;
					}
					else {

						potent[GRADXVX][iS4] = potent[GRADXVX][iP];
						potent[GRADYVX][iS4] = potent[GRADYVX][iP];
						potent[GRADZVX][iS4] = potent[GRADZVX][iP];

						potent[GRADXVY][iS4] = potent[GRADXVY][iP];
						potent[GRADYVY][iS4] = potent[GRADYVY][iP];
						potent[GRADZVY][iS4] = potent[GRADZVY][iP];

						potent[GRADXVZ][iS4] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iS4] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iS4] = potent[GRADZVZ][iP];
					}
				}
			}


			if (iT2 > -1) {
				if (bT2) {

					doublereal dspeed = sqrt((potent[VXCOR][iT2]) * (potent[VXCOR][iT2]) + (potent[VYCOR][iT2]) * (potent[VYCOR][iT2]) + (potent[VZCOR][iT2]) * (potent[VZCOR][iT2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iT2] = 0.0;
						potent[GRADYVX][iT2] = 0.0;
						potent[GRADZVX][iT2] = 0.0;

						potent[GRADXVY][iT2] = 0.0;
						potent[GRADYVY][iT2] = 0.0;
						potent[GRADZVY][iT2] = 0.0;

						potent[GRADXVZ][iT2] = 0.0;
						potent[GRADYVZ][iT2] = 0.0;
						potent[GRADZVZ][iT2] = 0.0;
					}
					else {

						potent[GRADXVX][iT2] = potent[GRADXVX][iP];
						potent[GRADYVX][iT2] = potent[GRADYVX][iP];
						potent[GRADZVX][iT2] = potent[GRADZVX][iP];

						potent[GRADXVY][iT2] = potent[GRADXVY][iP];
						potent[GRADYVY][iT2] = potent[GRADYVY][iP];
						potent[GRADZVY][iT2] = potent[GRADZVY][iP];

						potent[GRADXVZ][iT2] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iT2] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iT2] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					doublereal dspeed = sqrt((potent[VXCOR][iT3]) * (potent[VXCOR][iT3]) + (potent[VYCOR][iT3]) * (potent[VYCOR][iT3]) + (potent[VZCOR][iT3]) * (potent[VZCOR][iT3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iT3] = 0.0;
						potent[GRADYVX][iT3] = 0.0;
						potent[GRADZVX][iT3] = 0.0;

						potent[GRADXVY][iT3] = 0.0;
						potent[GRADYVY][iT3] = 0.0;
						potent[GRADZVY][iT3] = 0.0;

						potent[GRADXVZ][iT3] = 0.0;
						potent[GRADYVZ][iT3] = 0.0;
						potent[GRADZVZ][iT3] = 0.0;
					}
					else {

						potent[GRADXVX][iT3] = potent[GRADXVX][iP];
						potent[GRADYVX][iT3] = potent[GRADYVX][iP];
						potent[GRADZVX][iT3] = potent[GRADZVX][iP];

						potent[GRADXVY][iT3] = potent[GRADXVY][iP];
						potent[GRADYVY][iT3] = potent[GRADYVY][iP];
						potent[GRADZVY][iT3] = potent[GRADZVY][iP];

						potent[GRADXVZ][iT3] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iT3] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iT3] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					doublereal dspeed = sqrt((potent[VXCOR][iT4]) * (potent[VXCOR][iT4]) + (potent[VYCOR][iT4]) * (potent[VYCOR][iT4]) + (potent[VZCOR][iT4]) * (potent[VZCOR][iT4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iT4] = 0.0;
						potent[GRADYVX][iT4] = 0.0;
						potent[GRADZVX][iT4] = 0.0;

						potent[GRADXVY][iT4] = 0.0;
						potent[GRADYVY][iT4] = 0.0;
						potent[GRADZVY][iT4] = 0.0;

						potent[GRADXVZ][iT4] = 0.0;
						potent[GRADYVZ][iT4] = 0.0;
						potent[GRADZVZ][iT4] = 0.0;
					}
					else {

						potent[GRADXVX][iT4] = potent[GRADXVX][iP];
						potent[GRADYVX][iT4] = potent[GRADYVX][iP];
						potent[GRADZVX][iT4] = potent[GRADZVX][iP];

						potent[GRADXVY][iT4] = potent[GRADXVY][iP];
						potent[GRADYVY][iT4] = potent[GRADYVY][iP];
						potent[GRADZVY][iT4] = potent[GRADZVY][iP];

						potent[GRADXVZ][iT4] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iT4] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iT4] = potent[GRADZVZ][iP];
					}
				}
			}


			if (iB2 > -1) {
				if (bB2) {
					doublereal dspeed = sqrt((potent[VXCOR][iB2]) * (potent[VXCOR][iB2]) + (potent[VYCOR][iB2]) * (potent[VYCOR][iB2]) + (potent[VZCOR][iB2]) * (potent[VZCOR][iB2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iB2] = 0.0;
						potent[GRADYVX][iB2] = 0.0;
						potent[GRADZVX][iB2] = 0.0;

						potent[GRADXVY][iB2] = 0.0;
						potent[GRADYVY][iB2] = 0.0;
						potent[GRADZVY][iB2] = 0.0;

						potent[GRADXVZ][iB2] = 0.0;
						potent[GRADYVZ][iB2] = 0.0;
						potent[GRADZVZ][iB2] = 0.0;
					}
					else {

						potent[GRADXVX][iB2] = potent[GRADXVX][iP];
						potent[GRADYVX][iB2] = potent[GRADYVX][iP];
						potent[GRADZVX][iB2] = potent[GRADZVX][iP];

						potent[GRADXVY][iB2] = potent[GRADXVY][iP];
						potent[GRADYVY][iB2] = potent[GRADYVY][iP];
						potent[GRADZVY][iB2] = potent[GRADZVY][iP];

						potent[GRADXVZ][iB2] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iB2] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iB2] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					doublereal dspeed = sqrt((potent[VXCOR][iB3]) * (potent[VXCOR][iB3]) + (potent[VYCOR][iB3]) * (potent[VYCOR][iB3]) + (potent[VZCOR][iB3]) * (potent[VZCOR][iB3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iB3] = 0.0;
						potent[GRADYVX][iB3] = 0.0;
						potent[GRADZVX][iB3] = 0.0;

						potent[GRADXVY][iB3] = 0.0;
						potent[GRADYVY][iB3] = 0.0;
						potent[GRADZVY][iB3] = 0.0;

						potent[GRADXVZ][iB3] = 0.0;
						potent[GRADYVZ][iB3] = 0.0;
						potent[GRADZVZ][iB3] = 0.0;
					}
					else {

						potent[GRADXVX][iB3] = potent[GRADXVX][iP];
						potent[GRADYVX][iB3] = potent[GRADYVX][iP];
						potent[GRADZVX][iB3] = potent[GRADZVX][iP];

						potent[GRADXVY][iB3] = potent[GRADXVY][iP];
						potent[GRADYVY][iB3] = potent[GRADYVY][iP];
						potent[GRADZVY][iB3] = potent[GRADZVY][iP];

						potent[GRADXVZ][iB3] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iB3] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iB3] = potent[GRADZVZ][iP];
					}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					doublereal dspeed = sqrt((potent[VXCOR][iB4]) * (potent[VXCOR][iB4]) + (potent[VYCOR][iB4]) * (potent[VYCOR][iB4]) + (potent[VZCOR][iB4]) * (potent[VZCOR][iB4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXVX][iB4] = 0.0;
						potent[GRADYVX][iB4] = 0.0;
						potent[GRADZVX][iB4] = 0.0;

						potent[GRADXVY][iB4] = 0.0;
						potent[GRADYVY][iB4] = 0.0;
						potent[GRADZVY][iB4] = 0.0;

						potent[GRADXVZ][iB4] = 0.0;
						potent[GRADYVZ][iB4] = 0.0;
						potent[GRADZVZ][iB4] = 0.0;
					}
					else {

						potent[GRADXVX][iB4] = potent[GRADXVX][iP];
						potent[GRADYVX][iB4] = potent[GRADYVX][iP];
						potent[GRADZVX][iB4] = potent[GRADZVX][iP];

						potent[GRADXVY][iB4] = potent[GRADXVY][iP];
						potent[GRADYVY][iB4] = potent[GRADYVY][iP];
						potent[GRADZVY][iB4] = potent[GRADZVY][iP];

						potent[GRADXVZ][iB4] = potent[GRADXVZ][iP];
						potent[GRADYVZ][iB4] = potent[GRADYVZ][iP];
						potent[GRADZVZ][iB4] = potent[GRADZVZ][iP];
					}
				}
			}
		    
			}
		}
		else
		{

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

		if (bE) {
			potent[GRADXVX][iE]=potent[GRADXVX][iP]+(dxe/dxw)*(potent[GRADXVX][iP]-potent[GRADXVX][iW]);
			potent[GRADYVX][iE]=potent[GRADYVX][iP]+(dxe/dxw)*(potent[GRADYVX][iP]-potent[GRADYVX][iW]);
			potent[GRADZVX][iE]=potent[GRADZVX][iP]+(dxe/dxw)*(potent[GRADZVX][iP]-potent[GRADZVX][iW]);

			potent[GRADXVY][iE]=potent[GRADXVY][iP]+(dxe/dxw)*(potent[GRADXVY][iP]-potent[GRADXVY][iW]);
			potent[GRADYVY][iE]=potent[GRADYVY][iP]+(dxe/dxw)*(potent[GRADYVY][iP]-potent[GRADYVY][iW]);
			potent[GRADZVY][iE]=potent[GRADZVY][iP]+(dxe/dxw)*(potent[GRADZVY][iP]-potent[GRADZVY][iW]);

			potent[GRADXVZ][iE]=potent[GRADXVZ][iP]+(dxe/dxw)*(potent[GRADXVZ][iP]-potent[GRADXVZ][iW]);
			potent[GRADYVZ][iE]=potent[GRADYVZ][iP]+(dxe/dxw)*(potent[GRADYVZ][iP]-potent[GRADYVZ][iW]);
			potent[GRADZVZ][iE]=potent[GRADZVZ][iP]+(dxe/dxw)*(potent[GRADZVZ][iP]-potent[GRADZVZ][iW]);
		}

		if (bW) {
			potent[GRADXVX][iW]=potent[GRADXVX][iP]+(dxw/dxe)*(potent[GRADXVX][iP]-potent[GRADXVX][iE]);
			potent[GRADYVX][iW]=potent[GRADYVX][iP]+(dxw/dxe)*(potent[GRADYVX][iP]-potent[GRADYVX][iE]);
			potent[GRADZVX][iW]=potent[GRADZVX][iP]+(dxw/dxe)*(potent[GRADZVX][iP]-potent[GRADZVX][iE]);

			potent[GRADXVY][iW]=potent[GRADXVY][iP]+(dxw/dxe)*(potent[GRADXVY][iP]-potent[GRADXVY][iE]);
			potent[GRADYVY][iW]=potent[GRADYVY][iP]+(dxw/dxe)*(potent[GRADYVY][iP]-potent[GRADYVY][iE]);
			potent[GRADZVY][iW]=potent[GRADZVY][iP]+(dxw/dxe)*(potent[GRADZVY][iP]-potent[GRADZVY][iE]);

			potent[GRADXVZ][iW]=potent[GRADXVZ][iP]+(dxw/dxe)*(potent[GRADXVZ][iP]-potent[GRADXVZ][iE]);
			potent[GRADYVZ][iW]=potent[GRADYVZ][iP]+(dxw/dxe)*(potent[GRADYVZ][iP]-potent[GRADYVZ][iE]);
			potent[GRADZVZ][iW]=potent[GRADZVZ][iP]+(dxw/dxe)*(potent[GRADZVZ][iP]-potent[GRADZVZ][iE]);
		}

		if (bN) {
			potent[GRADXVX][iN]=potent[GRADXVX][iP]+(dyn/dys)*(potent[GRADXVX][iP]-potent[GRADXVX][iS]);
			potent[GRADYVX][iN]=potent[GRADYVX][iP]+(dyn/dys)*(potent[GRADYVX][iP]-potent[GRADYVX][iS]);
			potent[GRADZVX][iN]=potent[GRADZVX][iP]+(dyn/dys)*(potent[GRADZVX][iP]-potent[GRADZVX][iS]);

			potent[GRADXVY][iN]=potent[GRADXVY][iP]+(dyn/dys)*(potent[GRADXVY][iP]-potent[GRADXVY][iS]);
			potent[GRADYVY][iN]=potent[GRADYVY][iP]+(dyn/dys)*(potent[GRADYVY][iP]-potent[GRADYVY][iS]);
			potent[GRADZVY][iN]=potent[GRADZVY][iP]+(dyn/dys)*(potent[GRADZVY][iP]-potent[GRADZVY][iS]);

			potent[GRADXVZ][iN]=potent[GRADXVZ][iP]+(dyn/dys)*(potent[GRADXVZ][iP]-potent[GRADXVZ][iS]);
			potent[GRADYVZ][iN]=potent[GRADYVZ][iP]+(dyn/dys)*(potent[GRADYVZ][iP]-potent[GRADYVZ][iS]);
			potent[GRADZVZ][iN]=potent[GRADZVZ][iP]+(dyn/dys)*(potent[GRADZVZ][iP]-potent[GRADZVZ][iS]);
		}

		if (bS) {
			potent[GRADXVX][iS]=potent[GRADXVX][iP]+(dys/dyn)*(potent[GRADXVX][iP]-potent[GRADXVX][iN]);
			potent[GRADYVX][iS]=potent[GRADYVX][iP]+(dys/dyn)*(potent[GRADYVX][iP]-potent[GRADYVX][iN]);
			potent[GRADZVX][iS]=potent[GRADZVX][iP]+(dys/dyn)*(potent[GRADZVX][iP]-potent[GRADZVX][iN]);

			potent[GRADXVY][iS]=potent[GRADXVY][iP]+(dys/dyn)*(potent[GRADXVY][iP]-potent[GRADXVY][iN]);
			potent[GRADYVY][iS]=potent[GRADYVY][iP]+(dys/dyn)*(potent[GRADYVY][iP]-potent[GRADYVY][iN]);
			potent[GRADZVY][iS]=potent[GRADZVY][iP]+(dys/dyn)*(potent[GRADZVY][iP]-potent[GRADZVY][iN]);

			potent[GRADXVZ][iS]=potent[GRADXVZ][iP]+(dys/dyn)*(potent[GRADXVZ][iP]-potent[GRADXVZ][iN]);
			potent[GRADYVZ][iS]=potent[GRADYVZ][iP]+(dys/dyn)*(potent[GRADYVZ][iP]-potent[GRADYVZ][iN]);
			potent[GRADZVZ][iS]=potent[GRADZVZ][iP]+(dys/dyn)*(potent[GRADZVZ][iP]-potent[GRADZVZ][iN]);
		}

		if (bT) {
			potent[GRADXVX][iT]=potent[GRADXVX][iP]+(dzt/dzb)*(potent[GRADXVX][iP]-potent[GRADXVX][iB]);
			potent[GRADYVX][iT]=potent[GRADYVX][iP]+(dzt/dzb)*(potent[GRADYVX][iP]-potent[GRADYVX][iB]);
			potent[GRADZVX][iT]=potent[GRADZVX][iP]+(dzt/dzb)*(potent[GRADZVX][iP]-potent[GRADZVX][iB]);

			potent[GRADXVY][iT]=potent[GRADXVY][iP]+(dzt/dzb)*(potent[GRADXVY][iP]-potent[GRADXVY][iB]);
			potent[GRADYVY][iT]=potent[GRADYVY][iP]+(dzt/dzb)*(potent[GRADYVY][iP]-potent[GRADYVY][iB]);
			potent[GRADZVY][iT]=potent[GRADZVY][iP]+(dzt/dzb)*(potent[GRADZVY][iP]-potent[GRADZVY][iB]);

			potent[GRADXVZ][iT]=potent[GRADXVZ][iP]+(dzt/dzb)*(potent[GRADXVZ][iP]-potent[GRADXVZ][iB]);
			potent[GRADYVZ][iT]=potent[GRADYVZ][iP]+(dzt/dzb)*(potent[GRADYVZ][iP]-potent[GRADYVZ][iB]);
			potent[GRADZVZ][iT]=potent[GRADZVZ][iP]+(dzt/dzb)*(potent[GRADZVZ][iP]-potent[GRADZVZ][iB]);
		}

		if (bB) {
			potent[GRADXVX][iB]=potent[GRADXVX][iP]+(dzb/dzt)*(potent[GRADXVX][iP]-potent[GRADXVX][iT]);
			potent[GRADYVX][iB]=potent[GRADYVX][iP]+(dzb/dzt)*(potent[GRADYVX][iP]-potent[GRADYVX][iT]);
			potent[GRADZVX][iB]=potent[GRADZVX][iP]+(dzb/dzt)*(potent[GRADZVX][iP]-potent[GRADZVX][iT]);

			potent[GRADXVY][iB]=potent[GRADXVY][iP]+(dzb/dzt)*(potent[GRADXVY][iP]-potent[GRADXVY][iT]);
			potent[GRADYVY][iB]=potent[GRADYVY][iP]+(dzb/dzt)*(potent[GRADYVY][iP]-potent[GRADYVY][iT]);
			potent[GRADZVY][iB]=potent[GRADZVY][iP]+(dzb/dzt)*(potent[GRADZVY][iP]-potent[GRADZVY][iT]);

			potent[GRADXVZ][iB]=potent[GRADXVZ][iP]+(dzb/dzt)*(potent[GRADXVZ][iP]-potent[GRADXVZ][iT]);
			potent[GRADYVZ][iB]=potent[GRADYVZ][iP]+(dzb/dzt)*(potent[GRADYVZ][iP]-potent[GRADYVZ][iT]);
			potent[GRADZVZ][iB]=potent[GRADZVZ][iP]+(dzb/dzt)*(potent[GRADZVZ][iP]-potent[GRADZVZ][iT]);
		}

		}
	}

} // green_gaussO1

// 14.04.2019; 29.09.2019. Работает на АЛИС сетке.
// Вычисление градиентов модифицированной кинематической турбулентной вязкости
// в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
// Поскольку интерполяция линейная то точность данной формулы O(h). 
// По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
// что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
// совпадает со взвешенным методом наименьших квадратов.
void green_gauss_SpallartAllmares(integer iP,
	doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND* &border_neighbor, integer *ilevel_alice, TOCHKA*& volume_loc) {


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу 1983г.
	//bool bRCh = false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	//volume3D(iP, nvtx, pa, dx, dy, dz);
	//dx = fabs(dx);
	//dy = fabs(dy);
	//dz = fabs(dz);

	TOCHKA point_loc = volume_loc[iP];
	dx = point_loc.x;
	dy = point_loc.y;
	dz = point_loc.z;

	doublereal dxe = 0.5*dx, dxw = 0.5*dx, dyn = 0.5*dy, dys = 0.5*dy, dzt = 0.5*dz, dzb = 0.5*dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5*dx / dxe;
	fwplus = 0.5*dx / dxw;
	// y-direction
	fnplus = 0.5*dy / dyn;
	fsplus = 0.5*dy / dys;
	// z-direction
	ftplus = 0.5*dz / dzt;
	fbplus = 0.5*dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE4];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW4];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN4];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS4];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT4];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB4];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}





	// линейная интерполяция скорости VX на грань КО.
	doublereal fe = 0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// VX
		if (iE > -1) {
			if (!bE) {
				fe = feplus * potent[NUSHA][iE] + (1.0 - feplus)*potent[NUSHA][iP];
			}
			else fe = potent[NUSHA][iE];
		}
		if (iW > -1) {
			if (!bW) {
				fw = fwplus * potent[NUSHA][iW] + (1.0 - fwplus)*potent[NUSHA][iP];
			}
			else fw = potent[NUSHA][iW];
		}

		if (iE2 > -1) {
			if (!bE2) {
				fe2 = feplus2 * potent[NUSHA][iE2] + (1.0 - feplus2)*potent[NUSHA][iP];
			}
			else fe2 = potent[NUSHA][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) {
				fw2 = fwplus2 * potent[NUSHA][iW2] + (1.0 - fwplus2)*potent[NUSHA][iP];
			}
			else fw2 = potent[NUSHA][iW2];
		}

		if (iE3 > -1) {
			if (!bE3) {
				fe3 = feplus3 * potent[NUSHA][iE3] + (1.0 - feplus3)*potent[NUSHA][iP];
			}
			else fe3 = potent[NUSHA][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) {
				fw3 = fwplus3 * potent[NUSHA][iW3] + (1.0 - fwplus3)*potent[NUSHA][iP];
			}
			else fw3 = potent[NUSHA][iW3];
		}

		if (iE4 > -1) {
			if (!bE4) {
				fe4 = feplus4 * potent[NUSHA][iE4] + (1.0 - feplus4)*potent[NUSHA][iP];
			}
			else fe4 = potent[NUSHA][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) {
				fw4 = fwplus4 * potent[NUSHA][iW4] + (1.0 - fwplus4)*potent[NUSHA][iP];
			}
			else fw4 = potent[NUSHA][iW4];
		}
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus * potent[NUSHA][iN] + (1.0 - fnplus)*potent[NUSHA][iP]; else fn = potent[NUSHA][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus * potent[NUSHA][iS] + (1.0 - fsplus)*potent[NUSHA][iP]; else fs = potent[NUSHA][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus * potent[NUSHA][iT] + (1.0 - ftplus)*potent[NUSHA][iP]; else ft = potent[NUSHA][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus * potent[NUSHA][iB] + (1.0 - fbplus)*potent[NUSHA][iP]; else fb = potent[NUSHA][iB];
		}

		if (iN2 > -1) {
			if (!bN2) fn2 = fnplus2 * potent[NUSHA][iN2] + (1.0 - fnplus2)*potent[NUSHA][iP]; else fn2 = potent[NUSHA][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) fs2 = fsplus2 * potent[NUSHA][iS2] + (1.0 - fsplus2)*potent[NUSHA][iP]; else fs2 = potent[NUSHA][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) ft2 = ftplus2 * potent[NUSHA][iT2] + (1.0 - ftplus2)*potent[NUSHA][iP]; else ft2 = potent[NUSHA][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) fb2 = fbplus2 * potent[NUSHA][iB2] + (1.0 - fbplus2)*potent[NUSHA][iP]; else fb2 = potent[NUSHA][iB2];
		}

		if (iN3 > -1) {
			if (!bN3) fn3 = fnplus3 * potent[NUSHA][iN3] + (1.0 - fnplus3)*potent[NUSHA][iP]; else fn3 = potent[NUSHA][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) fs3 = fsplus3 * potent[NUSHA][iS3] + (1.0 - fsplus3)*potent[NUSHA][iP]; else fs3 = potent[NUSHA][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) ft3 = ftplus3 * potent[NUSHA][iT3] + (1.0 - ftplus3)*potent[NUSHA][iP]; else ft3 = potent[NUSHA][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) fb3 = fbplus3 * potent[NUSHA][iB3] + (1.0 - fbplus3)*potent[NUSHA][iP]; else fb3 = potent[NUSHA][iB3];
		}

		if (iN4 > -1) {
			if (!bN4) fn4 = fnplus4 * potent[NUSHA][iN4] + (1.0 - fnplus4)*potent[NUSHA][iP]; else fn4 = potent[NUSHA][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) fs4 = fsplus4 * potent[NUSHA][iS4] + (1.0 - fsplus4)*potent[NUSHA][iP]; else fs4 = potent[NUSHA][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) ft4 = ftplus4 * potent[NUSHA][iT4] + (1.0 - ftplus4)*potent[NUSHA][iP]; else ft4 = potent[NUSHA][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) fb4 = fbplus4 * potent[NUSHA][iB4] + (1.0 - fbplus4)*potent[NUSHA][iP]; else fb4 = potent[NUSHA][iB4];
		}
		// градиент NUSHA
		//potent[GRADXNUSHA][iP]=(fe-fw)/dx;
		//potent[GRADYNUSHA][iP]=(fn-fs)/dy;
		//potent[GRADZNUSHA][iP]=(ft-fb)/dz;
		potent[GRADXNUSHA][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYNUSHA][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZNUSHA][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;


	}
	else {

		if (1) {



			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE])*(potent[VXCOR][iE]) + (potent[VYCOR][iE])*(potent[VYCOR][iE]) + (potent[VZCOR][iE])*(potent[VZCOR][iE]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iE] = 0.0;
						potent[GRADYNUSHA][iE] = 0.0;
						potent[GRADZNUSHA][iE] = 0.0;



						/*
						potent[GRADXNUSHA][iP] = 0.0;
						potent[GRADYNUSHA][iP] = 0.0;
						potent[GRADZNUSHA][iP] = 0.0;


						*/
					}
					else {

						potent[GRADXNUSHA][iE] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iE] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iE] = potent[GRADZNUSHA][iP];


					}
				}
			}



			if (iE2 > -1) {
				if (bE2) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE2])*(potent[VXCOR][iE2]) + (potent[VYCOR][iE2])*(potent[VYCOR][iE2]) + (potent[VZCOR][iE2])*(potent[VZCOR][iE2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iE2] = 0.0;
						potent[GRADYNUSHA][iE2] = 0.0;
						potent[GRADZNUSHA][iE2] = 0.0;




					}
					else {

						potent[GRADXNUSHA][iE2] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iE2] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iE2] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iE3 > -1) {
				if (bE3) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE3])*(potent[VXCOR][iE3]) + (potent[VYCOR][iE3])*(potent[VYCOR][iE3]) + (potent[VZCOR][iE3])*(potent[VZCOR][iE3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iE3] = 0.0;
						potent[GRADYNUSHA][iE3] = 0.0;
						potent[GRADZNUSHA][iE3] = 0.0;




					}
					else {

						potent[GRADXNUSHA][iE3] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iE3] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iE3] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iE4 > -1) {
				if (bE4) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE4])*(potent[VXCOR][iE4]) + (potent[VYCOR][iE4])*(potent[VYCOR][iE4]) + (potent[VZCOR][iE4])*(potent[VZCOR][iE4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iE4] = 0.0;
						potent[GRADYNUSHA][iE4] = 0.0;
						potent[GRADZNUSHA][iE4] = 0.0;




					}
					else {

						potent[GRADXNUSHA][iE4] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iE4] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iE4] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iW > -1) {
				if (bW) {

					doublereal dspeed = sqrt((potent[VXCOR][iW])*(potent[VXCOR][iW]) + (potent[VYCOR][iW])*(potent[VYCOR][iW]) + (potent[VZCOR][iW])*(potent[VZCOR][iW]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iW] = 0.0;
						potent[GRADYNUSHA][iW] = 0.0;
						potent[GRADZNUSHA][iW] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iW] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iW] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iW] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iW2 > -1) {
				if (bW2) {

					doublereal dspeed = sqrt((potent[VXCOR][iW2])*(potent[VXCOR][iW2]) + (potent[VYCOR][iW2])*(potent[VYCOR][iW2]) + (potent[VZCOR][iW2])*(potent[VZCOR][iW2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iW2] = 0.0;
						potent[GRADYNUSHA][iW2] = 0.0;
						potent[GRADZNUSHA][iW2] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iW2] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iW2] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iW2] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					doublereal dspeed = sqrt((potent[VXCOR][iW3])*(potent[VXCOR][iW3]) + (potent[VYCOR][iW3])*(potent[VYCOR][iW3]) + (potent[VZCOR][iW3])*(potent[VZCOR][iW3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iW3] = 0.0;
						potent[GRADYNUSHA][iW3] = 0.0;
						potent[GRADZNUSHA][iW3] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iW3] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iW3] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iW3] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					doublereal dspeed = sqrt((potent[VXCOR][iW4])*(potent[VXCOR][iW4]) + (potent[VYCOR][iW4])*(potent[VYCOR][iW4]) + (potent[VZCOR][iW4])*(potent[VZCOR][iW4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iW4] = 0.0;
						potent[GRADYNUSHA][iW4] = 0.0;
						potent[GRADZNUSHA][iW4] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iW4] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iW4] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iW4] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iN > -1) {
				if (bN) {

					doublereal dspeed = sqrt((potent[VXCOR][iN])*(potent[VXCOR][iN]) + (potent[VYCOR][iN])*(potent[VYCOR][iN]) + (potent[VZCOR][iN])*(potent[VZCOR][iN]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iN] = 0.0;
						potent[GRADYNUSHA][iN] = 0.0;
						potent[GRADZNUSHA][iN] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iN] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iN] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iN] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iN2 > -1) {
				if (bN2) {

					doublereal dspeed = sqrt((potent[VXCOR][iN2])*(potent[VXCOR][iN2]) + (potent[VYCOR][iN2])*(potent[VYCOR][iN2]) + (potent[VZCOR][iN2])*(potent[VZCOR][iN2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iN2] = 0.0;
						potent[GRADYNUSHA][iN2] = 0.0;
						potent[GRADZNUSHA][iN2] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iN2] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iN2] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iN2] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					doublereal dspeed = sqrt((potent[VXCOR][iN3])*(potent[VXCOR][iN3]) + (potent[VYCOR][iN3])*(potent[VYCOR][iN3]) + (potent[VZCOR][iN3])*(potent[VZCOR][iN3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iN3] = 0.0;
						potent[GRADYNUSHA][iN3] = 0.0;
						potent[GRADZNUSHA][iN3] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iN3] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iN3] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iN3] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					doublereal dspeed = sqrt((potent[VXCOR][iN4])*(potent[VXCOR][iN4]) + (potent[VYCOR][iN4])*(potent[VYCOR][iN4]) + (potent[VZCOR][iN4])*(potent[VZCOR][iN4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iN4] = 0.0;
						potent[GRADYNUSHA][iN4] = 0.0;
						potent[GRADZNUSHA][iN4] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iN4] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iN4] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iN4] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iS > -1) {
				if (bS) {

					doublereal dspeed = sqrt((potent[VXCOR][iS])*(potent[VXCOR][iS]) + (potent[VYCOR][iS])*(potent[VYCOR][iS]) + (potent[VZCOR][iS])*(potent[VZCOR][iS]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iS] = 0.0;
						potent[GRADYNUSHA][iS] = 0.0;
						potent[GRADZNUSHA][iS] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iS] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iS] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iS] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iS2 > -1) {
				if (bS2) {

					doublereal dspeed = sqrt((potent[VXCOR][iS2])*(potent[VXCOR][iS2]) + (potent[VYCOR][iS2])*(potent[VYCOR][iS2]) + (potent[VZCOR][iS2])*(potent[VZCOR][iS2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iS2] = 0.0;
						potent[GRADYNUSHA][iS2] = 0.0;
						potent[GRADZNUSHA][iS2] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iS2] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iS2] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iS2] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					doublereal dspeed = sqrt((potent[VXCOR][iS3])*(potent[VXCOR][iS3]) + (potent[VYCOR][iS3])*(potent[VYCOR][iS3]) + (potent[VZCOR][iS3])*(potent[VZCOR][iS3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iS3] = 0.0;
						potent[GRADYNUSHA][iS3] = 0.0;
						potent[GRADZNUSHA][iS3] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iS3] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iS3] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iS3] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					doublereal dspeed = sqrt((potent[VXCOR][iS4])*(potent[VXCOR][iS4]) + (potent[VYCOR][iS4])*(potent[VYCOR][iS4]) + (potent[VZCOR][iS4])*(potent[VZCOR][iS4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iS4] = 0.0;
						potent[GRADYNUSHA][iS4] = 0.0;
						potent[GRADZNUSHA][iS4] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iS4] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iS4] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iS4] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iT > -1) {
				if (bT) {

					doublereal dspeed = sqrt((potent[VXCOR][iT])*(potent[VXCOR][iT]) + (potent[VYCOR][iT])*(potent[VYCOR][iT]) + (potent[VZCOR][iT])*(potent[VZCOR][iT]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iT] = 0.0;
						potent[GRADYNUSHA][iT] = 0.0;
						potent[GRADZNUSHA][iT] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iT] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iT] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iT] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iT2 > -1) {
				if (bT2) {

					doublereal dspeed = sqrt((potent[VXCOR][iT2])*(potent[VXCOR][iT2]) + (potent[VYCOR][iT2])*(potent[VYCOR][iT2]) + (potent[VZCOR][iT2])*(potent[VZCOR][iT2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iT2] = 0.0;
						potent[GRADYNUSHA][iT2] = 0.0;
						potent[GRADZNUSHA][iT2] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iT2] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iT2] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iT2] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					doublereal dspeed = sqrt((potent[VXCOR][iT3])*(potent[VXCOR][iT3]) + (potent[VYCOR][iT3])*(potent[VYCOR][iT3]) + (potent[VZCOR][iT3])*(potent[VZCOR][iT3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iT3] = 0.0;
						potent[GRADYNUSHA][iT3] = 0.0;
						potent[GRADZNUSHA][iT3] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iT3] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iT3] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iT3] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					doublereal dspeed = sqrt((potent[VXCOR][iT4])*(potent[VXCOR][iT4]) + (potent[VYCOR][iT4])*(potent[VYCOR][iT4]) + (potent[VZCOR][iT4])*(potent[VZCOR][iT4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iT4] = 0.0;
						potent[GRADYNUSHA][iT4] = 0.0;
						potent[GRADZNUSHA][iT4] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iT4] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iT4] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iT4] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iB > -1) {
				if (bB) {
					doublereal dspeed = sqrt((potent[VXCOR][iB])*(potent[VXCOR][iB]) + (potent[VYCOR][iB])*(potent[VYCOR][iB]) + (potent[VZCOR][iB])*(potent[VZCOR][iB]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iB] = 0.0;
						potent[GRADYNUSHA][iB] = 0.0;
						potent[GRADZNUSHA][iB] = 0.0;


					}
					else {

						potent[GRADXNUSHA][iB] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iB] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iB] = potent[GRADZNUSHA][iP];


					}
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					doublereal dspeed = sqrt((potent[VXCOR][iB2])*(potent[VXCOR][iB2]) + (potent[VYCOR][iB2])*(potent[VYCOR][iB2]) + (potent[VZCOR][iB2])*(potent[VZCOR][iB2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iB2] = 0.0;
						potent[GRADYNUSHA][iB2] = 0.0;
						potent[GRADZNUSHA][iB2] = 0.0;

					}
					else {

						potent[GRADXNUSHA][iB2] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iB2] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iB2] = potent[GRADZNUSHA][iP];

					}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					doublereal dspeed = sqrt((potent[VXCOR][iB3])*(potent[VXCOR][iB3]) + (potent[VYCOR][iB3])*(potent[VYCOR][iB3]) + (potent[VZCOR][iB3])*(potent[VZCOR][iB3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iB3] = 0.0;
						potent[GRADYNUSHA][iB3] = 0.0;
						potent[GRADZNUSHA][iB3] = 0.0;
					}
					else {

						potent[GRADXNUSHA][iB3] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iB3] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iB3] = potent[GRADZNUSHA][iP];
					}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					doublereal dspeed = sqrt((potent[VXCOR][iB4])*(potent[VXCOR][iB4]) + (potent[VYCOR][iB4])*(potent[VYCOR][iB4]) + (potent[VZCOR][iB4])*(potent[VZCOR][iB4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXNUSHA][iB4] = 0.0;
						potent[GRADYNUSHA][iB4] = 0.0;
						potent[GRADZNUSHA][iB4] = 0.0;
					}
					else {

						potent[GRADXNUSHA][iB4] = potent[GRADXNUSHA][iP];
						potent[GRADYNUSHA][iB4] = potent[GRADYNUSHA][iP];
						potent[GRADZNUSHA][iB4] = potent[GRADZNUSHA][iP];

					}
				}
			}

		}
		else
		{

			// граничные узлы.
			// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

			if (bE) {
				potent[GRADXNUSHA][iE] = potent[GRADXNUSHA][iP] + (dxe / dxw)*(potent[GRADXNUSHA][iP] - potent[GRADXNUSHA][iW]);
				potent[GRADYNUSHA][iE] = potent[GRADYNUSHA][iP] + (dxe / dxw)*(potent[GRADYNUSHA][iP] - potent[GRADYNUSHA][iW]);
				potent[GRADZNUSHA][iE] = potent[GRADZNUSHA][iP] + (dxe / dxw)*(potent[GRADZNUSHA][iP] - potent[GRADZNUSHA][iW]);

			}

			if (bW) {
				potent[GRADXNUSHA][iW] = potent[GRADXNUSHA][iP] + (dxw / dxe)*(potent[GRADXNUSHA][iP] - potent[GRADXNUSHA][iE]);
				potent[GRADYNUSHA][iW] = potent[GRADYNUSHA][iP] + (dxw / dxe)*(potent[GRADYNUSHA][iP] - potent[GRADYNUSHA][iE]);
				potent[GRADZNUSHA][iW] = potent[GRADZNUSHA][iP] + (dxw / dxe)*(potent[GRADZNUSHA][iP] - potent[GRADZNUSHA][iE]);

			}

			if (bN) {
				potent[GRADXNUSHA][iN] = potent[GRADXNUSHA][iP] + (dyn / dys)*(potent[GRADXNUSHA][iP] - potent[GRADXNUSHA][iS]);
				potent[GRADYNUSHA][iN] = potent[GRADYNUSHA][iP] + (dyn / dys)*(potent[GRADYNUSHA][iP] - potent[GRADYNUSHA][iS]);
				potent[GRADZNUSHA][iN] = potent[GRADZNUSHA][iP] + (dyn / dys)*(potent[GRADZNUSHA][iP] - potent[GRADZNUSHA][iS]);

			}

			if (bS) {
				potent[GRADXNUSHA][iS] = potent[GRADXNUSHA][iP] + (dys / dyn)*(potent[GRADXNUSHA][iP] - potent[GRADXNUSHA][iN]);
				potent[GRADYNUSHA][iS] = potent[GRADYNUSHA][iP] + (dys / dyn)*(potent[GRADYNUSHA][iP] - potent[GRADYNUSHA][iN]);
				potent[GRADZNUSHA][iS] = potent[GRADZNUSHA][iP] + (dys / dyn)*(potent[GRADZNUSHA][iP] - potent[GRADZNUSHA][iN]);

			}

			if (bT) {
				potent[GRADXNUSHA][iT] = potent[GRADXNUSHA][iP] + (dzt / dzb)*(potent[GRADXNUSHA][iP] - potent[GRADXNUSHA][iB]);
				potent[GRADYNUSHA][iT] = potent[GRADYNUSHA][iP] + (dzt / dzb)*(potent[GRADYNUSHA][iP] - potent[GRADYNUSHA][iB]);
				potent[GRADZNUSHA][iT] = potent[GRADZNUSHA][iP] + (dzt / dzb)*(potent[GRADZNUSHA][iP] - potent[GRADZNUSHA][iB]);

			}

			if (bB) {
				potent[GRADXNUSHA][iB] = potent[GRADXNUSHA][iP] + (dzb / dzt)*(potent[GRADXNUSHA][iP] - potent[GRADXNUSHA][iT]);
				potent[GRADYNUSHA][iB] = potent[GRADYNUSHA][iP] + (dzb / dzt)*(potent[GRADYNUSHA][iP] - potent[GRADYNUSHA][iT]);
				potent[GRADZNUSHA][iB] = potent[GRADZNUSHA][iP] + (dzb / dzt)*(potent[GRADZNUSHA][iP] - potent[GRADZNUSHA][iT]);

			}

		}
	}

} // green_gauss_SpallartAllmares


// Производные от k на твердой неподвижной стенке ?
// Надо попробовать оба варианта и остановиться на каком -то одном.
// один вариант 0 другой вариант скопировать из ближайшего внутреннего узла.
// Это касается всех моделей турбулентности RANS. 03.10.2019.

// 14.04.2019; 03.10.2019. Работает на АЛИС сетке.
// Вычисление градиентов кинетической энергии турбулентных пульсаций
// в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
// Поскольку интерполяция линейная то точность данной формулы O(h). 
// По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
// что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
// совпадает со взвешенным методом наименьших квадратов.
void green_gauss_turbulent_kinetik_energy_MenterSST(integer iP,
	doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND* &border_neighbor, integer *ilevel_alice)
{


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу.
	//bool bRCh = false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов


	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}


	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe = 0.5*dx, dxw = 0.5*dx, dyn = 0.5*dy, dys = 0.5*dy, dzt = 0.5*dz, dzb = 0.5*dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5*dx / dxe;
	fwplus = 0.5*dx / dxw;
	// y-direction
	fnplus = 0.5*dy / dyn;
	fsplus = 0.5*dy / dys;
	// z-direction
	ftplus = 0.5*dz / dzt;
	fbplus = 0.5*dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}





	// линейная интерполяция скорости VX на грань КО.
	doublereal fe = 0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// TURBULENT_KINETIK_ENERGY
		if (iE > -1) {
			if (!bE) {
				fe = feplus * potent[TURBULENT_KINETIK_ENERGY][iE] + (1.0 - feplus)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fe = potent[TURBULENT_KINETIK_ENERGY][iE];
		}
		if (iW > -1) {
			if (!bW) {
				fw = fwplus * potent[TURBULENT_KINETIK_ENERGY][iW] + (1.0 - fwplus)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fw = potent[TURBULENT_KINETIK_ENERGY][iW];
		}

		if (iE2 > -1) {
			if (!bE2) {
				fe2 = feplus2 * potent[TURBULENT_KINETIK_ENERGY][iE2] + (1.0 - feplus2)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fe2 = potent[TURBULENT_KINETIK_ENERGY][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) {
				fw2 = fwplus2 * potent[TURBULENT_KINETIK_ENERGY][iW2] + (1.0 - fwplus2)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fw2 = potent[TURBULENT_KINETIK_ENERGY][iW2];
		}

		if (iE3 > -1) {
			if (!bE3) {
				fe3 = feplus3 * potent[TURBULENT_KINETIK_ENERGY][iE3] + (1.0 - feplus3)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fe3 = potent[TURBULENT_KINETIK_ENERGY][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) {
				fw3 = fwplus3 * potent[TURBULENT_KINETIK_ENERGY][iW3] + (1.0 - fwplus3)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fw3 = potent[TURBULENT_KINETIK_ENERGY][iW3];
		}

		if (iE4 > -1) {
			if (!bE4) {
				fe4 = feplus4 * potent[TURBULENT_KINETIK_ENERGY][iE4] + (1.0 - feplus4)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fe4 = potent[TURBULENT_KINETIK_ENERGY][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) {
				fw4 = fwplus4 * potent[TURBULENT_KINETIK_ENERGY][iW4] + (1.0 - fwplus4)*potent[TURBULENT_KINETIK_ENERGY][iP];
			}
			else fw4 = potent[TURBULENT_KINETIK_ENERGY][iW4];
		}
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus * potent[TURBULENT_KINETIK_ENERGY][iN] + (1.0 - fnplus)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fn = potent[TURBULENT_KINETIK_ENERGY][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus * potent[TURBULENT_KINETIK_ENERGY][iS] + (1.0 - fsplus)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fs = potent[TURBULENT_KINETIK_ENERGY][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus * potent[TURBULENT_KINETIK_ENERGY][iT] + (1.0 - ftplus)*potent[TURBULENT_KINETIK_ENERGY][iP]; else ft = potent[TURBULENT_KINETIK_ENERGY][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus * potent[TURBULENT_KINETIK_ENERGY][iB] + (1.0 - fbplus)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fb = potent[TURBULENT_KINETIK_ENERGY][iB];
		}

		if (iN2 > -1) {
			if (!bN2) fn2 = fnplus2 * potent[TURBULENT_KINETIK_ENERGY][iN2] + (1.0 - fnplus2)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fn2 = potent[TURBULENT_KINETIK_ENERGY][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) fs2 = fsplus2 * potent[TURBULENT_KINETIK_ENERGY][iS2] + (1.0 - fsplus2)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fs2 = potent[TURBULENT_KINETIK_ENERGY][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) ft2 = ftplus2 * potent[TURBULENT_KINETIK_ENERGY][iT2] + (1.0 - ftplus2)*potent[TURBULENT_KINETIK_ENERGY][iP]; else ft2 = potent[TURBULENT_KINETIK_ENERGY][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) fb2 = fbplus2 * potent[TURBULENT_KINETIK_ENERGY][iB2] + (1.0 - fbplus2)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fb2 = potent[TURBULENT_KINETIK_ENERGY][iB2];
		}

		if (iN3 > -1) {
			if (!bN3) fn3 = fnplus3 * potent[TURBULENT_KINETIK_ENERGY][iN3] + (1.0 - fnplus3)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fn3 = potent[TURBULENT_KINETIK_ENERGY][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) fs3 = fsplus3 * potent[TURBULENT_KINETIK_ENERGY][iS3] + (1.0 - fsplus3)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fs3 = potent[TURBULENT_KINETIK_ENERGY][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) ft3 = ftplus3 * potent[TURBULENT_KINETIK_ENERGY][iT3] + (1.0 - ftplus3)*potent[TURBULENT_KINETIK_ENERGY][iP]; else ft3 = potent[TURBULENT_KINETIK_ENERGY][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) fb3 = fbplus3 * potent[TURBULENT_KINETIK_ENERGY][iB3] + (1.0 - fbplus3)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fb3 = potent[TURBULENT_KINETIK_ENERGY][iB3];
		}

		if (iN4 > -1) {
			if (!bN4) fn4 = fnplus4 * potent[TURBULENT_KINETIK_ENERGY][iN4] + (1.0 - fnplus4)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fn4 = potent[TURBULENT_KINETIK_ENERGY][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) fs4 = fsplus4 * potent[TURBULENT_KINETIK_ENERGY][iS4] + (1.0 - fsplus4)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fs4 = potent[TURBULENT_KINETIK_ENERGY][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) ft4 = ftplus4 * potent[TURBULENT_KINETIK_ENERGY][iT4] + (1.0 - ftplus4)*potent[TURBULENT_KINETIK_ENERGY][iP]; else ft4 = potent[TURBULENT_KINETIK_ENERGY][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) fb4 = fbplus4 * potent[TURBULENT_KINETIK_ENERGY][iB4] + (1.0 - fbplus4)*potent[TURBULENT_KINETIK_ENERGY][iP]; else fb4 = potent[TURBULENT_KINETIK_ENERGY][iB4];
		}
		// градиент NUSHA
		//potent[GRADXTURBULENT_KINETIK_ENERGY][iP]=(fe-fw)/dx;
		//potent[GRADYTURBULENT_KINETIK_ENERGY][iP]=(fn-fs)/dy;
		//potent[GRADZTURBULENT_KINETIK_ENERGY][iP]=(ft-fb)/dz;
		potent[GRADXTURBULENT_KINETIK_ENERGY][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYTURBULENT_KINETIK_ENERGY][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZTURBULENT_KINETIK_ENERGY][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;


	}
	else {

		if (1) {



			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE])*(potent[VXCOR][iE]) + (potent[VYCOR][iE])*(potent[VYCOR][iE]) + (potent[VZCOR][iE])*(potent[VZCOR][iE]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iE] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE] = 0.0;



						/*
						potent[GRADXTURBULENT_KINETIK_ENERGY][iP] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iP] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iP] = 0.0;


						*/
					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iE] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}



			if (iE2 > -1) {
				if (bE2) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE2])*(potent[VXCOR][iE2]) + (potent[VYCOR][iE2])*(potent[VYCOR][iE2]) + (potent[VZCOR][iE2])*(potent[VZCOR][iE2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iE2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE2] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iE2] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE2] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE2] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iE3 > -1) {
				if (bE3) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE3])*(potent[VXCOR][iE3]) + (potent[VYCOR][iE3])*(potent[VYCOR][iE3]) + (potent[VZCOR][iE3])*(potent[VZCOR][iE3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iE3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE3] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iE3] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE3] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE3] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iE4 > -1) {
				if (bE4) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE4])*(potent[VXCOR][iE4]) + (potent[VYCOR][iE4])*(potent[VYCOR][iE4]) + (potent[VZCOR][iE4])*(potent[VZCOR][iE4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iE4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE4] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iE4] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iE4] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iE4] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iW > -1) {
				if (bW) {

					doublereal dspeed = sqrt((potent[VXCOR][iW])*(potent[VXCOR][iW]) + (potent[VYCOR][iW])*(potent[VYCOR][iW]) + (potent[VZCOR][iW])*(potent[VZCOR][iW]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iW] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iW] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iW2 > -1) {
				if (bW2) {

					doublereal dspeed = sqrt((potent[VXCOR][iW2])*(potent[VXCOR][iW2]) + (potent[VYCOR][iW2])*(potent[VYCOR][iW2]) + (potent[VZCOR][iW2])*(potent[VZCOR][iW2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iW2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iW2] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW2] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW2] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					doublereal dspeed = sqrt((potent[VXCOR][iW3])*(potent[VXCOR][iW3]) + (potent[VYCOR][iW3])*(potent[VYCOR][iW3]) + (potent[VZCOR][iW3])*(potent[VZCOR][iW3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iW3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iW3] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW3] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW3] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					doublereal dspeed = sqrt((potent[VXCOR][iW4])*(potent[VXCOR][iW4]) + (potent[VYCOR][iW4])*(potent[VYCOR][iW4]) + (potent[VZCOR][iW4])*(potent[VZCOR][iW4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iW4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iW4] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iW4] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iW4] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iN > -1) {
				if (bN) {

					doublereal dspeed = sqrt((potent[VXCOR][iN])*(potent[VXCOR][iN]) + (potent[VYCOR][iN])*(potent[VYCOR][iN]) + (potent[VZCOR][iN])*(potent[VZCOR][iN]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iN] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iN] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iN2 > -1) {
				if (bN2) {

					doublereal dspeed = sqrt((potent[VXCOR][iN2])*(potent[VXCOR][iN2]) + (potent[VYCOR][iN2])*(potent[VYCOR][iN2]) + (potent[VZCOR][iN2])*(potent[VZCOR][iN2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iN2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iN2] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN2] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN2] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					doublereal dspeed = sqrt((potent[VXCOR][iN3])*(potent[VXCOR][iN3]) + (potent[VYCOR][iN3])*(potent[VYCOR][iN3]) + (potent[VZCOR][iN3])*(potent[VZCOR][iN3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iN3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iN3] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN3] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN3] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					doublereal dspeed = sqrt((potent[VXCOR][iN4])*(potent[VXCOR][iN4]) + (potent[VYCOR][iN4])*(potent[VYCOR][iN4]) + (potent[VZCOR][iN4])*(potent[VZCOR][iN4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iN4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iN4] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iN4] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iN4] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iS > -1) {
				if (bS) {

					doublereal dspeed = sqrt((potent[VXCOR][iS])*(potent[VXCOR][iS]) + (potent[VYCOR][iS])*(potent[VYCOR][iS]) + (potent[VZCOR][iS])*(potent[VZCOR][iS]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iS] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iS] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iS2 > -1) {
				if (bS2) {

					doublereal dspeed = sqrt((potent[VXCOR][iS2])*(potent[VXCOR][iS2]) + (potent[VYCOR][iS2])*(potent[VYCOR][iS2]) + (potent[VZCOR][iS2])*(potent[VZCOR][iS2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iS2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iS2] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS2] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS2] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					doublereal dspeed = sqrt((potent[VXCOR][iS3])*(potent[VXCOR][iS3]) + (potent[VYCOR][iS3])*(potent[VYCOR][iS3]) + (potent[VZCOR][iS3])*(potent[VZCOR][iS3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iS3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iS3] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS3] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS3] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					doublereal dspeed = sqrt((potent[VXCOR][iS4])*(potent[VXCOR][iS4]) + (potent[VYCOR][iS4])*(potent[VYCOR][iS4]) + (potent[VZCOR][iS4])*(potent[VZCOR][iS4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iS4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iS4] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iS4] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iS4] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iT > -1) {
				if (bT) {

					doublereal dspeed = sqrt((potent[VXCOR][iT])*(potent[VXCOR][iT]) + (potent[VYCOR][iT])*(potent[VYCOR][iT]) + (potent[VZCOR][iT])*(potent[VZCOR][iT]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iT] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iT] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iT2 > -1) {
				if (bT2) {

					doublereal dspeed = sqrt((potent[VXCOR][iT2])*(potent[VXCOR][iT2]) + (potent[VYCOR][iT2])*(potent[VYCOR][iT2]) + (potent[VZCOR][iT2])*(potent[VZCOR][iT2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iT2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iT2] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT2] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT2] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					doublereal dspeed = sqrt((potent[VXCOR][iT3])*(potent[VXCOR][iT3]) + (potent[VYCOR][iT3])*(potent[VYCOR][iT3]) + (potent[VZCOR][iT3])*(potent[VZCOR][iT3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iT3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iT3] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT3] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT3] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					doublereal dspeed = sqrt((potent[VXCOR][iT4])*(potent[VXCOR][iT4]) + (potent[VYCOR][iT4])*(potent[VYCOR][iT4]) + (potent[VZCOR][iT4])*(potent[VZCOR][iT4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iT4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iT4] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iT4] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iT4] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iB > -1) {
				if (bB) {
					doublereal dspeed = sqrt((potent[VXCOR][iB])*(potent[VXCOR][iB]) + (potent[VYCOR][iB])*(potent[VYCOR][iB]) + (potent[VZCOR][iB])*(potent[VZCOR][iB]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iB] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iB] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];


					}
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					doublereal dspeed = sqrt((potent[VXCOR][iB2])*(potent[VXCOR][iB2]) + (potent[VYCOR][iB2])*(potent[VYCOR][iB2]) + (potent[VZCOR][iB2])*(potent[VZCOR][iB2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iB2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB2] = 0.0;

					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iB2] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB2] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB2] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];

					}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					doublereal dspeed = sqrt((potent[VXCOR][iB3])*(potent[VXCOR][iB3]) + (potent[VYCOR][iB3])*(potent[VYCOR][iB3]) + (potent[VZCOR][iB3])*(potent[VZCOR][iB3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iB3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB3] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iB3] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB3] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB3] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];
					}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					doublereal dspeed = sqrt((potent[VXCOR][iB4])*(potent[VXCOR][iB4]) + (potent[VYCOR][iB4])*(potent[VYCOR][iB4]) + (potent[VZCOR][iB4])*(potent[VZCOR][iB4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY][iB4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB4] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY][iB4] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY][iB4] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY][iB4] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP];

					}
				}
			}

		}
		else
		{

			// граничные узлы.
			// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

			if (bE) {
				potent[GRADXTURBULENT_KINETIK_ENERGY][iE] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP] + (dxe / dxw)*(potent[GRADXTURBULENT_KINETIK_ENERGY][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY][iW]);
				potent[GRADYTURBULENT_KINETIK_ENERGY][iE] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP] + (dxe / dxw)*(potent[GRADYTURBULENT_KINETIK_ENERGY][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY][iW]);
				potent[GRADZTURBULENT_KINETIK_ENERGY][iE] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP] + (dxe / dxw)*(potent[GRADZTURBULENT_KINETIK_ENERGY][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY][iW]);

			}

			if (bW) {
				potent[GRADXTURBULENT_KINETIK_ENERGY][iW] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP] + (dxw / dxe)*(potent[GRADXTURBULENT_KINETIK_ENERGY][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY][iE]);
				potent[GRADYTURBULENT_KINETIK_ENERGY][iW] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP] + (dxw / dxe)*(potent[GRADYTURBULENT_KINETIK_ENERGY][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY][iE]);
				potent[GRADZTURBULENT_KINETIK_ENERGY][iW] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP] + (dxw / dxe)*(potent[GRADZTURBULENT_KINETIK_ENERGY][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY][iE]);

			}

			if (bN) {
				potent[GRADXTURBULENT_KINETIK_ENERGY][iN] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP] + (dyn / dys)*(potent[GRADXTURBULENT_KINETIK_ENERGY][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY][iS]);
				potent[GRADYTURBULENT_KINETIK_ENERGY][iN] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP] + (dyn / dys)*(potent[GRADYTURBULENT_KINETIK_ENERGY][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY][iS]);
				potent[GRADZTURBULENT_KINETIK_ENERGY][iN] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP] + (dyn / dys)*(potent[GRADZTURBULENT_KINETIK_ENERGY][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY][iS]);

			}

			if (bS) {
				potent[GRADXTURBULENT_KINETIK_ENERGY][iS] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP] + (dys / dyn)*(potent[GRADXTURBULENT_KINETIK_ENERGY][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY][iN]);
				potent[GRADYTURBULENT_KINETIK_ENERGY][iS] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP] + (dys / dyn)*(potent[GRADYTURBULENT_KINETIK_ENERGY][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY][iN]);
				potent[GRADZTURBULENT_KINETIK_ENERGY][iS] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP] + (dys / dyn)*(potent[GRADZTURBULENT_KINETIK_ENERGY][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY][iN]);

			}

			if (bT) {
				potent[GRADXTURBULENT_KINETIK_ENERGY][iT] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP] + (dzt / dzb)*(potent[GRADXTURBULENT_KINETIK_ENERGY][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY][iB]);
				potent[GRADYTURBULENT_KINETIK_ENERGY][iT] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP] + (dzt / dzb)*(potent[GRADYTURBULENT_KINETIK_ENERGY][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY][iB]);
				potent[GRADZTURBULENT_KINETIK_ENERGY][iT] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP] + (dzt / dzb)*(potent[GRADZTURBULENT_KINETIK_ENERGY][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY][iB]);

			}

			if (bB) {
				potent[GRADXTURBULENT_KINETIK_ENERGY][iB] = potent[GRADXTURBULENT_KINETIK_ENERGY][iP] + (dzb / dzt)*(potent[GRADXTURBULENT_KINETIK_ENERGY][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY][iT]);
				potent[GRADYTURBULENT_KINETIK_ENERGY][iB] = potent[GRADYTURBULENT_KINETIK_ENERGY][iP] + (dzb / dzt)*(potent[GRADYTURBULENT_KINETIK_ENERGY][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY][iT]);
				potent[GRADZTURBULENT_KINETIK_ENERGY][iB] = potent[GRADZTURBULENT_KINETIK_ENERGY][iP] + (dzb / dzt)*(potent[GRADZTURBULENT_KINETIK_ENERGY][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY][iT]);

			}

		}
	}

} // green_gauss_turbulent_kinetik_energy_MenterSST

  // 14.04.2019; 03.10.2019; 24.10.2019. Работает на АЛИС сетке.
  // Вычисление градиентов кинетической энергии турбулентных пульсаций
  // в центрах внутренних КО
  // и на границах с помощью линейной интерполяции.
  // Поскольку интерполяция линейная то точность данной формулы O(h). 
  // По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
  // что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
  // совпадает со взвешенным методом наименьших квадратов.
void green_gauss_turbulent_kinetik_energy_standart_k_epsilon(integer iP,
	doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND* &border_neighbor, integer *ilevel_alice)
{


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу.
	//bool bRCh = false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
				
		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
				
		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe = 0.5*dx, dxw = 0.5*dx, dyn = 0.5*dy, dys = 0.5*dy, dzt = 0.5*dz, dzb = 0.5*dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5*dx / dxe;
	fwplus = 0.5*dx / dxw;
	// y-direction
	fnplus = 0.5*dy / dyn;
	fsplus = 0.5*dy / dys;
	// z-direction
	ftplus = 0.5*dz / dzt;
	fbplus = 0.5*dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}





	// линейная интерполяция скорости VX на грань КО.
	doublereal fe = 0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// TURBULENT_KINETIK_ENERGY
		if (iE > -1) {
			if (!bE) {
				fe = feplus * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] + (1.0 - feplus)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fe = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE];
		}
		if (iW > -1) {
			if (!bW) {
				fw = fwplus * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] + (1.0 - fwplus)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fw = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW];
		}

		if (iE2 > -1) {
			if (!bE2) {
				fe2 = feplus2 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] + (1.0 - feplus2)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fe2 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) {
				fw2 = fwplus2 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] + (1.0 - fwplus2)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fw2 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2];
		}

		if (iE3 > -1) {
			if (!bE3) {
				fe3 = feplus3 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] + (1.0 - feplus3)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fe3 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) {
				fw3 = fwplus3 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] + (1.0 - fwplus3)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fw3 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3];
		}

		if (iE4 > -1) {
			if (!bE4) {
				fe4 = feplus4 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] + (1.0 - feplus4)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fe4 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) {
				fw4 = fwplus4 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] + (1.0 - fwplus4)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
			}
			else fw4 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4];
		}
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] + (1.0 - fnplus)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fn = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] + (1.0 - fsplus)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fs = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] + (1.0 - ftplus)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else ft = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] + (1.0 - fbplus)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fb = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB];
		}

		if (iN2 > -1) {
			if (!bN2) fn2 = fnplus2 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] + (1.0 - fnplus2)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fn2 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) fs2 = fsplus2 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] + (1.0 - fsplus2)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fs2 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) ft2 = ftplus2 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] + (1.0 - ftplus2)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else ft2 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) fb2 = fbplus2 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] + (1.0 - fbplus2)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fb2 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2];
		}

		if (iN3 > -1) {
			if (!bN3) fn3 = fnplus3 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] + (1.0 - fnplus3)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fn3 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) fs3 = fsplus3 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] + (1.0 - fsplus3)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fs3 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) ft3 = ftplus3 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] + (1.0 - ftplus3)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else ft3 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) fb3 = fbplus3 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] + (1.0 - fbplus3)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fb3 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3];
		}

		if (iN4 > -1) {
			if (!bN4) fn4 = fnplus4 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] + (1.0 - fnplus4)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fn4 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) fs4 = fsplus4 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] + (1.0 - fsplus4)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fs4 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) ft4 = ftplus4 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] + (1.0 - ftplus4)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else ft4 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) fb4 = fbplus4 * potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] + (1.0 - fbplus4)*potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]; else fb4 = potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4];
		}
		// градиент NUSHA
		//potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]=(fe-fw)/dx;
		//potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]=(fn-fs)/dy;
		//potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP]=(ft-fb)/dz;
		potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;


	}
	else {

		if (1) {



			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE])*(potent[VXCOR][iE]) + (potent[VYCOR][iE])*(potent[VYCOR][iE]) + (potent[VZCOR][iE])*(potent[VZCOR][iE]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = 0.0;



						/*
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] = 0.0;


						*/
					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}



			if (iE2 > -1) {
				if (bE2) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE2])*(potent[VXCOR][iE2]) + (potent[VYCOR][iE2])*(potent[VYCOR][iE2]) + (potent[VZCOR][iE2])*(potent[VZCOR][iE2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE2] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iE3 > -1) {
				if (bE3) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE3])*(potent[VXCOR][iE3]) + (potent[VYCOR][iE3])*(potent[VYCOR][iE3]) + (potent[VZCOR][iE3])*(potent[VZCOR][iE3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE3] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iE4 > -1) {
				if (bE4) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE4])*(potent[VXCOR][iE4]) + (potent[VYCOR][iE4])*(potent[VYCOR][iE4]) + (potent[VZCOR][iE4])*(potent[VZCOR][iE4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE4] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iW > -1) {
				if (bW) {

					doublereal dspeed = sqrt((potent[VXCOR][iW])*(potent[VXCOR][iW]) + (potent[VYCOR][iW])*(potent[VYCOR][iW]) + (potent[VZCOR][iW])*(potent[VZCOR][iW]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iW2 > -1) {
				if (bW2) {

					doublereal dspeed = sqrt((potent[VXCOR][iW2])*(potent[VXCOR][iW2]) + (potent[VYCOR][iW2])*(potent[VYCOR][iW2]) + (potent[VZCOR][iW2])*(potent[VZCOR][iW2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW2] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					doublereal dspeed = sqrt((potent[VXCOR][iW3])*(potent[VXCOR][iW3]) + (potent[VYCOR][iW3])*(potent[VYCOR][iW3]) + (potent[VZCOR][iW3])*(potent[VZCOR][iW3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW3] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					doublereal dspeed = sqrt((potent[VXCOR][iW4])*(potent[VXCOR][iW4]) + (potent[VYCOR][iW4])*(potent[VYCOR][iW4]) + (potent[VZCOR][iW4])*(potent[VZCOR][iW4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW4] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iN > -1) {
				if (bN) {

					doublereal dspeed = sqrt((potent[VXCOR][iN])*(potent[VXCOR][iN]) + (potent[VYCOR][iN])*(potent[VYCOR][iN]) + (potent[VZCOR][iN])*(potent[VZCOR][iN]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iN2 > -1) {
				if (bN2) {

					doublereal dspeed = sqrt((potent[VXCOR][iN2])*(potent[VXCOR][iN2]) + (potent[VYCOR][iN2])*(potent[VYCOR][iN2]) + (potent[VZCOR][iN2])*(potent[VZCOR][iN2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN2] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					doublereal dspeed = sqrt((potent[VXCOR][iN3])*(potent[VXCOR][iN3]) + (potent[VYCOR][iN3])*(potent[VYCOR][iN3]) + (potent[VZCOR][iN3])*(potent[VZCOR][iN3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN3] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					doublereal dspeed = sqrt((potent[VXCOR][iN4])*(potent[VXCOR][iN4]) + (potent[VYCOR][iN4])*(potent[VYCOR][iN4]) + (potent[VZCOR][iN4])*(potent[VZCOR][iN4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN4] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iS > -1) {
				if (bS) {

					doublereal dspeed = sqrt((potent[VXCOR][iS])*(potent[VXCOR][iS]) + (potent[VYCOR][iS])*(potent[VYCOR][iS]) + (potent[VZCOR][iS])*(potent[VZCOR][iS]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iS2 > -1) {
				if (bS2) {

					doublereal dspeed = sqrt((potent[VXCOR][iS2])*(potent[VXCOR][iS2]) + (potent[VYCOR][iS2])*(potent[VYCOR][iS2]) + (potent[VZCOR][iS2])*(potent[VZCOR][iS2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS2] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					doublereal dspeed = sqrt((potent[VXCOR][iS3])*(potent[VXCOR][iS3]) + (potent[VYCOR][iS3])*(potent[VYCOR][iS3]) + (potent[VZCOR][iS3])*(potent[VZCOR][iS3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS3] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					doublereal dspeed = sqrt((potent[VXCOR][iS4])*(potent[VXCOR][iS4]) + (potent[VYCOR][iS4])*(potent[VYCOR][iS4]) + (potent[VZCOR][iS4])*(potent[VZCOR][iS4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS4] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iT > -1) {
				if (bT) {

					doublereal dspeed = sqrt((potent[VXCOR][iT])*(potent[VXCOR][iT]) + (potent[VYCOR][iT])*(potent[VYCOR][iT]) + (potent[VZCOR][iT])*(potent[VZCOR][iT]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iT2 > -1) {
				if (bT2) {

					doublereal dspeed = sqrt((potent[VXCOR][iT2])*(potent[VXCOR][iT2]) + (potent[VYCOR][iT2])*(potent[VYCOR][iT2]) + (potent[VZCOR][iT2])*(potent[VZCOR][iT2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT2] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					doublereal dspeed = sqrt((potent[VXCOR][iT3])*(potent[VXCOR][iT3]) + (potent[VYCOR][iT3])*(potent[VYCOR][iT3]) + (potent[VZCOR][iT3])*(potent[VZCOR][iT3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT3] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					doublereal dspeed = sqrt((potent[VXCOR][iT4])*(potent[VXCOR][iT4]) + (potent[VYCOR][iT4])*(potent[VYCOR][iT4]) + (potent[VZCOR][iT4])*(potent[VZCOR][iT4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT4] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iB > -1) {
				if (bB) {
					doublereal dspeed = sqrt((potent[VXCOR][iB])*(potent[VXCOR][iB]) + (potent[VYCOR][iB])*(potent[VYCOR][iB]) + (potent[VZCOR][iB])*(potent[VZCOR][iB]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];


					}
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					doublereal dspeed = sqrt((potent[VXCOR][iB2])*(potent[VXCOR][iB2]) + (potent[VYCOR][iB2])*(potent[VYCOR][iB2]) + (potent[VZCOR][iB2])*(potent[VZCOR][iB2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] = 0.0;

					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB2] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];

					}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					doublereal dspeed = sqrt((potent[VXCOR][iB3])*(potent[VXCOR][iB3]) + (potent[VYCOR][iB3])*(potent[VYCOR][iB3]) + (potent[VZCOR][iB3])*(potent[VZCOR][iB3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB3] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
					}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					doublereal dspeed = sqrt((potent[VXCOR][iB4])*(potent[VXCOR][iB4]) + (potent[VYCOR][iB4])*(potent[VYCOR][iB4]) + (potent[VZCOR][iB4])*(potent[VZCOR][iB4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] = 0.0;
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] = 0.0;
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];
						potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB4] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP];

					}
				}
			}

		}
		else
		{

			// граничные узлы.
			// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

			if (bE) {
				potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dxe / dxw)*(potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW]);
				potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dxe / dxw)*(potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW]);
				potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dxe / dxw)*(potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW]);

			}

			if (bW) {
				potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dxw / dxe)*(potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE]);
				potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dxw / dxe)*(potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE]);
				potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iW] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dxw / dxe)*(potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iE]);

			}

			if (bN) {
				potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dyn / dys)*(potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS]);
				potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dyn / dys)*(potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS]);
				potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dyn / dys)*(potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS]);

			}

			if (bS) {
				potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dys / dyn)*(potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN]);
				potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dys / dyn)*(potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN]);
				potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iS] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dys / dyn)*(potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iN]);

			}

			if (bT) {
				potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dzt / dzb)*(potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB]);
				potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dzt / dzb)*(potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB]);
				potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dzt / dzb)*(potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB]);

			}

			if (bB) {
				potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dzb / dzt)*(potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT]);
				potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dzb / dzt)*(potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT]);
				potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iB] = potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] + (dzb / dzt)*(potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iP] - potent[GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS][iT]);

			}

		}
	}

} // green_gauss_turbulent_kinetik_energy_standart_k_epsilon

  // 14.04.2019; 03.10.2019; 24.10.2019. Работает на АЛИС сетке.
  // Вычисление градиентов скорости диссипации кинетической энергии турбулентных пульсаций
  // в центрах внутренних КО
  // и на границах с помощью линейной интерполяции.
  // Поскольку интерполяция линейная то точность данной формулы O(h). 
  // По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
  // что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
  // совпадает со взвешенным методом наименьших квадратов.
void green_gauss_turbulent_dissipation_rate_epsilon_standart_k_epsilon(integer iP,
	doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND* &border_neighbor, integer *ilevel_alice)
{


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу.
	//bool bRCh = false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {

		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe = 0.5*dx, dxw = 0.5*dx, dyn = 0.5*dy, dys = 0.5*dy, dzt = 0.5*dz, dzb = 0.5*dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5*dx / dxe;
	fwplus = 0.5*dx / dxw;
	// y-direction
	fnplus = 0.5*dy / dyn;
	fsplus = 0.5*dy / dys;
	// z-direction
	ftplus = 0.5*dz / dzt;
	fbplus = 0.5*dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}





	// линейная интерполяция скорости VX на грань КО.
	doublereal fe = 0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// TURBULENT_KINETIK_ENERGY
		if (iE > -1) {
			if (!bE) {
				fe = feplus * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] + (1.0 - feplus)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fe = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE];
		}
		if (iW > -1) {
			if (!bW) {
				fw = fwplus * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] + (1.0 - fwplus)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fw = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW];
		}

		if (iE2 > -1) {
			if (!bE2) {
				fe2 = feplus2 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] + (1.0 - feplus2)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fe2 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) {
				fw2 = fwplus2 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] + (1.0 - fwplus2)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fw2 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2];
		}

		if (iE3 > -1) {
			if (!bE3) {
				fe3 = feplus3 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] + (1.0 - feplus3)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fe3 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) {
				fw3 = fwplus3 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] + (1.0 - fwplus3)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fw3 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3];
		}

		if (iE4 > -1) {
			if (!bE4) {
				fe4 = feplus4 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] + (1.0 - feplus4)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fe4 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) {
				fw4 = fwplus4 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] + (1.0 - fwplus4)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
			}
			else fw4 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4];
		}
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] + (1.0 - fnplus)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fn = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] + (1.0 - fsplus)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fs = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] + (1.0 - ftplus)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else ft = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] + (1.0 - fbplus)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fb = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB];
		}

		if (iN2 > -1) {
			if (!bN2) fn2 = fnplus2 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] + (1.0 - fnplus2)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fn2 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) fs2 = fsplus2 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] + (1.0 - fsplus2)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fs2 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) ft2 = ftplus2 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] + (1.0 - ftplus2)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else ft2 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) fb2 = fbplus2 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] + (1.0 - fbplus2)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fb2 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2];
		}

		if (iN3 > -1) {
			if (!bN3) fn3 = fnplus3 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] + (1.0 - fnplus3)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fn3 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) fs3 = fsplus3 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] + (1.0 - fsplus3)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fs3 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) ft3 = ftplus3 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] + (1.0 - ftplus3)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else ft3 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) fb3 = fbplus3 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] + (1.0 - fbplus3)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fb3 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3];
		}

		if (iN4 > -1) {
			if (!bN4) fn4 = fnplus4 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] + (1.0 - fnplus4)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fn4 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) fs4 = fsplus4 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] + (1.0 - fsplus4)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fs4 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) ft4 = ftplus4 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] + (1.0 - ftplus4)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else ft4 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) fb4 = fbplus4 * potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] + (1.0 - fbplus4)*potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]; else fb4 = potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4];
		}
		// градиент NUSHA
		//potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]=(fe-fw)/dx;
		//potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]=(fn-fs)/dy;
		//potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP]=(ft-fb)/dz;
		potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;


	}
	else {

		if (1) {



			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE])*(potent[VXCOR][iE]) + (potent[VYCOR][iE])*(potent[VYCOR][iE]) + (potent[VZCOR][iE])*(potent[VZCOR][iE]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = 0.0;



						/*
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] = 0.0;


						*/
					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}



			if (iE2 > -1) {
				if (bE2) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE2])*(potent[VXCOR][iE2]) + (potent[VYCOR][iE2])*(potent[VYCOR][iE2]) + (potent[VZCOR][iE2])*(potent[VZCOR][iE2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE2] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iE3 > -1) {
				if (bE3) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE3])*(potent[VXCOR][iE3]) + (potent[VYCOR][iE3])*(potent[VYCOR][iE3]) + (potent[VZCOR][iE3])*(potent[VZCOR][iE3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE3] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iE4 > -1) {
				if (bE4) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE4])*(potent[VXCOR][iE4]) + (potent[VYCOR][iE4])*(potent[VYCOR][iE4]) + (potent[VZCOR][iE4])*(potent[VZCOR][iE4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE4] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iW > -1) {
				if (bW) {

					doublereal dspeed = sqrt((potent[VXCOR][iW])*(potent[VXCOR][iW]) + (potent[VYCOR][iW])*(potent[VYCOR][iW]) + (potent[VZCOR][iW])*(potent[VZCOR][iW]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iW2 > -1) {
				if (bW2) {

					doublereal dspeed = sqrt((potent[VXCOR][iW2])*(potent[VXCOR][iW2]) + (potent[VYCOR][iW2])*(potent[VYCOR][iW2]) + (potent[VZCOR][iW2])*(potent[VZCOR][iW2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW2] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					doublereal dspeed = sqrt((potent[VXCOR][iW3])*(potent[VXCOR][iW3]) + (potent[VYCOR][iW3])*(potent[VYCOR][iW3]) + (potent[VZCOR][iW3])*(potent[VZCOR][iW3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW3] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					doublereal dspeed = sqrt((potent[VXCOR][iW4])*(potent[VXCOR][iW4]) + (potent[VYCOR][iW4])*(potent[VYCOR][iW4]) + (potent[VZCOR][iW4])*(potent[VZCOR][iW4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW4] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iN > -1) {
				if (bN) {

					doublereal dspeed = sqrt((potent[VXCOR][iN])*(potent[VXCOR][iN]) + (potent[VYCOR][iN])*(potent[VYCOR][iN]) + (potent[VZCOR][iN])*(potent[VZCOR][iN]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iN2 > -1) {
				if (bN2) {

					doublereal dspeed = sqrt((potent[VXCOR][iN2])*(potent[VXCOR][iN2]) + (potent[VYCOR][iN2])*(potent[VYCOR][iN2]) + (potent[VZCOR][iN2])*(potent[VZCOR][iN2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN2] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					doublereal dspeed = sqrt((potent[VXCOR][iN3])*(potent[VXCOR][iN3]) + (potent[VYCOR][iN3])*(potent[VYCOR][iN3]) + (potent[VZCOR][iN3])*(potent[VZCOR][iN3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN3] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					doublereal dspeed = sqrt((potent[VXCOR][iN4])*(potent[VXCOR][iN4]) + (potent[VYCOR][iN4])*(potent[VYCOR][iN4]) + (potent[VZCOR][iN4])*(potent[VZCOR][iN4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN4] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iS > -1) {
				if (bS) {

					doublereal dspeed = sqrt((potent[VXCOR][iS])*(potent[VXCOR][iS]) + (potent[VYCOR][iS])*(potent[VYCOR][iS]) + (potent[VZCOR][iS])*(potent[VZCOR][iS]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iS2 > -1) {
				if (bS2) {

					doublereal dspeed = sqrt((potent[VXCOR][iS2])*(potent[VXCOR][iS2]) + (potent[VYCOR][iS2])*(potent[VYCOR][iS2]) + (potent[VZCOR][iS2])*(potent[VZCOR][iS2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS2] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					doublereal dspeed = sqrt((potent[VXCOR][iS3])*(potent[VXCOR][iS3]) + (potent[VYCOR][iS3])*(potent[VYCOR][iS3]) + (potent[VZCOR][iS3])*(potent[VZCOR][iS3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS3] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					doublereal dspeed = sqrt((potent[VXCOR][iS4])*(potent[VXCOR][iS4]) + (potent[VYCOR][iS4])*(potent[VYCOR][iS4]) + (potent[VZCOR][iS4])*(potent[VZCOR][iS4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS4] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iT > -1) {
				if (bT) {

					doublereal dspeed = sqrt((potent[VXCOR][iT])*(potent[VXCOR][iT]) + (potent[VYCOR][iT])*(potent[VYCOR][iT]) + (potent[VZCOR][iT])*(potent[VZCOR][iT]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iT2 > -1) {
				if (bT2) {

					doublereal dspeed = sqrt((potent[VXCOR][iT2])*(potent[VXCOR][iT2]) + (potent[VYCOR][iT2])*(potent[VYCOR][iT2]) + (potent[VZCOR][iT2])*(potent[VZCOR][iT2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT2] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					doublereal dspeed = sqrt((potent[VXCOR][iT3])*(potent[VXCOR][iT3]) + (potent[VYCOR][iT3])*(potent[VYCOR][iT3]) + (potent[VZCOR][iT3])*(potent[VZCOR][iT3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT3] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					doublereal dspeed = sqrt((potent[VXCOR][iT4])*(potent[VXCOR][iT4]) + (potent[VYCOR][iT4])*(potent[VYCOR][iT4]) + (potent[VZCOR][iT4])*(potent[VZCOR][iT4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT4] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iB > -1) {
				if (bB) {
					doublereal dspeed = sqrt((potent[VXCOR][iB])*(potent[VXCOR][iB]) + (potent[VYCOR][iB])*(potent[VYCOR][iB]) + (potent[VZCOR][iB])*(potent[VZCOR][iB]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];


					}
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					doublereal dspeed = sqrt((potent[VXCOR][iB2])*(potent[VXCOR][iB2]) + (potent[VYCOR][iB2])*(potent[VYCOR][iB2]) + (potent[VZCOR][iB2])*(potent[VZCOR][iB2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] = 0.0;

					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB2] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];

					}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					doublereal dspeed = sqrt((potent[VXCOR][iB3])*(potent[VXCOR][iB3]) + (potent[VYCOR][iB3])*(potent[VYCOR][iB3]) + (potent[VZCOR][iB3])*(potent[VZCOR][iB3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB3] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
					}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					doublereal dspeed = sqrt((potent[VXCOR][iB4])*(potent[VXCOR][iB4]) + (potent[VYCOR][iB4])*(potent[VYCOR][iB4]) + (potent[VZCOR][iB4])*(potent[VZCOR][iB4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] = 0.0;
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] = 0.0;
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];
						potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB4] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP];

					}
				}
			}

		}
		else
		{

			// граничные узлы.
			// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

			if (bE) {
				potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dxe / dxw)*(potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW]);
				potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dxe / dxw)*(potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW]);
				potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dxe / dxw)*(potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW]);

			}

			if (bW) {
				potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dxw / dxe)*(potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE]);
				potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dxw / dxe)*(potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE]);
				potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iW] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dxw / dxe)*(potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iE]);

			}

			if (bN) {
				potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dyn / dys)*(potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS]);
				potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dyn / dys)*(potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS]);
				potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dyn / dys)*(potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS]);

			}

			if (bS) {
				potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dys / dyn)*(potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN]);
				potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dys / dyn)*(potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN]);
				potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iS] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dys / dyn)*(potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iN]);

			}

			if (bT) {
				potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dzt / dzb)*(potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB]);
				potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dzt / dzb)*(potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB]);
				potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dzt / dzb)*(potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB]);

			}

			if (bB) {
				potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dzb / dzt)*(potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT]);
				potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dzb / dzt)*(potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT]);
				potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iB] = potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] + (dzb / dzt)*(potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iP] - potent[GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][iT]);

			}

		}
	}

} // green_gauss_turbulent_dissipation_rate_epsilon_standart_k_epsilon


// 14.04.2019; 03.10.2019. Работает на АЛИС сетке.
// Сделано по полному образцу кинетической энергии турбулентных пульсаций и 
// модели Спаларта Аллмареса.
// Вычисление градиентов удельной скорости диссипации кинетической энергии турбулентных пульсаций
// в центрах внутренних КО в модели SST Ментера.
// и на границах с помощью линейной интерполяции.
// Поскольку интерполяция линейная то точность данной формулы O(h). 
// По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
// что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
// совпадает со взвешенным методом наименьших квадратов.
void green_gauss_specific_dissipation_rate_omega_MenterSST(integer iP,
	doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND* &border_neighbor, integer *ilevel_alice) {


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу.
	//bool bRCh = false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}
	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe = 0.5*dx, dxw = 0.5*dx, dyn = 0.5*dy, dys = 0.5*dy, dzt = 0.5*dz, dzb = 0.5*dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus = 0.5*dx / dxe;
	fwplus = 0.5*dx / dxw;
	// y-direction
	fnplus = 0.5*dy / dyn;
	fsplus = 0.5*dy / dys;
	// z-direction
	ftplus = 0.5*dz / dzt;
	fbplus = 0.5*dz / dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}





	// линейная интерполяция скорости VX на грань КО.
	doublereal fe = 0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// VX
		if (iE > -1) {
			if (!bE) {
				fe = feplus * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] + (1.0 - feplus)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fe = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE];
		}
		if (iW > -1) {
			if (!bW) {
				fw = fwplus * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] + (1.0 - fwplus)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fw = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW];
		}

		if (iE2 > -1) {
			if (!bE2) {
				fe2 = feplus2 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] + (1.0 - feplus2)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fe2 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) {
				fw2 = fwplus2 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] + (1.0 - fwplus2)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fw2 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2];
		}

		if (iE3 > -1) {
			if (!bE3) {
				fe3 = feplus3 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] + (1.0 - feplus3)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fe3 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) {
				fw3 = fwplus3 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] + (1.0 - fwplus3)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fw3 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3];
		}

		if (iE4 > -1) {
			if (!bE4) {
				fe4 = feplus4 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] + (1.0 - feplus4)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fe4 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) {
				fw4 = fwplus4 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] + (1.0 - fwplus4)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
			}
			else fw4 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4];
		}
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] + (1.0 - fnplus)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fn = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] + (1.0 - fsplus)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fs = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] + (1.0 - ftplus)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else ft = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] + (1.0 - fbplus)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fb = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB];
		}

		if (iN2 > -1) {
			if (!bN2) fn2 = fnplus2 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] + (1.0 - fnplus2)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fn2 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) fs2 = fsplus2 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] + (1.0 - fsplus2)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fs2 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) ft2 = ftplus2 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] + (1.0 - ftplus2)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else ft2 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) fb2 = fbplus2 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] + (1.0 - fbplus2)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fb2 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2];
		}

		if (iN3 > -1) {
			if (!bN3) fn3 = fnplus3 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] + (1.0 - fnplus3)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fn3 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) fs3 = fsplus3 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] + (1.0 - fsplus3)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fs3 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) ft3 = ftplus3 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] + (1.0 - ftplus3)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else ft3 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) fb3 = fbplus3 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] + (1.0 - fbplus3)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fb3 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3];
		}

		if (iN4 > -1) {
			if (!bN4) fn4 = fnplus4 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] + (1.0 - fnplus4)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fn4 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) fs4 = fsplus4 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] + (1.0 - fsplus4)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fs4 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) ft4 = ftplus4 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] + (1.0 - ftplus4)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else ft4 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) fb4 = fbplus4 * potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] + (1.0 - fbplus4)*potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]; else fb4 = potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4];
		}
		// градиент NUSHA
		//potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]=(fe-fw)/dx;
		//potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]=(fn-fs)/dy;
		//potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP]=(ft-fb)/dz;
		potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] = (fe*dSqe / (dy*dz) + fe2 * dSqe2 / (dy*dz) + fe3 * dSqe3 / (dy*dz) + fe4 * dSqe4 / (dy*dz) - (fw*dSqw / (dy*dz) + fw2 * dSqw2 / (dy*dz) + fw3 * dSqw3 / (dy*dz) + fw4 * dSqw4 / (dy*dz))) / dx;
		potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] = (fn*dSqn / (dx*dz) + fn2 * dSqn2 / (dx*dz) + fn3 * dSqn3 / (dx*dz) + fn4 * dSqn4 / (dx*dz) - (fs*dSqs / (dx*dz) + fs2 * dSqs2 / (dx*dz) + fs3 * dSqs3 / (dx*dz) + fs4 * dSqs4 / (dx*dz))) / dy;
		potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] = (ft*dSqt / (dx*dy) + ft2 * dSqt2 / (dx*dy) + ft3 * dSqt3 / (dx*dy) + ft4 * dSqt4 / (dx*dy) - (fb*dSqb / (dx*dy) + fb2 * dSqb2 / (dx*dy) + fb3 * dSqb3 / (dx*dy) + fb4 * dSqb4 / (dx*dy))) / dz;


	}
	else {

		if (1) {



			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE])*(potent[VXCOR][iE]) + (potent[VYCOR][iE])*(potent[VYCOR][iE]) + (potent[VZCOR][iE])*(potent[VZCOR][iE]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = 0.0;



						/*
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] = 0.0;


						*/
					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}



			if (iE2 > -1) {
				if (bE2) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE2])*(potent[VXCOR][iE2]) + (potent[VYCOR][iE2])*(potent[VYCOR][iE2]) + (potent[VZCOR][iE2])*(potent[VZCOR][iE2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE2] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iE3 > -1) {
				if (bE3) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE3])*(potent[VXCOR][iE3]) + (potent[VYCOR][iE3])*(potent[VYCOR][iE3]) + (potent[VZCOR][iE3])*(potent[VZCOR][iE3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE3] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iE4 > -1) {
				if (bE4) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					doublereal dspeed = sqrt((potent[VXCOR][iE4])*(potent[VXCOR][iE4]) + (potent[VYCOR][iE4])*(potent[VYCOR][iE4]) + (potent[VZCOR][iE4])*(potent[VZCOR][iE4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] = 0.0;




					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE4] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iW > -1) {
				if (bW) {

					doublereal dspeed = sqrt((potent[VXCOR][iW])*(potent[VXCOR][iW]) + (potent[VYCOR][iW])*(potent[VYCOR][iW]) + (potent[VZCOR][iW])*(potent[VZCOR][iW]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iW2 > -1) {
				if (bW2) {

					doublereal dspeed = sqrt((potent[VXCOR][iW2])*(potent[VXCOR][iW2]) + (potent[VYCOR][iW2])*(potent[VYCOR][iW2]) + (potent[VZCOR][iW2])*(potent[VZCOR][iW2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW2] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					doublereal dspeed = sqrt((potent[VXCOR][iW3])*(potent[VXCOR][iW3]) + (potent[VYCOR][iW3])*(potent[VYCOR][iW3]) + (potent[VZCOR][iW3])*(potent[VZCOR][iW3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW3] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					doublereal dspeed = sqrt((potent[VXCOR][iW4])*(potent[VXCOR][iW4]) + (potent[VYCOR][iW4])*(potent[VYCOR][iW4]) + (potent[VZCOR][iW4])*(potent[VZCOR][iW4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW4] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iN > -1) {
				if (bN) {

					doublereal dspeed = sqrt((potent[VXCOR][iN])*(potent[VXCOR][iN]) + (potent[VYCOR][iN])*(potent[VYCOR][iN]) + (potent[VZCOR][iN])*(potent[VZCOR][iN]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iN2 > -1) {
				if (bN2) {

					doublereal dspeed = sqrt((potent[VXCOR][iN2])*(potent[VXCOR][iN2]) + (potent[VYCOR][iN2])*(potent[VYCOR][iN2]) + (potent[VZCOR][iN2])*(potent[VZCOR][iN2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN2] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					doublereal dspeed = sqrt((potent[VXCOR][iN3])*(potent[VXCOR][iN3]) + (potent[VYCOR][iN3])*(potent[VYCOR][iN3]) + (potent[VZCOR][iN3])*(potent[VZCOR][iN3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN3] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					doublereal dspeed = sqrt((potent[VXCOR][iN4])*(potent[VXCOR][iN4]) + (potent[VYCOR][iN4])*(potent[VYCOR][iN4]) + (potent[VZCOR][iN4])*(potent[VZCOR][iN4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN4] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iS > -1) {
				if (bS) {

					doublereal dspeed = sqrt((potent[VXCOR][iS])*(potent[VXCOR][iS]) + (potent[VYCOR][iS])*(potent[VYCOR][iS]) + (potent[VZCOR][iS])*(potent[VZCOR][iS]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iS2 > -1) {
				if (bS2) {

					doublereal dspeed = sqrt((potent[VXCOR][iS2])*(potent[VXCOR][iS2]) + (potent[VYCOR][iS2])*(potent[VYCOR][iS2]) + (potent[VZCOR][iS2])*(potent[VZCOR][iS2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS2] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					doublereal dspeed = sqrt((potent[VXCOR][iS3])*(potent[VXCOR][iS3]) + (potent[VYCOR][iS3])*(potent[VYCOR][iS3]) + (potent[VZCOR][iS3])*(potent[VZCOR][iS3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS3] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					doublereal dspeed = sqrt((potent[VXCOR][iS4])*(potent[VXCOR][iS4]) + (potent[VYCOR][iS4])*(potent[VYCOR][iS4]) + (potent[VZCOR][iS4])*(potent[VZCOR][iS4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS4] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iT > -1) {
				if (bT) {

					doublereal dspeed = sqrt((potent[VXCOR][iT])*(potent[VXCOR][iT]) + (potent[VYCOR][iT])*(potent[VYCOR][iT]) + (potent[VZCOR][iT])*(potent[VZCOR][iT]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iT2 > -1) {
				if (bT2) {

					doublereal dspeed = sqrt((potent[VXCOR][iT2])*(potent[VXCOR][iT2]) + (potent[VYCOR][iT2])*(potent[VYCOR][iT2]) + (potent[VZCOR][iT2])*(potent[VZCOR][iT2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT2] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					doublereal dspeed = sqrt((potent[VXCOR][iT3])*(potent[VXCOR][iT3]) + (potent[VYCOR][iT3])*(potent[VYCOR][iT3]) + (potent[VZCOR][iT3])*(potent[VZCOR][iT3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT3] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					doublereal dspeed = sqrt((potent[VXCOR][iT4])*(potent[VXCOR][iT4]) + (potent[VYCOR][iT4])*(potent[VYCOR][iT4]) + (potent[VZCOR][iT4])*(potent[VZCOR][iT4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT4] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iB > -1) {
				if (bB) {
					doublereal dspeed = sqrt((potent[VXCOR][iB])*(potent[VXCOR][iB]) + (potent[VYCOR][iB])*(potent[VYCOR][iB]) + (potent[VZCOR][iB])*(potent[VZCOR][iB]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = 0.0;


					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					}
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					doublereal dspeed = sqrt((potent[VXCOR][iB2])*(potent[VXCOR][iB2]) + (potent[VYCOR][iB2])*(potent[VYCOR][iB2]) + (potent[VZCOR][iB2])*(potent[VZCOR][iB2]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] = 0.0;

					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB2] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];

					}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					doublereal dspeed = sqrt((potent[VXCOR][iB3])*(potent[VXCOR][iB3]) + (potent[VYCOR][iB3])*(potent[VYCOR][iB3]) + (potent[VZCOR][iB3])*(potent[VZCOR][iB3]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB3] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
					}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					doublereal dspeed = sqrt((potent[VXCOR][iB4])*(potent[VXCOR][iB4]) + (potent[VYCOR][iB4])*(potent[VYCOR][iB4]) + (potent[VZCOR][iB4])*(potent[VZCOR][iB4]));

					if (dspeed < 1.0e-10) {
						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] = 0.0;
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] = 0.0;
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] = 0.0;
					}
					else {

						potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
						potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB4] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];

					}
				}
			}

		}
		else
		{

			// граничные узлы.
			// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

			if (bE) {
				potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dxe / dxw)*(potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW]);
				potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dxe / dxw)*(potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW]);
				potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dxe / dxw)*(potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW]);

			}

			if (bW) {
				potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dxw / dxe)*(potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE]);
				potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dxw / dxe)*(potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE]);
				potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iW] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dxw / dxe)*(potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iE]);

			}

			if (bN) {
				potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dyn / dys)*(potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS]);
				potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dyn / dys)*(potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS]);
				potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dyn / dys)*(potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS]);

			}

			if (bS) {
				potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dys / dyn)*(potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN]);
				potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dys / dyn)*(potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN]);
				potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iS] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dys / dyn)*(potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iN]);

			}

			if (bT) {
				potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dzt / dzb)*(potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB]);
				potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dzt / dzb)*(potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB]);
				potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dzt / dzb)*(potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB]);

			}

			if (bB) {
				potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dzb / dzt)*(potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT]);
				potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dzb / dzt)*(potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT]);
				potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iB] = potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] + (dzb / dzt)*(potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP] - potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iT]);

			}

		}
	}

} // green_gauss_specific_dissipation_rate_omega_MenterSST


// вычисление градиентов поправки давления с помощью теоремы Грина-Гаусса. 
void green_gaussPAM(integer iP, doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	                int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
					BOUND* &border_neighbor, integer ls, integer lw, WALL* &w, bool bLRfree,
	                integer *ilevel_alice, int* ptr, TOCHKA*& volume_loc) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты поправки давления для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {

		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];


		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];


		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;


	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	//volume3D(iP, nvtx, pa, dx, dy, dz);
	//dx = fabs(dx);
	//dy = fabs(dy);
	//dz = fabs(dz);

	TOCHKA point_loc = volume_loc[iP];
	dx = point_loc.x;
	dy = point_loc.y;
	dz = point_loc.z;

	doublereal dxe = 0.5*dx, dxw = 0.5*dx, dyn = 0.5*dy, dys = 0.5*dy, dzt = 0.5*dz, dzb = 0.5*dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	if (b_on_adaptive_local_refinement_mesh) {

		// т.к. известна нумерация вершин куба, то здесь она используется
		// x - direction
		if (iE2 > -1) {
			if (!bE2) dxe2 = 0.5 * (pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
			if (!bE2) dxe2 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		}
		if (iW2 > -1) {
			if (!bW2) dxw2 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
			if (!bW2) dxw2 -= 0.5 * (pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
		}
		// y - direction
		if (iN2 > -1) {
			if (!bN2) dyn2 = 0.5 * (pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
			if (!bN2) dyn2 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		}
		if (iS2 > -1) {
			if (!bS2) dys2 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
			if (!bS2) dys2 -= 0.5 * (pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
		}
		// z - direction
		if (iT2 > -1) {
			if (!bT2) dzt2 = 0.5 * (pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
			if (!bT2) dzt2 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		}
		if (iB2 > -1) {
			if (!bB2) dzb2 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
			if (!bB2) dzb2 -= 0.5 * (pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
		}

		// т.к. известна нумерация вершин куба, то здесь она используется
		// x - direction
		if (iE3 > -1) {
			if (!bE3) dxe3 = 0.5 * (pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
			if (!bE3) dxe3 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		}
		if (iW3 > -1) {
			if (!bW3) dxw3 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
			if (!bW3) dxw3 -= 0.5 * (pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
		}
		// y - direction
		if (iN3 > -1) {
			if (!bN3) dyn3 = 0.5 * (pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
			if (!bN3) dyn3 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		}
		if (iS3 > -1) {
			if (!bS3) dys3 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
			if (!bS3) dys3 -= 0.5 * (pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
		}
		// z - direction
		if (iT3 > -1) {
			if (!bT3) dzt3 = 0.5 * (pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
			if (!bT3) dzt3 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		}
		if (iB3 > -1) {
			if (!bB3) dzb3 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
			if (!bB3) dzb3 -= 0.5 * (pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
		}

		// т.к. известна нумерация вершин куба, то здесь она используется
		// x - direction
		if (iE4 > -1) {
			if (!bE4) dxe4 = 0.5 * (pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
			if (!bE4) dxe4 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		}
		if (iW4 > -1) {
			if (!bW4) dxw4 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
			if (!bW4) dxw4 -= 0.5 * (pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
		}
		// y - direction
		if (iN4 > -1) {
			if (!bN4) dyn4 = 0.5 * (pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
			if (!bN4) dyn4 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		}
		if (iS4 > -1) {
			if (!bS4) dys4 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
			if (!bS4) dys4 -= 0.5 * (pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
		}
		// z - direction
		if (iT4 > -1) {
			if (!bT4) dzt4 = 0.5 * (pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
			if (!bT4) dzt4 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		}
		if (iB4 > -1) {
			if (!bB4) dzb4 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
			if (!bB4) dzb4 -= 0.5 * (pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
		}
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus=0.5*dx/dxe;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;


	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.

	


	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe = dy_loc * dz_loc;
			}
		}


		
	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw = dy_loc * dz_loc;
			}
		}

		
	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn = dx_loc * dz_loc;
			}
		}

		
	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs = dx_loc * dz_loc;
			}
		}

		
	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;


				dSqt = dx_loc * dy_loc;
			}
		}

		
	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb = dx_loc * dy_loc;
			}
		}

		
	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.


	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 > -1) {

			dSqe2 = dy * dz;

			if (bE2) {
				// граничный узел.
				dSqe2 = border_neighbor[iE2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
					dSqe2 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE2];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqe2 = dy_loc * dz_loc;
				}
			}


		}


		if (iW2 > -1) {
			dSqw2 = dy * dz;

			if (bW2) {
				// граничный узел.
				dSqw2 = border_neighbor[iW2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW2]]) {
					dSqw2 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW2];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqw2 = dy_loc * dz_loc;
				}
			}


		}


		if (iN2 > -1) {

			dSqn2 = dx * dz;

			if (bN2) {
				// граничный узел.
				dSqn2 = border_neighbor[iN2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
					dSqn2 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN2];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqn2 = dx_loc * dz_loc;
				}
			}


		}


		if (iS2 > -1) {

			dSqs2 = dx * dz;

			if (bS2) {
				// граничный узел.
				dSqs2 = border_neighbor[iS2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
					dSqs2 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS2];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqs2 = dx_loc * dz_loc;
				}
			}


		}


		if (iT2 > -1) {

			dSqt2 = dx * dy;

			if (bT2) {
				// граничный узел.
				dSqt2 = border_neighbor[iT2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
					dSqt2 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT2];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqt2 = dx_loc * dy_loc;
				}
			}


		}


		if (iB2 > -1) {

			dSqb2 = dx * dy;

			if (bB2) {
				// граничный узел.
				dSqb2 = border_neighbor[iB2 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
					dSqb2 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB2];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqb2 = dx_loc * dy_loc;
				}
			}


		}






		if (iE3 > -1) {

			dSqe3 = dy * dz;

			if (bE3) {
				// граничный узел.
				dSqe3 = border_neighbor[iE3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
					dSqe3 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqe3 = dy_loc * dz_loc;
				}
			}


		}


		if (iW3 > -1) {

			dSqw3 = dy * dz;

			if (bW3) {
				// граничный узел.
				dSqw3 = border_neighbor[iW3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
					dSqw3 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqw3 = dy_loc * dz_loc;
				}
			}


		}


		if (iN3 > -1) {

			dSqn3 = dx * dz;

			if (bN3) {
				// граничный узел.
				dSqn3 = border_neighbor[iN3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
					dSqn3 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqn3 = dx_loc * dz_loc;
				}
			}


		}


		if (iS3 > -1) {

			dSqs3 = dx * dz;

			if (bS3) {
				// граничный узел.
				dSqs3 = border_neighbor[iS3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
					dSqs3 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqs3 = dx_loc * dz_loc;
				}
			}


		}


		if (iT3 > -1) {

			dSqt3 = dx * dy;

			if (bT3) {
				// граничный узел.
				dSqt3 = border_neighbor[iT3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
					dSqt3 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqt3 = dx_loc * dy_loc;
				}
			}


		}


		if (iB3 > -1) {

			dSqb3 = dx * dy;

			if (bB3) {
				// граничный узел.
				dSqb3 = border_neighbor[iB3 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
					dSqb3 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB3];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqb3 = dx_loc * dy_loc;
				}
			}


		}




		if (iE4 > -1) {

			dSqe4 = dy * dz;

			if (bE4) {
				// граничный узел.
				dSqe4 = border_neighbor[iE4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
					dSqe4 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqe4 = dy_loc * dz_loc;
				}
			}


		}


		if (iW4 > -1) {

			dSqw4 = dy * dz;

			if (bW4) {
				// граничный узел.
				dSqw4 = border_neighbor[iW4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
					dSqw4 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqw4 = dy_loc * dz_loc;
				}
			}


		}


		if (iN4 > -1) {

			dSqn4 = dx * dz;

			if (bN4) {
				// граничный узел.
				dSqn4 = border_neighbor[iN4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
					dSqn4 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqn4 = dx_loc * dz_loc;
				}
			}


		}


		if (iS4 > -1) {

			dSqs4 = dx * dz;

			if (bS4) {
				// граничный узел.
				dSqs4 = border_neighbor[iS4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
					dSqs4 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqs4 = dx_loc * dz_loc;
				}
			}


		}


		if (iT4 > -1) {

			dSqt4 = dx * dy;

			if (bT4) {
				// граничный узел.
				dSqt4 = border_neighbor[iT4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
					dSqt4 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqt4 = dx_loc * dy_loc;
				}
			}


		}


		if (iB4 > -1) {

			dSqb4 = dx * dy;

			if (bB4) {
				// граничный узел.
				dSqb4 = border_neighbor[iB4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
					dSqb4 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqb4 = dx_loc * dy_loc;
				}
			}


		}
	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}


	doublereal PAMe = 0.0, PAMw = 0.0, PAMn = 0.0, PAMs = 0.0, PAMt = 0.0, PAMb = 0.0;
	doublereal PAMe2 = 0.0, PAMw2 = 0.0, PAMn2 = 0.0, PAMs2 = 0.0, PAMt2 = 0.0, PAMb2 = 0.0;
	doublereal PAMe3 = 0.0, PAMw3 = 0.0, PAMn3 = 0.0, PAMs3 = 0.0, PAMt3 = 0.0, PAMb3 = 0.0;
	doublereal PAMe4 = 0.0, PAMw4 = 0.0, PAMn4 = 0.0, PAMs4 = 0.0, PAMt4 = 0.0, PAMb4 = 0.0;

    if (!bbond ) {

		if (bLRfree && ((!b_on_adaptive_local_refinement_mesh))) {

			// не работает на АЛИС.
			//if (b_on_adaptive_local_refinement_mesh) {
				//printf("function green_gaussPAM in module greengauss.c bLRfree not worked in ALICE mesh...\n ");
				//getchar();
				//exit(1);
			//}

			// Так гораздо хуже на opening тесте.

			// В случае bLRFree мы будем линейно интерполировать давление в граничные узлы чтобы сохранить значение градиента.
			// градиент давления очень важен для естественно конвективных течний, т.к. они протекают под действием этого градиента
			// и это для них особенно важно вблизи твёрдых стенок.
			
			if (iE > -1) {
				if (!bE) PAMe = feplus * potent[PAM][iE] + (1.0 - feplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iW > -1) {
						PAMe = potent[PAM][iP] + (dxe / dxw)*(potent[PAM][iP] - potent[PAM][iW]);
					}
					else if (iW2 > -1) {
						PAMe = potent[PAM][iP] + (dxe / dxw2)*(potent[PAM][iP] - potent[PAM][iW2]);
					}
					else if (iW3 > -1) {
						PAMe = potent[PAM][iP] + (dxe / dxw3)*(potent[PAM][iP] - potent[PAM][iW3]);
					}
					else if (iW4 > -1) {
						PAMe = potent[PAM][iP] + (dxe / dxw4)*(potent[PAM][iP] - potent[PAM][iW4]);
					}
				}
			}

			if (iE2 > -1) {
				if (!bE2) PAMe2 = feplus2 * potent[PAM][iE2] + (1.0 - feplus2)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iW > -1) {
						PAMe2 = potent[PAM][iP] + (dxe2 / dxw)*(potent[PAM][iP] - potent[PAM][iW]);
					}
					else if (iW2 > -1) {
						PAMe2 = potent[PAM][iP] + (dxe2 / dxw2)*(potent[PAM][iP] - potent[PAM][iW2]);
					}
					else if (iW3 > -1) {
						PAMe2 = potent[PAM][iP] + (dxe2 / dxw3)*(potent[PAM][iP] - potent[PAM][iW3]);
					}
					else if (iW4 > -1) {
						PAMe2 = potent[PAM][iP] + (dxe2 / dxw4)*(potent[PAM][iP] - potent[PAM][iW4]);
					}
				}
			}

			if (iE3 > -1) {
				if (!bE3) PAMe3 = feplus3 * potent[PAM][iE3] + (1.0 - feplus3)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iW > -1) {
						PAMe3 = potent[PAM][iP] + (dxe3 / dxw)*(potent[PAM][iP] - potent[PAM][iW]);
					}
					else if (iW2 > -1) {
						PAMe3 = potent[PAM][iP] + (dxe3 / dxw2)*(potent[PAM][iP] - potent[PAM][iW2]);
					}
					else if (iW3 > -1) {
						PAMe3 = potent[PAM][iP] + (dxe3 / dxw3)*(potent[PAM][iP] - potent[PAM][iW3]);
					}
					else if (iW4 > -1) {
						PAMe3 = potent[PAM][iP] + (dxe3 / dxw4)*(potent[PAM][iP] - potent[PAM][iW4]);
					}
				}
			}

			if (iE4 > -1) {
				if (!bE4) PAMe4 = feplus4 * potent[PAM][iE4] + (1.0 - feplus4)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iW > -1) {
						PAMe4 = potent[PAM][iP] + (dxe4 / dxw)*(potent[PAM][iP] - potent[PAM][iW]);
					}
					else if (iW2 > -1) {
						PAMe4 = potent[PAM][iP] + (dxe4 / dxw2)*(potent[PAM][iP] - potent[PAM][iW2]);
					}
					else if (iW3 > -1) {
						PAMe4 = potent[PAM][iP] + (dxe4 / dxw3)*(potent[PAM][iP] - potent[PAM][iW3]);
					}
					else if (iW4 > -1) {
						PAMe4 = potent[PAM][iP] + (dxe4 / dxw4)*(potent[PAM][iP] - potent[PAM][iW4]);
					}
				}
			}

			if (iW > -1) {
				if (!bW) PAMw = fwplus * potent[PAM][iW] + (1.0 - fwplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iE > -1) {
						PAMw = potent[PAM][iP] + (dxw / dxe)*(potent[PAM][iP] - potent[PAM][iE]);
					} 
					else if (iE2 > -1) {
						PAMw = potent[PAM][iP] + (dxw / dxe2)*(potent[PAM][iP] - potent[PAM][iE2]);
					}
					else if (iE3 > -1) {
						PAMw = potent[PAM][iP] + (dxw / dxe3)*(potent[PAM][iP] - potent[PAM][iE3]);
					}
					else if (iE4 > -1) {
						PAMw = potent[PAM][iP] + (dxw / dxe4)*(potent[PAM][iP] - potent[PAM][iE4]);
					}
					
				}
			}

			if (iW2 > -1) {
				if (!bW2) PAMw2 = fwplus2 * potent[PAM][iW2] + (1.0 - fwplus2)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iE > -1) {
						PAMw2 = potent[PAM][iP] + (dxw2 / dxe)*(potent[PAM][iP] - potent[PAM][iE]);
					}
					else if (iE2 > -1) {
						PAMw2 = potent[PAM][iP] + (dxw2 / dxe2)*(potent[PAM][iP] - potent[PAM][iE2]);
					}
					else if (iE3 > -1) {
						PAMw2 = potent[PAM][iP] + (dxw2 / dxe3)*(potent[PAM][iP] - potent[PAM][iE3]);
					}
					else if (iE4 > -1) {
						PAMw2 = potent[PAM][iP] + (dxw2 / dxe4)*(potent[PAM][iP] - potent[PAM][iE4]);
					}

				}
			}

			if (iW3 > -1) {
				if (!bW3) PAMw3 = fwplus3 * potent[PAM][iW3] + (1.0 - fwplus3)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iE > -1) {
						PAMw3 = potent[PAM][iP] + (dxw3 / dxe)*(potent[PAM][iP] - potent[PAM][iE]);
					}
					else if (iE2 > -1) {
						PAMw3 = potent[PAM][iP] + (dxw3 / dxe2)*(potent[PAM][iP] - potent[PAM][iE2]);
					}
					else if (iE3 > -1) {
						PAMw3 = potent[PAM][iP] + (dxw3 / dxe3)*(potent[PAM][iP] - potent[PAM][iE3]);
					}
					else if (iE4 > -1) {
						PAMw3 = potent[PAM][iP] + (dxw3 / dxe4)*(potent[PAM][iP] - potent[PAM][iE4]);
					}

				}
			}

			if (iW4 > -1) {
				if (!bW4) PAMw4 = fwplus4 * potent[PAM][iW4] + (1.0 - fwplus4)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iE > -1) {
						PAMw4 = potent[PAM][iP] + (dxw4 / dxe)*(potent[PAM][iP] - potent[PAM][iE]);
					}
					else if (iE2 > -1) {
						PAMw4 = potent[PAM][iP] + (dxw4 / dxe2)*(potent[PAM][iP] - potent[PAM][iE2]);
					}
					else if (iE3 > -1) {
						PAMw4 = potent[PAM][iP] + (dxw4 / dxe3)*(potent[PAM][iP] - potent[PAM][iE3]);
					}
					else if (iE4 > -1) {
						PAMw4 = potent[PAM][iP] + (dxw4 / dxe4)*(potent[PAM][iP] - potent[PAM][iE4]);
					}

				}
			}

			if (iN > -1) {
				if (!bN) PAMn = fnplus * potent[PAM][iN] + (1.0 - fnplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iS > -1) {
						PAMn = potent[PAM][iP] + (dyn / dys)*(potent[PAM][iP] - potent[PAM][iS]);
					}
					else if (iS2 > -1) {
						PAMn = potent[PAM][iP] + (dyn / dys2)*(potent[PAM][iP] - potent[PAM][iS2]);
					}
					else if (iS3 > -1) {
						PAMn = potent[PAM][iP] + (dyn / dys3)*(potent[PAM][iP] - potent[PAM][iS3]);
					}
					else if (iS4 > -1) {
						PAMn = potent[PAM][iP] + (dyn / dys4)*(potent[PAM][iP] - potent[PAM][iS4]);
					}
				
				}
			}

			if (iN2 > -1) {
				if (!bN2) PAMn2 = fnplus2 * potent[PAM][iN2] + (1.0 - fnplus2)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iS > -1) {
						PAMn2 = potent[PAM][iP] + (dyn2 / dys)*(potent[PAM][iP] - potent[PAM][iS]);
					}
					else if (iS2 > -1) {
						PAMn2 = potent[PAM][iP] + (dyn2 / dys2)*(potent[PAM][iP] - potent[PAM][iS2]);
					}
					else if (iS3 > -1) {
						PAMn2 = potent[PAM][iP] + (dyn2 / dys3)*(potent[PAM][iP] - potent[PAM][iS3]);
					}
					else if (iS4 > -1) {
						PAMn2 = potent[PAM][iP] + (dyn2 / dys4)*(potent[PAM][iP] - potent[PAM][iS4]);
					}

				}
			}

			if (iN3 > -1) {
				if (!bN3) PAMn3 = fnplus3 * potent[PAM][iN3] + (1.0 - fnplus3)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iS > -1) {
						PAMn3 = potent[PAM][iP] + (dyn3 / dys)*(potent[PAM][iP] - potent[PAM][iS]);
					}
					else if (iS2 > -1) {
						PAMn3 = potent[PAM][iP] + (dyn3 / dys2)*(potent[PAM][iP] - potent[PAM][iS2]);
					}
					else if (iS3 > -1) {
						PAMn3 = potent[PAM][iP] + (dyn3 / dys3)*(potent[PAM][iP] - potent[PAM][iS3]);
					}
					else if (iS4 > -1) {
						PAMn3 = potent[PAM][iP] + (dyn3 / dys4)*(potent[PAM][iP] - potent[PAM][iS4]);
					}

				}
			}

			if (iN4 > -1) {
				if (!bN4) PAMn4 = fnplus4 * potent[PAM][iN4] + (1.0 - fnplus4)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iS > -1) {
						PAMn4 = potent[PAM][iP] + (dyn4 / dys)*(potent[PAM][iP] - potent[PAM][iS]);
					}
					else if (iS2 > -1) {
						PAMn4 = potent[PAM][iP] + (dyn4 / dys2)*(potent[PAM][iP] - potent[PAM][iS2]);
					}
					else if (iS3 > -1) {
						PAMn4 = potent[PAM][iP] + (dyn4 / dys3)*(potent[PAM][iP] - potent[PAM][iS3]);
					}
					else if (iS4 > -1) {
						PAMn4 = potent[PAM][iP] + (dyn4 / dys4)*(potent[PAM][iP] - potent[PAM][iS4]);
					}

				}
			}

			if (iS > -1) {
				if (!bS) PAMs = fsplus * potent[PAM][iS] + (1.0 - fsplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iN > -1) {
						PAMs = potent[PAM][iP] + (dys / dyn)*(potent[PAM][iP] - potent[PAM][iN]);
					}
					else if (iN2 > -1) {
						PAMs = potent[PAM][iP] + (dys / dyn2)*(potent[PAM][iP] - potent[PAM][iN2]);
					}
					else if (iN3 > -1) {
						PAMs = potent[PAM][iP] + (dys / dyn3)*(potent[PAM][iP] - potent[PAM][iN3]);
					}
					else if (iN4 > -1) {
						PAMs = potent[PAM][iP] + (dys / dyn4)*(potent[PAM][iP] - potent[PAM][iN4]);
					}
					
				}
			}

			if (iS2 > -1) {
				if (!bS2) PAMs2 = fsplus2 * potent[PAM][iS2] + (1.0 - fsplus2)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iN > -1) {
						PAMs2 = potent[PAM][iP] + (dys2 / dyn)*(potent[PAM][iP] - potent[PAM][iN]);
					}
					else if (iN2 > -1) {
						PAMs2 = potent[PAM][iP] + (dys2 / dyn2)*(potent[PAM][iP] - potent[PAM][iN2]);
					}
					else if (iN3 > -1) {
						PAMs2 = potent[PAM][iP] + (dys2 / dyn3)*(potent[PAM][iP] - potent[PAM][iN3]);
					}
					else if (iN4 > -1) {
						PAMs2 = potent[PAM][iP] + (dys2 / dyn4)*(potent[PAM][iP] - potent[PAM][iN4]);
					}

				}
			}

			if (iS3 > -1) {
				if (!bS3) PAMs3 = fsplus3 * potent[PAM][iS3] + (1.0 - fsplus3)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iN > -1) {
						PAMs3 = potent[PAM][iP] + (dys3 / dyn)*(potent[PAM][iP] - potent[PAM][iN]);
					}
					else if (iN2 > -1) {
						PAMs3 = potent[PAM][iP] + (dys3 / dyn2)*(potent[PAM][iP] - potent[PAM][iN2]);
					}
					else if (iN3 > -1) {
						PAMs3 = potent[PAM][iP] + (dys3 / dyn3)*(potent[PAM][iP] - potent[PAM][iN3]);
					}
					else if (iN4 > -1) {
						PAMs3 = potent[PAM][iP] + (dys3 / dyn4)*(potent[PAM][iP] - potent[PAM][iN4]);
					}

				}
			}

			if (iS4 > -1) {
				if (!bS4) PAMs4 = fsplus4 * potent[PAM][iS4] + (1.0 - fsplus4)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iN > -1) {
						PAMs4 = potent[PAM][iP] + (dys4 / dyn)*(potent[PAM][iP] - potent[PAM][iN]);
					}
					else if (iN2 > -1) {
						PAMs4 = potent[PAM][iP] + (dys4 / dyn2)*(potent[PAM][iP] - potent[PAM][iN2]);
					}
					else if (iN3 > -1) {
						PAMs4 = potent[PAM][iP] + (dys4 / dyn3)*(potent[PAM][iP] - potent[PAM][iN3]);
					}
					else if (iN4 > -1) {
						PAMs4 = potent[PAM][iP] + (dys4 / dyn4)*(potent[PAM][iP] - potent[PAM][iN4]);
					}

				}
			}

			if (iT > -1) {
				if (!bT) PAMt = ftplus * potent[PAM][iT] + (1.0 - ftplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iB > -1) {
						PAMt = potent[PAM][iP] + (dzt / dzb)*(potent[PAM][iP] - potent[PAM][iB]);
					} 
					else if (iB2 > -1) {
						PAMt = potent[PAM][iP] + (dzt / dzb2)*(potent[PAM][iP] - potent[PAM][iB2]);
					}
					else if (iB3 > -1) {
						PAMt = potent[PAM][iP] + (dzt / dzb3)*(potent[PAM][iP] - potent[PAM][iB3]);
					}
					else if (iB4 > -1) {
						PAMt = potent[PAM][iP] + (dzt / dzb4)*(potent[PAM][iP] - potent[PAM][iB4]);
					}
					
				}
			}


			if (iT2 > -1) {
				if (!bT2) PAMt2 = ftplus2 * potent[PAM][iT2] + (1.0 - ftplus2)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iB > -1) {
						PAMt2 = potent[PAM][iP] + (dzt2 / dzb)*(potent[PAM][iP] - potent[PAM][iB]);
					}
					else if (iB2 > -1) {
						PAMt2 = potent[PAM][iP] + (dzt2 / dzb2)*(potent[PAM][iP] - potent[PAM][iB2]);
					}
					else if (iB3 > -1) {
						PAMt2 = potent[PAM][iP] + (dzt2 / dzb3)*(potent[PAM][iP] - potent[PAM][iB3]);
					}
					else if (iB4 > -1) {
						PAMt2 = potent[PAM][iP] + (dzt2 / dzb4)*(potent[PAM][iP] - potent[PAM][iB4]);
					}

				}
			}

			if (iT3 > -1) {
				if (!bT3) PAMt3 = ftplus3 * potent[PAM][iT3] + (1.0 - ftplus3)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iB > -1) {
						PAMt3 = potent[PAM][iP] + (dzt3 / dzb)*(potent[PAM][iP] - potent[PAM][iB]);
					}
					else if (iB2 > -1) {
						PAMt3 = potent[PAM][iP] + (dzt3 / dzb2)*(potent[PAM][iP] - potent[PAM][iB2]);
					}
					else if (iB3 > -1) {
						PAMt3 = potent[PAM][iP] + (dzt3 / dzb3)*(potent[PAM][iP] - potent[PAM][iB3]);
					}
					else if (iB4 > -1) {
						PAMt3 = potent[PAM][iP] + (dzt3 / dzb4)*(potent[PAM][iP] - potent[PAM][iB4]);
					}

				}
			}

			if (iT4 > -1) {
				if (!bT4) PAMt4 = ftplus4 * potent[PAM][iT4] + (1.0 - ftplus4)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iB > -1) {
						PAMt4 = potent[PAM][iP] + (dzt4 / dzb)*(potent[PAM][iP] - potent[PAM][iB]);
					}
					else if (iB2 > -1) {
						PAMt4 = potent[PAM][iP] + (dzt4 / dzb2)*(potent[PAM][iP] - potent[PAM][iB2]);
					}
					else if (iB3 > -1) {
						PAMt4 = potent[PAM][iP] + (dzt4 / dzb3)*(potent[PAM][iP] - potent[PAM][iB3]);
					}
					else if (iB4 > -1) {
						PAMt4 = potent[PAM][iP] + (dzt4 / dzb4)*(potent[PAM][iP] - potent[PAM][iB4]);
					}

				}
			}

			if (iB > -1) {
				if (!bB) PAMb = fbplus * potent[PAM][iB] + (1.0 - fbplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iT > -1) {
						PAMb = potent[PAM][iP] + (dzb / dzt)*(potent[PAM][iP] - potent[PAM][iT]);
					}
					else if (iT2 > -1) {
						PAMb = potent[PAM][iP] + (dzb / dzt2)*(potent[PAM][iP] - potent[PAM][iT2]);
					}
					else  if (iT3 > -1) {
						PAMb = potent[PAM][iP] + (dzb / dzt3)*(potent[PAM][iP] - potent[PAM][iT3]);
					}
					else if (iT4 > -1) {
						PAMb = potent[PAM][iP] + (dzb / dzt4)*(potent[PAM][iP] - potent[PAM][iT4]);
					}					
				}
			}

			if (iB2 > -1) {
				if (!bB2) PAMb2 = fbplus2 * potent[PAM][iB2] + (1.0 - fbplus2)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iT > -1) {
						PAMb2 = potent[PAM][iP] + (dzb2 / dzt)*(potent[PAM][iP] - potent[PAM][iT]);
					}
					else if (iT2 > -1) {
						PAMb2 = potent[PAM][iP] + (dzb2 / dzt2)*(potent[PAM][iP] - potent[PAM][iT2]);
					}
					else  if (iT3 > -1) {
						PAMb2 = potent[PAM][iP] + (dzb2 / dzt3)*(potent[PAM][iP] - potent[PAM][iT3]);
					}
					else if (iT4 > -1) {
						PAMb2 = potent[PAM][iP] + (dzb2 / dzt4)*(potent[PAM][iP] - potent[PAM][iT4]);
					}
				}
			}

			if (iB3 > -1) {
				if (!bB3) PAMb3 = fbplus3 * potent[PAM][iB3] + (1.0 - fbplus3)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iT > -1) {
						PAMb3 = potent[PAM][iP] + (dzb3 / dzt)*(potent[PAM][iP] - potent[PAM][iT]);
					}
					else if (iT2 > -1) {
						PAMb3 = potent[PAM][iP] + (dzb3 / dzt2)*(potent[PAM][iP] - potent[PAM][iT2]);
					}
					else  if (iT3 > -1) {
						PAMb3 = potent[PAM][iP] + (dzb3 / dzt3)*(potent[PAM][iP] - potent[PAM][iT3]);
					}
					else if (iT4 > -1) {
						PAMb3 = potent[PAM][iP] + (dzb3 / dzt4)*(potent[PAM][iP] - potent[PAM][iT4]);
					}
				}
			}

			if (iB4 > -1) {
				if (!bB4) PAMb4 = fbplus4 * potent[PAM][iB4] + (1.0 - fbplus4)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					if (iT > -1) {
						PAMb4 = potent[PAM][iP] + (dzb4 / dzt)*(potent[PAM][iP] - potent[PAM][iT]);
					}
					else if (iT2 > -1) {
						PAMb4 = potent[PAM][iP] + (dzb4 / dzt2)*(potent[PAM][iP] - potent[PAM][iT2]);
					}
					else  if (iT3 > -1) {
						PAMb4 = potent[PAM][iP] + (dzb4 / dzt3)*(potent[PAM][iP] - potent[PAM][iT3]);
					}
					else if (iT4 > -1) {
						PAMb4 = potent[PAM][iP] + (dzb4 / dzt4)*(potent[PAM][iP] - potent[PAM][iT4]);
					}
				}
			}


			
             // градиент Давления.
	         //potent[GRADXPAM][iP]=(PAMe-PAMw)/dx;
	         //potent[GRADYPAM][iP]=(PAMn-PAMs)/dy;
	         //potent[GRADZPAM][iP]=(PAMt-PAMb)/dz;
			 potent[GRADXPAM][iP] = (PAMe*dSqe / (dy*dz) + PAMe2 * dSqe2 / (dy*dz) + PAMe3 * dSqe3 / (dy*dz) + PAMe4 * dSqe4 / (dy*dz) - (PAMw*dSqw / (dy*dz) + PAMw2 * dSqw2 / (dy*dz) + PAMw3 * dSqw3 / (dy*dz) + PAMw4 * dSqw4 / (dy*dz))) / dx;
			 potent[GRADYPAM][iP] = (PAMn*dSqn / (dx*dz) + PAMn2 * dSqn2 / (dx*dz) + PAMn3 * dSqn3 / (dx*dz) + PAMn4 * dSqn4 / (dx*dz) - (PAMs*dSqs / (dx*dz) + PAMs2 * dSqs2 / (dx*dz) + PAMs3 * dSqs3 / (dx*dz) + PAMs4 * dSqs4 / (dx*dz))) / dy;
			 potent[GRADZPAM][iP] = (PAMt*dSqt / (dx*dy) + PAMt2 * dSqt2 / (dx*dy) + PAMt3 * dSqt3 / (dx*dy) + PAMt4 * dSqt4 / (dx*dy) - (PAMb*dSqb / (dx*dy) + PAMb2 * dSqb2 / (dx*dy) + PAMb3 * dSqb3 / (dx*dy) + PAMb4 * dSqb4 / (dx*dy))) / dz;
			 if (potent[GRADXPAM][iP] != potent[GRADXPAM][iP]) {
				 printf("Error: potent[GRADXPAM][%lld]=%e is NAN or INF\n",iP, potent[GRADXPAM][iP]);
				 printf("PAMe=%e PAMe2=%e PAMe3=%e PAMe4=%e\n", PAMe, PAMe2, PAMe3, PAMe4);
				 printf("PAMw=%e PAMw2=%e PAMw3=%e PAMw4=%e\n", PAMw, PAMw2, PAMw3, PAMw4);
				 printf("dx=%e dy=%e dz=%e\n",dx,dy,dz);
				 printf("dSqe=%e dSqe2=%e dSqe3=%e dSqe4=%e\n", dSqe, dSqe2, dSqe3, dSqe4);
				 printf("dSqw=%e dSqw2=%e dSqw3=%e dSqw4=%e\n", dSqw, dSqw2, dSqw3, dSqw4);
				 system("pause");
			 }
			 if (potent[GRADYPAM][iP] != potent[GRADYPAM][iP]) {
				 printf("Error: potent[GRADYPAM][%lld]=%e is NAN or INF\n", iP, potent[GRADYPAM][iP]);
				 printf("PAMn=%e PAMn2=%e PAMn3=%e PAMn4=%e\n", PAMn, PAMn2, PAMn3, PAMn4);
				 printf("PAMs=%e PAMs2=%e PAMs3=%e PAMs4=%e\n", PAMs, PAMs2, PAMs3, PAMs4);
				 printf("dx=%e dy=%e dz=%e\n", dx, dy, dz);
				 printf("dSqn=%e dSqn2=%e dSqn3=%e dSqn4=%e\n", dSqn, dSqn2, dSqn3, dSqn4);
				 printf("dSqs=%e dSqs2=%e dSqs3=%e dSqs4=%e\n", dSqs, dSqs2, dSqs3, dSqs4);
				 system("pause");
			 }
			 if (potent[GRADZPAM][iP] != potent[GRADZPAM][iP]) {
				 printf("Error: potent[GRADZPAM][%lld]=%e is NAN or INF\n", iP, potent[GRADZPAM][iP]);
				 printf("PAMt=%e PAMt2=%e PAMt3=%e PAMt4=%e\n", PAMt, PAMt2, PAMt3, PAMt4);
				 printf("PAMb=%e PAMb2=%e PAMb3=%e PAMb4=%e\n", PAMb, PAMb2, PAMb3, PAMb4);
				 printf("dx=%e dy=%e dz=%e\n", dx, dy, dz);
				 printf("dSqt=%e dSqt2=%e dSqt3=%e dSqt4=%e\n", dSqt, dSqt2, dSqt3, dSqt4);
				 printf("dSqb=%e dSqb2=%e dSqb3=%e dSqb4=%e\n", dSqb, dSqb2, dSqb3, dSqb4);
				 system("pause");
			 }

		}
		else {


		    // внутренние КО.

		    // Линейно интерполируем поправку давления на грань контрольного объёма,
		    // а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 
			if (iE > -1) {
				if (!bE) PAMe = feplus * potent[PAM][iE] + (1.0 - feplus)*potent[PAM][iP]; else PAMe = potent[PAM][iE];
			}
			if (iW > -1) {
				if (!bW) PAMw = fwplus * potent[PAM][iW] + (1.0 - fwplus)*potent[PAM][iP]; else PAMw = potent[PAM][iW];
			}
			if (iN > -1) {
				if (!bN) PAMn = fnplus * potent[PAM][iN] + (1.0 - fnplus)*potent[PAM][iP]; else PAMn = potent[PAM][iN];
			}
			if (iS > -1) {
				if (!bS) PAMs = fsplus * potent[PAM][iS] + (1.0 - fsplus)*potent[PAM][iP]; else PAMs = potent[PAM][iS];
			}
			if (iT > -1) {
				if (!bT) PAMt = ftplus * potent[PAM][iT] + (1.0 - ftplus)*potent[PAM][iP]; else PAMt = potent[PAM][iT];
			}
			if (iB > -1) {
				if (!bB) PAMb = fbplus * potent[PAM][iB] + (1.0 - fbplus)*potent[PAM][iP]; else PAMb = potent[PAM][iB];
			}

			if (iE2 > -1) {
				if (!bE2) PAMe2 = feplus2 * potent[PAM][iE2] + (1.0 - feplus2)*potent[PAM][iP]; else PAMe2 = potent[PAM][iE2];
			}
			if (iW2 > -1) {
				if (!bW2) PAMw2 = fwplus2 * potent[PAM][iW2] + (1.0 - fwplus2)*potent[PAM][iP]; else PAMw2 = potent[PAM][iW2];
			}
			if (iN2 > -1) {
				if (!bN2) PAMn2 = fnplus2 * potent[PAM][iN2] + (1.0 - fnplus2)*potent[PAM][iP]; else PAMn2 = potent[PAM][iN2];
			}
			if (iS2 > -1) {
				if (!bS2) PAMs2 = fsplus2 * potent[PAM][iS2] + (1.0 - fsplus2)*potent[PAM][iP]; else PAMs2 = potent[PAM][iS2];
			}
			if (iT2 > -1) {
				if (!bT2) PAMt2 = ftplus2 * potent[PAM][iT2] + (1.0 - ftplus2)*potent[PAM][iP]; else PAMt2 = potent[PAM][iT2];
			}
			if (iB2 > -1) {
				if (!bB2) PAMb2 = fbplus2 * potent[PAM][iB2] + (1.0 - fbplus2)*potent[PAM][iP]; else PAMb2 = potent[PAM][iB2];
			}

			if (iE3 > -1) {
				if (!bE3) PAMe3 = feplus3 * potent[PAM][iE3] + (1.0 - feplus3)*potent[PAM][iP]; else PAMe3 = potent[PAM][iE3];
			}
			if (iW3 > -1) {
				if (!bW3) PAMw3 = fwplus3 * potent[PAM][iW3] + (1.0 - fwplus3)*potent[PAM][iP]; else PAMw3 = potent[PAM][iW3];
			}
			if (iN3 > -1) {
				if (!bN3) PAMn3 = fnplus3 * potent[PAM][iN3] + (1.0 - fnplus3)*potent[PAM][iP]; else PAMn3 = potent[PAM][iN3];
			}
			if (iS3 > -1) {
				if (!bS3) PAMs3 = fsplus3 * potent[PAM][iS3] + (1.0 - fsplus3)*potent[PAM][iP]; else PAMs3 = potent[PAM][iS3];
			}
			if (iT3 > -1) {
				if (!bT3) PAMt3 = ftplus3 * potent[PAM][iT3] + (1.0 - ftplus3)*potent[PAM][iP]; else PAMt3 = potent[PAM][iT3];
			}
			if (iB3 > -1) {
				if (!bB3) PAMb3 = fbplus3 * potent[PAM][iB3] + (1.0 - fbplus3)*potent[PAM][iP]; else PAMb3 = potent[PAM][iB3];
			}

			if (iE4 > -1) {
				if (!bE4) PAMe4 = feplus4 * potent[PAM][iE4] + (1.0 - feplus4)*potent[PAM][iP]; else PAMe4 = potent[PAM][iE4];
			}
			if (iW4 > -1) {
				if (!bW4) PAMw4 = fwplus4 * potent[PAM][iW4] + (1.0 - fwplus4)*potent[PAM][iP]; else PAMw4 = potent[PAM][iW4];
			}
			if (iN4 > -1) {
				if (!bN4) PAMn4 = fnplus4 * potent[PAM][iN4] + (1.0 - fnplus4)*potent[PAM][iP]; else PAMn4 = potent[PAM][iN4];
			}
			if (iS4 > -1) {
				if (!bS4) PAMs4 = fsplus4 * potent[PAM][iS4] + (1.0 - fsplus4)*potent[PAM][iP]; else PAMs4 = potent[PAM][iS4];
			}
			if (iT4 > -1) {
				if (!bT4) PAMt4 = ftplus4 * potent[PAM][iT4] + (1.0 - ftplus4)*potent[PAM][iP]; else PAMt4 = potent[PAM][iT4];
			}
			if (iB4 > -1) {
				if (!bB4) PAMb4 = fbplus4 * potent[PAM][iB4] + (1.0 - fbplus4)*potent[PAM][iP]; else PAMb4 = potent[PAM][iB4];
			}

			// градиент Поправки давления.
	        //potent[GRADXPAM][iP]=(PAMe-PAMw)/dx;
	        //potent[GRADYPAM][iP]=(PAMn-PAMs)/dy;
	        //potent[GRADZPAM][iP]=(PAMt-PAMb)/dz;
			potent[GRADXPAM][iP] = (PAMe*dSqe/(dy*dz)+ PAMe2 * dSqe2 / (dy*dz)+ PAMe3 * dSqe3 / (dy*dz)+ PAMe4 * dSqe4 / (dy*dz)-(PAMw*dSqw / (dy*dz) + PAMw2 * dSqw2 / (dy*dz) + PAMw3 * dSqw3 / (dy*dz) + PAMw4 * dSqw4 / (dy*dz))) / dx;
			potent[GRADYPAM][iP] = (PAMn*dSqn/(dx*dz)+ PAMn2 * dSqn2 / (dx*dz)+ PAMn3 * dSqn3 / (dx*dz)+ PAMn4 * dSqn4 / (dx*dz) -( PAMs*dSqs/(dx*dz)+ PAMs2 * dSqs2 / (dx*dz)+ PAMs3 * dSqs3 / (dx*dz) + PAMs4 * dSqs4 / (dx*dz))) / dy;
			potent[GRADZPAM][iP] = (PAMt*dSqt/(dx*dy)+ PAMt2 * dSqt2 / (dx*dy)+ PAMt3 * dSqt3 / (dx*dy)+ PAMt4 * dSqt4 / (dx*dy) - (PAMb*dSqb/(dx*dy)+ PAMb2 * dSqb2 / (dx*dy)+ PAMb3 * dSqb3 / (dx*dy)+ PAMb4 * dSqb4 / (dx*dy))) / dz;

		}

	}
   else {
	    const integer interpol=0; // 0 или 1 при линейной интерполяции.

#if (interpol==0)
		{

			if (1) {

				// Используется по умолчанию и даёт отличные результаты на большом наборе расчётных задач.
				// лучший выбор 6.05.2017.

				if (iE > -1) {
					if (bE) {
						potent[GRADXPAM][iE] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iE] = potent[GRADYPAM][iP];
						potent[GRADZPAM][iE] = potent[GRADZPAM][iP];
					}
				}

				if (iW > -1) {
					if (bW) {
						potent[GRADXPAM][iW] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iW] = potent[GRADYPAM][iP];
						potent[GRADZPAM][iW] = potent[GRADZPAM][iP];
					}
				}

				if (iN > -1) {
					if (bN) {
						potent[GRADYPAM][iN] = potent[GRADYPAM][iP];
						potent[GRADXPAM][iN] = potent[GRADXPAM][iP];
						potent[GRADZPAM][iN] = potent[GRADZPAM][iP];
					}
				}

				if (iS > -1) {
					if (bS) {
						// до 10.02.2017 был косяк
						potent[GRADYPAM][iS] = potent[GRADYPAM][iP];
						potent[GRADXPAM][iS] = potent[GRADXPAM][iP];
						potent[GRADZPAM][iS] = potent[GRADZPAM][iP];
					}
				}

				if (iT > -1) {
					if (bT) {
						potent[GRADZPAM][iT] = potent[GRADZPAM][iP];
						potent[GRADXPAM][iT] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iT] = potent[GRADYPAM][iP];
					}
				}

				if (iB > -1) {
					if (bB) {
						potent[GRADZPAM][iB] = potent[GRADZPAM][iP];
						potent[GRADXPAM][iB] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iB] = potent[GRADYPAM][iP];
					}
				}

				if (b_on_adaptive_local_refinement_mesh) {

					if (iE2 > -1) {
						if (bE2) {
							potent[GRADXPAM][iE2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iE2] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE2] = potent[GRADZPAM][iP];
						}
					}

					if (iW2 > -1) {
						if (bW2) {
							potent[GRADXPAM][iW2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iW2] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW2] = potent[GRADZPAM][iP];
						}
					}

					if (iN2 > -1) {
						if (bN2) {
							potent[GRADYPAM][iN2] = potent[GRADYPAM][iP];
							potent[GRADXPAM][iN2] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN2] = potent[GRADZPAM][iP];
						}
					}

					if (iS2 > -1) {
						if (bS2) {
							// до 10.02.2017 был косяк
							potent[GRADYPAM][iS2] = potent[GRADYPAM][iP];
							potent[GRADXPAM][iS2] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS2] = potent[GRADZPAM][iP];
						}
					}

					if (iT2 > -1) {
						if (bT2) {
							potent[GRADZPAM][iT2] = potent[GRADZPAM][iP];
							potent[GRADXPAM][iT2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT2] = potent[GRADYPAM][iP];
						}
					}

					if (iB2 > -1) {
						if (bB2) {
							potent[GRADZPAM][iB2] = potent[GRADZPAM][iP];
							potent[GRADXPAM][iB2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB2] = potent[GRADYPAM][iP];
						}
					}

					if (iE3 > -1) {
						if (bE3) {
							potent[GRADXPAM][iE3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iE3] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE3] = potent[GRADZPAM][iP];
						}
					}

					if (iW3 > -1) {
						if (bW3) {
							potent[GRADXPAM][iW3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iW3] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW3] = potent[GRADZPAM][iP];
						}
					}

					if (iN3 > -1) {
						if (bN3) {
							potent[GRADYPAM][iN3] = potent[GRADYPAM][iP];
							potent[GRADXPAM][iN3] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN3] = potent[GRADZPAM][iP];
						}
					}

					if (iS3 > -1) {
						if (bS3) {
							// до 10.02.2017 был косяк
							potent[GRADYPAM][iS3] = potent[GRADYPAM][iP];
							potent[GRADXPAM][iS3] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS3] = potent[GRADZPAM][iP];
						}
					}

					if (iT3 > -1) {
						if (bT3) {
							potent[GRADZPAM][iT3] = potent[GRADZPAM][iP];
							potent[GRADXPAM][iT3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT3] = potent[GRADYPAM][iP];
						}
					}

					if (iB3 > -1) {
						if (bB3) {
							potent[GRADZPAM][iB3] = potent[GRADZPAM][iP];
							potent[GRADXPAM][iB3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB3] = potent[GRADYPAM][iP];
						}
					}

					if (iE4 > -1) {
						if (bE4) {
							potent[GRADXPAM][iE4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iE4] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE4] = potent[GRADZPAM][iP];
						}
					}

					if (iW4 > -1) {
						if (bW4) {
							potent[GRADXPAM][iW4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iW4] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW4] = potent[GRADZPAM][iP];
						}
					}

					if (iN4 > -1) {
						if (bN4) {
							potent[GRADYPAM][iN4] = potent[GRADYPAM][iP];
							potent[GRADXPAM][iN4] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN4] = potent[GRADZPAM][iP];
						}
					}

					if (iS4 > -1) {
						if (bS4) {
							// до 10.02.2017 был косяк
							potent[GRADYPAM][iS4] = potent[GRADYPAM][iP];
							potent[GRADXPAM][iS4] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS4] = potent[GRADZPAM][iP];
						}
					}

					if (iT4 > -1) {
						if (bT4) {
							potent[GRADZPAM][iT4] = potent[GRADZPAM][iP];
							potent[GRADXPAM][iT4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT4] = potent[GRADYPAM][iP];
						}
					}

					if (iB4 > -1) {
						if (bB4) {
							potent[GRADZPAM][iB4] = potent[GRADZPAM][iP];
							potent[GRADXPAM][iB4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB4] = potent[GRADYPAM][iP];
						}
					}
				}

			}
			else if (0) {
				// Это хуже!!!.

				// Рассмотрим границу расчётной области.
				// Вблизи неё есть нормальные и тангенсальные ориентации.
				// Градиент давления вдоль тангенсальных (касательных) направлений к границе может присутствовать это бесспорно.
				// Поэтому на границе он берётся из ближайшего внутреннего узла.
				// В нормальном к границе направлении если на границе задана нормальная компонента скорости, то
				// в нормальном к границе направлении на самой границе градиент давления равен нулю в случае opening границы
				// расчётной области.
				// Это экспериментальная модификация алгоритма расчёта гидродинамических характеристик 6.05.2017.				
				
				if (iE > -1) {
					if (bE) {
						integer inumber = iE - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADXPAM][iE] = 0.0;
						}
						else {
							potent[GRADXPAM][iE] = potent[GRADXPAM][iP];
						}

						// Тангенциальные компоненты.
						potent[GRADYPAM][iE] = potent[GRADYPAM][iP];
						potent[GRADZPAM][iE] = potent[GRADZPAM][iP];
					}
				}

				if (iW > -1) {
					if (bW) {
						integer inumber = iW - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADXPAM][iW] = 0.0;
						}
						else {
							potent[GRADXPAM][iW] = potent[GRADXPAM][iP];
						}

						// Тангенциальные компоненты.					
						potent[GRADYPAM][iW] = potent[GRADYPAM][iP];
						potent[GRADZPAM][iW] = potent[GRADZPAM][iP];
					}
				}

				if (iN > -1) {
					if (bN) {
						integer inumber = iN - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADYPAM][iN] = 0.0;
						}
						else {
							potent[GRADYPAM][iN] = potent[GRADYPAM][iP];
						}

						// Тангенциальные компоненты.					
						potent[GRADXPAM][iN] = potent[GRADXPAM][iP];
						potent[GRADZPAM][iN] = potent[GRADZPAM][iP];
					}
				}

				if (iS > -1) {
					if (bS) {
						integer inumber = iS - maxelm;

						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADYPAM][iS] = 0.0;
						}
						else {
							potent[GRADYPAM][iS] = potent[GRADYPAM][iP];
						}

						// Тангенциальные компоненты.
						// до 10.02.2017 был косяк
						potent[GRADXPAM][iS] = potent[GRADXPAM][iP];
						potent[GRADZPAM][iS] = potent[GRADZPAM][iP];
					}
				}

				if (iT > -1) {
					if (bT) {
						integer inumber = iT - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADZPAM][iT] = 0.0;
						}
						else {
							potent[GRADZPAM][iT] = potent[GRADZPAM][iP];
						}

						// Тангенциальные компоненты.						
						potent[GRADXPAM][iT] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iT] = potent[GRADYPAM][iP];
					}
				}

				if (iB > -1) {
					if (bB) {
						integer inumber = iB - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADZPAM][iB] = 0.0;
						}
						else {
							potent[GRADZPAM][iB] = potent[GRADZPAM][iP];
						}

						// Тангенциальные компоненты.					
						potent[GRADXPAM][iB] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iB] = potent[GRADYPAM][iP];
					}
				}

				if (b_on_adaptive_local_refinement_mesh) {

					if (iE2 > -1) {
						if (bE2) {
							integer inumber = iE2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iE2] = 0.0;
							}
							else {
								potent[GRADXPAM][iE2] = potent[GRADXPAM][iP];
							}

							// Тангенциальные компоненты.
							potent[GRADYPAM][iE2] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE2] = potent[GRADZPAM][iP];
						}
					}

					if (iW2 > -1) {
						if (bW2) {
							integer inumber = iW2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iW2] = 0.0;
							}
							else {
								potent[GRADXPAM][iW2] = potent[GRADXPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADYPAM][iW2] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW2] = potent[GRADZPAM][iP];
						}
					}

					if (iN2 > -1) {
						if (bN2) {
							integer inumber = iN2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iN2] = 0.0;
							}
							else {
								potent[GRADYPAM][iN2] = potent[GRADYPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADXPAM][iN2] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN2] = potent[GRADZPAM][iP];
						}
					}

					if (iS2 > -1) {
						if (bS2) {
							integer inumber = iS2 - maxelm;

							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iS2] = 0.0;
							}
							else {
								potent[GRADYPAM][iS2] = potent[GRADYPAM][iP];
							}

							// Тангенциальные компоненты.
							// до 10.02.2017 был косяк
							potent[GRADXPAM][iS2] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS2] = potent[GRADZPAM][iP];
						}
					}

					if (iT2 > -1) {
						if (bT2) {
							integer inumber = iT2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iT2] = 0.0;
							}
							else {
								potent[GRADZPAM][iT2] = potent[GRADZPAM][iP];
							}

							// Тангенциальные компоненты.						
							potent[GRADXPAM][iT2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT2] = potent[GRADYPAM][iP];
						}
					}

					if (iB2 > -1) {
						if (bB2) {
							integer inumber = iB2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iB2] = 0.0;
							}
							else {
								potent[GRADZPAM][iB2] = potent[GRADZPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADXPAM][iB2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB2] = potent[GRADYPAM][iP];
						}
					}

					if (iE3 > -1) {
						if (bE3) {
							integer inumber = iE3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iE3] = 0.0;
							}
							else {
								potent[GRADXPAM][iE3] = potent[GRADXPAM][iP];
							}

							// Тангенциальные компоненты.
							potent[GRADYPAM][iE3] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE3] = potent[GRADZPAM][iP];
						}
					}

					if (iW3 > -1) {
						if (bW3) {
							integer inumber = iW3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iW3] = 0.0;
							}
							else {
								potent[GRADXPAM][iW3] = potent[GRADXPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADYPAM][iW3] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW3] = potent[GRADZPAM][iP];
						}
					}

					if (iN3 > -1) {
						if (bN3) {
							integer inumber = iN3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iN3] = 0.0;
							}
							else {
								potent[GRADYPAM][iN3] = potent[GRADYPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADXPAM][iN3] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN3] = potent[GRADZPAM][iP];
						}
					}

					if (iS3 > -1) {
						if (bS3) {
							integer inumber = iS3 - maxelm;

							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iS3] = 0.0;
							}
							else {
								potent[GRADYPAM][iS3] = potent[GRADYPAM][iP];
							}

							// Тангенциальные компоненты.
							// до 10.02.2017 был косяк
							potent[GRADXPAM][iS3] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS3] = potent[GRADZPAM][iP];
						}
					}

					if (iT3 > -1) {
						if (bT3) {
							integer inumber = iT3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iT3] = 0.0;
							}
							else {
								potent[GRADZPAM][iT3] = potent[GRADZPAM][iP];
							}

							// Тангенциальные компоненты.						
							potent[GRADXPAM][iT3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT3] = potent[GRADYPAM][iP];
						}
					}

					if (iB3 > -1) {
						if (bB3) {
							integer inumber = iB3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iB3] = 0.0;
							}
							else {
								potent[GRADZPAM][iB3] = potent[GRADZPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADXPAM][iB3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB3] = potent[GRADYPAM][iP];
						}
					}

					if (iE4 > -1) {
						if (bE4) {
							integer inumber = iE4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iE4] = 0.0;
							}
							else {
								potent[GRADXPAM][iE4] = potent[GRADXPAM][iP];
							}

							// Тангенциальные компоненты.
							potent[GRADYPAM][iE4] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE4] = potent[GRADZPAM][iP];
						}
					}

					if (iW4 > -1) {
						if (bW4) {
							integer inumber = iW4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iW4] = 0.0;
							}
							else {
								potent[GRADXPAM][iW4] = potent[GRADXPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADYPAM][iW4] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW4] = potent[GRADZPAM][iP];
						}
					}

					if (iN4 > -1) {
						if (bN4) {
							integer inumber = iN4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iN4] = 0.0;
							}
							else {
								potent[GRADYPAM][iN4] = potent[GRADYPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADXPAM][iN4] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN4] = potent[GRADZPAM][iP];
						}
					}

					if (iS4 > -1) {
						if (bS4) {
							integer inumber = iS4 - maxelm;

							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iS4] = 0.0;
							}
							else {
								potent[GRADYPAM][iS4] = potent[GRADYPAM][iP];
							}

							// Тангенциальные компоненты.
							// до 10.02.2017 был косяк
							potent[GRADXPAM][iS4] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS4] = potent[GRADZPAM][iP];
						}
					}

					if (iT4 > -1) {
						if (bT4) {
							integer inumber = iT4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iT4] = 0.0;
							}
							else {
								potent[GRADZPAM][iT4] = potent[GRADZPAM][iP];
							}

							// Тангенциальные компоненты.						
							potent[GRADXPAM][iT4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT4] = potent[GRADYPAM][iP];
						}
					}

					if (iB4 > -1) {
						if (bB4) {
							integer inumber = iB4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iB4] = 0.0;
							}
							else {
								potent[GRADZPAM][iB4] = potent[GRADZPAM][iP];
							}

							// Тангенциальные компоненты.					
							potent[GRADXPAM][iB4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB4] = potent[GRADYPAM][iP];
						}
					}
				}

			}
			else {

				if (iE > -1) {
					if (bE) {
						integer inumber = iE - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iE] = potent[GRADXPAM][iP];
						}
						else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADXPAM][iE] = potent[GRADXPAM][iP];
						}
						else {
							potent[GRADXPAM][iE] = 0.0;
						}

						//potent[GRADXPAM][iE] = potent[GRADXPAM][iP];

						potent[GRADYPAM][iE] = potent[GRADYPAM][iP];
						potent[GRADZPAM][iE] = potent[GRADZPAM][iP];
					}
				}

				if (iW > -1) {
					if (bW) {
						integer inumber = iW - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iW] = potent[GRADXPAM][iP];
						}
						else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADXPAM][iW] = potent[GRADXPAM][iP];
						}
						else {
							// производная от давления по нормали равна нулю, таковы граничные условия.
							potent[GRADXPAM][iW] = 0.0;
						}

						//potent[GRADXPAM][iW] = potent[GRADXPAM][iP];

						potent[GRADYPAM][iW] = potent[GRADYPAM][iP];
						potent[GRADZPAM][iW] = potent[GRADZPAM][iP];
					}
				}

				if (iN > -1) {
					if (bN) {

						integer inumber = iN - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iN] = potent[GRADYPAM][iP];
						}
						else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADYPAM][iN] = potent[GRADYPAM][iP];
						}
						else {
							// производная от давления по нормали равна нулю, таковы граничные условия.
							potent[GRADYPAM][iN] = 0.0;
						}

						//potent[GRADYPAM][iN] = potent[GRADYPAM][iP];

						potent[GRADXPAM][iN] = potent[GRADXPAM][iP];
						potent[GRADZPAM][iN] = potent[GRADZPAM][iP];
					}
				}

				if (iS > -1) {
					if (bS) {

						integer inumber = iS - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iS] = potent[GRADYPAM][iP];
						}
						else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADYPAM][iS] = potent[GRADYPAM][iP];
						}
						else {
							// производная от давления по нормали равна нулю, таковы граничные условия.
							potent[GRADYPAM][iS] = 0.0;
						}

						//potent[GRADYPAM][iS] = potent[GRADYPAM][iP];

						potent[GRADXPAM][iS] = potent[GRADXPAM][iP];
						potent[GRADZPAM][iS] = potent[GRADZPAM][iP];
					}
				}

				if (iT > -1) {
					if (bT) {

						integer inumber = iT - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iT] = potent[GRADZPAM][iP];
						}
						else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADZPAM][iT] = potent[GRADZPAM][iP];
						}
						else {
							// производная от давления по нормали равна нулю, таковы граничные условия.
							potent[GRADZPAM][iT] = 0.0;
						}

						//potent[GRADZPAM][iT] = potent[GRADZPAM][iP];

						potent[GRADXPAM][iT] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iT] = potent[GRADYPAM][iP];
					}
				}

				if (iB > -1) {
					if (bB) {

						integer inumber = iB - maxelm;
						if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iB] = potent[GRADZPAM][iP];
						}
						else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
							potent[GRADZPAM][iB] = potent[GRADZPAM][iP];
						}
						else {
							// производная от давления по нормали равна нулю, таковы граничные условия.
							potent[GRADZPAM][iB] = 0.0;
						}

						//potent[GRADZPAM][iB] = potent[GRADZPAM][iP];

						potent[GRADXPAM][iB] = potent[GRADXPAM][iP];
						potent[GRADYPAM][iB] = potent[GRADYPAM][iP];
					}
				}

				if (b_on_adaptive_local_refinement_mesh) {

					if (iE2 > -1) {
						if (bE2) {
							integer inumber = iE2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADXPAM][iE2] = potent[GRADXPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iE2] = potent[GRADXPAM][iP];
							}
							else {
								potent[GRADXPAM][iE2] = 0.0;
							}

							//potent[GRADXPAM][iE2] = potent[GRADXPAM][iP];

							potent[GRADYPAM][iE2] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE2] = potent[GRADZPAM][iP];
						}
					}

					if (iW2 > -1) {
						if (bW2) {
							integer inumber = iW2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADXPAM][iW2] = potent[GRADXPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iW2] = potent[GRADXPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADXPAM][iW2] = 0.0;
							}

							//potent[GRADXPAM][iW2] = potent[GRADXPAM][iP];

							potent[GRADYPAM][iW2] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW2] = potent[GRADZPAM][iP];
						}
					}

					if (iN2 > -1) {
						if (bN2) {

							integer inumber = iN2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADYPAM][iN2] = potent[GRADYPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iN2] = potent[GRADYPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADYPAM][iN2] = 0.0;
							}

							//potent[GRADYPAM][iN] = potent[GRADYPAM][iP];

							potent[GRADXPAM][iN2] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN2] = potent[GRADZPAM][iP];
						}
					}

					if (iS2 > -1) {
						if (bS2) {

							integer inumber = iS2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADYPAM][iS2] = potent[GRADYPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iS2] = potent[GRADYPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADYPAM][iS2] = 0.0;
							}

							//potent[GRADYPAM][iS2] = potent[GRADYPAM][iP];

							potent[GRADXPAM][iS2] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS2] = potent[GRADZPAM][iP];
						}
					}

					if (iT2 > -1) {
						if (bT2) {

							integer inumber = iT2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADZPAM][iT2] = potent[GRADZPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iT2] = potent[GRADZPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADZPAM][iT2] = 0.0;
							}

							//potent[GRADZPAM][iT2] = potent[GRADZPAM][iP];

							potent[GRADXPAM][iT2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT2] = potent[GRADYPAM][iP];
						}
					}

					if (iB2 > -1) {
						if (bB2) {

							integer inumber = iB2 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADZPAM][iB2] = potent[GRADZPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iB2] = potent[GRADZPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADZPAM][iB2] = 0.0;
							}

							//potent[GRADZPAM][iB2] = potent[GRADZPAM][iP];

							potent[GRADXPAM][iB2] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB2] = potent[GRADYPAM][iP];
						}
					}


					if (iE3 > -1) {
						if (bE3) {
							integer inumber = iE3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADXPAM][iE3] = potent[GRADXPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iE3] = potent[GRADXPAM][iP];
							}
							else {
								potent[GRADXPAM][iE3] = 0.0;
							}

							//potent[GRADXPAM][iE3] = potent[GRADXPAM][iP];

							potent[GRADYPAM][iE3] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE3] = potent[GRADZPAM][iP];
						}
					}

					if (iW3 > -1) {
						if (bW3) {
							integer inumber = iW3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADXPAM][iW3] = potent[GRADXPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iW3] = potent[GRADXPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADXPAM][iW3] = 0.0;
							}

							//potent[GRADXPAM][iW3] = potent[GRADXPAM][iP];

							potent[GRADYPAM][iW3] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW3] = potent[GRADZPAM][iP];
						}
					}

					if (iN3 > -1) {
						if (bN3) {

							integer inumber = iN3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADYPAM][iN3] = potent[GRADYPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iN3] = potent[GRADYPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADYPAM][iN3] = 0.0;
							}

							//potent[GRADYPAM][iN3] = potent[GRADYPAM][iP];

							potent[GRADXPAM][iN3] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN3] = potent[GRADZPAM][iP];
						}
					}

					if (iS3 > -1) {
						if (bS3) {

							integer inumber = iS3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADYPAM][iS3] = potent[GRADYPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iS3] = potent[GRADYPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADYPAM][iS3] = 0.0;
							}

							//potent[GRADYPAM][iS3] = potent[GRADYPAM][iP];

							potent[GRADXPAM][iS3] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS3] = potent[GRADZPAM][iP];
						}
					}

					if (iT3 > -1) {
						if (bT3) {

							integer inumber = iT3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADZPAM][iT3] = potent[GRADZPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iT3] = potent[GRADZPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADZPAM][iT3] = 0.0;
							}

							//potent[GRADZPAM][iT3] = potent[GRADZPAM][iP];

							potent[GRADXPAM][iT3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT3] = potent[GRADYPAM][iP];
						}
					}

					if (iB3 > -1) {
						if (bB3) {

							integer inumber = iB3 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADZPAM][iB3] = potent[GRADZPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iB3] = potent[GRADZPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADZPAM][iB3] = 0.0;
							}

							//potent[GRADZPAM][iB3] = potent[GRADZPAM][iP];

							potent[GRADXPAM][iB3] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB3] = potent[GRADYPAM][iP];
						}
					}

					if (iE4 > -1) {
						if (bE4) {
							integer inumber = iE4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADXPAM][iE4] = potent[GRADXPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iE4] = potent[GRADXPAM][iP];
							}
							else {
								potent[GRADXPAM][iE4] = 0.0;
							}

							//potent[GRADXPAM][iE4] = potent[GRADXPAM][iP];

							potent[GRADYPAM][iE4] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iE4] = potent[GRADZPAM][iP];
						}
					}

					if (iW4 > -1) {
						if (bW4) {
							integer inumber = iW4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADXPAM][iW4] = potent[GRADXPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADXPAM][iW4] = potent[GRADXPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADXPAM][iW4] = 0.0;
							}

							//potent[GRADXPAM][iW4] = potent[GRADXPAM][iP];

							potent[GRADYPAM][iW4] = potent[GRADYPAM][iP];
							potent[GRADZPAM][iW4] = potent[GRADZPAM][iP];
						}
					}

					if (iN4 > -1) {
						if (bN4) {

							integer inumber = iN4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADYPAM][iN4] = potent[GRADYPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iN4] = potent[GRADYPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADYPAM][iN4] = 0.0;
							}

							//potent[GRADYPAM][iN4] = potent[GRADYPAM][iP];

							potent[GRADXPAM][iN4] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iN4] = potent[GRADZPAM][iP];
						}
					}

					if (iS4 > -1) {
						if (bS4) {

							integer inumber = iS4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADYPAM][iS4] = potent[GRADYPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADYPAM][iS4] = potent[GRADYPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADYPAM][iS4] = 0.0;
							}

							//potent[GRADYPAM][iS4] = potent[GRADYPAM][iP];

							potent[GRADXPAM][iS4] = potent[GRADXPAM][iP];
							potent[GRADZPAM][iS4] = potent[GRADZPAM][iP];
						}
					}

					if (iT4 > -1) {
						if (bT4) {

							integer inumber = iT4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADZPAM][iT4] = potent[GRADZPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iT4] = potent[GRADZPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADZPAM][iT4] = 0.0;
							}

							//potent[GRADZPAM][iT4] = potent[GRADZPAM][iP];

							potent[GRADXPAM][iT4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iT4] = potent[GRADYPAM][iP];
						}
					}

					if (iB4 > -1) {
						if (bB4) {

							integer inumber = iB4 - maxelm;
							if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
								potent[GRADZPAM][iB4] = potent[GRADZPAM][iP];
							}
							else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
								potent[GRADZPAM][iB4] = potent[GRADZPAM][iP];
							}
							else {
								// производная от давления по нормали равна нулю, таковы граничные условия.
								potent[GRADZPAM][iB4] = 0.0;
							}

							//potent[GRADZPAM][iB4] = potent[GRADZPAM][iP];

							potent[GRADXPAM][iB4] = potent[GRADXPAM][iP];
							potent[GRADYPAM][iB4] = potent[GRADYPAM][iP];
						}
					}
				}

			}

		}
#endif

#if (interpol==1) 
{

			// не работает на АЛИС.
			if (b_on_adaptive_local_refinement_mesh) {
				printf("function green_gaussPAM in module greengauss.c else if (interpol==1) { not worked in ALICE mesh...\n ");
				system("pause");
				exit(1);
			}

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

		// Если строка с пометкой <-- раскомментирована то градиент поправки давления линейно интерполируется на границу 
		// расчётной области изнутри расчётной области.

		if (bE) {

			integer inumber=iE-maxelm;
			if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
				//potent[GRADXPAM][iE]=potent[GRADXPAM][iP]+(dxe/dxw)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iW]);
				potent[GRADXPAM][iE]=0.0;
			}
			else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
				//potent[GRADXPAM][iE]=potent[GRADXPAM][iP]+(dxe/dxw)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iW]);
				potent[GRADXPAM][iE] = 0.0;
			}
			else {
				// производная от давления по нормали равна нулю, таковы граничные условия.
				potent[GRADXPAM][iE]=0.0;
			}

			potent[GRADYPAM][iE]=potent[GRADYPAM][iP]+(dxe/dxw)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iW]);
			potent[GRADZPAM][iE]=potent[GRADZPAM][iP]+(dxe/dxw)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iW]);
			potent[GRADXPAM][iE]=potent[GRADXPAM][iP]+(dxe/dxw)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iW]); // <--

		}

		if (bW) {

			integer inumber=iW-maxelm;
			if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
				//potent[GRADXPAM][iW]=potent[GRADXPAM][iP]+(dxw/dxe)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iE]);
				potent[GRADXPAM][iW]=0.0;
			}
			else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
				//potent[GRADXPAM][iW]=potent[GRADXPAM][iP]+(dxw/dxe)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iE]);
				potent[GRADXPAM][iW] = 0.0;
			}
			else {
				// производная от давления по нормали равна нулю, таковы граничные условия.
				potent[GRADXPAM][iW]=0.0;
			}

			potent[GRADYPAM][iW]=potent[GRADYPAM][iP]+(dxw/dxe)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iE]);
			potent[GRADZPAM][iW]=potent[GRADZPAM][iP]+(dxw/dxe)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iE]);
			potent[GRADXPAM][iW]=potent[GRADXPAM][iP]+(dxw/dxe)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iE]); // <--

		}

		if (bN) {

			integer inumber=iN-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADYPAM][iN]=potent[GRADYPAM][iP]+(dyn/dys)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iS]);
					potent[GRADYPAM][iN]=0.0;
				}
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
					//potent[GRADYPAM][iN]=potent[GRADYPAM][iP]+(dyn/dys)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iS]);
					potent[GRADYPAM][iN] = 0.0;
				}
				else {
					// производная от давления по нормали равна нулю, таковы граничные условия.
					potent[GRADYPAM][iN]=0.0;
				}

			potent[GRADXPAM][iN]=potent[GRADXPAM][iP]+(dyn/dys)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iS]);
			potent[GRADZPAM][iN]=potent[GRADZPAM][iP]+(dyn/dys)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iS]);
			potent[GRADYPAM][iN]=potent[GRADYPAM][iP]+(dyn/dys)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iS]); // <--
		}

		if (bS) {

			integer inumber=iS-maxelm;
			if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
				//potent[GRADYPAM][iS]=potent[GRADYPAM][iP]+(dys/dyn)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iN]);
				potent[GRADYPAM][iS]=0.0;
			}
			else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
				//potent[GRADYPAM][iS]=potent[GRADYPAM][iP]+(dys/dyn)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iN]);
				potent[GRADYPAM][iS] = 0.0;
			}
			else {
			    // производная от давления по нормали равна нулю, таковы граничные условия.
				potent[GRADYPAM][iS]=0.0;
			}


			potent[GRADXPAM][iS]=potent[GRADXPAM][iP]+(dys/dyn)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iN]);
			potent[GRADZPAM][iS]=potent[GRADZPAM][iP]+(dys/dyn)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iN]);
			potent[GRADYPAM][iS]=potent[GRADYPAM][iP]+(dys/dyn)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iN]); //<--

		}

		if (bT) {

            integer inumber=iT-maxelm;
			if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
				//potent[GRADZPAM][iT]=potent[GRADZPAM][iP]+(dzt/dzb)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iB]);
				potent[GRADZPAM][iT]=0.0;
			}
			else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
				//potent[GRADZPAM][iT]=potent[GRADZPAM][iP]+(dzt/dzb)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iB]);
				potent[GRADZPAM][iT] = 0.0;
			}
			else {
				// производная от давления по нормали равна нулю, таковы граничные условия.
				potent[GRADZPAM][iT]=0.0;
			}

			potent[GRADXPAM][iT]=potent[GRADXPAM][iP]+(dzt/dzb)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iB]);
			potent[GRADYPAM][iT]=potent[GRADYPAM][iP]+(dzt/dzb)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iB]);
			potent[GRADZPAM][iT]=potent[GRADZPAM][iP]+(dzt/dzb)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iB]); // <--
		}

		if (bB) {

			integer inumber=iB-maxelm;
			if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
				//potent[GRADZPAM][iB]=potent[GRADZPAM][iP]+(dzb/dzt)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iT]);
				potent[GRADZPAM][iB]=0.0;
			}
			else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
				//potent[GRADZPAM][iB]=potent[GRADZPAM][iP]+(dzb/dzt)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iT]);
				potent[GRADZPAM][iB] = 0.0;
			}
			else {
				// производная от давления по нормали равна нулю, таковы граничные условия.
				potent[GRADZPAM][iB]=0.0;
			}

			potent[GRADXPAM][iB]=potent[GRADXPAM][iP]+(dzb/dzt)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iT]);
			potent[GRADYPAM][iB]=potent[GRADYPAM][iP]+(dzb/dzt)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iT]);
			potent[GRADZPAM][iB]=potent[GRADZPAM][iP]+(dzb/dzt)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iT]); // <--
		}
		}
#endif
	}


} // green_gaussPAM
  
// 13 апреля 2015 года.
// вычисление градиентов механических величин с помощью теоремы Грина-Гаусса. 
// Механические величины не определены на границах расчетной области являющихся также границами ко.
// Сносим значение из центра ко на его границу.
void green_gaussMechanical(integer iP, doublereal*& potent, int**& nvtx, TOCHKA*& pa,
	int***& neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND*& border_neighbor, doublereal*& Tx, doublereal*& Ty, doublereal*& Tz,
	integer* ilevel_alice) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты поправки давления для внутренних КО.
	// если bbond   то будут вычислены значения только в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// Внимание ! память под Tx, Ty, Tz предполагается выделенной заранее.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}
	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe = 0.5 * dx, dxw = 0.5 * dx, dyn = 0.5 * dy, dys = 0.5 * dy, dzt = 0.5 * dz, dzb = 0.5 * dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5 * (pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5 * (pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5 * (pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5 * (pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5 * (pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5 * (pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}

	doublereal dxe2 = 0.5 * dx, dxw2 = 0.5 * dx, dyn2 = 0.5 * dy, dys2 = 0.5 * dy, dzt2 = 0.5 * dz, dzb2 = 0.5 * dz;
	doublereal dxe3 = 0.5 * dx, dxw3 = 0.5 * dx, dyn3 = 0.5 * dy, dys3 = 0.5 * dy, dzt3 = 0.5 * dz, dzb3 = 0.5 * dz;
	doublereal dxe4 = 0.5 * dx, dxw4 = 0.5 * dx, dyn4 = 0.5 * dy, dys4 = 0.5 * dy, dzt4 = 0.5 * dz, dzb4 = 0.5 * dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5 * (pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5 * (pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5 * (pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5 * (pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5 * (pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5 * (pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5 * (pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5 * (pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5 * (pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5 * (pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5 * (pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5 * (pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5 * (pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5 * (pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5 * (pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5 * (pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5 * (pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5 * (pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
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

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy * dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy * dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx * dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx * dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx * dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx * dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}

	doublereal Te = 0.0, Tw = 0.0, Tn = 0.0, Ts = 0.0, Tt = 0.0, Tb = 0.0;
	doublereal Te2 = 0.0, Tw2 = 0.0, Tn2 = 0.0, Ts2 = 0.0, Tt2 = 0.0, Tb2 = 0.0;
	doublereal Te3 = 0.0, Tw3 = 0.0, Tn3 = 0.0, Ts3 = 0.0, Tt3 = 0.0, Tb3 = 0.0;
	doublereal Te4 = 0.0, Tw4 = 0.0, Tn4 = 0.0, Ts4 = 0.0, Tt4 = 0.0, Tb4 = 0.0;

	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем поправку давления на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		if (iE > -1) {
			if (!bE) Te = feplus * potent[iE] + (1.0 - feplus) * potent[iP]; else Te = potent[iP];
		}
		if (iW > -1) {
			if (!bW) Tw = fwplus * potent[iW] + (1.0 - fwplus) * potent[iP]; else Tw = potent[iP];
		}
		if (iN > -1) {
			if (!bN) Tn = fnplus * potent[iN] + (1.0 - fnplus) * potent[iP]; else Tn = potent[iP];
		}
		if (iS > -1) {
			if (!bS) Ts = fsplus * potent[iS] + (1.0 - fsplus) * potent[iP]; else Ts = potent[iP];
		}
		if (iT > -1) {
			if (!bT) Tt = ftplus * potent[iT] + (1.0 - ftplus) * potent[iP]; else Tt = potent[iP];
		}
		if (iB > -1) {
			if (!bB) Tb = fbplus * potent[iB] + (1.0 - fbplus) * potent[iP]; else Tb = potent[iP];
		}

		if (iE2 > -1) {
			if (!bE2) Te2 = feplus2 * potent[iE2] + (1.0 - feplus2) * potent[iP]; else Te2 = potent[iP];
		}
		if (iW2 > -1) {
			if (!bW2) Tw2 = fwplus2 * potent[iW2] + (1.0 - fwplus2) * potent[iP]; else Tw2 = potent[iP];
		}
		if (iN2 > -1) {
			if (!bN2) Tn2 = fnplus2 * potent[iN2] + (1.0 - fnplus2) * potent[iP]; else Tn2 = potent[iP];
		}
		if (iS2 > -1) {
			if (!bS2) Ts2 = fsplus2 * potent[iS2] + (1.0 - fsplus2) * potent[iP]; else Ts2 = potent[iP];
		}
		if (iT2 > -1) {
			if (!bT2) Tt2 = ftplus2 * potent[iT2] + (1.0 - ftplus2) * potent[iP]; else Tt2 = potent[iP];
		}
		if (iB2 > -1) {
			if (!bB2) Tb2 = fbplus2 * potent[iB2] + (1.0 - fbplus2) * potent[iP]; else Tb2 = potent[iP];
		}

		if (iE3 > -1) {
			if (!bE3) Te3 = feplus3 * potent[iE3] + (1.0 - feplus3) * potent[iP]; else Te3 = potent[iP];
		}
		if (iW3 > -1) {
			if (!bW3) Tw3 = fwplus3 * potent[iW3] + (1.0 - fwplus3) * potent[iP]; else Tw3 = potent[iP];
		}
		if (iN3 > -1) {
			if (!bN3) Tn3 = fnplus3 * potent[iN3] + (1.0 - fnplus3) * potent[iP]; else Tn3 = potent[iP];
		}
		if (iS3 > -1) {
			if (!bS3) Ts3 = fsplus3 * potent[iS3] + (1.0 - fsplus3) * potent[iP]; else Ts3 = potent[iP];
		}
		if (iT3 > -1) {
			if (!bT3) Tt3 = ftplus3 * potent[iT3] + (1.0 - ftplus3) * potent[iP]; else Tt3 = potent[iP];
		}
		if (iB3 > -1) {
			if (!bB3) Tb3 = fbplus3 * potent[iB3] + (1.0 - fbplus3) * potent[iP]; else Tb3 = potent[iP];
		}

		if (iE4 > -1) {
			if (!bE4) Te4 = feplus4 * potent[iE4] + (1.0 - feplus4) * potent[iP]; else Te4 = potent[iP];
		}
		if (iW4 > -1) {
			if (!bW4) Tw4 = fwplus4 * potent[iW4] + (1.0 - fwplus4) * potent[iP]; else Tw4 = potent[iP];
		}
		if (iN4 > -1) {
			if (!bN4) Tn4 = fnplus4 * potent[iN4] + (1.0 - fnplus4) * potent[iP]; else Tn4 = potent[iP];
		}
		if (iS4 > -1) {
			if (!bS4) Ts4 = fsplus4 * potent[iS4] + (1.0 - fsplus4) * potent[iP]; else Ts4 = potent[iP];
		}
		if (iT4 > -1) {
			if (!bT4) Tt4 = ftplus4 * potent[iT4] + (1.0 - ftplus4) * potent[iP]; else Tt4 = potent[iP];
		}
		if (iB4 > -1) {
			if (!bB4) Tb4 = fbplus4 * potent[iB4] + (1.0 - fbplus4) * potent[iP]; else Tb4 = potent[iP];
		}

		// градиент Температуры. 20.03.2019
		//Tx[iP]=(Te-Tw)/dx;
		//Ty[iP]=(Tn-Ts)/dy;
		//Tz[iP]=(Tt-Tb)/dz;
		Tx[iP] = (Te * dSqe / (dy * dz) + Te2 * dSqe2 / (dy * dz) + Te3 * dSqe3 / (dy * dz) + Te4 * dSqe4 / (dy * dz) - (Tw * dSqw / (dy * dz) + Tw2 * dSqw2 / (dy * dz) + Tw3 * dSqw3 / (dy * dz) + Tw4 * dSqw4 / (dy * dz))) / dx;
		Ty[iP] = (Tn * dSqn / (dx * dz) + Tn2 * dSqn2 / (dx * dz) + Tn3 * dSqn3 / (dx * dz) + Tn4 * dSqn4 / (dx * dz) - (Ts * dSqs / (dx * dz) + Ts2 * dSqs2 / (dx * dz) + Ts3 * dSqs3 / (dx * dz) + Ts4 * dSqs4 / (dx * dz))) / dy;
		Tz[iP] = (Tt * dSqt / (dx * dy) + Tt2 * dSqt2 / (dx * dy) + Tt3 * dSqt3 / (dx * dy) + Tt4 * dSqt4 / (dx * dy) - (Tb * dSqb / (dx * dy) + Tb2 * dSqb2 / (dx * dy) + Tb3 * dSqb3 / (dx * dy) + Tb4 * dSqb4 / (dx * dy))) / dz;


	}
	else {
		// На АЛИС сетках работает только значение interpol==0.
		const integer interpol = 0; // 0 для переноса из центра на грань или 1 при линейной интерполяции.


#if (interpol==0) 
		{

			if (iE > -1) {
				if (bE) {
					Tx[iE] = Tx[iP];
					Ty[iE] = Ty[iP];
					Tz[iE] = Tz[iP];
				}
			}

			if (iW > -1) {
				if (bW) {
					Tx[iW] = Tx[iP];
					Ty[iW] = Ty[iP];
					Tz[iW] = Tz[iP];
				}
			}

			if (iN > -1) {
				if (bN) {
					Tx[iN] = Tx[iP];
					Ty[iN] = Ty[iP];
					Tz[iN] = Tz[iP];
				}
			}

			if (iS > -1) {
				if (bS) {
					Tx[iS] = Tx[iP];
					Ty[iS] = Ty[iP];
					Tz[iS] = Tz[iP];
				}
			}

			if (iT > -1) {
				if (bT) {
					Tx[iT] = Tx[iP];
					Ty[iT] = Ty[iP];
					Tz[iT] = Tz[iP];
				}
			}

			if (iB > -1) {
				if (bB) {
					Tx[iB] = Tx[iP];
					Ty[iB] = Ty[iP];
					Tz[iB] = Tz[iP];
				}
			}

			if (iE2 > -1) {
				if (bE2) {
					Tx[iE2] = Tx[iP];
					Ty[iE2] = Ty[iP];
					Tz[iE2] = Tz[iP];
				}
			}

			if (iW2 > -1) {
				if (bW2) {
					Tx[iW2] = Tx[iP];
					Ty[iW2] = Ty[iP];
					Tz[iW2] = Tz[iP];
				}
			}

			if (iN2 > -1) {
				if (bN2) {
					Tx[iN2] = Tx[iP];
					Ty[iN2] = Ty[iP];
					Tz[iN2] = Tz[iP];
				}
			}

			if (iS2 > -1) {
				if (bS2) {
					Tx[iS2] = Tx[iP];
					Ty[iS2] = Ty[iP];
					Tz[iS2] = Tz[iP];
				}
			}

			if (iT2 > -1) {
				if (bT2) {
					Tx[iT2] = Tx[iP];
					Ty[iT2] = Ty[iP];
					Tz[iT2] = Tz[iP];
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					Tx[iB2] = Tx[iP];
					Ty[iB2] = Ty[iP];
					Tz[iB2] = Tz[iP];
				}
			}

			if (iE3 > -1) {
				if (bE3) {
					Tx[iE3] = Tx[iP];
					Ty[iE3] = Ty[iP];
					Tz[iE3] = Tz[iP];
				}
			}

			if (iW3 > -1) {
				if (bW3) {
					Tx[iW3] = Tx[iP];
					Ty[iW3] = Ty[iP];
					Tz[iW3] = Tz[iP];
				}
			}

			if (iN3 > -1) {
				if (bN3) {
					Tx[iN3] = Tx[iP];
					Ty[iN3] = Ty[iP];
					Tz[iN3] = Tz[iP];
				}
			}

			if (iS3 > -1) {
				if (bS3) {
					Tx[iS3] = Tx[iP];
					Ty[iS3] = Ty[iP];
					Tz[iS3] = Tz[iP];
				}
			}

			if (iT3 > -1) {
				if (bT3) {
					Tx[iT3] = Tx[iP];
					Ty[iT3] = Ty[iP];
					Tz[iT3] = Tz[iP];
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					Tx[iB3] = Tx[iP];
					Ty[iB3] = Ty[iP];
					Tz[iB3] = Tz[iP];
				}
			}

			if (iE4 > -1) {
				if (bE4) {
					Tx[iE4] = Tx[iP];
					Ty[iE4] = Ty[iP];
					Tz[iE4] = Tz[iP];
				}
			}

			if (iW4 > -1) {
				if (bW4) {
					Tx[iW4] = Tx[iP];
					Ty[iW4] = Ty[iP];
					Tz[iW4] = Tz[iP];
				}
			}

			if (iN4 > -1) {
				if (bN4) {
					Tx[iN4] = Tx[iP];
					Ty[iN4] = Ty[iP];
					Tz[iN4] = Tz[iP];
				}
			}

			if (iS4 > -1) {
				if (bS4) {
					Tx[iS4] = Tx[iP];
					Ty[iS4] = Ty[iP];
					Tz[iS4] = Tz[iP];
				}
			}

			if (iT4 > -1) {
				if (bT4) {
					Tx[iT4] = Tx[iP];
					Ty[iT4] = Ty[iP];
					Tz[iT4] = Tz[iP];
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					Tx[iB4] = Tx[iP];
					Ty[iB4] = Ty[iP];
					Tz[iB4] = Tz[iP];
				}
			}


		}
#endif

#if (interpol == 1)
		{
			if (b_on_adaptive_local_refinement_mesh) {
				printf("Linear interpolation not work on adaptive local refinement mesh. !!!\n");
				printf("LOCATION: function green_gaussTemperature in module greengauss.c\n");
				system("PAUSE");
				exit(1);
			}
			else {
				// граничные узлы.
				// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

				// Если строка с пометкой <-- раскомментирована то градиент поправки давления линейно интерполлируется на границу 
				// расчётной области изнутри расчётной области.

				if (bE) {


					Ty[iE] = Ty[iP] + (dxe / dxw) * (Ty[iP] - Ty[iW]);
					Tz[iE] = Tz[iP] + (dxe / dxw) * (Tz[iP] - Tz[iW]);
					Tx[iE] = Tx[iP] + (dxe / dxw) * (Tx[iP] - Tx[iW]); // <--

				}

				if (bW) {


					Ty[iW] = Ty[iP] + (dxw / dxe) * (Ty[iP] - Ty[iE]);
					Tz[iW] = Tz[iP] + (dxw / dxe) * (Tz[iP] - Tz[iE]);
					Tx[iW] = Tx[iP] + (dxw / dxe) * (Tx[iP] - Tx[iE]); // <--

				}

				if (bN) {



					Tx[iN] = Tx[iP] + (dyn / dys) * (Tx[iP] - Tx[iS]);
					Tz[iN] = Tz[iP] + (dyn / dys) * (Tz[iP] - Tz[iS]);
					Ty[iN] = Ty[iP] + (dyn / dys) * (Ty[iP] - Ty[iS]); // <--
				}

				if (bS) {



					Tx[iS] = Tx[iP] + (dys / dyn) * (Tx[iP] - Tx[iN]);
					Tz[iS] = Tz[iP] + (dys / dyn) * (Tz[iP] - Tz[iN]);
					Ty[iS] = Ty[iP] + (dys / dyn) * (Ty[iP] - Ty[iN]); //<--

				}

				if (bT) {



					Tx[iT] = Tx[iP] + (dzt / dzb) * (Tx[iP] - Tx[iB]);
					Ty[iT] = Ty[iP] + (dzt / dzb) * (Ty[iP] - Ty[iB]);
					Tz[iT] = Tz[iP] + (dzt / dzb) * (Tz[iP] - Tz[iB]); // <--
				}

				if (bB) {



					Tx[iB] = Tx[iP] + (dzb / dzt) * (Tx[iP] - Tx[iT]);
					Ty[iB] = Ty[iP] + (dzb / dzt) * (Ty[iP] - Ty[iT]);
					Tz[iB] = Tz[iP] + (dzb / dzt) * (Tz[iP] - Tz[iT]); // <--
				}
			}
		}
#endif
	}


} // green_gaussMechanical

// 13 апреля 2015 года.
// вычисление градиентов Температуры с помощью теоремы Грина-Гаусса. 
void green_gaussTemperature(integer iP, doublereal* &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
					BOUND* &border_neighbor, doublereal* &Tx, doublereal* &Ty, doublereal* &Tz,
	integer *ilevel_alice) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты поправки давления для внутренних КО.
	// если bbond   то будут вычислены значения только в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// Внимание ! память под Tx, Ty, Tz предполагается выделенной заранее.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP]; 
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
		
		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
		
		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}
	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}

	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus=0.5*dx/dxe;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}

	doublereal Te = 0.0, Tw = 0.0, Tn = 0.0, Ts = 0.0, Tt = 0.0, Tb = 0.0;
	doublereal Te2 = 0.0, Tw2 = 0.0, Tn2 = 0.0, Ts2 = 0.0, Tt2 = 0.0, Tb2 = 0.0;
	doublereal Te3 = 0.0, Tw3 = 0.0, Tn3 = 0.0, Ts3 = 0.0, Tt3 = 0.0, Tb3 = 0.0;
	doublereal Te4 = 0.0, Tw4 = 0.0, Tn4 = 0.0, Ts4 = 0.0, Tt4 = 0.0, Tb4 = 0.0;

    if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем поправку давления на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		if (iE > -1) {
			if (!bE) Te = feplus*potent[iE] + (1.0 - feplus)*potent[iP]; else Te = potent[iE];
		}
		if (iW > -1) {
			if (!bW) Tw = fwplus*potent[iW] + (1.0 - fwplus)*potent[iP]; else Tw = potent[iW];
		}
		if (iN > -1) {
			if (!bN) Tn = fnplus*potent[iN] + (1.0 - fnplus)*potent[iP]; else Tn = potent[iN];
		}
		if (iS > -1) {
			if (!bS) Ts = fsplus*potent[iS] + (1.0 - fsplus)*potent[iP]; else Ts = potent[iS];
		}
		if (iT > -1) {
			if (!bT) Tt = ftplus*potent[iT] + (1.0 - ftplus)*potent[iP]; else Tt = potent[iT];
		}
		if (iB > -1) {
			if (!bB) Tb = fbplus*potent[iB] + (1.0 - fbplus)*potent[iP]; else Tb = potent[iB];
		}

		if (iE2 > -1) {
			if (!bE2) Te2 = feplus2*potent[iE2] + (1.0 - feplus2)*potent[iP]; else Te2 = potent[iE2];
		}
		if (iW2 > -1) {
			if (!bW2) Tw2 = fwplus2*potent[iW2] + (1.0 - fwplus2)*potent[iP]; else Tw2 = potent[iW2];
		}
		if (iN2 > -1) {
			if (!bN2) Tn2 = fnplus2*potent[iN2] + (1.0 - fnplus2)*potent[iP]; else Tn2 = potent[iN2];
		}
		if (iS2 > -1) {
			if (!bS2) Ts2 = fsplus2*potent[iS2] + (1.0 - fsplus2)*potent[iP]; else Ts2 = potent[iS2];
		}
		if (iT2 > -1) {
			if (!bT2) Tt2 = ftplus2*potent[iT2] + (1.0 - ftplus2)*potent[iP]; else Tt2 = potent[iT2];
		}
		if (iB2 > -1) {
			if (!bB2) Tb2 = fbplus2*potent[iB2] + (1.0 - fbplus2)*potent[iP]; else Tb2 = potent[iB2];
		}

		if (iE3 > -1) {
			if (!bE3) Te3 = feplus3*potent[iE3] + (1.0 - feplus3)*potent[iP]; else Te3 = potent[iE3];
		}
		if (iW3 > -1) {
			if (!bW3) Tw3 = fwplus3*potent[iW3] + (1.0 - fwplus3)*potent[iP]; else Tw3 = potent[iW3];
		}
		if (iN3 > -1) {
			if (!bN3) Tn3 = fnplus3*potent[iN3] + (1.0 - fnplus3)*potent[iP]; else Tn3 = potent[iN3];
		}
		if (iS3 > -1) {
			if (!bS3) Ts3 = fsplus3*potent[iS3] + (1.0 - fsplus3)*potent[iP]; else Ts3 = potent[iS3];
		}
		if (iT3 > -1) {
			if (!bT3) Tt3 = ftplus3*potent[iT3] + (1.0 - ftplus3)*potent[iP]; else Tt3 = potent[iT3];
		}
		if (iB3 > -1) {
			if (!bB3) Tb3 = fbplus3*potent[iB3] + (1.0 - fbplus3)*potent[iP]; else Tb3 = potent[iB3];
		}

		if (iE4 > -1) {
			if (!bE4) Te4 = feplus4*potent[iE4] + (1.0 - feplus4)*potent[iP]; else Te4 = potent[iE4];
		}
		if (iW4 > -1) {
			if (!bW4) Tw4 = fwplus4*potent[iW4] + (1.0 - fwplus4)*potent[iP]; else Tw4 = potent[iW4];
		}
		if (iN4 > -1) {
			if (!bN4) Tn4 = fnplus4*potent[iN4] + (1.0 - fnplus4)*potent[iP]; else Tn4 = potent[iN4];
		}
		if (iS4 > -1) {
			if (!bS4) Ts4 = fsplus4*potent[iS4] + (1.0 - fsplus4)*potent[iP]; else Ts4 = potent[iS4];
		}
		if (iT4 > -1) {
			if (!bT4) Tt4 = ftplus4*potent[iT4] + (1.0 - ftplus4)*potent[iP]; else Tt4 = potent[iT4];
		}
		if (iB4 > -1) {
			if (!bB4) Tb4 = fbplus4*potent[iB4] + (1.0 - fbplus4)*potent[iP]; else Tb4 = potent[iB4];
		}

        // градиент Температуры. 20.03.2019
	    //Tx[iP]=(Te-Tw)/dx;
	    //Ty[iP]=(Tn-Ts)/dy;
	    //Tz[iP]=(Tt-Tb)/dz;
		Tx[iP] = (Te*dSqe / (dy*dz) + Te2 * dSqe2 / (dy*dz) + Te3 * dSqe3 / (dy*dz) + Te4 * dSqe4 / (dy*dz) - (Tw*dSqw / (dy*dz) + Tw2 * dSqw2 / (dy*dz) + Tw3 * dSqw3 / (dy*dz) + Tw4 * dSqw4 / (dy*dz))) / dx;
		Ty[iP] = (Tn*dSqn / (dx*dz) + Tn2 * dSqn2 / (dx*dz) + Tn3 * dSqn3 / (dx*dz) + Tn4 * dSqn4 / (dx*dz) - (Ts*dSqs / (dx*dz) + Ts2 * dSqs2 / (dx*dz) + Ts3 * dSqs3 / (dx*dz) + Ts4 * dSqs4 / (dx*dz))) / dy;
		Tz[iP] = (Tt*dSqt / (dx*dy) + Tt2 * dSqt2 / (dx*dy) + Tt3 * dSqt3 / (dx*dy) + Tt4 * dSqt4 / (dx*dy) - (Tb*dSqb / (dx*dy) + Tb2 * dSqb2 / (dx*dy) + Tb3 * dSqb3 / (dx*dy) + Tb4 * dSqb4 / (dx*dy))) / dz;


   }
   else {
	   // На АЛИС сетках работает только значение interpol==0.
	    const integer interpol=0; // 0 для переноса из центра на грань или 1 при линейной интерполяции.


#if (interpol==0) 
		{

			if (iE > -1) {
				if (bE) {
					Tx[iE] = Tx[iP];
					Ty[iE] = Ty[iP];
					Tz[iE] = Tz[iP];
				}
			}

			if (iW > -1) {
				if (bW) {
					Tx[iW] = Tx[iP];
					Ty[iW] = Ty[iP];
					Tz[iW] = Tz[iP];
				}
			}

			if (iN > -1) {
				if (bN) {
					Tx[iN] = Tx[iP];
					Ty[iN] = Ty[iP];
					Tz[iN] = Tz[iP];
				}
			}

			if (iS > -1) {
				if (bS) {
					Tx[iS] = Tx[iP];
					Ty[iS] = Ty[iP];
					Tz[iS] = Tz[iP];
				}
			}

			if (iT > -1) {
				if (bT) {
					Tx[iT] = Tx[iP];
					Ty[iT] = Ty[iP];
					Tz[iT] = Tz[iP];
				}
			}

			if (iB > -1) {
				if (bB) {
					Tx[iB] = Tx[iP];
					Ty[iB] = Ty[iP];
					Tz[iB] = Tz[iP];
				}
			}

			if (iE2 > -1) {
				if (bE2) {
					Tx[iE2] = Tx[iP];
					Ty[iE2] = Ty[iP];
					Tz[iE2] = Tz[iP];
				}
			}

			if (iW2 > -1) {
				if (bW2) {
					Tx[iW2] = Tx[iP];
					Ty[iW2] = Ty[iP];
					Tz[iW2] = Tz[iP];
				}
			}

			if (iN2 > -1) {
				if (bN2) {
					Tx[iN2] = Tx[iP];
					Ty[iN2] = Ty[iP];
					Tz[iN2] = Tz[iP];
				}
			}

			if (iS2 > -1) {
				if (bS2) {
					Tx[iS2] = Tx[iP];
					Ty[iS2] = Ty[iP];
					Tz[iS2] = Tz[iP];
				}
			}

			if (iT2 > -1) {
				if (bT2) {
					Tx[iT2] = Tx[iP];
					Ty[iT2] = Ty[iP];
					Tz[iT2] = Tz[iP];
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					Tx[iB2] = Tx[iP];
					Ty[iB2] = Ty[iP];
					Tz[iB2] = Tz[iP];
				}
			}

			if (iE3 > -1) {
				if (bE3) {
					Tx[iE3] = Tx[iP];
					Ty[iE3] = Ty[iP];
					Tz[iE3] = Tz[iP];
				}
			}

			if (iW3 > -1) {
				if (bW3) {
					Tx[iW3] = Tx[iP];
					Ty[iW3] = Ty[iP];
					Tz[iW3] = Tz[iP];
				}
			}

			if (iN3 > -1) {
				if (bN3) {
					Tx[iN3] = Tx[iP];
					Ty[iN3] = Ty[iP];
					Tz[iN3] = Tz[iP];
				}
			}

			if (iS3 > -1) {
				if (bS3) {
					Tx[iS3] = Tx[iP];
					Ty[iS3] = Ty[iP];
					Tz[iS3] = Tz[iP];
				}
			}

			if (iT3 > -1) {
				if (bT3) {
					Tx[iT3] = Tx[iP];
					Ty[iT3] = Ty[iP];
					Tz[iT3] = Tz[iP];
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					Tx[iB3] = Tx[iP];
					Ty[iB3] = Ty[iP];
					Tz[iB3] = Tz[iP];
				}
			}

			if (iE4 > -1) {
				if (bE4) {
					Tx[iE4] = Tx[iP];
					Ty[iE4] = Ty[iP];
					Tz[iE4] = Tz[iP];
				}
			}

			if (iW4 > -1) {
				if (bW4) {
					Tx[iW4] = Tx[iP];
					Ty[iW4] = Ty[iP];
					Tz[iW4] = Tz[iP];
				}
			}

			if (iN4 > -1) {
				if (bN4) {
					Tx[iN4] = Tx[iP];
					Ty[iN4] = Ty[iP];
					Tz[iN4] = Tz[iP];
				}
			}

			if (iS4 > -1) {
				if (bS4) {
					Tx[iS4] = Tx[iP];
					Ty[iS4] = Ty[iP];
					Tz[iS4] = Tz[iP];
				}
			}

			if (iT4 > -1) {
				if (bT4) {
					Tx[iT4] = Tx[iP];
					Ty[iT4] = Ty[iP];
					Tz[iT4] = Tz[iP];
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					Tx[iB4] = Tx[iP];
					Ty[iB4] = Ty[iP];
					Tz[iB4] = Tz[iP];
				}
			}


		}
#endif
		
#if (interpol == 1)
 {
			if (b_on_adaptive_local_refinement_mesh) {
				printf("Linear interpolation not work on adaptive local refinement mesh. !!!\n");
				printf("LOCATION: function green_gaussTemperature in module greengauss.c\n");
				system("PAUSE");
				exit(1);
			}
			else {
				// граничные узлы.
				// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

				// Если строка с пометкой <-- раскомментирована то градиент поправки давления линейно интерполлируется на границу 
				// расчётной области изнутри расчётной области.

				if (bE) {


					Ty[iE] = Ty[iP] + (dxe / dxw)*(Ty[iP] - Ty[iW]);
					Tz[iE] = Tz[iP] + (dxe / dxw)*(Tz[iP] - Tz[iW]);
					Tx[iE] = Tx[iP] + (dxe / dxw)*(Tx[iP] - Tx[iW]); // <--

				}

				if (bW) {


					Ty[iW] = Ty[iP] + (dxw / dxe)*(Ty[iP] - Ty[iE]);
					Tz[iW] = Tz[iP] + (dxw / dxe)*(Tz[iP] - Tz[iE]);
					Tx[iW] = Tx[iP] + (dxw / dxe)*(Tx[iP] - Tx[iE]); // <--

				}

				if (bN) {



					Tx[iN] = Tx[iP] + (dyn / dys)*(Tx[iP] - Tx[iS]);
					Tz[iN] = Tz[iP] + (dyn / dys)*(Tz[iP] - Tz[iS]);
					Ty[iN] = Ty[iP] + (dyn / dys)*(Ty[iP] - Ty[iS]); // <--
				}

				if (bS) {



					Tx[iS] = Tx[iP] + (dys / dyn)*(Tx[iP] - Tx[iN]);
					Tz[iS] = Tz[iP] + (dys / dyn)*(Tz[iP] - Tz[iN]);
					Ty[iS] = Ty[iP] + (dys / dyn)*(Ty[iP] - Ty[iN]); //<--

				}

				if (bT) {



					Tx[iT] = Tx[iP] + (dzt / dzb)*(Tx[iP] - Tx[iB]);
					Ty[iT] = Ty[iP] + (dzt / dzb)*(Ty[iP] - Ty[iB]);
					Tz[iT] = Tz[iP] + (dzt / dzb)*(Tz[iP] - Tz[iB]); // <--
				}

				if (bB) {



					Tx[iB] = Tx[iP] + (dzb / dzt)*(Tx[iP] - Tx[iT]);
					Ty[iB] = Ty[iP] + (dzb / dzt)*(Ty[iP] - Ty[iT]);
					Tz[iB] = Tz[iP] + (dzb / dzt)*(Tz[iP] - Tz[iT]); // <--
				}
			}
		}
#endif
	}


} // green_gaussTemperature



// Минимаксное сглаживание градиента поправки давления.
// по моему это весьма неудачная идея.
void green_gaussPAMminmax(doublereal** &potent, int*** &neighbors_for_the_internal_node, integer maxelm, integer maxbound) {
	

	doublereal** minmaxlimitergrad=new doublereal*[3];
	for (integer i=0; i<3; i++) {
		minmaxlimitergrad[i]=new doublereal[maxelm+maxbound];
	}

	for (integer iP=0; iP<maxelm; iP++) {
	   // iP - номер внутреннего контрольного объёма
	   // iP изменяется от 0 до maxelm-1.
	   integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	   iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP]; iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	   // минимаксное ограничение против появления ложных максимумов.
	   minmaxlimitergrad[VELOCITY_X_COMPONENT][iP]=fmax(fmin(fmax(potent[GRADXPAM][iE],potent[GRADXPAM][iW]),potent[GRADXPAM][iP]),fmin(potent[GRADXPAM][iE],potent[GRADXPAM][iW]));
	   minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP]=fmax(fmin(fmax(potent[GRADYPAM][iN],potent[GRADYPAM][iS]),potent[GRADYPAM][iP]),fmin(potent[GRADYPAM][iN],potent[GRADYPAM][iS]));
	   minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP]=fmax(fmin(fmax(potent[GRADZPAM][iT],potent[GRADZPAM][iB]),potent[GRADZPAM][iP]),fmin(potent[GRADZPAM][iT],potent[GRADZPAM][iB]));
	}

	for (integer iP=0; iP<maxelm; iP++) {

		// iP - номер внутреннего контрольного объёма
	    // iP изменяется от 0 до maxelm-1.
	    integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

		// Если с одной из сторон стоит граница расчётной области
	    // то соответствующая переменная равна true
	    bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	    if (iE>=maxelm) bE=true;
	    if (iN>=maxelm) bN=true;
	    if (iT>=maxelm) bT=true;
        if (iW>=maxelm) bW=true;
	    if (iS>=maxelm) bS=true;
	    if (iB>=maxelm) bB=true;

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются простым снесением из ближайшего внутреннего узла.

		if (bE) {
			minmaxlimitergrad[VELOCITY_X_COMPONENT][iE]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Y_COMPONENT][iE]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Z_COMPONENT][iE]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
		}

		if (bW) {
			minmaxlimitergrad[VELOCITY_X_COMPONENT][iW]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Y_COMPONENT][iW]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Z_COMPONENT][iW]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
		}

		if (bN) {
			minmaxlimitergrad[VELOCITY_X_COMPONENT][iN]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Y_COMPONENT][iN]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Z_COMPONENT][iN]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
		}

		if (bS) {
			minmaxlimitergrad[VELOCITY_X_COMPONENT][iS]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Y_COMPONENT][iS]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Z_COMPONENT][iS]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
		}

		if (bT) {
			minmaxlimitergrad[VELOCITY_X_COMPONENT][iT]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Y_COMPONENT][iT]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Z_COMPONENT][iT]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
		}

		if (bB) {
			minmaxlimitergrad[VELOCITY_X_COMPONENT][iB]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Y_COMPONENT][iB]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
			minmaxlimitergrad[VELOCITY_Z_COMPONENT][iB]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
		}
	}

	// обратное копирование.
	for (integer iP=0; iP<maxelm+maxbound; iP++) {
		potent[GRADXPAM][iP]=minmaxlimitergrad[VELOCITY_X_COMPONENT][iP];
		potent[GRADYPAM][iP]=minmaxlimitergrad[VELOCITY_Y_COMPONENT][iP];
		potent[GRADZPAM][iP]=minmaxlimitergrad[VELOCITY_Z_COMPONENT][iP];
	}

	// Освобождение памяти.
	if (minmaxlimitergrad != NULL) {
		for (integer i = 0; i < 3; i++) {
			if (minmaxlimitergrad[i] != NULL) {
				delete[] minmaxlimitergrad[i];
				minmaxlimitergrad[i] = NULL;
			}
		}
		delete[] minmaxlimitergrad;
		minmaxlimitergrad = NULL;
	}

} // gaussPAMminmax

// вычисление градиентов давления с помощью теоремы Грина-Гаусса. 
// begin 21 июня 2012 года.
void green_gaussPRESS(integer iP, doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond,
					BOUND* &border_neighbor, integer ls, integer lw, WALL* &w, bool bLRfree,
	integer *ilevel_alice, int* ptr, TOCHKA* & volume_loc) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты давления для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {

		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	//volume3D(iP, nvtx, pa, dx, dy, dz);
	//dx = fabs(dx);
	//dy = fabs(dy);
	//dz = fabs(dz);

	TOCHKA point_loc = volume_loc[iP];
	dx = point_loc.x;
	dy = point_loc.y;
	dz = point_loc.z;

	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5*(pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5*(pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5*(pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5*(pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5*(pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5*(pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5*dx, dxw2 = 0.5*dx, dyn2 = 0.5*dy, dys2 = 0.5*dy, dzt2 = 0.5*dz, dzb2 = 0.5*dz;
	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5*(pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5*(pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5*(pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5*(pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5*(pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5*(pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5*(pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5*(pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5*(pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5*(pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5*(pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5*(pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5*(pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5*(pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5*(pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5*(pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5*(pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5*(pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);


	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus=0.5*dx/dxe;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;

	doublereal feplus2, fwplus2, fnplus2, fsplus2, ftplus2, fbplus2;
	// x-direction
	feplus2 = 0.5*dx / dxe2;
	fwplus2 = 0.5*dx / dxw2;
	// y-direction
	fnplus2 = 0.5*dy / dyn2;
	fsplus2 = 0.5*dy / dys2;
	// z-direction
	ftplus2 = 0.5*dz / dzt2;
	fbplus2 = 0.5*dz / dzb2;

	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	// x-direction
	feplus3 = 0.5*dx / dxe3;
	fwplus3 = 0.5*dx / dxw3;
	// y-direction
	fnplus3 = 0.5*dy / dyn3;
	fsplus3 = 0.5*dy / dys3;
	// z-direction
	ftplus3 = 0.5*dz / dzt3;
	fbplus3 = 0.5*dz / dzb3;

	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;
	// x-direction
	feplus4 = 0.5*dx / dxe4;
	fwplus4 = 0.5*dx / dxw4;
	// y-direction
	fnplus4 = 0.5*dy / dyn4;
	fsplus4 = 0.5*dy / dys4;
	// z-direction
	ftplus4 = 0.5*dz / dzt4;
	fbplus4 = 0.5*dz / dzb4;


	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.

	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.

	if (b_on_adaptive_local_refinement_mesh) {

	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW2]]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB2];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	

	

	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iE3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iW3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iN3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iS3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iT3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				//volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				point_loc = volume_loc[iB3];
				dx_loc = point_loc.x;
				dy_loc = point_loc.y;
				dz_loc = point_loc.z;

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	

	

		if (iE4 > -1) {

			dSqe4 = dy * dz;

			if (bE4) {
				// граничный узел.
				dSqe4 = border_neighbor[iE4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
					dSqe4 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iE4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqe4 = dy_loc * dz_loc;
				}
			}


		}


		if (iW4 > -1) {

			dSqw4 = dy * dz;

			if (bW4) {
				// граничный узел.
				dSqw4 = border_neighbor[iW4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
					dSqw4 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iW4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqw4 = dy_loc * dz_loc;
				}
			}


		}


		if (iN4 > -1) {

			dSqn4 = dx * dz;

			if (bN4) {
				// граничный узел.
				dSqn4 = border_neighbor[iN4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
					dSqn4 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iN4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqn4 = dx_loc * dz_loc;
				}
			}


		}


		if (iS4 > -1) {

			dSqs4 = dx * dz;

			if (bS4) {
				// граничный узел.
				dSqs4 = border_neighbor[iS4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
					dSqs4 = dx * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iS4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqs4 = dx_loc * dz_loc;
				}
			}


		}


		if (iT4 > -1) {

			dSqt4 = dx * dy;

			if (bT4) {
				// граничный узел.
				dSqt4 = border_neighbor[iT4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
					dSqt4 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iT4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqt4 = dx_loc * dy_loc;
				}
			}


		}


		if (iB4 > -1) {

			dSqb4 = dx * dy;

			if (bB4) {
				// граничный узел.
				dSqb4 = border_neighbor[iB4 - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
					dSqb4 = dx * dy;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					//volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

					point_loc = volume_loc[iB4];
					dx_loc = point_loc.x;
					dy_loc = point_loc.y;
					dz_loc = point_loc.z;

					dSqb4 = dx_loc * dy_loc;
				}
			}


		}
	}
	
	// 28.04.2019
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy*dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy*dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx*dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx*dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx*dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx*dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}
	

	doublereal PRESSe=0.0, PRESSw=0.0, PRESSn=0.0, PRESSs=0.0, PRESSt=0.0, PRESSb=0.0;
	doublereal PRESSe2 = 0.0, PRESSw2 = 0.0, PRESSn2 = 0.0, PRESSs2 = 0.0, PRESSt2 = 0.0, PRESSb2 = 0.0;
	doublereal PRESSe3 = 0.0, PRESSw3 = 0.0, PRESSn3 = 0.0, PRESSs3 = 0.0, PRESSt3 = 0.0, PRESSb3 = 0.0;
	doublereal PRESSe4 = 0.0, PRESSw4 = 0.0, PRESSn4 = 0.0, PRESSs4 = 0.0, PRESSt4 = 0.0, PRESSb4 = 0.0;

    if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		if (bLRfree&&((!b_on_adaptive_local_refinement_mesh))) {

			// не работает на АЛИС.
			//if (b_on_adaptive_local_refinement_mesh) {
				//printf("function green_gaussPRESS in module greengauss.c bLRfree not worked in ALICE mesh...\n ");
				//getchar();
				//exit(1);
			//}

			if (1) {
				// В случае bLRFree мы будем линейно интерполировать давление в граничные узлы чтобы сохранить значение градиента.
				// градиент давления очень важен для естественно конвективных течний, т.к. они протекают под действием этого градиента
				// и это для них особенно важно вблизи твёрдых стенок.
				/*
				if (!bE) PRESSe=feplus*potent[PRESS][iE]+(1.0-feplus)*potent[PRESS][iP]; else {
					//PRESSe=potent[PRESS][iE];
					//potent[GRADXPRESS][iE]=potent[GRADXPRESS][iP]+(dxe/dxw)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iW]);
					// линейная интерполяция давления на граничный узел !!!
					PRESSe=potent[PRESS][iP]+(dxe/dxw)*(potent[PRESS][iP]-potent[PRESS][iW]);
				}
				 if (!bW) PRESSw=fwplus*potent[PRESS][iW]+(1.0-fwplus)*potent[PRESS][iP]; else {
					 //PRESSw=potent[PRESS][iW];
					 //potent[GRADXPRESS][iW]=potent[GRADXPRESS][iP]+(dxw/dxe)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iE]);
					 // линейная интерполяция давления на граничный узел !!!
					 PRESSw=potent[PRESS][iP]+(dxw/dxe)*(potent[PRESS][iP]-potent[PRESS][iE]);
				 }
				 if (!bN) PRESSn=fnplus*potent[PRESS][iN]+(1.0-fnplus)*potent[PRESS][iP]; else {
					 //PRESSn=potent[PRESS][iN];
					 //potent[GRADYPRESS][iN]=potent[GRADYPRESS][iP]+(dyn/dys)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iS]);
					 // линейная интерполяция давления на граничный узел !!!
					 PRESSn=potent[PRESS][iP]+(dyn/dys)*(potent[PRESS][iP]-potent[PRESS][iS]);
				 }
				 if (!bS) PRESSs=fsplus*potent[PRESS][iS]+(1.0-fsplus)*potent[PRESS][iP]; else {
					 //PRESSs=potent[PRESS][iS];
					 //potent[GRADYPRESS][iS]=potent[GRADYPRESS][iP]+(dys/dyn)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iN]);
					 // линейная интерполяция давления на граничный узел !!!
					 PRESSs=potent[PRESS][iP]+(dys/dyn)*(potent[PRESS][iP]-potent[PRESS][iN]);
				 }
				 if (!bT) PRESSt=ftplus*potent[PRESS][iT]+(1.0-ftplus)*potent[PRESS][iP]; else {
					 PRESSt=potent[PRESS][iT];
					 //potent[GRADZPRESS][iT]=potent[GRADZPRESS][iP]+(dzt/dzb)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iB]);
					 // линейная интерполяция давления на граничный узел !!!
					 PRESSt=potent[PRESS][iP]+(dzt/dzb)*(potent[PRESS][iP]-potent[PRESS][iB]);
				 }
				 if (!bB) PRESSb=fbplus*potent[PRESS][iB]+(1.0-fbplus)*potent[PRESS][iP]; else {
					 PRESSb=potent[PRESS][iB];
					 //potent[GRADZPRESS][iB]=potent[GRADZPRESS][iP]+(dzb/dzt)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iT]);
					 // линейная интерполяция давления на граничный узел !!!
					 PRESSb=potent[PRESS][iP]+(dzb/dzt)*(potent[PRESS][iP]-potent[PRESS][iT]);
				 }
				 */



				if (iE > -1) {
					if (!bE) PRESSe = feplus * potent[PRESS][iE] + (1.0 - feplus)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iW > -1) {
							PRESSe = potent[PRESS][iP] + (dxe / dxw)*(potent[PRESS][iP] - potent[PRESS][iW]);
						}
						else if (iW2 > -1) {
							PRESSe = potent[PRESS][iP] + (dxe / dxw2)*(potent[PRESS][iP] - potent[PRESS][iW2]);
						}
						else if (iW3 > -1) {
							PRESSe = potent[PRESS][iP] + (dxe / dxw3)*(potent[PRESS][iP] - potent[PRESS][iW3]);
						}
						else if (iW4 > -1) {
							PRESSe = potent[PRESS][iP] + (dxe / dxw4)*(potent[PRESS][iP] - potent[PRESS][iW4]);
						}
					}
				}

				if (iE2 > -1) {
					if (!bE2) PRESSe2 = feplus2 * potent[PRESS][iE2] + (1.0 - feplus2)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iW > -1) {
							PRESSe2 = potent[PRESS][iP] + (dxe2 / dxw)*(potent[PRESS][iP] - potent[PRESS][iW]);
						}
						else if (iW2 > -1) {
							PRESSe2 = potent[PRESS][iP] + (dxe2 / dxw2)*(potent[PRESS][iP] - potent[PRESS][iW2]);
						}
						else if (iW3 > -1) {
							PRESSe2 = potent[PRESS][iP] + (dxe2 / dxw3)*(potent[PRESS][iP] - potent[PRESS][iW3]);
						}
						else if (iW4 > -1) {
							PRESSe2 = potent[PRESS][iP] + (dxe2 / dxw4)*(potent[PRESS][iP] - potent[PRESS][iW4]);
						}
					}
				}

				if (iE3 > -1) {
					if (!bE3) PRESSe3 = feplus3 * potent[PRESS][iE3] + (1.0 - feplus3)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iW > -1) {
							PRESSe3 = potent[PRESS][iP] + (dxe3 / dxw)*(potent[PRESS][iP] - potent[PRESS][iW]);
						}
						else if (iW2 > -1) {
							PRESSe3 = potent[PRESS][iP] + (dxe3 / dxw2)*(potent[PRESS][iP] - potent[PRESS][iW2]);
						}
						else if (iW3 > -1) {
							PRESSe3 = potent[PRESS][iP] + (dxe3 / dxw3)*(potent[PRESS][iP] - potent[PRESS][iW3]);
						}
						else if (iW4 > -1) {
							PRESSe3 = potent[PRESS][iP] + (dxe3 / dxw4)*(potent[PRESS][iP] - potent[PRESS][iW4]);
						}
					}
				}

				if (iE4 > -1) {
					if (!bE4) PRESSe4 = feplus4 * potent[PRESS][iE4] + (1.0 - feplus4)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iW > -1) {
							PRESSe4 = potent[PRESS][iP] + (dxe4 / dxw)*(potent[PRESS][iP] - potent[PRESS][iW]);
						}
						else if (iW2 > -1) {
							PRESSe4 = potent[PRESS][iP] + (dxe4 / dxw2)*(potent[PRESS][iP] - potent[PRESS][iW2]);
						}
						else if (iW3 > -1) {
							PRESSe4 = potent[PRESS][iP] + (dxe4 / dxw3)*(potent[PRESS][iP] - potent[PRESS][iW3]);
						}
						else if (iW4 > -1) {
							PRESSe4 = potent[PRESS][iP] + (dxe4 / dxw4)*(potent[PRESS][iP] - potent[PRESS][iW4]);
						}
					}
				}

				if (iW > -1) {
					if (!bW) PRESSw = fwplus * potent[PRESS][iW] + (1.0 - fwplus)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iE > -1) {
							PRESSw = potent[PRESS][iP] + (dxw / dxe)*(potent[PRESS][iP] - potent[PRESS][iE]);
						}
						else if (iE2 > -1) {
							PRESSw = potent[PRESS][iP] + (dxw / dxe2)*(potent[PRESS][iP] - potent[PRESS][iE2]);
						}
						else if (iE3 > -1) {
							PRESSw = potent[PRESS][iP] + (dxw / dxe3)*(potent[PRESS][iP] - potent[PRESS][iE3]);
						}
						else if (iE4 > -1) {
							PRESSw = potent[PRESS][iP] + (dxw / dxe4)*(potent[PRESS][iP] - potent[PRESS][iE4]);
						}

					}
				}

				if (iW2 > -1) {
					if (!bW2) PRESSw2 = fwplus2 * potent[PRESS][iW2] + (1.0 - fwplus2)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iE > -1) {
							PRESSw2 = potent[PRESS][iP] + (dxw2 / dxe)*(potent[PRESS][iP] - potent[PRESS][iE]);
						}
						else if (iE2 > -1) {
							PRESSw2 = potent[PRESS][iP] + (dxw2 / dxe2)*(potent[PRESS][iP] - potent[PRESS][iE2]);
						}
						else if (iE3 > -1) {
							PRESSw2 = potent[PRESS][iP] + (dxw2 / dxe3)*(potent[PRESS][iP] - potent[PRESS][iE3]);
						}
						else if (iE4 > -1) {
							PRESSw2 = potent[PRESS][iP] + (dxw2 / dxe4)*(potent[PRESS][iP] - potent[PRESS][iE4]);
						}

					}
				}

				if (iW3 > -1) {
					if (!bW3) PRESSw3 = fwplus3 * potent[PRESS][iW3] + (1.0 - fwplus3)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iE > -1) {
							PRESSw3 = potent[PRESS][iP] + (dxw3 / dxe)*(potent[PRESS][iP] - potent[PRESS][iE]);
						}
						else if (iE2 > -1) {
							PRESSw3 = potent[PRESS][iP] + (dxw3 / dxe2)*(potent[PRESS][iP] - potent[PRESS][iE2]);
						}
						else if (iE3 > -1) {
							PRESSw3 = potent[PRESS][iP] + (dxw3 / dxe3)*(potent[PRESS][iP] - potent[PRESS][iE3]);
						}
						else if (iE4 > -1) {
							PRESSw3 = potent[PRESS][iP] + (dxw3 / dxe4)*(potent[PRESS][iP] - potent[PRESS][iE4]);
						}

					}
				}

				if (iW4 > -1) {
					if (!bW4) PRESSw4 = fwplus4 * potent[PRESS][iW4] + (1.0 - fwplus4)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iE > -1) {
							PRESSw4 = potent[PRESS][iP] + (dxw4 / dxe)*(potent[PRESS][iP] - potent[PRESS][iE]);
						}
						else if (iE2 > -1) {
							PRESSw4 = potent[PRESS][iP] + (dxw4 / dxe2)*(potent[PRESS][iP] - potent[PRESS][iE2]);
						}
						else if (iE3 > -1) {
							PRESSw4 = potent[PRESS][iP] + (dxw4 / dxe3)*(potent[PRESS][iP] - potent[PRESS][iE3]);
						}
						else if (iE4 > -1) {
							PRESSw4 = potent[PRESS][iP] + (dxw4 / dxe4)*(potent[PRESS][iP] - potent[PRESS][iE4]);
						}

					}
				}

				if (iN > -1) {
					if (!bN) PRESSn = fnplus * potent[PRESS][iN] + (1.0 - fnplus)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iS > -1) {
							PRESSn = potent[PRESS][iP] + (dyn / dys)*(potent[PRESS][iP] - potent[PRESS][iS]);
						}
						else if (iS2 > -1) {
							PRESSn = potent[PRESS][iP] + (dyn / dys2)*(potent[PRESS][iP] - potent[PRESS][iS2]);
						}
						else if (iS3 > -1) {
							PRESSn = potent[PRESS][iP] + (dyn / dys3)*(potent[PRESS][iP] - potent[PRESS][iS3]);
						}
						else if (iS4 > -1) {
							PRESSn = potent[PRESS][iP] + (dyn / dys4)*(potent[PRESS][iP] - potent[PRESS][iS4]);
						}

					}
				}

				if (iN2 > -1) {
					if (!bN2) PRESSn2 = fnplus2 * potent[PRESS][iN2] + (1.0 - fnplus2)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iS > -1) {
							PRESSn2 = potent[PRESS][iP] + (dyn2 / dys)*(potent[PRESS][iP] - potent[PRESS][iS]);
						}
						else if (iS2 > -1) {
							PRESSn2 = potent[PRESS][iP] + (dyn2 / dys2)*(potent[PRESS][iP] - potent[PRESS][iS2]);
						}
						else if (iS3 > -1) {
							PRESSn2 = potent[PRESS][iP] + (dyn2 / dys3)*(potent[PRESS][iP] - potent[PRESS][iS3]);
						}
						else if (iS4 > -1) {
							PRESSn2 = potent[PRESS][iP] + (dyn2 / dys4)*(potent[PRESS][iP] - potent[PRESS][iS4]);
						}

					}
				}

				if (iN3 > -1) {
					if (!bN3) PRESSn3 = fnplus3 * potent[PRESS][iN3] + (1.0 - fnplus3)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iS > -1) {
							PRESSn3 = potent[PRESS][iP] + (dyn3 / dys)*(potent[PRESS][iP] - potent[PRESS][iS]);
						}
						else if (iS2 > -1) {
							PRESSn3 = potent[PRESS][iP] + (dyn3 / dys2)*(potent[PRESS][iP] - potent[PRESS][iS2]);
						}
						else if (iS3 > -1) {
							PRESSn3 = potent[PRESS][iP] + (dyn3 / dys3)*(potent[PRESS][iP] - potent[PRESS][iS3]);
						}
						else if (iS4 > -1) {
							PRESSn3 = potent[PRESS][iP] + (dyn3 / dys4)*(potent[PRESS][iP] - potent[PRESS][iS4]);
						}

					}
				}

				if (iN4 > -1) {
					if (!bN4) PRESSn4 = fnplus4 * potent[PRESS][iN4] + (1.0 - fnplus4)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iS > -1) {
							PRESSn4 = potent[PRESS][iP] + (dyn4 / dys)*(potent[PRESS][iP] - potent[PRESS][iS]);
						}
						else if (iS2 > -1) {
							PRESSn4 = potent[PRESS][iP] + (dyn4 / dys2)*(potent[PRESS][iP] - potent[PRESS][iS2]);
						}
						else if (iS3 > -1) {
							PRESSn4 = potent[PRESS][iP] + (dyn4 / dys3)*(potent[PRESS][iP] - potent[PRESS][iS3]);
						}
						else if (iS4 > -1) {
							PRESSn4 = potent[PRESS][iP] + (dyn4 / dys4)*(potent[PRESS][iP] - potent[PRESS][iS4]);
						}

					}
				}

				if (iS > -1) {
					if (!bS) PRESSs = fsplus * potent[PRESS][iS] + (1.0 - fsplus)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iN > -1) {
							PRESSs = potent[PRESS][iP] + (dys / dyn)*(potent[PRESS][iP] - potent[PRESS][iN]);
						}
						else if (iN2 > -1) {
							PRESSs = potent[PRESS][iP] + (dys / dyn2)*(potent[PRESS][iP] - potent[PRESS][iN2]);
						}
						else if (iN3 > -1) {
							PRESSs = potent[PRESS][iP] + (dys / dyn3)*(potent[PRESS][iP] - potent[PRESS][iN3]);
						}
						else if (iN4 > -1) {
							PRESSs = potent[PRESS][iP] + (dys / dyn4)*(potent[PRESS][iP] - potent[PRESS][iN4]);
						}

					}
				}

				if (iS2 > -1) {
					if (!bS2) PRESSs2 = fsplus2 * potent[PRESS][iS2] + (1.0 - fsplus2)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iN > -1) {
							PRESSs2 = potent[PRESS][iP] + (dys2 / dyn)*(potent[PRESS][iP] - potent[PRESS][iN]);
						}
						else if (iN2 > -1) {
							PRESSs2 = potent[PRESS][iP] + (dys2 / dyn2)*(potent[PRESS][iP] - potent[PRESS][iN2]);
						}
						else if (iN3 > -1) {
							PRESSs2 = potent[PRESS][iP] + (dys2 / dyn3)*(potent[PRESS][iP] - potent[PRESS][iN3]);
						}
						else if (iN4 > -1) {
							PRESSs2 = potent[PRESS][iP] + (dys2 / dyn4)*(potent[PRESS][iP] - potent[PRESS][iN4]);
						}

					}
				}

				if (iS3 > -1) {
					if (!bS3) PRESSs3 = fsplus3 * potent[PRESS][iS3] + (1.0 - fsplus3)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iN > -1) {
							PRESSs3 = potent[PRESS][iP] + (dys3 / dyn)*(potent[PRESS][iP] - potent[PRESS][iN]);
						}
						else if (iN2 > -1) {
							PRESSs3 = potent[PRESS][iP] + (dys3 / dyn2)*(potent[PRESS][iP] - potent[PRESS][iN2]);
						}
						else if (iN3 > -1) {
							PRESSs3 = potent[PRESS][iP] + (dys3 / dyn3)*(potent[PRESS][iP] - potent[PRESS][iN3]);
						}
						else if (iN4 > -1) {
							PRESSs3 = potent[PRESS][iP] + (dys3 / dyn4)*(potent[PRESS][iP] - potent[PRESS][iN4]);
						}

					}
				}

				if (iS4 > -1) {
					if (!bS4) PRESSs4 = fsplus4 * potent[PRESS][iS4] + (1.0 - fsplus4)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iN > -1) {
							PRESSs4 = potent[PRESS][iP] + (dys4 / dyn)*(potent[PRESS][iP] - potent[PRESS][iN]);
						}
						else if (iN2 > -1) {
							PRESSs4 = potent[PRESS][iP] + (dys4 / dyn2)*(potent[PRESS][iP] - potent[PRESS][iN2]);
						}
						else if (iN3 > -1) {
							PRESSs4 = potent[PRESS][iP] + (dys4 / dyn3)*(potent[PRESS][iP] - potent[PRESS][iN3]);
						}
						else if (iN4 > -1) {
							PRESSs4 = potent[PRESS][iP] + (dys4 / dyn4)*(potent[PRESS][iP] - potent[PRESS][iN4]);
						}

					}
				}

				if (iT > -1) {
					if (!bT) PRESSt = ftplus * potent[PRESS][iT] + (1.0 - ftplus)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iB > -1) {
							PRESSt = potent[PRESS][iP] + (dzt / dzb)*(potent[PRESS][iP] - potent[PRESS][iB]);
						}
						else if (iB2 > -1) {
							PRESSt = potent[PRESS][iP] + (dzt / dzb2)*(potent[PRESS][iP] - potent[PRESS][iB2]);
						}
						else if (iB3 > -1) {
							PRESSt = potent[PRESS][iP] + (dzt / dzb3)*(potent[PRESS][iP] - potent[PRESS][iB3]);
						}
						else if (iB4 > -1) {
							PRESSt = potent[PRESS][iP] + (dzt / dzb4)*(potent[PRESS][iP] - potent[PRESS][iB4]);
						}

					}
				}


				if (iT2 > -1) {
					if (!bT2) PRESSt2 = ftplus2 * potent[PRESS][iT2] + (1.0 - ftplus2)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iB > -1) {
							PRESSt2 = potent[PRESS][iP] + (dzt2 / dzb)*(potent[PRESS][iP] - potent[PRESS][iB]);
						}
						else if (iB2 > -1) {
							PRESSt2 = potent[PRESS][iP] + (dzt2 / dzb2)*(potent[PRESS][iP] - potent[PRESS][iB2]);
						}
						else if (iB3 > -1) {
							PRESSt2 = potent[PRESS][iP] + (dzt2 / dzb3)*(potent[PRESS][iP] - potent[PRESS][iB3]);
						}
						else if (iB4 > -1) {
							PRESSt2 = potent[PRESS][iP] + (dzt2 / dzb4)*(potent[PRESS][iP] - potent[PRESS][iB4]);
						}

					}
				}

				if (iT3 > -1) {
					if (!bT3) PRESSt3 = ftplus3 * potent[PRESS][iT3] + (1.0 - ftplus3)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iB > -1) {
							PRESSt3 = potent[PRESS][iP] + (dzt3 / dzb)*(potent[PRESS][iP] - potent[PRESS][iB]);
						}
						else if (iB2 > -1) {
							PRESSt3 = potent[PRESS][iP] + (dzt3 / dzb2)*(potent[PRESS][iP] - potent[PRESS][iB2]);
						}
						else if (iB3 > -1) {
							PRESSt3 = potent[PRESS][iP] + (dzt3 / dzb3)*(potent[PRESS][iP] - potent[PRESS][iB3]);
						}
						else if (iB4 > -1) {
							PRESSt3 = potent[PRESS][iP] + (dzt3 / dzb4)*(potent[PRESS][iP] - potent[PRESS][iB4]);
						}

					}
				}

				if (iT4 > -1) {
					if (!bT4) PRESSt4 = ftplus4 * potent[PRESS][iT4] + (1.0 - ftplus4)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iB > -1) {
							PRESSt4 = potent[PRESS][iP] + (dzt4 / dzb)*(potent[PRESS][iP] - potent[PRESS][iB]);
						}
						else if (iB2 > -1) {
							PRESSt4 = potent[PRESS][iP] + (dzt4 / dzb2)*(potent[PRESS][iP] - potent[PRESS][iB2]);
						}
						else if (iB3 > -1) {
							PRESSt4 = potent[PRESS][iP] + (dzt4 / dzb3)*(potent[PRESS][iP] - potent[PRESS][iB3]);
						}
						else if (iB4 > -1) {
							PRESSt4 = potent[PRESS][iP] + (dzt4 / dzb4)*(potent[PRESS][iP] - potent[PRESS][iB4]);
						}

					}
				}

				if (iB > -1) {
					if (!bB) PRESSb = fbplus * potent[PRESS][iB] + (1.0 - fbplus)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iT > -1) {
							PRESSb = potent[PRESS][iP] + (dzb / dzt)*(potent[PRESS][iP] - potent[PRESS][iT]);
						}
						else if (iT2 > -1) {
							PRESSb = potent[PRESS][iP] + (dzb / dzt2)*(potent[PRESS][iP] - potent[PRESS][iT2]);
						}
						else  if (iT3 > -1) {
							PRESSb = potent[PRESS][iP] + (dzb / dzt3)*(potent[PRESS][iP] - potent[PRESS][iT3]);
						}
						else if (iT4 > -1) {
							PRESSb = potent[PRESS][iP] + (dzb / dzt4)*(potent[PRESS][iP] - potent[PRESS][iT4]);
						}
					}
				}

				if (iB2 > -1) {
					if (!bB2) PRESSb2 = fbplus2 * potent[PRESS][iB2] + (1.0 - fbplus2)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iT > -1) {
							PRESSb2 = potent[PRESS][iP] + (dzb2 / dzt)*(potent[PRESS][iP] - potent[PRESS][iT]);
						}
						else if (iT2 > -1) {
							PRESSb2 = potent[PRESS][iP] + (dzb2 / dzt2)*(potent[PRESS][iP] - potent[PRESS][iT2]);
						}
						else  if (iT3 > -1) {
							PRESSb2 = potent[PRESS][iP] + (dzb2 / dzt3)*(potent[PRESS][iP] - potent[PRESS][iT3]);
						}
						else if (iT4 > -1) {
							PRESSb2 = potent[PRESS][iP] + (dzb2 / dzt4)*(potent[PRESS][iP] - potent[PRESS][iT4]);
						}
					}
				}

				if (iB3 > -1) {
					if (!bB3) PRESSb3 = fbplus3 * potent[PRESS][iB3] + (1.0 - fbplus3)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iT > -1) {
							PRESSb3 = potent[PRESS][iP] + (dzb3 / dzt)*(potent[PRESS][iP] - potent[PRESS][iT]);
						}
						else if (iT2 > -1) {
							PRESSb3 = potent[PRESS][iP] + (dzb3 / dzt2)*(potent[PRESS][iP] - potent[PRESS][iT2]);
						}
						else  if (iT3 > -1) {
							PRESSb3 = potent[PRESS][iP] + (dzb3 / dzt3)*(potent[PRESS][iP] - potent[PRESS][iT3]);
						}
						else if (iT4 > -1) {
							PRESSb3 = potent[PRESS][iP] + (dzb3 / dzt4)*(potent[PRESS][iP] - potent[PRESS][iT4]);
						}
					}
				}

				if (iB4 > -1) {
					if (!bB4) PRESSb4 = fbplus4 * potent[PRESS][iB4] + (1.0 - fbplus4)*potent[PRESS][iP]; else {
						// линейная интерполяция давления на граничный узел !!!
						if (iT > -1) {
							PRESSb4 = potent[PRESS][iP] + (dzb4 / dzt)*(potent[PRESS][iP] - potent[PRESS][iT]);
						}
						else if (iT2 > -1) {
							PRESSb4 = potent[PRESS][iP] + (dzb4 / dzt2)*(potent[PRESS][iP] - potent[PRESS][iT2]);
						}
						else  if (iT3 > -1) {
							PRESSb4 = potent[PRESS][iP] + (dzb4 / dzt3)*(potent[PRESS][iP] - potent[PRESS][iT3]);
						}
						else if (iT4 > -1) {
							PRESSb4 = potent[PRESS][iP] + (dzb4 / dzt4)*(potent[PRESS][iP] - potent[PRESS][iT4]);
						}
					}
				}
			}


			 // градиент Давления.
			 //potent[GRADXPRESS][iP]=(PRESSe-PRESSw)/dx;
			 //potent[GRADYPRESS][iP]=(PRESSn-PRESSs)/dy;
			 //potent[GRADZPRESS][iP]=(PRESSt-PRESSb)/dz;
			 potent[GRADXPRESS][iP] = (PRESSe*dSqe / (dy*dz) + PRESSe2 * dSqe2 / (dy*dz) + PRESSe3 * dSqe3 / (dy*dz) + PRESSe4 * dSqe4 / (dy*dz) - (PRESSw*dSqw / (dy*dz) + PRESSw2 * dSqw2 / (dy*dz) + PRESSw3 * dSqw3 / (dy*dz) + PRESSw4 * dSqw4 / (dy*dz))) / dx;
			 potent[GRADYPRESS][iP] = (PRESSn*dSqn / (dx*dz) + PRESSn2 * dSqn2 / (dx*dz) + PRESSn3 * dSqn3 / (dx*dz) + PRESSn4 * dSqn4 / (dx*dz) - (PRESSs*dSqs / (dx*dz) + PRESSs2 * dSqs2 / (dx*dz) + PRESSs3 * dSqs3 / (dx*dz) + PRESSs4 * dSqs4 / (dx*dz))) / dy;
			 potent[GRADZPRESS][iP] = (PRESSt*dSqt / (dx*dy) + PRESSt2 * dSqt2 / (dx*dy) + PRESSt3 * dSqt3 / (dx*dy) + PRESSt4 * dSqt4 / (dx*dy) - (PRESSb*dSqb / (dx*dy) + PRESSb2 * dSqb2 / (dx*dy) + PRESSb3 * dSqb3 / (dx*dy) + PRESSb4 * dSqb4 / (dx*dy))) / dz;



		}
		else {

			if (iE > -1) {
				if (!bE) PRESSe = feplus * potent[PRESS][iE] + (1.0 - feplus)*potent[PRESS][iP]; else PRESSe = potent[PRESS][iE]; // проверено !
			}
			if (iW > -1) {
				if (!bW) PRESSw = fwplus * potent[PRESS][iW] + (1.0 - fwplus)*potent[PRESS][iP]; else PRESSw = potent[PRESS][iW];
			}
			if (iN > -1) {
				if (!bN) PRESSn = fnplus * potent[PRESS][iN] + (1.0 - fnplus)*potent[PRESS][iP]; else PRESSn = potent[PRESS][iN];
			}
			if (iS > -1) {
				if (!bS) PRESSs = fsplus * potent[PRESS][iS] + (1.0 - fsplus)*potent[PRESS][iP]; else PRESSs = potent[PRESS][iS];
			}
			if (iT > -1) {
				if (!bT) PRESSt = ftplus * potent[PRESS][iT] + (1.0 - ftplus)*potent[PRESS][iP]; else PRESSt = potent[PRESS][iT];
			}
			if (iB > -1) {
				if (!bB) PRESSb = fbplus * potent[PRESS][iB] + (1.0 - fbplus)*potent[PRESS][iP]; else PRESSb = potent[PRESS][iB];
			}

			if (iE2 > -1) {
				if (!bE2) PRESSe2 = feplus2 * potent[PRESS][iE2] + (1.0 - feplus2)*potent[PRESS][iP]; else PRESSe2 = potent[PRESS][iE2]; // проверено !
			}
			if (iW2 > -1) {
				if (!bW2) PRESSw2 = fwplus2 * potent[PRESS][iW2] + (1.0 - fwplus2)*potent[PRESS][iP]; else PRESSw2 = potent[PRESS][iW2];
			}
			if (iN2 > -1) {
				if (!bN2) PRESSn2 = fnplus2 * potent[PRESS][iN2] + (1.0 - fnplus2)*potent[PRESS][iP]; else PRESSn2 = potent[PRESS][iN2];
			}
			if (iS2 > -1) {
				if (!bS2) PRESSs2 = fsplus2 * potent[PRESS][iS2] + (1.0 - fsplus2)*potent[PRESS][iP]; else PRESSs2 = potent[PRESS][iS2];
			}
			if (iT2 > -1) {
				if (!bT2) PRESSt2 = ftplus2 * potent[PRESS][iT2] + (1.0 - ftplus2)*potent[PRESS][iP]; else PRESSt2 = potent[PRESS][iT2];
			}
			if (iB2 > -1) {
				if (!bB2) PRESSb2 = fbplus2 * potent[PRESS][iB2] + (1.0 - fbplus2)*potent[PRESS][iP]; else PRESSb2 = potent[PRESS][iB2];
			}

			if (iE3 > -1) {
				if (!bE3) PRESSe3 = feplus3 * potent[PRESS][iE3] + (1.0 - feplus3)*potent[PRESS][iP]; else PRESSe3 = potent[PRESS][iE3]; // проверено !
			}
			if (iW3 > -1) {
				if (!bW3) PRESSw3 = fwplus3 * potent[PRESS][iW3] + (1.0 - fwplus3)*potent[PRESS][iP]; else PRESSw3 = potent[PRESS][iW3];
			}
			if (iN3 > -1) {
				if (!bN3) PRESSn3 = fnplus3 * potent[PRESS][iN3] + (1.0 - fnplus3)*potent[PRESS][iP]; else PRESSn3 = potent[PRESS][iN3];
			}
			if (iS3 > -1) {
				if (!bS3) PRESSs3 = fsplus3 * potent[PRESS][iS3] + (1.0 - fsplus3)*potent[PRESS][iP]; else PRESSs3 = potent[PRESS][iS3];
			}
			if (iT3 > -1) {
				if (!bT3) PRESSt3 = ftplus3 * potent[PRESS][iT3] + (1.0 - ftplus3)*potent[PRESS][iP]; else PRESSt3 = potent[PRESS][iT3];
			}
			if (iB3 > -1) {
				if (!bB3) PRESSb3 = fbplus3 * potent[PRESS][iB3] + (1.0 - fbplus3)*potent[PRESS][iP]; else PRESSb3 = potent[PRESS][iB3];
			}

			if (iE4 > -1) {
				if (!bE4) PRESSe4 = feplus4 * potent[PRESS][iE4] + (1.0 - feplus4)*potent[PRESS][iP]; else PRESSe4 = potent[PRESS][iE4]; // проверено !
			}
			if (iW4 > -1) {
				if (!bW4) PRESSw4 = fwplus4 * potent[PRESS][iW4] + (1.0 - fwplus4)*potent[PRESS][iP]; else PRESSw4 = potent[PRESS][iW4];
			}
			if (iN4 > -1) {
				if (!bN4) PRESSn4 = fnplus4 * potent[PRESS][iN4] + (1.0 - fnplus4)*potent[PRESS][iP]; else PRESSn4 = potent[PRESS][iN4];
			}
			if (iS4 > -1) {
				if (!bS4) PRESSs4 = fsplus4 * potent[PRESS][iS4] + (1.0 - fsplus4)*potent[PRESS][iP]; else PRESSs4 = potent[PRESS][iS4];
			}
			if (iT4 > -1) {
				if (!bT4) PRESSt4 = ftplus4 * potent[PRESS][iT4] + (1.0 - ftplus4)*potent[PRESS][iP]; else PRESSt4 = potent[PRESS][iT4];
			}
			if (iB4 > -1) {
				if (!bB4) PRESSb4 = fbplus4 * potent[PRESS][iB4] + (1.0 - fbplus4)*potent[PRESS][iP]; else PRESSb4 = potent[PRESS][iB4];
			}

             // градиент Давления.
	         //potent[GRADXPRESS][iP]=(PRESSe-PRESSw)/dx;
	         //potent[GRADYPRESS][iP]=(PRESSn-PRESSs)/dy;
	         //potent[GRADZPRESS][iP]=(PRESSt-PRESSb)/dz;
			 potent[GRADXPRESS][iP] = (PRESSe*dSqe / (dy*dz) + PRESSe2 * dSqe2 / (dy*dz) + PRESSe3 * dSqe3 / (dy*dz) + PRESSe4 * dSqe4 / (dy*dz) - (PRESSw*dSqw / (dy*dz) + PRESSw2 * dSqw2 / (dy*dz) + PRESSw3 * dSqw3 / (dy*dz) + PRESSw4 * dSqw4 / (dy*dz))) / dx;
			 potent[GRADYPRESS][iP] = (PRESSn*dSqn / (dx*dz) + PRESSn2 * dSqn2 / (dx*dz) + PRESSn3 * dSqn3 / (dx*dz) + PRESSn4 * dSqn4 / (dx*dz) - (PRESSs*dSqs / (dx*dz) + PRESSs2 * dSqs2 / (dx*dz) + PRESSs3 * dSqs3 / (dx*dz) + PRESSs4 * dSqs4 / (dx*dz))) / dy;
			 potent[GRADZPRESS][iP] = (PRESSt*dSqt / (dx*dy) + PRESSt2 * dSqt2 / (dx*dy) + PRESSt3 * dSqt3 / (dx*dy) + PRESSt4 * dSqt4 / (dx*dy) - (PRESSb*dSqb / (dx*dy) + PRESSb2 * dSqb2 / (dx*dy) + PRESSb3 * dSqb3 / (dx*dy) + PRESSb4 * dSqb4 / (dx*dy))) / dz;



		}

	}
   else {

	   const integer interpol=0; // 0 или 1 при линейной интерполяции.

#if (interpol==0)
	   {

		   if (1) {

			   if (iE > -1) {
				   if (bE) {
					   potent[GRADXPRESS][iE] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iE] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW > -1) {
				   if (bW) {
					   potent[GRADXPRESS][iW] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iW] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN > -1) {
				   if (bN) {
					   potent[GRADYPRESS][iN] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iN] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS > -1) {
				   if (bS) {
					   potent[GRADYPRESS][iS] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iS] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT > -1) {
				   if (bT) {
					   potent[GRADZPRESS][iT] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iT] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB > -1) {
				   if (bB) {
					   potent[GRADZPRESS][iB] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iB] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iE2 > -1) {
				   if (bE2) {
					   potent[GRADXPRESS][iE2] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iE2] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW2 > -1) {
				   if (bW2) {
					   potent[GRADXPRESS][iW2] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iW2] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN2 > -1) {
				   if (bN2) {
					   potent[GRADYPRESS][iN2] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iN2] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS2 > -1) {
				   if (bS2) {
					   potent[GRADYPRESS][iS2] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iS2] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT2 > -1) {
				   if (bT2) {
					   potent[GRADZPRESS][iT2] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iT2] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT2] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB2 > -1) {
				   if (bB2) {
					   potent[GRADZPRESS][iB2] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iB2] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB2] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iE3 > -1) {
				   if (bE3) {
					   potent[GRADXPRESS][iE3] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iE3] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW3 > -1) {
				   if (bW3) {
					   potent[GRADXPRESS][iW3] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iW3] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN3 > -1) {
				   if (bN3) {
					   potent[GRADYPRESS][iN3] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iN3] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS3 > -1) {
				   if (bS3) {
					   potent[GRADYPRESS][iS3] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iS3] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT3 > -1) {
				   if (bT3) {
					   potent[GRADZPRESS][iT3] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iT3] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT3] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB3 > -1) {
				   if (bB3) {
					   potent[GRADZPRESS][iB3] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iB3] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB3] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iE4 > -1) {
				   if (bE4) {
					   potent[GRADXPRESS][iE4] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iE4] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW4 > -1) {
				   if (bW4) {
					   potent[GRADXPRESS][iW4] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iW4] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN4 > -1) {
				   if (bN4) {
					   potent[GRADYPRESS][iN4] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iN4] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS4 > -1) {
				   if (bS4) {
					   potent[GRADYPRESS][iS4] = potent[GRADYPRESS][iP];
					   potent[GRADXPRESS][iS4] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT4 > -1) {
				   if (bT4) {
					   potent[GRADZPRESS][iT4] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iT4] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT4] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB4 > -1) {
				   if (bB4) {
					   potent[GRADZPRESS][iB4] = potent[GRADZPRESS][iP];
					   potent[GRADXPRESS][iB4] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB4] = potent[GRADYPRESS][iP];
				   }
			   }


		   }
		   else {

			   if (iE > -1) {
				   if (bE) {

					   integer inumber = iE - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iE] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iE] = 0.0;
					   }

					   potent[GRADYPRESS][iE] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW > -1) {
				   if (bW) {

					   integer inumber = iW - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iW] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iW] = 0.0;
					   }

					   potent[GRADYPRESS][iW] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN > -1) {
				   if (bN) {

					   integer inumber = iN - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iN] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iN] = 0.0;
					   }

					   potent[GRADXPRESS][iN] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS > -1) {
				   if (bS) {

					   integer inumber = iS - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iS] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iS] = 0.0;
					   }

					   potent[GRADXPRESS][iS] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT > -1) {
				   if (bT) {

					   integer inumber = iT - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iT] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iT] = 0.0;
					   }

					   potent[GRADXPRESS][iT] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB > -1) {
				   if (bB) {

					   integer inumber = iB - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iB] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iB] = 0.0;
					   }

					   potent[GRADXPRESS][iB] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iE2 > -1) {
				   if (bE2) {

					   integer inumber = iE2 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE2] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iE2] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iE2] = 0.0;
					   }

					   potent[GRADYPRESS][iE2] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW2 > -1) {
				   if (bW2) {

					   integer inumber = iW2 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW2] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iW2] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iW2] = 0.0;
					   }

					   potent[GRADYPRESS][iW2] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN2 > -1) {
				   if (bN2) {

					   integer inumber = iN2 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN2] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iN2] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iN2] = 0.0;
					   }

					   potent[GRADXPRESS][iN2] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS2 > -1) {
				   if (bS2) {

					   integer inumber = iS2 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS2] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iS2] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iS2] = 0.0;
					   }

					   potent[GRADXPRESS][iS2] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS2] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT2 > -1) {
				   if (bT2) {

					   integer inumber = iT2 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT2] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iT2] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iT2] = 0.0;
					   }

					   potent[GRADXPRESS][iT2] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT2] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB2 > -1) {
				   if (bB2) {

					   integer inumber = iB2 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB2] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iB2] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iB2] = 0.0;
					   }

					   potent[GRADXPRESS][iB2] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB2] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iE3 > -1) {
				   if (bE3) {

					   integer inumber = iE3 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE3] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iE3] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iE3] = 0.0;
					   }

					   potent[GRADYPRESS][iE3] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW3 > -1) {
				   if (bW3) {

					   integer inumber = iW3 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW3] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iW3] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iW3] = 0.0;
					   }

					   potent[GRADYPRESS][iW3] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN3 > -1) {
				   if (bN3) {

					   integer inumber = iN3 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN3] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iN3] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iN3] = 0.0;
					   }

					   potent[GRADXPRESS][iN3] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS3 > -1) {
				   if (bS3) {

					   integer inumber = iS3 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS3] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iS3] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iS3] = 0.0;
					   }

					   potent[GRADXPRESS][iS3] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS3] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT3 > -1) {
				   if (bT3) {

					   integer inumber = iT3 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT3] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iT3] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iT3] = 0.0;
					   }

					   potent[GRADXPRESS][iT3] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT3] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB3 > -1) {
				   if (bB3) {

					   integer inumber = iB3 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB3] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iB3] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iB3] = 0.0;
					   }

					   potent[GRADXPRESS][iB3] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB3] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iE4 > -1) {
				   if (bE4) {

					   integer inumber = iE4 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE4] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iE4] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iE4] = 0.0;
					   }

					   potent[GRADYPRESS][iE4] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iE4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iW4 > -1) {
				   if (bW4) {

					   integer inumber = iW4 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW4] = potent[GRADXPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADXPRESS][iW4] = potent[GRADXPRESS][iP];
					   }
					   else {
						   potent[GRADXPRESS][iW4] = 0.0;
					   }

					   potent[GRADYPRESS][iW4] = potent[GRADYPRESS][iP];
					   potent[GRADZPRESS][iW4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iN4 > -1) {
				   if (bN4) {

					   integer inumber = iN4 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN4] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iN4] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iN4] = 0.0;
					   }

					   potent[GRADXPRESS][iN4] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iN4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iS4 > -1) {
				   if (bS4) {

					   integer inumber = iS4 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS4] = potent[GRADYPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADYPRESS][iS4] = potent[GRADYPRESS][iP];
					   }
					   else {
						   potent[GRADYPRESS][iS4] = 0.0;
					   }

					   potent[GRADXPRESS][iS4] = potent[GRADXPRESS][iP];
					   potent[GRADZPRESS][iS4] = potent[GRADZPRESS][iP];
				   }
			   }

			   if (iT4 > -1) {
				   if (bT4) {

					   integer inumber = iT4 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT4] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iT4] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iT4] = 0.0;
					   }

					   potent[GRADXPRESS][iT4] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iT4] = potent[GRADYPRESS][iP];
				   }
			   }

			   if (iB4 > -1) {
				   if (bB4) {

					   integer inumber = iB4 - maxelm;
					   if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB4] = potent[GRADZPRESS][iP];
					   }
					   else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
						   potent[GRADZPRESS][iB4] = potent[GRADZPRESS][iP];
					   }
					   else {
						   potent[GRADZPRESS][iB4] = 0.0;
					   }

					   potent[GRADXPRESS][iB4] = potent[GRADXPRESS][iP];
					   potent[GRADYPRESS][iB4] = potent[GRADYPRESS][iP];
				   }
			   }

		   }
	   }
	   
#endif
#if (interpol==1)
	   {

        
		   // не работает на АЛИС.
		   if (b_on_adaptive_local_refinement_mesh) {
			   printf("function green_gaussPRESS in module greengauss.c else if (interpol==1) { not worked in ALICE mesh...\n ");
			   system("pause");
			   exit(1);
		   }

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

		if (bE) {

			 integer inumber=iE-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADXPRESS][iE]=potent[GRADXPRESS][iP]+(dxe/dxw)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iW]);
					potent[GRADXPRESS][iE]=0.0;
				}
				else {
					potent[GRADXPRESS][iE]=0.0;
				}

			
			potent[GRADYPRESS][iE]=potent[GRADYPRESS][iP]+(dxe/dxw)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iW]);
			potent[GRADZPRESS][iE]=potent[GRADZPRESS][iP]+(dxe/dxw)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iW]);
			//potent[GRADXPRESS][iE]=potent[GRADXPRESS][iP]+(dxe/dxw)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iW]); // <--
		}

		if (bW) {

                integer inumber=iW-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADXPRESS][iW]=potent[GRADXPRESS][iP]+(dxw/dxe)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iE]);
					potent[GRADXPRESS][iW]=0.0;
				}
				else {
					potent[GRADXPRESS][iW]=0.0;
				}

			
			potent[GRADYPRESS][iW]=potent[GRADYPRESS][iP]+(dxw/dxe)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iE]);
			potent[GRADZPRESS][iW]=potent[GRADZPRESS][iP]+(dxw/dxe)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iE]);
            //potent[GRADXPRESS][iW]=potent[GRADXPRESS][iP]+(dxw/dxe)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iE]); // <--
		}

		if (bN) {

                integer inumber=iN-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADYPRESS][iN]=potent[GRADYPRESS][iP]+(dyn/dys)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iS]);
					potent[GRADYPRESS][iN]=0.0;
				}
				else {
					potent[GRADYPRESS][iN]=0.0;
				}

			potent[GRADXPRESS][iN]=potent[GRADXPRESS][iP]+(dyn/dys)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iS]);
			potent[GRADZPRESS][iN]=potent[GRADZPRESS][iP]+(dyn/dys)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iS]);
			//potent[GRADYPRESS][iN]=potent[GRADYPRESS][iP]+(dyn/dys)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iS]); // <--
		}

		if (bS) {

			integer inumber=iS-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADYPRESS][iS]=potent[GRADYPRESS][iP]+(dys/dyn)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iN]);
					potent[GRADYPRESS][iS]=0.0;
				}
				else {
					potent[GRADYPRESS][iS]=0.0;
				}

			potent[GRADXPRESS][iS]=potent[GRADXPRESS][iP]+(dys/dyn)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iN]);
			potent[GRADZPRESS][iS]=potent[GRADZPRESS][iP]+(dys/dyn)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iN]);
			//potent[GRADYPRESS][iS]=potent[GRADYPRESS][iP]+(dys/dyn)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iN]); // <--
		}

		if (bT) {

			integer inumber=iT-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADZPRESS][iT]=potent[GRADZPRESS][iP]+(dzt/dzb)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iB]);
					potent[GRADZPRESS][iT]=0.0;
				}
				else {
					potent[GRADZPRESS][iT]=0.0;
				}

			potent[GRADXPRESS][iT]=potent[GRADXPRESS][iP]+(dzt/dzb)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iB]);
			potent[GRADYPRESS][iT]=potent[GRADYPRESS][iP]+(dzt/dzb)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iB]);
			//potent[GRADZPRESS][iT]=potent[GRADZPRESS][iP]+(dzt/dzb)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iB]);//<--
			
		}

		if (bB) {

                integer inumber=iB-maxelm;
				if (((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure)) {
					//potent[GRADZPRESS][iB]=potent[GRADZPRESS][iP]+(dzb/dzt)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iT]);
					potent[GRADZPRESS][iB]=0.0;
				}
				else {
					potent[GRADZPRESS][iB]=0.0;
				}

			potent[GRADXPRESS][iB]=potent[GRADXPRESS][iP]+(dzb/dzt)*(potent[GRADXPRESS][iP]-potent[GRADXPRESS][iT]);
			potent[GRADYPRESS][iB]=potent[GRADYPRESS][iP]+(dzb/dzt)*(potent[GRADYPRESS][iP]-potent[GRADYPRESS][iT]);
			//potent[GRADZPRESS][iB]=potent[GRADZPRESS][iP]+(dzb/dzt)*(potent[GRADZPRESS][iP]-potent[GRADZPRESS][iT]);//<--
			
		}

	   }
#endif	
}


} // green_gaussPRESS

// Вычисление градиентов скоростей в центрах внутренних КО
// и на границах с помощью квадратичной интерполяции.
// Поскольку интерполяция квадратичная то точность данной формулы O(h^2).
// данная функция реализована 15 мая 2012 года.
void green_gaussO2(integer iP, doublereal** &potent, int** &nvtx, TOCHKA* &pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond) {
	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP]; iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);


	if (!bbond) {
		// внутренние КО.

		// Воспользуемся формулой из книги Г.З. Гарбера обеспечивающую второй порядок точности.


		doublereal hxminus, hxplus, hyminus, hyplus, hzminus, hzplus;
		// компоненты скорости.
		doublereal VXP, VXW, VXE, VXN, VXS, VXT, VXB;
		doublereal VYP, VYW, VYE, VYN, VYS, VYT, VYB;
		doublereal VZP, VZW, VZE, VZN, VZS, VZT, VZB;

		VXP=potent[VELOCITY_X_COMPONENT][iP]; VYP=potent[VELOCITY_Y_COMPONENT][iP]; VZP=potent[VELOCITY_Z_COMPONENT][iP]; // значения скоростей в центральном контрольном объёме.
		VXW=potent[VELOCITY_X_COMPONENT][iW]; VYW=potent[VELOCITY_Y_COMPONENT][iW]; VZW=potent[VELOCITY_Z_COMPONENT][iW];
		VXE=potent[VELOCITY_X_COMPONENT][iE]; VYE=potent[VELOCITY_Y_COMPONENT][iE]; VZE=potent[VELOCITY_Z_COMPONENT][iE];
		VXS=potent[VELOCITY_X_COMPONENT][iS]; VYS=potent[VELOCITY_Y_COMPONENT][iS]; VZS=potent[VELOCITY_Z_COMPONENT][iS]; 
		VXN=potent[VELOCITY_X_COMPONENT][iN]; VYN=potent[VELOCITY_Y_COMPONENT][iN]; VZN=potent[VELOCITY_Z_COMPONENT][iN]; 
		VXB=potent[VELOCITY_X_COMPONENT][iB]; VYB=potent[VELOCITY_Y_COMPONENT][iB]; VZB=potent[VELOCITY_Z_COMPONENT][iB];  
		VXT=potent[VELOCITY_X_COMPONENT][iT]; VYT=potent[VELOCITY_Y_COMPONENT][iT]; VZT=potent[VELOCITY_Z_COMPONENT][iT];

		if (!bW) {
			TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iW, nvtx, pa, pb,W_SIDE);
		    hxminus=fabs(pp.x-pb.x);		    
	    } else {
		    // узел W граничный
			hxminus=0.5*dx;
	    }
	    if (!bE) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iE, nvtx, pa, pb,E_SIDE);
		    hxplus=fabs(pb.x-pp.x);		
	    } else {
            // узел E граничный
		    hxplus=0.5*dx;
       	}
	    if (!bS) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iS, nvtx, pa, pb,S_SIDE);
		    hyminus=fabs(pp.y-pb.y);
	    } else {
		    // узел S граничный
		    hyminus=0.5*dy;
		}
	    if (!bN) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iN, nvtx, pa, pb,N_SIDE);
		    hyplus=fabs(pb.y-pp.y);		
	    } else {
		    // узел N граничный
		    hyplus=0.5*dy;
	    } 
		if (!bB) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iB, nvtx, pa, pb,B_SIDE);
		    hzminus=fabs(pp.z-pb.z);
		} else {
		    // узел B граничный
		    hzminus=0.5*dz;
	    }
	    if (!bT) { 
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iT, nvtx, pa, pb,T_SIDE);
		    hzplus=fabs(pb.z-pp.z);
	    } else {
		    // узел T граничный
		    hzplus=0.5*dz;
	    }

		// VX
		potent[GRADXVX][iP]=rgradF(VXW, VXP, VXE, hxminus, hxplus); // второй порядок точности.
		potent[GRADYVX][iP]=rgradF(VXS, VXP, VXN, hyminus, hyplus);
	    potent[GRADZVX][iP]=rgradF(VXB, VXP, VXT, hzminus, hzplus);
		// VY
		potent[GRADXVY][iP]=rgradF(VYW, VYP, VYE, hxminus, hxplus);
	    potent[GRADYVY][iP]=rgradF(VYS, VYP, VYN, hyminus, hyplus);
	    potent[GRADZVY][iP]=rgradF(VYB, VYP, VYT, hzminus, hzplus);
		// VZ
		potent[GRADXVZ][iP]=rgradF(VZW, VZP, VZE, hxminus, hxplus);
	    potent[GRADYVZ][iP]=rgradF(VZS, VZP, VZN, hyminus, hyplus);
	    potent[GRADZVZ][iP]=rgradF(VZB, VZP, VZT, hzminus, hzplus);

	}
	else {
		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью квадратичной интерполяции.

		if (bE) {

			// квадратичная интерполяция.

			TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iW, nvtx, pa, pb,W_SIDE);
			center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb,WW_SIDE);
					
			potent[GRADXVX][iE] = my_quadratic_interpolation('+', potent[GRADXVX][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADXVX][iW], potent[GRADXVX][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADYVX][iE] = my_quadratic_interpolation('+', potent[GRADYVX][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADYVX][iW], potent[GRADYVX][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADZVX][iE] = my_quadratic_interpolation('+', potent[GRADZVX][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADZVX][iW], potent[GRADZVX][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);

			potent[GRADXVY][iE] = my_quadratic_interpolation('+', potent[GRADXVY][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADXVY][iW], potent[GRADXVY][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADYVY][iE] = my_quadratic_interpolation('+', potent[GRADYVY][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADYVY][iW], potent[GRADYVY][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADZVY][iE] = my_quadratic_interpolation('+', potent[GRADZVY][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADZVY][iW], potent[GRADZVY][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);

			potent[GRADXVZ][iE] = my_quadratic_interpolation('+', potent[GRADXVZ][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADXVZ][iW], potent[GRADXVZ][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADYVZ][iE] = my_quadratic_interpolation('+', potent[GRADYVZ][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADYVZ][iW], potent[GRADYVZ][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADZVZ][iE] = my_quadratic_interpolation('+', potent[GRADZVZ][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[GRADZVZ][iW], potent[GRADZVZ][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
		}

		if (bW) {

			// квадратичная интерполяция.

			TOCHKA pp, pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iE, nvtx, pa, pb,E_SIDE);
			center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb,EE_SIDE);

			potent[GRADXVX][iW] = my_quadratic_interpolation('-', potent[GRADXVX][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADXVX][iE], potent[GRADXVX][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADYVX][iW] = my_quadratic_interpolation('-', potent[GRADYVX][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADYVX][iE], potent[GRADYVX][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADZVX][iW] = my_quadratic_interpolation('-', potent[GRADZVX][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADZVX][iE], potent[GRADZVX][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);

			potent[GRADXVY][iW] = my_quadratic_interpolation('-', potent[GRADXVY][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADXVY][iE], potent[GRADXVY][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADYVY][iW] = my_quadratic_interpolation('-', potent[GRADYVY][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADYVY][iE], potent[GRADYVY][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADZVY][iW] = my_quadratic_interpolation('-', potent[GRADZVY][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADZVY][iE], potent[GRADZVY][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);

			potent[GRADXVZ][iW] = my_quadratic_interpolation('-', potent[GRADXVZ][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADXVZ][iE], potent[GRADXVZ][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADYVZ][iW] = my_quadratic_interpolation('-', potent[GRADYVZ][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADYVZ][iE], potent[GRADYVZ][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADZVZ][iW] = my_quadratic_interpolation('-', potent[GRADZVZ][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[GRADZVZ][iE], potent[GRADZVZ][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
		}

		if (bN) {

			// квадратичная интерполяция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iS, nvtx, pa, pb,S_SIDE);
			center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb,SS_SIDE);

			potent[GRADXVX][iN] = my_quadratic_interpolation('+', potent[GRADXVX][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADXVX][iS], potent[GRADXVX][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADYVX][iN] = my_quadratic_interpolation('+', potent[GRADYVX][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADYVX][iS], potent[GRADYVX][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADZVX][iN] = my_quadratic_interpolation('+', potent[GRADZVX][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADZVX][iS], potent[GRADZVX][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);

			potent[GRADXVY][iN] = my_quadratic_interpolation('+', potent[GRADXVY][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADXVY][iS], potent[GRADXVY][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADYVY][iN] = my_quadratic_interpolation('+', potent[GRADYVY][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADYVY][iS], potent[GRADYVY][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADZVY][iN] = my_quadratic_interpolation('+', potent[GRADZVY][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADZVY][iS], potent[GRADZVY][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);

			potent[GRADXVZ][iN] = my_quadratic_interpolation('+', potent[GRADXVZ][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADXVZ][iS], potent[GRADXVZ][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADYVZ][iN] = my_quadratic_interpolation('+', potent[GRADYVZ][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADYVZ][iS], potent[GRADYVZ][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADZVZ][iN] = my_quadratic_interpolation('+', potent[GRADZVZ][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[GRADZVZ][iS], potent[GRADZVZ][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
		}

		if (bS) {

            // квадратичная интерполяция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp, 100);
		    center_cord3D(iN, nvtx, pa, pb, N_SIDE);
			center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

			potent[GRADXVX][iS] = my_quadratic_interpolation('-', potent[GRADXVX][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADXVX][iN], potent[GRADXVX][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADYVX][iS] = my_quadratic_interpolation('-', potent[GRADYVX][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADYVX][iN], potent[GRADYVX][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADZVX][iS] = my_quadratic_interpolation('-', potent[GRADZVX][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADZVX][iN], potent[GRADZVX][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);

			potent[GRADXVY][iS] = my_quadratic_interpolation('-', potent[GRADXVY][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADXVY][iN], potent[GRADXVY][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADYVY][iS] = my_quadratic_interpolation('-', potent[GRADYVY][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADYVY][iN], potent[GRADYVY][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADZVY][iS] = my_quadratic_interpolation('-', potent[GRADZVY][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADZVY][iN], potent[GRADZVY][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);

			potent[GRADXVZ][iS] = my_quadratic_interpolation('-', potent[GRADXVZ][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADXVZ][iN], potent[GRADXVZ][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADYVZ][iS] = my_quadratic_interpolation('-', potent[GRADYVZ][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADYVZ][iN], potent[GRADYVZ][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADZVZ][iS] = my_quadratic_interpolation('-', potent[GRADZVZ][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[GRADZVZ][iN], potent[GRADZVZ][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
		}

		if (bT) {

			// квадратичная интерполяция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iB, nvtx, pa, pb,B_SIDE);
			center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb,BB_SIDE);
					
			potent[GRADXVX][iT] = my_quadratic_interpolation('+', potent[GRADXVX][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADXVX][iB], potent[GRADXVX][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADYVX][iT] = my_quadratic_interpolation('+', potent[GRADYVX][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADYVX][iB], potent[GRADYVX][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADZVX][iT] = my_quadratic_interpolation('+', potent[GRADZVX][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADZVX][iB], potent[GRADZVX][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);

			potent[GRADXVY][iT] = my_quadratic_interpolation('+', potent[GRADXVY][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADXVY][iB], potent[GRADXVY][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADYVY][iT] = my_quadratic_interpolation('+', potent[GRADYVY][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADYVY][iB], potent[GRADYVY][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADZVY][iT] = my_quadratic_interpolation('+', potent[GRADZVY][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADZVY][iB], potent[GRADZVY][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);

			potent[GRADXVZ][iT] = my_quadratic_interpolation('+', potent[GRADXVZ][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADXVZ][iB], potent[GRADXVZ][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADYVZ][iT] = my_quadratic_interpolation('+', potent[GRADYVZ][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADYVZ][iB], potent[GRADYVZ][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADZVZ][iT] = my_quadratic_interpolation('+', potent[GRADZVZ][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[GRADZVZ][iB], potent[GRADZVZ][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
		}

		if (bB) {

			// квадратичная интерполяция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iT, nvtx, pa, pb,T_SIDE);
			center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);


			potent[GRADXVX][iB] = my_quadratic_interpolation('-', potent[GRADXVX][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADXVX][iT], potent[GRADXVX][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADYVX][iB] = my_quadratic_interpolation('-', potent[GRADYVX][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADYVX][iT], potent[GRADYVX][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADZVX][iB] = my_quadratic_interpolation('-', potent[GRADZVX][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADZVX][iT], potent[GRADZVX][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);

			potent[GRADXVY][iB] = my_quadratic_interpolation('-', potent[GRADXVY][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADXVY][iT], potent[GRADXVY][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADYVY][iB] = my_quadratic_interpolation('-', potent[GRADYVY][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADYVY][iT], potent[GRADYVY][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADZVY][iB] = my_quadratic_interpolation('-', potent[GRADZVY][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADZVY][iT], potent[GRADZVY][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);

			potent[GRADXVZ][iB] = my_quadratic_interpolation('-', potent[GRADXVZ][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADXVZ][iT], potent[GRADXVZ][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADYVZ][iB] = my_quadratic_interpolation('-', potent[GRADYVZ][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADYVZ][iT], potent[GRADYVZ][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADZVZ][iB] = my_quadratic_interpolation('-', potent[GRADZVZ][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[GRADZVZ][iT], potent[GRADZVZ][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
		}
	}

} // green_gaussO2

// нахождение производных от скорости первого или второго порядка точности.
void green_gauss(integer iP, doublereal** &potent, int** &nvtx, TOCHKA* pa,
	int*** &neighbors_for_the_internal_node, integer maxelm, bool bbond, FLOW &f,
	BOUND* &border_neighbor, integer *ilevel_alice) {

	// если bsecondorder==true то производные будут вычисляться со вторым порядком точности.
	bool bsecondorder=false; // если false то внутри может применяться монотонизирующая поправка Рхи-Чоу.

	if (bsecondorder) {
		// второй порядок точности.
		green_gaussO2(iP, potent, nvtx, pa, neighbors_for_the_internal_node, maxelm, bbond);
	}
	else {
		// первый порядок точности.
		green_gaussO1(iP, potent, nvtx, pa, neighbors_for_the_internal_node, maxelm, bbond, f.mf[iP], f.prop[RHO], f.prop_b[RHO], border_neighbor, ilevel_alice, f.volume);
	}

} // green_gauss

// 14.04.2019; 03.10.2019. Работает на АЛИС сетке.
// Вычисление градиентов деформаций
// в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
// Поскольку интерполяция линейная то точность данной формулы O(h). 
// По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
// что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
// совпадает со взвешенным методом наименьших квадратов.
void green_gauss_Stress(integer iP,
	doublereal**& potent, int**& nvtx, TOCHKA*& pa,
	int***& neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND*& border_neighbor, integer* ilevel_alice, integer iDATA, integer iTARGET, LINE_DIRECTIONAL iDIRECTIONAL)
{

	// iDATA - обрабатываемые данные,
	// iTARGET -куда пишем результат обработки.
	// iDIRECTIONAL X_LINE_DIRECTIONAL-X, Y_LINE_DIRECTIONAL - Y, Z_LINE_DIRECTIONAL - Z.

	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу.
	//bool bRCh = false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
	// Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

	if (b_on_adaptive_local_refinement_mesh) {

		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;

	if (iE2 >= maxelm) bE2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iB2 >= maxelm) bB2 = true;

	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;

	if (iE3 >= maxelm) bE3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iB3 >= maxelm) bB3 = true;

	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE4 >= maxelm) bE4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);
	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	doublereal dxe = 0.5 * dx, dxw = 0.5 * dx, dyn = 0.5 * dy, dys = 0.5 * dy, dzt = 0.5 * dz, dzb = 0.5 * dz;
	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE > -1) {
		if (!bE) dxe = 0.5 * (pa[nvtx[1][iE] - 1].x + pa[nvtx[0][iE] - 1].x);
		if (!bE) dxe -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW > -1) {
		if (!bW) dxw = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW) dxw -= 0.5 * (pa[nvtx[1][iW] - 1].x + pa[nvtx[0][iW] - 1].x);
	}
	// y - direction
	if (iN > -1) {
		if (!bN) dyn = 0.5 * (pa[nvtx[2][iN] - 1].y + pa[nvtx[0][iN] - 1].y);
		if (!bN) dyn -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS > -1) {
		if (!bS) dys = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS) dys -= 0.5 * (pa[nvtx[2][iS] - 1].y + pa[nvtx[0][iS] - 1].y);
	}
	// z - direction
	if (iT > -1) {
		if (!bT) dzt = 0.5 * (pa[nvtx[4][iT] - 1].z + pa[nvtx[0][iT] - 1].z);
		if (!bT) dzt -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB > -1) {
		if (!bB) dzb = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB) dzb -= 0.5 * (pa[nvtx[4][iB] - 1].z + pa[nvtx[0][iB] - 1].z);
	}


	doublereal dxe2 = 0.5 * dx, dxw2 = 0.5 * dx, dyn2 = 0.5 * dy, dys2 = 0.5 * dy, dzt2 = 0.5 * dz, dzb2 = 0.5 * dz;
	doublereal dxe3 = 0.5 * dx, dxw3 = 0.5 * dx, dyn3 = 0.5 * dy, dys3 = 0.5 * dy, dzt3 = 0.5 * dz, dzb3 = 0.5 * dz;
	doublereal dxe4 = 0.5 * dx, dxw4 = 0.5 * dx, dyn4 = 0.5 * dy, dys4 = 0.5 * dy, dzt4 = 0.5 * dz, dzb4 = 0.5 * dz;

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE2 > -1) {
		if (!bE2) dxe2 = 0.5 * (pa[nvtx[1][iE2] - 1].x + pa[nvtx[0][iE2] - 1].x);
		if (!bE2) dxe2 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW2 > -1) {
		if (!bW2) dxw2 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW2) dxw2 -= 0.5 * (pa[nvtx[1][iW2] - 1].x + pa[nvtx[0][iW2] - 1].x);
	}
	// y - direction
	if (iN2 > -1) {
		if (!bN2) dyn2 = 0.5 * (pa[nvtx[2][iN2] - 1].y + pa[nvtx[0][iN2] - 1].y);
		if (!bN2) dyn2 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS2 > -1) {
		if (!bS2) dys2 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS2) dys2 -= 0.5 * (pa[nvtx[2][iS2] - 1].y + pa[nvtx[0][iS2] - 1].y);
	}
	// z - direction
	if (iT2 > -1) {
		if (!bT2) dzt2 = 0.5 * (pa[nvtx[4][iT2] - 1].z + pa[nvtx[0][iT2] - 1].z);
		if (!bT2) dzt2 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB2 > -1) {
		if (!bB2) dzb2 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB2) dzb2 -= 0.5 * (pa[nvtx[4][iB2] - 1].z + pa[nvtx[0][iB2] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE3 > -1) {
		if (!bE3) dxe3 = 0.5 * (pa[nvtx[1][iE3] - 1].x + pa[nvtx[0][iE3] - 1].x);
		if (!bE3) dxe3 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW3 > -1) {
		if (!bW3) dxw3 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW3) dxw3 -= 0.5 * (pa[nvtx[1][iW3] - 1].x + pa[nvtx[0][iW3] - 1].x);
	}
	// y - direction
	if (iN3 > -1) {
		if (!bN3) dyn3 = 0.5 * (pa[nvtx[2][iN3] - 1].y + pa[nvtx[0][iN3] - 1].y);
		if (!bN3) dyn3 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS3 > -1) {
		if (!bS3) dys3 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS3) dys3 -= 0.5 * (pa[nvtx[2][iS3] - 1].y + pa[nvtx[0][iS3] - 1].y);
	}
	// z - direction
	if (iT3 > -1) {
		if (!bT3) dzt3 = 0.5 * (pa[nvtx[4][iT3] - 1].z + pa[nvtx[0][iT3] - 1].z);
		if (!bT3) dzt3 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB3 > -1) {
		if (!bB3) dzb3 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB3) dzb3 -= 0.5 * (pa[nvtx[4][iB3] - 1].z + pa[nvtx[0][iB3] - 1].z);
	}

	// т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
	if (iE4 > -1) {
		if (!bE4) dxe4 = 0.5 * (pa[nvtx[1][iE4] - 1].x + pa[nvtx[0][iE4] - 1].x);
		if (!bE4) dxe4 -= 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	}
	if (iW4 > -1) {
		if (!bW4) dxw4 = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		if (!bW4) dxw4 -= 0.5 * (pa[nvtx[1][iW4] - 1].x + pa[nvtx[0][iW4] - 1].x);
	}
	// y - direction
	if (iN4 > -1) {
		if (!bN4) dyn4 = 0.5 * (pa[nvtx[2][iN4] - 1].y + pa[nvtx[0][iN4] - 1].y);
		if (!bN4) dyn4 -= 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	}
	if (iS4 > -1) {
		if (!bS4) dys4 = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		if (!bS4) dys4 -= 0.5 * (pa[nvtx[2][iS4] - 1].y + pa[nvtx[0][iS4] - 1].y);
	}
	// z - direction
	if (iT4 > -1) {
		if (!bT4) dzt4 = 0.5 * (pa[nvtx[4][iT4] - 1].z + pa[nvtx[0][iT4] - 1].z);
		if (!bT4) dzt4 -= 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
	}
	if (iB4 > -1) {
		if (!bB4) dzb4 = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);
		if (!bB4) dzb4 -= 0.5 * (pa[nvtx[4][iB4] - 1].z + pa[nvtx[0][iB4] - 1].z);
	}

	dxe = fabs(dxe);
	dxe2 = fabs(dxe2);
	dxe3 = fabs(dxe3);
	dxe4 = fabs(dxe4);

	dxw = fabs(dxw);
	dxw2 = fabs(dxw2);
	dxw3 = fabs(dxw3);
	dxw4 = fabs(dxw4);

	dyn = fabs(dyn);
	dyn2 = fabs(dyn2);
	dyn3 = fabs(dyn3);
	dyn4 = fabs(dyn4);

	dys = fabs(dys);
	dys2 = fabs(dys2);
	dys3 = fabs(dys3);
	dys4 = fabs(dys4);

	dzt = fabs(dzt);
	dzt2 = fabs(dzt2);
	dzt3 = fabs(dzt3);
	dzt4 = fabs(dzt4);

	dzb = fabs(dzb);
	dzb2 = fabs(dzb2);
	dzb3 = fabs(dzb3);
	dzb4 = fabs(dzb4);

	// Учёт неравномерности расчётной сетки:
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

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.




	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}


	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = border_neighbor[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = border_neighbor[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = border_neighbor[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = border_neighbor[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = border_neighbor[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.



	if (iE2 > -1) {

		dSqe2 = dy * dz;

		if (bE2) {
			// граничный узел.
			dSqe2 = border_neighbor[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE2]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW2) {
			// граничный узел.
			dSqw2 = border_neighbor[iW2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW2]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = border_neighbor[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN2]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = border_neighbor[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS2]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = border_neighbor[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT2]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = border_neighbor[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB2]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb2 = dx_loc * dy_loc;
			}
		}


	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.



	if (iE3 > -1) {

		dSqe3 = dy * dz;

		if (bE3) {
			// граничный узел.
			dSqe3 = border_neighbor[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE3]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = border_neighbor[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW3]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = border_neighbor[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN3]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = border_neighbor[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS3]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = border_neighbor[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT3]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = border_neighbor[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB3]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb3 = dx_loc * dy_loc;
			}
		}


	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.



	if (iE4 > -1) {

		dSqe4 = dy * dz;

		if (bE4) {
			// граничный узел.
			dSqe4 = border_neighbor[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iE4]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = border_neighbor[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iW4]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = border_neighbor[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iN4]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = border_neighbor[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iS4]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = border_neighbor[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iT4]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = border_neighbor[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[iP] >= ilevel_alice[iB4]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	// 28.04.2019	
	if (fabs(dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqe %e %e %e %e\n", dSqe, dSqe2, dSqe3, dSqe4);
		//printf("dSqw %e %e %e %e\n", dSqw, dSqw2, dSqw3, dSqw4);
		//printf("disbalanse: %e \n", dSqe + dSqe2 + dSqe3 + dSqe4 - dSqw - dSqw2 - dSqw3 - dSqw4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSE = dSqe + dSqe2 + dSqe3 + dSqe4;
		doublereal dSW = dSqw + dSqw2 + dSqw3 + dSqw4;
		doublereal km = (dy * dz) / dSE;
		dSqe *= km; dSqe2 *= km; dSqe3 *= km; dSqe4 *= km;
		km = (dy * dz) / dSW;
		dSqw *= km; dSqw2 *= km; dSqw3 *= km; dSqw4 *= km;
	}

	if (fabs(dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqn %e %e %e %e\n", dSqn, dSqn2, dSqn3, dSqn4);
		//printf("dSqs %e %e %e %e\n", dSqs, dSqs2, dSqs3, dSqs4);
		//printf("disbalanse: %e \n", dSqn + dSqn2 + dSqn3 + dSqn4 - dSqs - dSqs2 - dSqs3 - dSqs4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dSN = dSqn + dSqn2 + dSqn3 + dSqn4;
		doublereal dSS = dSqs + dSqs2 + dSqs3 + dSqs4;
		doublereal km = (dx * dz) / dSN;
		dSqn *= km; dSqn2 *= km; dSqn3 *= km; dSqn4 *= km;
		km = (dx * dz) / dSS;
		dSqs *= km; dSqs2 *= km; dSqs3 *= km; dSqs4 *= km;
	}

	if (fabs(dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4) > 1.0e-36) {
		// Небольшой дисбаланс присутствует.
		//printf("dSqt %e %e %e %e\n", dSqt, dSqt2, dSqt3, dSqt4);
		//printf("dSqb %e %e %e %e\n", dSqb, dSqb2, dSqb3, dSqb4);
		//printf("disbalanse: %e \n", dSqt + dSqt2 + dSqt3 + dSqt4 - dSqb - dSqb2 - dSqb3 - dSqb4);
		//system("PAUSE");
		// Вводим корректирующую поправку.
		doublereal dST = dSqt + dSqt2 + dSqt3 + dSqt4;
		doublereal dSB = dSqb + dSqb2 + dSqb3 + dSqb4;
		doublereal km = (dx * dy) / dST;
		dSqt *= km; dSqt2 *= km; dSqt3 *= km; dSqt4 *= km;
		km = (dx * dy) / dSB;
		dSqb *= km; dSqb2 *= km; dSqb3 *= km; dSqb4 *= km;
	}





	// линейная интерполяция скорости VX на грань КО.
	doublereal fe = 0.0, fw = 0.0, fn = 0.0, fs = 0.0, ft = 0.0, fb = 0.0;
	doublereal fe2 = 0.0, fw2 = 0.0, fn2 = 0.0, fs2 = 0.0, ft2 = 0.0, fb2 = 0.0;
	doublereal fe3 = 0.0, fw3 = 0.0, fn3 = 0.0, fs3 = 0.0, ft3 = 0.0, fb3 = 0.0;
	doublereal fe4 = 0.0, fw4 = 0.0, fn4 = 0.0, fs4 = 0.0, ft4 = 0.0, fb4 = 0.0;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// iDATA
		if (iE > -1) {
			if (!bE) {
				fe = feplus * potent[iDATA][iE] + (1.0 - feplus) * potent[iDATA][iP];
			}
			else {
				fe = potent[iDATA][iE];
				
			}
		}
		if (iW > -1) {
			if (!bW) {
				fw = fwplus * potent[iDATA][iW] + (1.0 - fwplus) * potent[iDATA][iP];
			}
			else fw = potent[iDATA][iW];
		}

		if (iE2 > -1) {
			if (!bE2) {
				fe2 = feplus2 * potent[iDATA][iE2] + (1.0 - feplus2) * potent[iDATA][iP];
			}
			else fe2 = potent[iDATA][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) {
				fw2 = fwplus2 * potent[iDATA][iW2] + (1.0 - fwplus2) * potent[iDATA][iP];
			}
			else fw2 = potent[iDATA][iW2];
		}

		if (iE3 > -1) {
			if (!bE3) {
				fe3 = feplus3 * potent[iDATA][iE3] + (1.0 - feplus3) * potent[iDATA][iP];
			}
			else fe3 = potent[iDATA][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) {
				fw3 = fwplus3 * potent[iDATA][iW3] + (1.0 - fwplus3) * potent[iDATA][iP];
			}
			else fw3 = potent[iDATA][iW3];
		}

		if (iE4 > -1) {
			if (!bE4) {
				fe4 = feplus4 * potent[iDATA][iE4] + (1.0 - feplus4) * potent[iDATA][iP];
			}
			else fe4 = potent[iDATA][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) {
				fw4 = fwplus4 * potent[iDATA][iW4] + (1.0 - fwplus4) * potent[iDATA][iP];
			}
			else fw4 = potent[iDATA][iW4];
		}
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнейшие компоненты скорости и производные аналогично для VX.
		if (iN > -1) {
			if (!bN) fn = fnplus * potent[iDATA][iN] + (1.0 - fnplus) * potent[iDATA][iP]; else fn = potent[iDATA][iN];
		}
		if (iS > -1) {
			if (!bS) fs = fsplus * potent[iDATA][iS] + (1.0 - fsplus) * potent[iDATA][iP]; else fs = potent[iDATA][iS];
		}
		if (iT > -1) {
			if (!bT) ft = ftplus * potent[iDATA][iT] + (1.0 - ftplus) * potent[iDATA][iP]; else ft = potent[iDATA][iT];
		}
		if (iB > -1) {
			if (!bB) fb = fbplus * potent[iDATA][iB] + (1.0 - fbplus) * potent[iDATA][iP]; else fb = potent[iDATA][iB];
		}

		if (iN2 > -1) {
			if (!bN2) fn2 = fnplus2 * potent[iDATA][iN2] + (1.0 - fnplus2) * potent[iDATA][iP]; else fn2 = potent[iDATA][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) fs2 = fsplus2 * potent[iDATA][iS2] + (1.0 - fsplus2) * potent[iDATA][iP]; else fs2 = potent[iDATA][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) ft2 = ftplus2 * potent[iDATA][iT2] + (1.0 - ftplus2) * potent[iDATA][iP]; else ft2 = potent[iDATA][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) fb2 = fbplus2 * potent[iDATA][iB2] + (1.0 - fbplus2) * potent[iDATA][iP]; else fb2 = potent[iDATA][iB2];
		}

		if (iN3 > -1) {
			if (!bN3) fn3 = fnplus3 * potent[iDATA][iN3] + (1.0 - fnplus3) * potent[iDATA][iP]; else fn3 = potent[iDATA][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) fs3 = fsplus3 * potent[iDATA][iS3] + (1.0 - fsplus3) * potent[iDATA][iP]; else fs3 = potent[iDATA][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) ft3 = ftplus3 * potent[iDATA][iT3] + (1.0 - ftplus3) * potent[iDATA][iP]; else ft3 = potent[iDATA][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) fb3 = fbplus3 * potent[iDATA][iB3] + (1.0 - fbplus3) * potent[iDATA][iP]; else fb3 = potent[iDATA][iB3];
		}

		if (iN4 > -1) {
			if (!bN4) fn4 = fnplus4 * potent[iDATA][iN4] + (1.0 - fnplus4) * potent[iDATA][iP]; else fn4 = potent[iDATA][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) fs4 = fsplus4 * potent[iDATA][iS4] + (1.0 - fsplus4) * potent[iDATA][iP]; else fs4 = potent[iDATA][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) ft4 = ftplus4 * potent[iDATA][iT4] + (1.0 - ftplus4) * potent[iDATA][iP]; else ft4 = potent[iDATA][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) fb4 = fbplus4 * potent[iDATA][iB4] + (1.0 - fbplus4) * potent[iDATA][iP]; else fb4 = potent[iDATA][iB4];
		}
		// градиент iDATA
		switch (iDIRECTIONAL) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // X
			//potent[iTARGET][iP]=(fe-fw)/dx;
			potent[iTARGET][iP] = (fe * dSqe / (dy * dz) + fe2 * dSqe2 / (dy * dz) + fe3 * dSqe3 / (dy * dz) + fe4 * dSqe4 / (dy * dz) - (fw * dSqw / (dy * dz) + fw2 * dSqw2 / (dy * dz) + fw3 * dSqw3 / (dy * dz) + fw4 * dSqw4 / (dy * dz))) / dx;
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: // Y
			//potent[iTARGET][iP]=(fn-fs)/dy;
			potent[iTARGET][iP] = (fn * dSqn / (dx * dz) + fn2 * dSqn2 / (dx * dz) + fn3 * dSqn3 / (dx * dz) + fn4 * dSqn4 / (dx * dz) - (fs * dSqs / (dx * dz) + fs2 * dSqs2 / (dx * dz) + fs3 * dSqs3 / (dx * dz) + fs4 * dSqs4 / (dx * dz))) / dy;
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // Z
			//potent[iTARGET][iP]=(ft-fb)/dz;
			potent[iTARGET][iP] = (ft * dSqt / (dx * dy) + ft2 * dSqt2 / (dx * dy) + ft3 * dSqt3 / (dx * dy) + ft4 * dSqt4 / (dx * dy) - (fb * dSqb / (dx * dy) + fb2 * dSqb2 / (dx * dy) + fb3 * dSqb3 / (dx * dy) + fb4 * dSqb4 / (dx * dy))) / dz;
			break;
			}

	}
	else {

		if (0) {



			// По простому: градиент на границе наследуем из ближайшего внутреннего узла.
			if (iE > -1) {
				if (bE) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					//doublereal dspeed = sqrt((potent[VXCOR][iE]) * (potent[VXCOR][iE]) + (potent[VYCOR][iE]) * (potent[VYCOR][iE]) + (potent[VZCOR][iE]) * (potent[VZCOR][iE]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iE] = 0.0;
						//potent[iTARGET][iE] = 0.0;
						//potent[iTARGET][iE] = 0.0;



						/*
						potent[iTARGET][iP] = 0.0;
						potent[iTARGET][iP] = 0.0;
						potent[iTARGET][iP] = 0.0;


						*/
					//}
					//else {

						potent[iTARGET][iE] = potent[iTARGET][iP];
						//potent[iTARGET][iE] = potent[iTARGET][iP];
						//potent[iTARGET][iE] = potent[iTARGET][iP];


					//}
				}
			}



			if (iE2 > -1) {
				if (bE2) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					//doublereal dspeed = sqrt((potent[VXCOR][iE2]) * (potent[VXCOR][iE2]) + (potent[VYCOR][iE2]) * (potent[VYCOR][iE2]) + (potent[VZCOR][iE2]) * (potent[VZCOR][iE2]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iE2] = 0.0;
						//potent[iTARGET][iE2] = 0.0;
						//potent[iTARGET][iE2] = 0.0;




					//}
					//else {

						potent[iTARGET][iE2] = potent[iTARGET][iP];
						//potent[iTARGET][iE2] = potent[iTARGET][iP];
						//potent[iTARGET][iE2] = potent[iTARGET][iP];


					//}
				}
			}

			if (iE3 > -1) {
				if (bE3) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					//doublereal dspeed = sqrt((potent[VXCOR][iE3]) * (potent[VXCOR][iE3]) + (potent[VYCOR][iE3]) * (potent[VYCOR][iE3]) + (potent[VZCOR][iE3]) * (potent[VZCOR][iE3]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iE3] = 0.0;
						//potent[iTARGET][iE3] = 0.0;
						//potent[iTARGET][iE3] = 0.0;




					//}
					//else {

						potent[iTARGET][iE3] = potent[iTARGET][iP];
						//potent[iTARGET][iE3] = potent[iTARGET][iP];
						//potent[iTARGET][iE3] = potent[iTARGET][iP];


					//}
				}
			}

			if (iE4 > -1) {
				if (bE4) {

					// 10.02.2017
					// Если на стенке выставлено условие прилипания то градиент скорости на стенке также тождественно равен нулю.

					//doublereal dspeed = sqrt((potent[VXCOR][iE4]) * (potent[VXCOR][iE4]) + (potent[VYCOR][iE4]) * (potent[VYCOR][iE4]) + (potent[VZCOR][iE4]) * (potent[VZCOR][iE4]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iE4] = 0.0;
						//potent[iTARGET][iE4] = 0.0;
						//potent[iTARGET][iE4] = 0.0;




					//}
					//else {

						potent[iTARGET][iE4] = potent[iTARGET][iP];
						//potent[iTARGET][iE4] = potent[iTARGET][iP];
						//potent[iTARGET][iE4] = potent[iTARGET][iP];


					//}
				}
			}

			if (iW > -1) {
				if (bW) {

					//doublereal dspeed = sqrt((potent[VXCOR][iW]) * (potent[VXCOR][iW]) + (potent[VYCOR][iW]) * (potent[VYCOR][iW]) + (potent[VZCOR][iW]) * (potent[VZCOR][iW]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iW] = 0.0;
						//potent[iTARGET][iW] = 0.0;
						//potent[iTARGET][iW] = 0.0;


					//}
					//else {

						potent[iTARGET][iW] = potent[iTARGET][iP];
						//potent[iTARGET][iW] = potent[iTARGET][iP];
						//potent[iTARGET][iW] = potent[iTARGET][iP];


					//}
				}
			}

			if (iW2 > -1) {
				if (bW2) {

					//doublereal dspeed = sqrt((potent[VXCOR][iW2]) * (potent[VXCOR][iW2]) + (potent[VYCOR][iW2]) * (potent[VYCOR][iW2]) + (potent[VZCOR][iW2]) * (potent[VZCOR][iW2]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iW2] = 0.0;
						//potent[iTARGET][iW2] = 0.0;
						//potent[iTARGET][iW2] = 0.0;


					//}
					//else {

						potent[iTARGET][iW2] = potent[iTARGET][iP];
						//potent[iTARGET][iW2] = potent[iTARGET][iP];
						//potent[iTARGET][iW2] = potent[iTARGET][iP];


					//}
				}
			}

			if (iW3 > -1) {
				if (bW3) {

					//doublereal dspeed = sqrt((potent[VXCOR][iW3]) * (potent[VXCOR][iW3]) + (potent[VYCOR][iW3]) * (potent[VYCOR][iW3]) + (potent[VZCOR][iW3]) * (potent[VZCOR][iW3]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iW3] = 0.0;
						//potent[iTARGET][iW3] = 0.0;
						//potent[iTARGET][iW3] = 0.0;


					//}
					//else {

						potent[iTARGET][iW3] = potent[iTARGET][iP];
						//potent[iTARGET][iW3] = potent[iTARGET][iP];
						//potent[iTARGET][iW3] = potent[iTARGET][iP];


					//}
				}
			}

			if (iW4 > -1) {
				if (bW4) {

					//doublereal dspeed = sqrt((potent[VXCOR][iW4]) * (potent[VXCOR][iW4]) + (potent[VYCOR][iW4]) * (potent[VYCOR][iW4]) + (potent[VZCOR][iW4]) * (potent[VZCOR][iW4]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iW4] = 0.0;
						//potent[iTARGET][iW4] = 0.0;
						//potent[iTARGET][iW4] = 0.0;


					//}
					//else {

						potent[iTARGET][iW4] = potent[iTARGET][iP];
						//potent[iTARGET][iW4] = potent[iTARGET][iP];
						//potent[iTARGET][iW4] = potent[iTARGET][iP];


					//}
				}
			}

			if (iN > -1) {
				if (bN) {

					//doublereal dspeed = sqrt((potent[VXCOR][iN]) * (potent[VXCOR][iN]) + (potent[VYCOR][iN]) * (potent[VYCOR][iN]) + (potent[VZCOR][iN]) * (potent[VZCOR][iN]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iN] = 0.0;
						//potent[iTARGET][iN] = 0.0;
						//potent[iTARGET][iN] = 0.0;


					//}
					//else {

						potent[iTARGET][iN] = potent[iTARGET][iP];
						//potent[iTARGET][iN] = potent[iTARGET][iP];
						//potent[iTARGET][iN] = potent[iTARGET][iP];


					//}
				}
			}

			if (iN2 > -1) {
				if (bN2) {

					//doublereal dspeed = sqrt((potent[VXCOR][iN2]) * (potent[VXCOR][iN2]) + (potent[VYCOR][iN2]) * (potent[VYCOR][iN2]) + (potent[VZCOR][iN2]) * (potent[VZCOR][iN2]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iN2] = 0.0;
						//potent[iTARGET][iN2] = 0.0;
						//potent[iTARGET][iN2] = 0.0;


					//}
					//else {

						potent[iTARGET][iN2] = potent[iTARGET][iP];
						//potent[iTARGET][iN2] = potent[iTARGET][iP];
						//potent[iTARGET][iN2] = potent[iTARGET][iP];


					//}
				}
			}

			if (iN3 > -1) {
				if (bN3) {

					//doublereal dspeed = sqrt((potent[VXCOR][iN3]) * (potent[VXCOR][iN3]) + (potent[VYCOR][iN3]) * (potent[VYCOR][iN3]) + (potent[VZCOR][iN3]) * (potent[VZCOR][iN3]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iN3] = 0.0;
						//potent[iTARGET][iN3] = 0.0;
						//potent[iTARGET][iN3] = 0.0;


					//}
					//else {

						potent[iTARGET][iN3] = potent[iTARGET][iP];
						//potent[iTARGET][iN3] = potent[iTARGET][iP];
						//potent[iTARGET][iN3] = potent[iTARGET][iP];


					//}
				}
			}

			if (iN4 > -1) {
				if (bN4) {

					//doublereal dspeed = sqrt((potent[VXCOR][iN4]) * (potent[VXCOR][iN4]) + (potent[VYCOR][iN4]) * (potent[VYCOR][iN4]) + (potent[VZCOR][iN4]) * (potent[VZCOR][iN4]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iN4] = 0.0;
						//potent[iTARGET][iN4] = 0.0;
						//potent[iTARGET][iN4] = 0.0;


					//}
					//else {

						potent[iTARGET][iN4] = potent[iTARGET][iP];
						//potent[iTARGET][iN4] = potent[iTARGET][iP];
						//potent[iTARGET][iN4] = potent[iTARGET][iP];


					//}
				}
			}

			if (iS > -1) {
				if (bS) {

					//doublereal dspeed = sqrt((potent[VXCOR][iS]) * (potent[VXCOR][iS]) + (potent[VYCOR][iS]) * (potent[VYCOR][iS]) + (potent[VZCOR][iS]) * (potent[VZCOR][iS]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iS] = 0.0;
						//potent[iTARGET][iS] = 0.0;
						//potent[iTARGET][iS] = 0.0;


					//}
					//else {

						potent[iTARGET][iS] = potent[iTARGET][iP];
						//potent[iTARGET][iS] = potent[iTARGET][iP];
						//potent[iTARGET][iS] = potent[iTARGET][iP];


					//}
				}
			}

			if (iS2 > -1) {
				if (bS2) {

					//doublereal dspeed = sqrt((potent[VXCOR][iS2]) * (potent[VXCOR][iS2]) + (potent[VYCOR][iS2]) * (potent[VYCOR][iS2]) + (potent[VZCOR][iS2]) * (potent[VZCOR][iS2]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iS2] = 0.0;
						//potent[iTARGET][iS2] = 0.0;
						//potent[iTARGET][iS2] = 0.0;


					//}
					//else {

						potent[iTARGET][iS2] = potent[iTARGET][iP];
						//potent[iTARGET][iS2] = potent[iTARGET][iP];
						//potent[iTARGET][iS2] = potent[iTARGET][iP];


					//}
				}
			}

			if (iS3 > -1) {
				if (bS3) {

					//doublereal dspeed = sqrt((potent[VXCOR][iS3]) * (potent[VXCOR][iS3]) + (potent[VYCOR][iS3]) * (potent[VYCOR][iS3]) + (potent[VZCOR][iS3]) * (potent[VZCOR][iS3]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iS3] = 0.0;
						//potent[iTARGET][iS3] = 0.0;
						//potent[iTARGET][iS3] = 0.0;


					//}
					//else {

						potent[iTARGET][iS3] = potent[iTARGET][iP];
						//potent[iTARGET][iS3] = potent[iTARGET][iP];
						//potent[iTARGET][iS3] = potent[iTARGET][iP];


					//}
				}
			}

			if (iS4 > -1) {
				if (bS4) {

					//doublereal dspeed = sqrt((potent[VXCOR][iS4]) * (potent[VXCOR][iS4]) + (potent[VYCOR][iS4]) * (potent[VYCOR][iS4]) + (potent[VZCOR][iS4]) * (potent[VZCOR][iS4]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iS4] = 0.0;
						//potent[iTARGET][iS4] = 0.0;
						//potent[iTARGET][iS4] = 0.0;


					//}
					//else {

						potent[iTARGET][iS4] = potent[iTARGET][iP];
						//potent[iTARGET][iS4] = potent[iTARGET][iP];
						//potent[iTARGET][iS4] = potent[iTARGET][iP];


					//}
				}
			}

			if (iT > -1) {
				if (bT) {

					//doublereal dspeed = sqrt((potent[VXCOR][iT]) * (potent[VXCOR][iT]) + (potent[VYCOR][iT]) * (potent[VYCOR][iT]) + (potent[VZCOR][iT]) * (potent[VZCOR][iT]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iT] = 0.0;
						//potent[iTARGET][iT] = 0.0;
						//potent[iTARGET][iT] = 0.0;


					//}
					//else {

						potent[iTARGET][iT] = potent[iTARGET][iP];
						//potent[iTARGET][iT] = potent[iTARGET][iP];
						//potent[iTARGET][iT] = potent[iTARGET][iP];


					//}
				}
			}

			if (iT2 > -1) {
				if (bT2) {

					//doublereal dspeed = sqrt((potent[VXCOR][iT2]) * (potent[VXCOR][iT2]) + (potent[VYCOR][iT2]) * (potent[VYCOR][iT2]) + (potent[VZCOR][iT2]) * (potent[VZCOR][iT2]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iT2] = 0.0;
						//potent[iTARGET][iT2] = 0.0;
						//potent[iTARGET][iT2] = 0.0;


					//}
					//else {

						potent[iTARGET][iT2] = potent[iTARGET][iP];
						//potent[iTARGET][iT2] = potent[iTARGET][iP];
						//potent[iTARGET][iT2] = potent[iTARGET][iP];


					//}
				}
			}

			if (iT3 > -1) {
				if (bT3) {

					//doublereal dspeed = sqrt((potent[VXCOR][iT3]) * (potent[VXCOR][iT3]) + (potent[VYCOR][iT3]) * (potent[VYCOR][iT3]) + (potent[VZCOR][iT3]) * (potent[VZCOR][iT3]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iT3] = 0.0;
						//potent[iTARGET][iT3] = 0.0;
						//potent[iTARGET][iT3] = 0.0;


					//}
					//else {

						potent[iTARGET][iT3] = potent[iTARGET][iP];
						//potent[iTARGET][iT3] = potent[iTARGET][iP];
						//potent[iTARGET][iT3] = potent[iTARGET][iP];


					//}
				}
			}

			if (iT4 > -1) {
				if (bT4) {

					//doublereal dspeed = sqrt((potent[VXCOR][iT4]) * (potent[VXCOR][iT4]) + (potent[VYCOR][iT4]) * (potent[VYCOR][iT4]) + (potent[VZCOR][iT4]) * (potent[VZCOR][iT4]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iT4] = 0.0;
						//potent[iTARGET][iT4] = 0.0;
						//potent[iTARGET][iT4] = 0.0;


					//}
					//else {

						potent[iTARGET][iT4] = potent[iTARGET][iP];
						//potent[iTARGET][iT4] = potent[iTARGET][iP];
						//potent[iTARGET][iT4] = potent[iTARGET][iP];


					//}
				}
			}

			if (iB > -1) {
				if (bB) {
					//doublereal dspeed = sqrt((potent[VXCOR][iB]) * (potent[VXCOR][iB]) + (potent[VYCOR][iB]) * (potent[VYCOR][iB]) + (potent[VZCOR][iB]) * (potent[VZCOR][iB]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iB] = 0.0;
						//potent[iTARGET][iB] = 0.0;
						//potent[iTARGET][iB] = 0.0;


					//}
					//else {

						potent[iTARGET][iB] = potent[iTARGET][iP];
						//potent[iTARGET][iB] = potent[iTARGET][iP];
						//potent[iTARGET][iB] = potent[iTARGET][iP];


					//}
				}
			}

			if (iB2 > -1) {
				if (bB2) {
					//doublereal dspeed = sqrt((potent[VXCOR][iB2]) * (potent[VXCOR][iB2]) + (potent[VYCOR][iB2]) * (potent[VYCOR][iB2]) + (potent[VZCOR][iB2]) * (potent[VZCOR][iB2]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iB2] = 0.0;
						//potent[iTARGET][iB2] = 0.0;
						//potent[iTARGET][iB2] = 0.0;

					//}
					//else {

						potent[iTARGET][iB2] = potent[iTARGET][iP];
						//potent[iTARGET][iB2] = potent[iTARGET][iP];
						//potent[iTARGET][iB2] = potent[iTARGET][iP];

					//}
				}
			}

			if (iB3 > -1) {
				if (bB3) {
					//doublereal dspeed = sqrt((potent[VXCOR][iB3]) * (potent[VXCOR][iB3]) + (potent[VYCOR][iB3]) * (potent[VYCOR][iB3]) + (potent[VZCOR][iB3]) * (potent[VZCOR][iB3]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iB3] = 0.0;
						//potent[iTARGET][iB3] = 0.0;
						//potent[iTARGET][iB3] = 0.0;
					//}
					//else {

						potent[iTARGET][iB3] = potent[iTARGET][iP];
						//potent[iTARGET][iB3] = potent[iTARGET][iP];
						//potent[iTARGET][iB3] = potent[iTARGET][iP];
					//}
				}
			}

			if (iB4 > -1) {
				if (bB4) {
					//doublereal dspeed = sqrt((potent[VXCOR][iB4]) * (potent[VXCOR][iB4]) + (potent[VYCOR][iB4]) * (potent[VYCOR][iB4]) + (potent[VZCOR][iB4]) * (potent[VZCOR][iB4]));

					//if (dspeed < 1.0e-10) {
						//potent[iTARGET][iB4] = 0.0;
						//potent[iTARGET][iB4] = 0.0;
						//potent[iTARGET][iB4] = 0.0;
					//}
					//else {

						potent[iTARGET][iB4] = potent[iTARGET][iP];
						//potent[iTARGET][iB4] = potent[iTARGET][iP];
						//potent[iTARGET][iB4] = potent[iTARGET][iP];

					//}
				}
			}

		}
		else
		{

			// граничные узлы.
			// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.
			
			if (bE) {
			//----->	potent[iTARGET][iE] = potent[iTARGET][iP] + (dxe / dxw) * (potent[iTARGET][iP] - potent[iTARGET][iW]);
				//potent[iTARGET][iE] = potent[iTARGET][iP] + (dxe / dxw) * (potent[iTARGET][iP] - potent[iTARGET][iW]);
				//potent[iTARGET][iE] = potent[iTARGET][iP] + (dxe / dxw) * (potent[iTARGET][iP] - potent[iTARGET][iW]);

				// Интерполяционный полином Лагранжа
				potent[iTARGET][iE] = potent[iTARGET][iP] * ((dxw + dxe) / dxw) - potent[iTARGET][iW] * (dxe / dxw);

			}

			if (bW) {
				//----->potent[iTARGET][iW] = potent[iTARGET][iP] + (dxw / dxe) * (potent[iTARGET][iP] - potent[iTARGET][iE]);
				//potent[iTARGET][iW] = potent[iTARGET][iP] + (dxw / dxe) * (potent[iTARGET][iP] - potent[iTARGET][iE]);
				//potent[iTARGET][iW] = potent[iTARGET][iP] + (dxw / dxe) * (potent[iTARGET][iP] - potent[iTARGET][iE]);

				// Интерполяционный полином Лагранжа
				potent[iTARGET][iW] = potent[iTARGET][iP] * ((dxw + dxe) / dxe) - potent[iTARGET][iE] * (dxw / dxe);

			}

			if (bN) {
				///------>potent[iTARGET][iN] = potent[iTARGET][iP] + (dyn / dys) * (potent[iTARGET][iP] - potent[iTARGET][iS]);
				//potent[iTARGET][iN] = potent[iTARGET][iP] + (dyn / dys) * (potent[iTARGET][iP] - potent[iTARGET][iS]);
				//potent[iTARGET][iN] = potent[iTARGET][iP] + (dyn / dys) * (potent[iTARGET][iP] - potent[iTARGET][iS]);

				// Интерполяционный полином Лагранжа
				potent[iTARGET][iN] = potent[iTARGET][iP] * ((dys + dyn) / dys) - potent[iTARGET][iS] * (dyn / dys);

			}

			if (bS) {
				//---->potent[iTARGET][iS] = potent[iTARGET][iP] + (dys / dyn) * (potent[iTARGET][iP] - potent[iTARGET][iN]);
				//potent[iTARGET][iS] = potent[iTARGET][iP] + (dys / dyn) * (potent[iTARGET][iP] - potent[iTARGET][iN]);
				//potent[iTARGET][iS] = potent[iTARGET][iP] + (dys / dyn) * (potent[iTARGET][iP] - potent[iTARGET][iN]);

				// Интерполяционный полином Лагранжа
				potent[iTARGET][iS] = potent[iTARGET][iP] * ((dys + dyn) / dyn) - potent[iTARGET][iN] * (dys / dyn);

			}

			if (bT) {
				//------->potent[iTARGET][iT] = potent[iTARGET][iP] + (dzt / dzb) * (potent[iTARGET][iP] - potent[iTARGET][iB]);
				//potent[iTARGET][iT] = potent[iTARGET][iP] + (dzt / dzb) * (potent[iTARGET][iP] - potent[iTARGET][iB]);
				//potent[iTARGET][iT] = potent[iTARGET][iP] + (dzt / dzb) * (potent[iTARGET][iP] - potent[iTARGET][iB]);

				// Интерполяционный полином Лагранжа
				potent[iTARGET][iT] = potent[iTARGET][iP] * ((dzb + dzt) / dzb) - potent[iTARGET][iB] * (dzt / dzb);


			}

			if (bB) {
				//--->potent[iTARGET][iB] = potent[iTARGET][iP] + (dzb / dzt) * (potent[iTARGET][iP] - potent[iTARGET][iT]);
				//potent[iTARGET][iB] = potent[iTARGET][iP] + (dzb / dzt) * (potent[iTARGET][iP] - potent[iTARGET][iT]);
				//potent[iTARGET][iB] = potent[iTARGET][iP] + (dzb / dzt) * (potent[iTARGET][iP] - potent[iTARGET][iT]);

				// Интерполяционный полином Лагранжа
				potent[iTARGET][iB] = potent[iTARGET][iP] * ((dzb + dzt) / dzt) - potent[iTARGET][iT] * (dzb / dzt);

			}
			

		}
	}

} // green_gauss_Stress

#endif
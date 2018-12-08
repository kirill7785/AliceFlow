// Файл greengauss.c вычисляет градиенты величин 
// в центрах контрольных объёмов.
// begin 20 октября 2011 г. 
// ссылки: Гаврилов Андрей опыт разработки cfd.
// Основано на теореме Грина-Гаусса.
// begin 15 мая 2012 года. Реализован второй порядок точности.

#ifndef GREEN_GAUSS_C
#define GREEN_GAUSS_C 1


// Вычисление градиентов скоростей в центрах внутренних КО
// и на границах с помощью линейной интерполляции.
// Поскольку интерполляция линейная то точность данной формулы O(h). 
// По поводу точности O(h) спорно, может быть и O(h^2). К тому же было выяснено
// что данный способ вычисления градиентов, для обычной прямоугольной неравномерной сетки
// совпадает со взвешенным методом наименьших квадратов.
void green_gaussO1(integer iP, doublereal** &potent, integer** &nvtx, TOCHKA* &pa,
	ALICE_PARTITION** &sosedi, integer maxelm, bool bbond,
								 doublereal* &mf, doublereal* &prop, doublereal* &prop_b) {


	// Рассчитывать ли скорость на грани с помощью поправки Рхи-Чоу.
	bool bRCh=false;

	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond == true то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

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
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);

	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
    if (!bE) dxe=0.5*(pa[nvtx[1][iE]-1].x+pa[nvtx[0][iE]-1].x);
	if (!bE) dxe-=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw-=0.5*(pa[nvtx[1][iW]-1].x+pa[nvtx[0][iW]-1].x);
    // y - direction
	if (!bN) dyn=0.5*(pa[nvtx[2][iN]-1].y+pa[nvtx[0][iN]-1].y);
	if (!bN) dyn-=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys-=0.5*(pa[nvtx[2][iS]-1].y+pa[nvtx[0][iS]-1].y);
    // z - direction
	if (!bT) dzt=0.5*(pa[nvtx[4][iT]-1].z+pa[nvtx[0][iT]-1].z);
	if (!bT) dzt-=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb-=0.5*(pa[nvtx[4][iB]-1].z+pa[nvtx[0][iB]-1].z);

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

	// плотность на грани КО аппроксимируется средним гармоническим
	doublereal rhoe=1.0, rhow=1.0, rhon=1.0, rhos=1.0, rhot=1.0, rhob=1.0;
	doublereal rP=1.0, rE=1.0, rN=1.0, rT=1.0, rW=1.0, rS=1.0, rB=1.0;
	if (bRCh) {
		rP=prop[iP];
	    if (!bE) rE=prop[iE]; else rE=prop_b[iE-maxelm];
        if (!bN) rN=prop[iN]; else rN=prop_b[iN-maxelm];
        if (!bT) rT=prop[iT]; else rT=prop_b[iT-maxelm];
	    if (!bW) rW=prop[iW]; else rW=prop_b[iW-maxelm];
        if (!bS) rS=prop[iS]; else rS=prop_b[iS-maxelm];
	    if (!bB) rB=prop[iB]; else rB=prop_b[iB-maxelm];
	    
	    rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP); // проверено.
	    rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP);
	    rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP);
	    rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP);
        rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP);
	    rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP);
	}

	// линейная интерполяция скорости VX на грань КО.
    doublereal fe, fw, fn, fs, ft, fb;
	if (!bbond) {
		// внутренние КО.

		// Линейно интерполлируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		// VX
	    if (!bE) {
			if (bRCh) {
				fe=mf[ESIDE]/(rhoe*dy*dz); // скорость на грани с учётом поправки Рхи-Чоу.
			}
			else {
				fe=feplus*potent[VXCOR][iE]+(1.0-feplus)*potent[VXCOR][iP]; 
			}
		}
		else fe=potent[VXCOR][iE];
        if (!bW) {
			if (bRCh) {
				fw=mf[WSIDE]/(rhow*dy*dz); // скорость на грани с учётом монотонизирующей поправки Рхи-Чоу.
			}
			else {
				fw=fwplus*potent[VXCOR][iW]+(1.0-fwplus)*potent[VXCOR][iP];
			}
		}
		else fw=potent[VXCOR][iW];
		// Эти компоненты скорости тоже по идее можно вычислять с помощью монотонизирующей поправки.
		// Вопрос о правомерности пока остаётся открытым. Дальнешие компоненты скорости и производные аналогично для VX.
	    if (!bN) fn=fnplus*potent[VXCOR][iN]+(1.0-fnplus)*potent[VXCOR][iP]; else fn=potent[VXCOR][iN];
        if (!bS) fs=fsplus*potent[VXCOR][iS]+(1.0-fsplus)*potent[VXCOR][iP]; else fs=potent[VXCOR][iS];
        if (!bT) ft=ftplus*potent[VXCOR][iT]+(1.0-ftplus)*potent[VXCOR][iP]; else ft=potent[VXCOR][iT];
        if (!bB) fb=fbplus*potent[VXCOR][iB]+(1.0-fbplus)*potent[VXCOR][iP]; else fb=potent[VXCOR][iB];
        // градиент VX
	    potent[GRADXVX][iP]=(fe-fw)/dx;
	    potent[GRADYVX][iP]=(fn-fs)/dy;
	    potent[GRADZVX][iP]=(ft-fb)/dz;

		// линейная интерполяция скорости VY на грань КО.
	    
	    if (!bE) fe=feplus*potent[VYCOR][iE]+(1.0-feplus)*potent[VYCOR][iP]; else fe=potent[VYCOR][iE];
        if (!bW) fw=fwplus*potent[VYCOR][iW]+(1.0-fwplus)*potent[VYCOR][iP]; else fw=potent[VYCOR][iW];
	    if (!bN) {
			if (bRCh) {
				fn=mf[NSIDE]/(rhon*dx*dz);
			}
			else {
				fn=fnplus*potent[VYCOR][iN]+(1.0-fnplus)*potent[VYCOR][iP]; 
			}
		} 
		else fn=potent[VYCOR][iN];
        if (!bS) {
			if (bRCh) {
				fs=mf[SSIDE]/(rhos*dx*dz);
			}
			else {
				fs=fsplus*potent[VYCOR][iS]+(1.0-fsplus)*potent[VYCOR][iP]; 
			}
		}
		else fs=potent[VYCOR][iS];
        if (!bT) ft=ftplus*potent[VYCOR][iT]+(1.0-ftplus)*potent[VYCOR][iP]; else ft=potent[VYCOR][iT];
        if (!bB) fb=fbplus*potent[VYCOR][iB]+(1.0-fbplus)*potent[VYCOR][iP]; else fb=potent[VYCOR][iB];
		
		// градиент VY

	    potent[GRADXVY][iP]=(fe-fw)/dx;
	    potent[GRADYVY][iP]=(fn-fs)/dy;
	    potent[GRADZVY][iP]=(ft-fb)/dz;

	    // линейная интерполяция скорости VZ на грань КО.
    
	    if (!bE) fe=feplus*potent[VZCOR][iE]+(1.0-feplus)*potent[VZCOR][iP]; else fe=potent[VZCOR][iE];
        if (!bW) fw=fwplus*potent[VZCOR][iW]+(1.0-fwplus)*potent[VZCOR][iP]; else fw=potent[VZCOR][iW];
	    if (!bN) fn=fnplus*potent[VZCOR][iN]+(1.0-fnplus)*potent[VZCOR][iP]; else fn=potent[VZCOR][iN];
        if (!bS) fs=fsplus*potent[VZCOR][iS]+(1.0-fsplus)*potent[VZCOR][iP]; else fs=potent[VZCOR][iS];
        if (!bT) {
			if (bRCh) {
				ft=mf[TSIDE]/(rhot*dx*dy); 
			}
			else {
               ft=ftplus*potent[VZCOR][iT]+(1.0-ftplus)*potent[VZCOR][iP]; 
			}
		} else ft=potent[VZCOR][iT];
        if (!bB) {
			if (bRCh) {
				fb=mf[BSIDE]/(rhob*dx*dy);
			}
			else {
				fb=fbplus*potent[VZCOR][iB]+(1.0-fbplus)*potent[VZCOR][iP];
			}
		} else fb=potent[VZCOR][iB];
		// градиент VZ
	    potent[GRADXVZ][iP]=(fe-fw)/dx;
	    potent[GRADYVZ][iP]=(fn-fs)/dy;
	    potent[GRADZVZ][iP]=(ft-fb)/dz;
	}
	else {

		if (1) {

			

			// По простому : градиент на границе наследуем из ближайшего внутреннего узла.

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

		if (bW) {

			doublereal dspeed = sqrt((potent[VXCOR][iW])*(potent[VXCOR][iW]) + (potent[VYCOR][iW])*(potent[VYCOR][iW]) + (potent[VZCOR][iW])*(potent[VZCOR][iW]));

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

		if (bN) {

			doublereal dspeed = sqrt((potent[VXCOR][iN])*(potent[VXCOR][iN]) + (potent[VYCOR][iN])*(potent[VYCOR][iN]) + (potent[VZCOR][iN])*(potent[VZCOR][iN]));

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

		if (bS) {

			doublereal dspeed = sqrt((potent[VXCOR][iS])*(potent[VXCOR][iS]) + (potent[VYCOR][iS])*(potent[VYCOR][iS]) + (potent[VZCOR][iS])*(potent[VZCOR][iS]));

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

		if (bT) {

			doublereal dspeed = sqrt((potent[VXCOR][iT])*(potent[VXCOR][iT]) + (potent[VYCOR][iT])*(potent[VYCOR][iT]) + (potent[VZCOR][iT])*(potent[VZCOR][iT]));

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

		if (bB) {
			doublereal dspeed = sqrt((potent[VXCOR][iB])*(potent[VXCOR][iB]) + (potent[VYCOR][iB])*(potent[VYCOR][iB]) + (potent[VZCOR][iB])*(potent[VZCOR][iB]));

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
		else
		{

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполляции.

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

// вычисление градиентов поправки давления с помощью теоремы Грина-Гаусса. 
void green_gaussPAM(integer iP, doublereal** &potent, integer** &nvtx, TOCHKA* &pa,
	                ALICE_PARTITION** &sosedi, integer maxelm, bool bbond,
					BOUND* &sosedb, integer ls, integer lw, WALL* &w, bool bLRfree,
	                integer *ilevel_alice, integer* ptr) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты поправки давления для внутренних КО.
	// если bbond == true то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1;
	iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

	integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
	iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
	iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

	integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
	iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
	iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

	integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
	iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
	iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;

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
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);

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
			dSqe = sosedb[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}


		
	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = sosedb[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}

		
	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = sosedb[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}

		
	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = sosedb[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}

		
	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = sosedb[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}

		
	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = sosedb[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			dSqe2 = sosedb[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}

		
	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw2 = sosedb[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}

		
	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = sosedb[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}

		
	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = sosedb[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}

		
	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = sosedb[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}

		
	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = sosedb[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			dSqe3 = sosedb[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}

		
	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = sosedb[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}

		
	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = sosedb[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}

		
	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = sosedb[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}

		
	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = sosedb[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}

		
	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = sosedb[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			dSqe4 = sosedb[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}

		
	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = sosedb[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}

		
	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = sosedb[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}

		
	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = sosedb[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}

		
	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = sosedb[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}

		
	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = sosedb[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}

		
	}

	doublereal PAMe = 0.0, PAMw = 0.0, PAMn = 0.0, PAMs = 0.0, PAMt = 0.0, PAMb = 0.0;
	doublereal PAMe2 = 0.0, PAMw2 = 0.0, PAMn2 = 0.0, PAMs2 = 0.0, PAMt2 = 0.0, PAMb2 = 0.0;
	doublereal PAMe3 = 0.0, PAMw3 = 0.0, PAMn3 = 0.0, PAMs3 = 0.0, PAMt3 = 0.0, PAMb3 = 0.0;
	doublereal PAMe4 = 0.0, PAMw4 = 0.0, PAMn4 = 0.0, PAMs4 = 0.0, PAMt4 = 0.0, PAMb4 = 0.0;

    if (!bbond) {

		if (bLRfree) {

			// не работает на АЛИС.
			if (b_on_adaptive_local_refinement_mesh) {
				printf("function green_gaussPAM in module greengauss.c bLRfree not worked in ALICE mesh...\n ");
				getchar();
				exit(1);
			}

			// Так гораздо хуже на opening тесте.

			// В случае bLRFree мы будем линейно интерполировать давление в граничные узлы чтобы сохранить значение градиента.
			// градиент давления очень важен для естественно конвективных течний, т.к. они протекают под действием этого градиента
			// и это для них особенно важно вблизи твёрдых стенок.
			
			if (iE > -1) {
				if (!bE) PAMe = feplus * potent[PAM][iE] + (1.0 - feplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					PAMe = potent[PAM][iP] + (dxe / dxw)*(potent[PAM][iP] - potent[PAM][iW]);
				}
			}

			if (iW > -1) {
				if (!bW) PAMw = fwplus * potent[PAM][iW] + (1.0 - fwplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					PAMw = potent[PAM][iP] + (dxw / dxe)*(potent[PAM][iP] - potent[PAM][iE]);
				}
			}

			if (iN > -1) {
				if (!bN) PAMn = fnplus * potent[PAM][iN] + (1.0 - fnplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					PAMn = potent[PAM][iP] + (dyn / dys)*(potent[PAM][iP] - potent[PAM][iS]);
				}
			}

			if (iS > -1) {
				if (!bS) PAMs = fsplus * potent[PAM][iS] + (1.0 - fsplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					PAMs = potent[PAM][iP] + (dys / dyn)*(potent[PAM][iP] - potent[PAM][iN]);
				}
			}

			if (iT > -1) {
				if (!bT) PAMt = ftplus * potent[PAM][iT] + (1.0 - ftplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					PAMt = potent[PAM][iP] + (dzt / dzb)*(potent[PAM][iP] - potent[PAM][iB]);
				}
			}

			if (iB > -1) {
				if (!bB) PAMb = fbplus * potent[PAM][iB] + (1.0 - fbplus)*potent[PAM][iP]; else {
					// линейная интерполяция давления на граничный узел !!!
					PAMb = potent[PAM][iP] + (dzb / dzt)*(potent[PAM][iP] - potent[PAM][iT]);
				}
			}

			

             // градиент Давления.
	         potent[GRADXPAM][iP]=(PAMe-PAMw)/dx;
	         potent[GRADYPAM][iP]=(PAMn-PAMs)/dy;
	         potent[GRADZPAM][iP]=(PAMt-PAMb)/dz;

		}
		else {


		    // внутренние КО.

		    // Линейно интерполлируем поправку давления на грань контрольного объёма,
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
	    integer interpol=0; // 0 или 1 при линейной интерполяции.

		if (interpol==0) {

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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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

						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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

				if (iE2 > -1) {
					if (bE2) {
						integer inumber = iE2 - maxelm;
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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

						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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

						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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

						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
			else {

				if (iE > -1) {
					if (bE) {
						integer inumber = iE - maxelm;
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iE] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iW] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iN] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iS] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iT] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iB] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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

				if (iE2 > -1) {
					if (bE2) {
						integer inumber = iE2 - maxelm;
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iE2] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iW2] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iN2] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iS2] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iT2] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iB2] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iE3] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iW3] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iN3] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iS3] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iT3] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iB3] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iE4] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADXPAM][iW4] = potent[GRADXPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iN4] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADYPAM][iS4] = potent[GRADYPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iT4] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
						if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
							potent[GRADZPAM][iB4] = potent[GRADZPAM][iP];
						}
						else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
		else if (interpol==1) {

			// не работает на АЛИС.
			if (b_on_adaptive_local_refinement_mesh) {
				printf("function green_gaussPAM in module greengauss.c else if (interpol==1) { not worked in ALICE mesh...\n ");
				getchar();
				exit(1);
			}

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполляции.

		// Если строка с пометкой <-- раскоментирована то градиент поправки давления линейно интерполлируется на границу 
		// расчётной области изнутри расчётной области.

		if (bE) {

			integer inumber=iE-maxelm;
			if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
				//potent[GRADXPAM][iE]=potent[GRADXPAM][iP]+(dxe/dxw)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iW]);
				potent[GRADXPAM][iE]=0.0;
			}
			else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
			if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
				//potent[GRADXPAM][iW]=potent[GRADXPAM][iP]+(dxw/dxe)*(potent[GRADXPAM][iP]-potent[GRADXPAM][iE]);
				potent[GRADXPAM][iW]=0.0;
			}
			else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
					//potent[GRADYPAM][iN]=potent[GRADYPAM][iP]+(dyn/dys)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iS]);
					potent[GRADYPAM][iN]=0.0;
				}
				else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
			if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
				//potent[GRADYPAM][iS]=potent[GRADYPAM][iP]+(dys/dyn)*(potent[GRADYPAM][iP]-potent[GRADYPAM][iN]);
				potent[GRADYPAM][iS]=0.0;
			}
			else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
			if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
				//potent[GRADZPAM][iT]=potent[GRADZPAM][iP]+(dzt/dzb)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iB]);
				potent[GRADZPAM][iT]=0.0;
			}
			else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
			if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
				//potent[GRADZPAM][iB]=potent[GRADZPAM][iP]+(dzb/dzt)*(potent[GRADZPAM][iP]-potent[GRADZPAM][iT]);
				potent[GRADZPAM][iB]=0.0;
			}
			else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
	}


} // green_gaussPAM

// 13 апреля 2015 года.
// вычисление градиентов Температуры с помощью теоремы Грина-Гаусса. 
void green_gaussTemperature(integer iP, doublereal* &potent, integer** &nvtx, TOCHKA* &pa,
	ALICE_PARTITION** &sosedi, integer maxelm, bool bbond,
					BOUND* &sosedb, doublereal* &Tx, doublereal* &Ty, doublereal* &Tz) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты поправки давления для внутренних КО.
	// если bbond == true то будут вычислены значения только в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// Внимание ! память под Tx, Ty, Tz предполагается выделенной заранее.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

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
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);

	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
    if (!bE) dxe=0.5*(pa[nvtx[1][iE]-1].x+pa[nvtx[0][iE]-1].x);
	if (!bE) dxe-=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw-=0.5*(pa[nvtx[1][iW]-1].x+pa[nvtx[0][iW]-1].x);
    // y - direction
	if (!bN) dyn=0.5*(pa[nvtx[2][iN]-1].y+pa[nvtx[0][iN]-1].y);
	if (!bN) dyn-=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys-=0.5*(pa[nvtx[2][iS]-1].y+pa[nvtx[0][iS]-1].y);
    // z - direction
	if (!bT) dzt=0.5*(pa[nvtx[4][iT]-1].z+pa[nvtx[0][iT]-1].z);
	if (!bT) dzt-=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb-=0.5*(pa[nvtx[4][iB]-1].z+pa[nvtx[0][iB]-1].z);

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

	doublereal Te, Tw, Tn, Ts, Tt, Tb;

    if (!bbond) {
		// внутренние КО.

		// Линейно интерполлируем поправку давления на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		if (!bE) Te=feplus*potent[iE]+(1.0-feplus)*potent[iP]; else Te=potent[iE];
        if (!bW) Tw=fwplus*potent[iW]+(1.0-fwplus)*potent[iP]; else Tw=potent[iW];
	    if (!bN) Tn=fnplus*potent[iN]+(1.0-fnplus)*potent[iP]; else Tn=potent[iN];
        if (!bS) Ts=fsplus*potent[iS]+(1.0-fsplus)*potent[iP]; else Ts=potent[iS];
        if (!bT) Tt=ftplus*potent[iT]+(1.0-ftplus)*potent[iP]; else Tt=potent[iT];
        if (!bB) Tb=fbplus*potent[iB]+(1.0-fbplus)*potent[iP]; else Tb=potent[iB];
        // градиент VX
	    Tx[iP]=(Te-Tw)/dx;
	    Ty[iP]=(Tn-Ts)/dy;
	    Tz[iP]=(Tt-Tb)/dz;

	}
   else {
	    integer interpol=0; // 0 или 1 при линейной интерполяции.

		if (interpol==0) {


				if (bE) {
				
			
				Tx[iE]=Tx[iP];
			    Ty[iE]=Ty[iP];
			    Tz[iE]=Tz[iP];
		    }

		if (bW) {
			

			Tx[iW]=Tx[iP];
			Ty[iW]=Ty[iP];
			Tz[iW]=Tz[iP];
		}

		if (bN) {

			Tx[iN]=Tx[iP];
			Ty[iN]=Ty[iP];
			Tz[iN]=Tz[iP];
		}

		if (bS) {
			Tx[iS]=Tx[iP];
			Ty[iS]=Ty[iP];
			Tz[iS]=Tz[iP];
		}

		if (bT) {

			

			Tx[iT]=Tx[iP];
			Ty[iT]=Ty[iP];
			Tz[iT]=Tz[iP];
		}

		if (bB) {

			

			Tx[iB]=Tx[iP];
			Ty[iB]=Ty[iP];
			Tz[iB]=Tz[iP];
		}


			

		}
		else if (interpol==1) {

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполляции.

		// Если строка с пометкой <-- раскоментирована то градиент поправки давления линейно интерполлируется на границу 
		// расчётной области изнутри расчётной области.

		if (bE) {
			

			Ty[iE]=Ty[iP]+(dxe/dxw)*(Ty[iP]-Ty[iW]);
			Tz[iE]=Tz[iP]+(dxe/dxw)*(Tz[iP]-Tz[iW]);
			Tx[iE]=Tx[iP]+(dxe/dxw)*(Tx[iP]-Tx[iW]); // <--

		}

		if (bW) {
			

			Ty[iW]=Ty[iP]+(dxw/dxe)*(Ty[iP]-Ty[iE]);
			Tz[iW]=Tz[iP]+(dxw/dxe)*(Tz[iP]-Tz[iE]);
			Tx[iW]=Tx[iP]+(dxw/dxe)*(Tx[iP]-Tx[iE]); // <--

		}

		if (bN) {

			

			Tx[iN]=Tx[iP]+(dyn/dys)*(Tx[iP]-Tx[iS]);
			Tz[iN]=Tz[iP]+(dyn/dys)*(Tz[iP]-Tz[iS]);
			Ty[iN]=Ty[iP]+(dyn/dys)*(Ty[iP]-Ty[iS]); // <--
		}

		if (bS) {

			

			Tx[iS]=Tx[iP]+(dys/dyn)*(Tx[iP]-Tx[iN]);
			Tz[iS]=Tz[iP]+(dys/dyn)*(Tz[iP]-Tz[iN]);
			Ty[iS]=Ty[iP]+(dys/dyn)*(Ty[iP]-Ty[iN]); //<--

		}

		if (bT) {

			

			Tx[iT]=Tx[iP]+(dzt/dzb)*(Tx[iP]-Tx[iB]);
			Ty[iT]=Ty[iP]+(dzt/dzb)*(Ty[iP]-Ty[iB]);
			Tz[iT]=Tz[iP]+(dzt/dzb)*(Tz[iP]-Tz[iB]); // <--
		}

		if (bB) {

			

			Tx[iB]=Tx[iP]+(dzb/dzt)*(Tx[iP]-Tx[iT]);
			Ty[iB]=Ty[iP]+(dzb/dzt)*(Ty[iP]-Ty[iT]);
			Tz[iB]=Tz[iP]+(dzb/dzt)*(Tz[iP]-Tz[iT]); // <--
		}
		}
	}


} // green_gaussTemperature



// Минимакстное сглаживание градиента поправки давления.
// помоему это весьма неудачная идея.
void green_gaussPAMminmax(doublereal** &potent, ALICE_PARTITION** &sosedi, integer maxelm, integer maxbound) {
	

	doublereal** minmaxlimitergrad=new doublereal*[3];
	for (integer i=0; i<3; i++) {
		minmaxlimitergrad[i]=new doublereal[maxelm+maxbound];
	}

	for (integer iP=0; iP<maxelm; iP++) {
	   // iP - номер внутреннего контрольного объёма
	   // iP изменяется от 0 до maxelm-1.
	   integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	   iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

	   // минимаксное ограничение против появления ложных максимумов.
	   minmaxlimitergrad[VX][iP]=fmax(fmin(fmax(potent[GRADXPAM][iE],potent[GRADXPAM][iW]),potent[GRADXPAM][iP]),fmin(potent[GRADXPAM][iE],potent[GRADXPAM][iW]));
	   minmaxlimitergrad[VY][iP]=fmax(fmin(fmax(potent[GRADYPAM][iN],potent[GRADYPAM][iS]),potent[GRADYPAM][iP]),fmin(potent[GRADYPAM][iN],potent[GRADYPAM][iS]));
	   minmaxlimitergrad[VZ][iP]=fmax(fmin(fmax(potent[GRADZPAM][iT],potent[GRADZPAM][iB]),potent[GRADZPAM][iP]),fmin(potent[GRADZPAM][iT],potent[GRADZPAM][iB]));
	}

	for (integer iP=0; iP<maxelm; iP++) {

		// iP - номер внутреннего контрольного объёма
	    // iP изменяется от 0 до maxelm-1.
	    integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1; iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;

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
			minmaxlimitergrad[VX][iE]=minmaxlimitergrad[VX][iP];
			minmaxlimitergrad[VY][iE]=minmaxlimitergrad[VY][iP];
			minmaxlimitergrad[VZ][iE]=minmaxlimitergrad[VZ][iP];
		}

		if (bW) {
			minmaxlimitergrad[VX][iW]=minmaxlimitergrad[VX][iP];
			minmaxlimitergrad[VY][iW]=minmaxlimitergrad[VY][iP];
			minmaxlimitergrad[VZ][iW]=minmaxlimitergrad[VZ][iP];
		}

		if (bN) {
			minmaxlimitergrad[VX][iN]=minmaxlimitergrad[VX][iP];
			minmaxlimitergrad[VY][iN]=minmaxlimitergrad[VY][iP];
			minmaxlimitergrad[VZ][iN]=minmaxlimitergrad[VZ][iP];
		}

		if (bS) {
			minmaxlimitergrad[VX][iS]=minmaxlimitergrad[VX][iP];
			minmaxlimitergrad[VY][iS]=minmaxlimitergrad[VY][iP];
			minmaxlimitergrad[VZ][iS]=minmaxlimitergrad[VZ][iP];
		}

		if (bT) {
			minmaxlimitergrad[VX][iT]=minmaxlimitergrad[VX][iP];
			minmaxlimitergrad[VY][iT]=minmaxlimitergrad[VY][iP];
			minmaxlimitergrad[VZ][iT]=minmaxlimitergrad[VZ][iP];
		}

		if (bB) {
			minmaxlimitergrad[VX][iB]=minmaxlimitergrad[VX][iP];
			minmaxlimitergrad[VY][iB]=minmaxlimitergrad[VY][iP];
			minmaxlimitergrad[VZ][iB]=minmaxlimitergrad[VZ][iP];
		}
	}

	// обратное копирование.
	for (integer iP=0; iP<maxelm+maxbound; iP++) {
		potent[GRADXPAM][iP]=minmaxlimitergrad[VX][iP];
		potent[GRADYPAM][iP]=minmaxlimitergrad[VY][iP];
		potent[GRADZPAM][iP]=minmaxlimitergrad[VZ][iP];
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
void green_gaussPRESS(integer iP, doublereal** &potent, integer** &nvtx, TOCHKA* &pa,
	ALICE_PARTITION** &sosedi, integer maxelm, bool bbond,
					BOUND* &sosedb, integer ls, integer lw, WALL* &w, bool bLRfree,
	integer *ilevel_alice, integer* ptr) {

	// maxelm - число внутренних КО.
	// Вычисляет градиенты давления для внутренних КО.
	// если bbond == true то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1;
	iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

	integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
	iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
	iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

	integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
	iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
	iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

	integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
	iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
	iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;

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
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);

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
			dSqe = sosedb[iE - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE]]) {
				dSqe = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe = dy_loc * dz_loc;
			}
		}



	}


	if (iW > -1) {

		dSqw = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw = sosedb[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				dSqw = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw = dy_loc * dz_loc;
			}
		}


	}


	if (iN > -1) {

		dSqn = dx * dz;

		if (bN) {
			// граничный узел.
			dSqn = sosedb[iN - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
				dSqn = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn = dx_loc * dz_loc;
			}
		}


	}


	if (iS > -1) {

		dSqs = dx * dz;

		if (bS) {
			// граничный узел.
			dSqs = sosedb[iS - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
				dSqs = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs = dx_loc * dz_loc;
			}
		}


	}


	if (iT > -1) {

		dSqt = dx * dy;

		if (bT) {
			// граничный узел.
			dSqt = sosedb[iT - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
				dSqt = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt = dx_loc * dy_loc;
			}
		}


	}


	if (iB > -1) {

		dSqb = dx * dy;

		if (bB) {
			// граничный узел.
			dSqb = sosedb[iB - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
				dSqb = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			dSqe2 = sosedb[iE2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
				dSqe2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe2 = dy_loc * dz_loc;
			}
		}


	}


	if (iW2 > -1) {
		dSqw2 = dy * dz;

		if (bW) {
			// граничный узел.
			dSqw2 = sosedb[iW - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
				dSqw2 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw2 = dy_loc * dz_loc;
			}
		}


	}


	if (iN2 > -1) {

		dSqn2 = dx * dz;

		if (bN2) {
			// граничный узел.
			dSqn2 = sosedb[iN2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN2]]) {
				dSqn2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn2 = dx_loc * dz_loc;
			}
		}


	}


	if (iS2 > -1) {

		dSqs2 = dx * dz;

		if (bS2) {
			// граничный узел.
			dSqs2 = sosedb[iS2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
				dSqs2 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs2 = dx_loc * dz_loc;
			}
		}


	}


	if (iT2 > -1) {

		dSqt2 = dx * dy;

		if (bT2) {
			// граничный узел.
			dSqt2 = sosedb[iT2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
				dSqt2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt2 = dx_loc * dy_loc;
			}
		}


	}


	if (iB2 > -1) {

		dSqb2 = dx * dy;

		if (bB2) {
			// граничный узел.
			dSqb2 = sosedb[iB2 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
				dSqb2 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			dSqe3 = sosedb[iE3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
				dSqe3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe3 = dy_loc * dz_loc;
			}
		}


	}


	if (iW3 > -1) {

		dSqw3 = dy * dz;

		if (bW3) {
			// граничный узел.
			dSqw3 = sosedb[iW3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
				dSqw3 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw3 = dy_loc * dz_loc;
			}
		}


	}


	if (iN3 > -1) {

		dSqn3 = dx * dz;

		if (bN3) {
			// граничный узел.
			dSqn3 = sosedb[iN3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
				dSqn3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn3 = dx_loc * dz_loc;
			}
		}


	}


	if (iS3 > -1) {

		dSqs3 = dx * dz;

		if (bS3) {
			// граничный узел.
			dSqs3 = sosedb[iS3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
				dSqs3 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs3 = dx_loc * dz_loc;
			}
		}


	}


	if (iT3 > -1) {

		dSqt3 = dx * dy;

		if (bT3) {
			// граничный узел.
			dSqt3 = sosedb[iT3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
				dSqt3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt3 = dx_loc * dy_loc;
			}
		}


	}


	if (iB3 > -1) {

		dSqb3 = dx * dy;

		if (bB3) {
			// граничный узел.
			dSqb3 = sosedb[iB3 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
				dSqb3 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			dSqe4 = sosedb[iE4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
				dSqe4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqe4 = dy_loc * dz_loc;
			}
		}


	}


	if (iW4 > -1) {

		dSqw4 = dy * dz;

		if (bW4) {
			// граничный узел.
			dSqw4 = sosedb[iW4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
				dSqw4 = dy * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqw4 = dy_loc * dz_loc;
			}
		}


	}


	if (iN4 > -1) {

		dSqn4 = dx * dz;

		if (bN4) {
			// граничный узел.
			dSqn4 = sosedb[iN4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
				dSqn4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqn4 = dx_loc * dz_loc;
			}
		}


	}


	if (iS4 > -1) {

		dSqs4 = dx * dz;

		if (bS4) {
			// граничный узел.
			dSqs4 = sosedb[iS4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
				dSqs4 = dx * dz;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqs4 = dx_loc * dz_loc;
			}
		}


	}


	if (iT4 > -1) {

		dSqt4 = dx * dy;

		if (bT4) {
			// граничный узел.
			dSqt4 = sosedb[iT4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
				dSqt4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqt4 = dx_loc * dy_loc;
			}
		}


	}


	if (iB4 > -1) {

		dSqb4 = dx * dy;

		if (bB4) {
			// граничный узел.
			dSqb4 = sosedb[iB4 - maxelm].dS;
		}
		else {
			if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
				dSqb4 = dx * dy;
			}
			else {
				// вычисление размеров соседнего контрольного объёма:
				doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
				volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

				dSqb4 = dx_loc * dy_loc;
			}
		}


	}

	doublereal PRESSe=0.0, PRESSw=0.0, PRESSn=0.0, PRESSs=0.0, PRESSt=0.0, PRESSb=0.0;
	doublereal PRESSe2 = 0.0, PRESSw2 = 0.0, PRESSn2 = 0.0, PRESSs2 = 0.0, PRESSt2 = 0.0, PRESSb2 = 0.0;
	doublereal PRESSe3 = 0.0, PRESSw3 = 0.0, PRESSn3 = 0.0, PRESSs3 = 0.0, PRESSt3 = 0.0, PRESSb3 = 0.0;
	doublereal PRESSe4 = 0.0, PRESSw4 = 0.0, PRESSn4 = 0.0, PRESSs4 = 0.0, PRESSt4 = 0.0, PRESSb4 = 0.0;

    if (!bbond) {
		// внутренние КО.

		// Линейно интерполлируем скорости на грань контрольного объёма,
		// а затем вычисляет производную в центре контрольного объёма по обычной конечно разностной формуле. 

		if (bLRfree) {

			// не работает на АЛИС.
			if (b_on_adaptive_local_refinement_mesh) {
				printf("function green_gaussPRESS in module greengauss.c bLRfree not worked in ALICE mesh...\n ");
				getchar();
				exit(1);
			}


			// В случае bLRFree мы будем линейно интерполировать давление в граничные узлы чтобы сохранить значение градиента.
			// градиент давления очень важен для естественно конвективных течний, т.к. они протекают под действием этого градиента
			// и это для них особенно важно вблизи твёрдых стенок.

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
             // градиент Давления.
	         potent[GRADXPRESS][iP]=(PRESSe-PRESSw)/dx;
	         potent[GRADYPRESS][iP]=(PRESSn-PRESSs)/dy;
	         potent[GRADZPRESS][iP]=(PRESSt-PRESSb)/dz;

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

	   integer interpol=0; // 0 или 1 при линейной интерполяции.

	   if (interpol==0) {

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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE2] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW2] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN2] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS2] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT2] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB2] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE3] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW3] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN3] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS3] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT3] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB3] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iE4] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADXPRESS][iW4] = potent[GRADXPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iN4] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADYPRESS][iS4] = potent[GRADYPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iT4] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
					   if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bpressure)) {
						   potent[GRADZPRESS][iB4] = potent[GRADZPRESS][iP];
					   }
					   else if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
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
	   else if (interpol==1) {

        
		   // не работает на АЛИС.
		   if (b_on_adaptive_local_refinement_mesh) {
			   printf("function green_gaussPRESS in module greengauss.c else if (interpol==1) { not worked in ALICE mesh...\n ");
			   getchar();
			   exit(1);
		   }

		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполляции.

		if (bE) {

			 integer inumber=iE-maxelm;
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
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
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
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
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
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
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
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
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
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
				if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure)) {
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
	}


} // green_gaussPRESS

// Вычисление градиентов скоростей в центрах внутренних КО
// и на границах с помощью квадратичной интерполляции.
// Поскольку интерполляция квадратичная то точность данной формулы O(h^2).
// данная функция реализована 15 мая 2012 года.
void green_gaussO2(integer iP, doublereal** &potent, integer** &nvtx, TOCHKA* &pa,
	ALICE_PARTITION** &sosedi, integer maxelm, bool bbond) {
	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond == true то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

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
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);


	if (!bbond) {
		// внутренние КО.

		// Воспользуемся формулой из книги Г.З. Гарбера обеспечивающую второй порядок точности.


		doublereal hxminus, hxplus, hyminus, hyplus, hzminus, hzplus;
		// компоненты сорости.
		doublereal VXP, VXW, VXE, VXN, VXS, VXT, VXB;
		doublereal VYP, VYW, VYE, VYN, VYS, VYT, VYB;
		doublereal VZP, VZW, VZE, VZN, VZS, VZT, VZB;

		VXP=potent[VX][iP]; VYP=potent[VY][iP]; VZP=potent[VZ][iP]; // значения скоростей в центральном контрольном объёме.
		VXW=potent[VX][iW]; VYW=potent[VY][iW]; VZW=potent[VZ][iW];
		VXE=potent[VX][iE]; VYE=potent[VY][iE]; VZE=potent[VZ][iE];
		VXS=potent[VX][iS]; VYS=potent[VY][iS]; VZS=potent[VZ][iS]; 
		VXN=potent[VX][iN]; VYN=potent[VY][iN]; VZN=potent[VZ][iN]; 
		VXB=potent[VX][iB]; VYB=potent[VY][iB]; VZB=potent[VZ][iB];  
		VXT=potent[VX][iT]; VYT=potent[VY][iT]; VZT=potent[VZ][iT];

		if (!bW) {
			TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iW, nvtx, pa, pb,WSIDE);
		    hxminus=fabs(pp.x-pb.x);		    
	    } else {
		    // узел W граничный
			hxminus=0.5*dx;
	    }
	    if (!bE) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iE, nvtx, pa, pb,ESIDE);
		    hxplus=fabs(pb.x-pp.x);		
	    } else {
            // узел E граничный
		    hxplus=0.5*dx;
       	}
	    if (!bS) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iS, nvtx, pa, pb,SSIDE);
		    hyminus=fabs(pp.y-pb.y);
	    } else {
		    // узел S граничный
		    hyminus=0.5*dy;
		}
	    if (!bN) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iN, nvtx, pa, pb,NSIDE);
		    hyplus=fabs(pb.y-pp.y);		
	    } else {
		    // узел N граничный
		    hyplus=0.5*dy;
	    } 
		if (!bB) {
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iB, nvtx, pa, pb,BSIDE);
		    hzminus=fabs(pp.z-pb.z);
		} else {
		    // узел B граничный
		    hzminus=0.5*dz;
	    }
	    if (!bT) { 
		    TOCHKA pp, pb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iT, nvtx, pa, pb,TSIDE);
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
		// градиенты в граничных узлах восстанавливаются с помощью квадратичной интерполляции.

		if (bE) {

			// квадратичная интерполляция.

			TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iW, nvtx, pa, pb,WSIDE);
			center_cord3D(sosedi[WSIDE][iW].iNODE1, nvtx, pa, pbb,WW);
					
			potent[GRADXVX][iE] = my_quadratic_interpolation('+', potent[GRADXVX][sosedi[WSIDE][iW].iNODE1], potent[GRADXVX][iW], potent[GRADXVX][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADYVX][iE] = my_quadratic_interpolation('+', potent[GRADYVX][sosedi[WSIDE][iW].iNODE1], potent[GRADYVX][iW], potent[GRADYVX][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADZVX][iE] = my_quadratic_interpolation('+', potent[GRADZVX][sosedi[WSIDE][iW].iNODE1], potent[GRADZVX][iW], potent[GRADZVX][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);

			potent[GRADXVY][iE] = my_quadratic_interpolation('+', potent[GRADXVY][sosedi[WSIDE][iW].iNODE1], potent[GRADXVY][iW], potent[GRADXVY][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADYVY][iE] = my_quadratic_interpolation('+', potent[GRADYVY][sosedi[WSIDE][iW].iNODE1], potent[GRADYVY][iW], potent[GRADYVY][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADZVY][iE] = my_quadratic_interpolation('+', potent[GRADZVY][sosedi[WSIDE][iW].iNODE1], potent[GRADZVY][iW], potent[GRADZVY][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);

			potent[GRADXVZ][iE] = my_quadratic_interpolation('+', potent[GRADXVZ][sosedi[WSIDE][iW].iNODE1], potent[GRADXVZ][iW], potent[GRADXVZ][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADYVZ][iE] = my_quadratic_interpolation('+', potent[GRADYVZ][sosedi[WSIDE][iW].iNODE1], potent[GRADYVZ][iW], potent[GRADYVZ][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			potent[GRADZVZ][iE] = my_quadratic_interpolation('+', potent[GRADZVZ][sosedi[WSIDE][iW].iNODE1], potent[GRADZVZ][iW], potent[GRADZVZ][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
		}

		if (bW) {

			// квадратичная интерполляция.

			TOCHKA pp, pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iE, nvtx, pa, pb,ESIDE);
			center_cord3D(sosedi[ESIDE][iE].iNODE1, nvtx, pa, pbb,EE);

			potent[GRADXVX][iW] = my_quadratic_interpolation('-', potent[GRADXVX][sosedi[ESIDE][iE].iNODE1], potent[GRADXVX][iE], potent[GRADXVX][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADYVX][iW] = my_quadratic_interpolation('-', potent[GRADYVX][sosedi[ESIDE][iE].iNODE1], potent[GRADYVX][iE], potent[GRADYVX][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADZVX][iW] = my_quadratic_interpolation('-', potent[GRADZVX][sosedi[ESIDE][iE].iNODE1], potent[GRADZVX][iE], potent[GRADZVX][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);

			potent[GRADXVY][iW] = my_quadratic_interpolation('-', potent[GRADXVY][sosedi[ESIDE][iE].iNODE1], potent[GRADXVY][iE], potent[GRADXVY][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADYVY][iW] = my_quadratic_interpolation('-', potent[GRADYVY][sosedi[ESIDE][iE].iNODE1], potent[GRADYVY][iE], potent[GRADYVY][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADZVY][iW] = my_quadratic_interpolation('-', potent[GRADZVY][sosedi[ESIDE][iE].iNODE1], potent[GRADZVY][iE], potent[GRADZVY][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);

			potent[GRADXVZ][iW] = my_quadratic_interpolation('-', potent[GRADXVZ][sosedi[ESIDE][iE].iNODE1], potent[GRADXVZ][iE], potent[GRADXVZ][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADYVZ][iW] = my_quadratic_interpolation('-', potent[GRADYVZ][sosedi[ESIDE][iE].iNODE1], potent[GRADYVZ][iE], potent[GRADYVZ][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
			potent[GRADZVZ][iW] = my_quadratic_interpolation('-', potent[GRADZVZ][sosedi[ESIDE][iE].iNODE1], potent[GRADZVZ][iE], potent[GRADZVZ][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
		}

		if (bN) {

			// квадратичная интерполляция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iS, nvtx, pa, pb,SSIDE);
			center_cord3D(sosedi[SSIDE][iS].iNODE1, nvtx, pa, pbb,SS);

			potent[GRADXVX][iN] = my_quadratic_interpolation('+', potent[GRADXVX][sosedi[SSIDE][iS].iNODE1], potent[GRADXVX][iS], potent[GRADXVX][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADYVX][iN] = my_quadratic_interpolation('+', potent[GRADYVX][sosedi[SSIDE][iS].iNODE1], potent[GRADYVX][iS], potent[GRADYVX][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADZVX][iN] = my_quadratic_interpolation('+', potent[GRADZVX][sosedi[SSIDE][iS].iNODE1], potent[GRADZVX][iS], potent[GRADZVX][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);

			potent[GRADXVY][iN] = my_quadratic_interpolation('+', potent[GRADXVY][sosedi[SSIDE][iS].iNODE1], potent[GRADXVY][iS], potent[GRADXVY][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADYVY][iN] = my_quadratic_interpolation('+', potent[GRADYVY][sosedi[SSIDE][iS].iNODE1], potent[GRADYVY][iS], potent[GRADYVY][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADZVY][iN] = my_quadratic_interpolation('+', potent[GRADZVY][sosedi[SSIDE][iS].iNODE1], potent[GRADZVY][iS], potent[GRADZVY][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);

			potent[GRADXVZ][iN] = my_quadratic_interpolation('+', potent[GRADXVZ][sosedi[SSIDE][iS].iNODE1], potent[GRADXVZ][iS], potent[GRADXVZ][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADYVZ][iN] = my_quadratic_interpolation('+', potent[GRADYVZ][sosedi[SSIDE][iS].iNODE1], potent[GRADYVZ][iS], potent[GRADYVZ][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			potent[GRADZVZ][iN] = my_quadratic_interpolation('+', potent[GRADZVZ][sosedi[SSIDE][iS].iNODE1], potent[GRADZVZ][iS], potent[GRADZVZ][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
		}

		if (bS) {

            // квадратичная интерполляция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp, 100);
		    center_cord3D(iN, nvtx, pa, pb, NSIDE);
			center_cord3D(sosedi[NSIDE][iN].iNODE1, nvtx, pa, pbb, NN);

			potent[GRADXVX][iS] = my_quadratic_interpolation('-', potent[GRADXVX][sosedi[NSIDE][iN].iNODE1], potent[GRADXVX][iN], potent[GRADXVX][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADYVX][iS] = my_quadratic_interpolation('-', potent[GRADYVX][sosedi[NSIDE][iN].iNODE1], potent[GRADYVX][iN], potent[GRADYVX][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADZVX][iS] = my_quadratic_interpolation('-', potent[GRADZVX][sosedi[NSIDE][iN].iNODE1], potent[GRADZVX][iN], potent[GRADZVX][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);

			potent[GRADXVY][iS] = my_quadratic_interpolation('-', potent[GRADXVY][sosedi[NSIDE][iN].iNODE1], potent[GRADXVY][iN], potent[GRADXVY][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADYVY][iS] = my_quadratic_interpolation('-', potent[GRADYVY][sosedi[NSIDE][iN].iNODE1], potent[GRADYVY][iN], potent[GRADYVY][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADZVY][iS] = my_quadratic_interpolation('-', potent[GRADZVY][sosedi[NSIDE][iN].iNODE1], potent[GRADZVY][iN], potent[GRADZVY][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);

			potent[GRADXVZ][iS] = my_quadratic_interpolation('-', potent[GRADXVZ][sosedi[NSIDE][iN].iNODE1], potent[GRADXVZ][iN], potent[GRADXVZ][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADYVZ][iS] = my_quadratic_interpolation('-', potent[GRADYVZ][sosedi[NSIDE][iN].iNODE1], potent[GRADYVZ][iN], potent[GRADYVZ][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
			potent[GRADZVZ][iS] = my_quadratic_interpolation('-', potent[GRADZVZ][sosedi[NSIDE][iN].iNODE1], potent[GRADZVZ][iN], potent[GRADZVZ][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
		}

		if (bT) {

			// квадратичная интерполляция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iB, nvtx, pa, pb,BSIDE);
			center_cord3D(sosedi[BSIDE][iB].iNODE1, nvtx, pa, pbb,BB);
					
			potent[GRADXVX][iT] = my_quadratic_interpolation('+', potent[GRADXVX][sosedi[BSIDE][iB].iNODE1], potent[GRADXVX][iB], potent[GRADXVX][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADYVX][iT] = my_quadratic_interpolation('+', potent[GRADYVX][sosedi[BSIDE][iB].iNODE1], potent[GRADYVX][iB], potent[GRADYVX][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADZVX][iT] = my_quadratic_interpolation('+', potent[GRADZVX][sosedi[BSIDE][iB].iNODE1], potent[GRADZVX][iB], potent[GRADZVX][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);

			potent[GRADXVY][iT] = my_quadratic_interpolation('+', potent[GRADXVY][sosedi[BSIDE][iB].iNODE1], potent[GRADXVY][iB], potent[GRADXVY][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADYVY][iT] = my_quadratic_interpolation('+', potent[GRADYVY][sosedi[BSIDE][iB].iNODE1], potent[GRADYVY][iB], potent[GRADYVY][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADZVY][iT] = my_quadratic_interpolation('+', potent[GRADZVY][sosedi[BSIDE][iB].iNODE1], potent[GRADZVY][iB], potent[GRADZVY][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);

			potent[GRADXVZ][iT] = my_quadratic_interpolation('+', potent[GRADXVZ][sosedi[BSIDE][iB].iNODE1], potent[GRADXVZ][iB], potent[GRADXVZ][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADYVZ][iT] = my_quadratic_interpolation('+', potent[GRADYVZ][sosedi[BSIDE][iB].iNODE1], potent[GRADYVZ][iB], potent[GRADYVZ][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			potent[GRADZVZ][iT] = my_quadratic_interpolation('+', potent[GRADZVZ][sosedi[BSIDE][iB].iNODE1], potent[GRADZVZ][iB], potent[GRADZVZ][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
		}

		if (bB) {

			// квадратичная интерполляция.

            TOCHKA pp,pb,pbb;
		    center_cord3D(iP, nvtx, pa, pp,100);
		    center_cord3D(iT, nvtx, pa, pb,TSIDE);
			center_cord3D(sosedi[TSIDE][iT].iNODE1, nvtx, pa, pbb, TTSIDE);


			potent[GRADXVX][iB] = my_quadratic_interpolation('-', potent[GRADXVX][sosedi[TSIDE][iT].iNODE1], potent[GRADXVX][iT], potent[GRADXVX][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADYVX][iB] = my_quadratic_interpolation('-', potent[GRADYVX][sosedi[TSIDE][iT].iNODE1], potent[GRADYVX][iT], potent[GRADYVX][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADZVX][iB] = my_quadratic_interpolation('-', potent[GRADZVX][sosedi[TSIDE][iT].iNODE1], potent[GRADZVX][iT], potent[GRADZVX][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);

			potent[GRADXVY][iB] = my_quadratic_interpolation('-', potent[GRADXVY][sosedi[TSIDE][iT].iNODE1], potent[GRADXVY][iT], potent[GRADXVY][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADYVY][iB] = my_quadratic_interpolation('-', potent[GRADYVY][sosedi[TSIDE][iT].iNODE1], potent[GRADYVY][iT], potent[GRADYVY][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADZVY][iB] = my_quadratic_interpolation('-', potent[GRADZVY][sosedi[TSIDE][iT].iNODE1], potent[GRADZVY][iT], potent[GRADZVY][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);

			potent[GRADXVZ][iB] = my_quadratic_interpolation('-', potent[GRADXVZ][sosedi[TSIDE][iT].iNODE1], potent[GRADXVZ][iT], potent[GRADXVZ][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADYVZ][iB] = my_quadratic_interpolation('-', potent[GRADYVZ][sosedi[TSIDE][iT].iNODE1], potent[GRADYVZ][iT], potent[GRADYVZ][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			potent[GRADZVZ][iB] = my_quadratic_interpolation('-', potent[GRADZVZ][sosedi[TSIDE][iT].iNODE1], potent[GRADZVZ][iT], potent[GRADZVZ][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
		}
	}

} // green_gaussO2

// нахождение производных от скорости первого или второго порядка точности.
void green_gauss(integer iP, doublereal** &potent, integer** &nvtx, TOCHKA* pa,
	ALICE_PARTITION** &sosedi, integer maxelm, bool bbond, FLOW &f) {

	// если bsecondorder==true то производные будут вычисляться со вторым порядком точности.
	bool bsecondorder=false; // если false то внутри может применяться монотонизирующая поправка Рхи-Чоу.

	if (bsecondorder) {
		// второй порядок точности.
		green_gaussO2(iP, potent, nvtx, pa, sosedi, maxelm, bbond);
	}
	else {
		// первый порядок точности.
		green_gaussO1(iP, potent, nvtx, pa, sosedi, maxelm, bbond, f.mf[iP], f.prop[RHO], f.prop_b[RHO]);
	}

} // green_gauss

#endif
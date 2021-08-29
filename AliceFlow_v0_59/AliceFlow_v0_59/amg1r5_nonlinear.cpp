// 11.01.2020 Специальная нелинейная версия алгоритма amg1r5 для решения
// нелинейного уравнения теплопроводности.

#pragma once
#ifndef AMG1R5_NONLINEAR_CPP
#define AMG1R5_NONLINEAR_CPP 1

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R5 SOLUTION-SUBROUTINES */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     SOLVE                                                SUBROUTINE */

/* ....................................................................... */
// non_linear_ version 28.05.2020

/* Subroutine */ integer solve_non_linear_(integer *madapt, integer *ncyc, integer *nrd,
	integer *nsolco, integer *nru, integer *iout, integer *ierr,
	doublereal *a, doublereal *u, doublereal *f, integer *ia, integer *ja,
	integer *iw, doublereal *eps, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *icg, integer *ifg, integer *nstcol,
	integer *iarr, /*real*/unsigned int *time, integer *ncyc0, integer *irow0, integer *
	levels, integer *nda, integer *ndja, integer *ndu, integer *ndf,
	integer *mda, integer *mdja, integer *mdu, integer *mdf, integer *iup,
	integer *ium, doublereal *resi, doublereal *res0, doublereal *res,
	BLOCK*& my_body, int& lb, integer maxelm_out, integer maxelm_plus_maxbound,
	int *& whot_is_block)
{

	
	/* Format strings */
	// static char fmt_9050[] = "(\002 *** ERROR IN SOLVE: NDU TOO SMALL ***"
	//    "\002)";
	//static char fmt_9060[] = "(\002 *** ERROR IN SOLVE: NDF TOO SMALL ***"
	//    "\002)";
	//static char fmt_9005[] = "(/\002 ************* CYCLING..... **********"
	//    "***\002/)";
	// static char fmt_9000[] = "(\002 CYCLE  0:\002,3x,\002RES=\002,d9.3)";
	//static char fmt_9040[] = "(/\002 CYCLING BETWEEN GRIDS 1 AND\002,i3"
	//  ",\002:\002/)";
	// static char fmt_9010[] = "(\002 CYCLE \002,i2,\002:   RESCG=\002,d9.3"
	//    ",\002   RES=\002,d9.3,\002   CFAC=\002,d9.3)";
	//static char fmt_9020[] = "(\002 CYCLE \002,i2,\002:   RES=\002,d9.3,\002"
	//    "   CFAC=\002,d9.3)";

	/* System generated locals */
	integer i__1 = 0, i__2 = 0;
	doublereal d__1 = 0.0, d__2 = 0.0, d__3 = 0.0;

	/* Builtin functions */
	// integer s_wsfe(cilist *), e_wsfe(void), i_sign(integer *, integer *), 
	//   do_fio(integer *, char *, ftnlen);
	integer i_sign(integer *, integer *);

	/* Local variables */
	integer i__ = 0, m = 0;
	integer /*l = 0,*/ n = 0;
	extern /* Subroutine */ integer cg_(integer *, integer *, integer *,
		doublereal *, doublereal *, doublereal *, integer *, integer *,
		integer *, integer *, integer *, integer *, integer *, /*real*/unsigned int *,
		integer *, integer *);
	doublereal fac = 0.0, ama = 0.0;
	extern /* Subroutine */ integer cyc_(integer *, integer *, integer *, integer
		*, integer *, integer *, integer *, integer *, integer *,
		doublereal *, doublereal *, doublereal *, integer *, integer *,
		integer *, integer *, integer *, integer *, integer *, integer *,
		integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *,
		integer *, integer *, integer *, integer *, integer *, integer *,
		integer *, integer *, integer *, doublereal *, doublereal *,
		integer *);
	integer nsc = 0;
	doublereal cfac = 0.0;
	extern /* Subroutine */ integer idec_(integer *, integer *, integer *,
		integer *);
	integer igam = 0, icgr = 0;
	doublereal fmax = 0.0, epsi = 0.0;
	integer msel = 0, iter = 0, nrcx = 0, nrdx = 0;
	doublereal umax = 0.0;
	integer nrux = 0;
	doublereal rescg = 0.0;
	extern /* Subroutine */ integer resid_(integer *, doublereal *, doublereal *,
		doublereal *, doublereal *, integer *, integer *, integer *,
		integer *, integer *, integer *);
	doublereal epsil = 0.0;
	integer iconv = 0;
	extern /* Subroutine */ integer usave_(integer *, integer *, doublereal *,
		integer *, integer *, integer *, integer *, /*real*/unsigned int *, integer *,
		integer *, integer *, integer *, integer *);
	integer ncycle = 0, ndigit = 0, nrdlen = 0;
	doublereal resold = 0.0;
	integer nrulen = 0, mfirst = 0, nrdtyp[10] = { 0 }, nrutyp[10] = { 0 };

	/* Fortran I/O blocks */
	// static cilist io___240 = { 0, 0, 0, fmt_9050, 0 };
	//static cilist io___241 = { 0, 0, 0, fmt_9060, 0 };
	// static cilist io___264 = { 0, 0, 0, fmt_9005, 0 };
	// static cilist io___265 = { 0, 0, 0, fmt_9000, 0 };
	//  static cilist io___269 = { 0, 0, 0, fmt_9040, 0 };
	// static cilist io___270 = { 0, 0, 0, fmt_9040, 0 };
	// static cilist io___273 = { 0, 0, 0, fmt_9010, 0 };
	//static cilist io___274 = { 0, 0, 0, fmt_9020, 0 };



	/*     SOLUTION PHASE OF AMG1R5 */


	/* ===> TEST OF AVAILABLE STORAGE */

	/* Parameter adjustments */
	--resi;
	--time;
	--iarr;
	--nstcol;
	--ifg;
	--icg;
	--imaxw;
	--iminw;
	--imax;
	--imin;
	--iw;
	--ja;
	--ia;
	--f;
	--u;
	--a;

	/* Function Body */
	if (*ndu < imax[*levels]) {
		//io___240.ciunit = *ium;
		//s_wsfe(&io___240);
		//e_wsfe();
		std::cout << "( *** ERROR IN SOLVE: NDU TOO SMALL ***)\n";
		*ierr = 4;
		return 0;
	}
	if (*ndf < imax[*levels]) {
		//io___241.ciunit = *ium;
		//s_wsfe(&io___241);
		//e_wsfe();
		std::cout << "( *** ERROR IN SOLVE: NDF TOO SMALL ***)\n";
		*ierr = 5;
		return 0;
	}

	m = *levels;
	*ncyc0 = 0;
	for (n = 11; n <= 20; ++n) {
		//time[n] = 0.f;
		time[n] = 0;
		/* L5: */
	}
	//if (*eps != 0.) {
	if (fabs(*eps) > MY_DOPUSK_AMG1R5) {
		epsi = *eps;
	}
	else {
		epsi = 1e-12;
	}

	/* ===> DECOMPOSE MADAPT */

	if (*madapt != 0) {
		idec_(madapt, &c__2, &ndigit, &iarr[1]);
		msel = iarr[1];
		if (msel == 2) {
			if (iarr[2] != 0) {
				fac = (doublereal)iarr[2];
				for (i__ = 1; i__ <= 100; ++i__) {
					fac /= 10.;
					if (fac <= 1.) {
						goto L9;
					}
					/* L8: */
				}
			}
			else {
				fac = .7;
			}
		}
	}
	else {
		msel = 2;
		fac = .7;
	}

	/* ===> DECOMPOSE NCYC */

L9:
	if (*ncyc != 0) {
		i__1 = abs(*ncyc);
		idec_(&i__1, &c__4, &ndigit, &iarr[1]);
		igam = i_sign(&iarr[1], ncyc);
		icgr = iarr[2];
		iconv = iarr[3];
		ncycle = iarr[4];
		if (ncycle == 0) {
			return 0;
		}
	}
	else {
		igam = 1;
		icgr = 0;
		iconv = 1;
		ncycle = 10;
	}

	/* ===> SET EPSI ACCORDING TO CONVERGENCE CRITERION GIVEN BY ICONV */

	if (iconv != 3) {
		if (iconv == 4) {
			ama = 0.;
			i__1 = imax[1];
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* Computing MAX */
				d__1 = ama, d__2 = a[ia[i__]];
				ama = myr_max(d__1, d__2);
				/* L6: */
			}
			epsi *= ama;
		}
	}
	else {
		fmax = 0.;
		i__1 = imax[1];
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing MAX */
			d__2 = fmax, d__3 = (d__1 = f[i__], fabs(d__1));
			fmax = myr_max(d__2, d__3);
			/* L7: */
		}
		epsi *= fmax;
	}

	/* ===> DECOMPOSE NRD */

	if (*nrd != 0) {
		idec_(nrd, &c__9, &ndigit, nrdtyp);
		nrdx = nrdtyp[1];
		nrdlen = ndigit - 2;
		i__1 = nrdlen;
		for (i__ = 1; i__ <= i__1; ++i__) {
			nrdtyp[i__ - 1] = nrdtyp[i__ + 1];
			/* L10: */
		}
	}
	else {
		nrdx = 1;
		nrdlen = 2;
		nrdtyp[0] = 3;
		nrdtyp[1] = 1;
	}

	/* ===> DECOMPOSE NRU */

	if (*nru != 0) {
		idec_(nru, &c__9, &ndigit, nrutyp);
		nrux = nrutyp[1];
		nrulen = ndigit - 2;
		i__1 = nrulen;
		for (i__ = 1; i__ <= i__1; ++i__) {
			nrutyp[i__ - 1] = nrutyp[i__ + 1];
			/* L40: */
		}
	}
	else {
		nrux = 1;
		nrulen = 2;
		nrutyp[0] = 3;
		nrutyp[1] = 1;
	}

	/* ===> DECOMPOSE NSOLCO */

	if (*nsolco != 0) {
		idec_(nsolco, &c__2, &ndigit, &iarr[1]);
		nsc = iarr[1];
		nrcx = iarr[2];

		/* ===> IN CASE OF YALE-SMP COARSE GRID SOLUTION, DON'T USE COARSEST */
		/*     GRID WITH LESS THAN 10 POINTS */

		if (nsc == 2) {
			for (i__ = m; i__ >= 1; --i__) {
				//l = i__;
				if (imax[i__] - imin[i__] >= 9) {
					goto L60;
				}
				/* L50: */
			}
		L60:
			m = i__;
			*levels = i__;
		}
	}
	else {
		nsc = 1;
		nrcx = 0;
	}

	/* ===> CYCLING */

	/* L100: */
	if (*iout != 0) {
		resid_(&c__1, res0, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
			imin[1], &imax[1], &iminw[1]);
		if (*iout == 3) {
			if (yes_print_amg) {
				std::cout << "( ************* CYCLING..... *************)\n";
				std::cout << "( CYCLE  0:,3x,RES=,"<< *res0 <<")\n";
			}

		}

		resold = *res0;
	}

	integer iadd1 = 1;//0
	integer n_a = maxelm_plus_maxbound;// размерность правой части.
	doublereal* u_old = nullptr;
	//x_temper = new doublerealT[n_a[0] + 1];
	u_old = (doublereal*)malloc(((integer)(n_a)+1) * sizeof(doublereal));

	doublereal* rthdsd_no_radiosity_patch_1 = nullptr;
	rthdsd_no_radiosity_patch_1 = (doublereal*)malloc(((integer)(n_a)+1) * sizeof(doublereal));
	for (integer i23 = 0; i23 < n_a; i23++) {
		rthdsd_no_radiosity_patch_1[i23] = f[i23 + iadd1];
	}

	doublereal* x_temper = nullptr;
	//x_temper = new doublerealT[n_a[0] + 1];
	x_temper = (doublereal*)malloc(((integer)(n_a)+1) * sizeof(doublereal));

	if (x_temper == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for x_temper my_agregat_amg.cpp..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}

	doublereal* rthdsd_loc123 = nullptr;
	//rthdsd_loc123 = new doublerealT[n_a[0] + 1];
	rthdsd_loc123 = (doublereal*)malloc(((integer)(n_a)+1) * sizeof(doublereal));
	if (rthdsd_loc123 == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for rthdsd_loc123 my_agregat_amg.cpp..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}
	

	//ncycle = 1300;
	//doublereal tmaxold = -1.0e30;
	
	doublereal rms = 1.0, rms1=1.0;
	bool brms1 = true;
	iter = 0;

	doublereal tmax = -1.0e30;
	for (integer i23 = 0; i23 < n_a; i23++) {
		
		if (u[i23 + iadd1] > tmax) {
			tmax = u[i23 + iadd1];
		}
	}

	

	i__1 = ncycle;
	while (/*(*res>0.1306)||(fabs(tmax- tmaxold)>2.0e-5)||*/(rms1/rms<15.0)&&(iter<300)) {

		iter++;
		
		for (integer i23 = 0; i23 < n_a; i23++) {
			u_old[i23] = u[i23 + iadd1];
			//if (u[i23 + iadd1] > tmax) {
				//tmax = u[i23 + iadd1];
			//}

			
		}
		
		//getchar();
		//tmaxold = tmax;


		usave_(&c__1, &icgr, &u[1], &imin[1], &imax[1], ndu, &m, &time[1],
			ierr, ium, mdu, ndf, mdf);
		cyc_(&c__1, &nrdx, nrdtyp, &nrdlen, &nrcx, &nrux, nrutyp, &nrulen, &
			igam, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &
			imax[1], &iminw[1], &imaxw[1], &ifg[1], &icg[1], &nstcol[1], &
			iarr[1], &time[1], irow0, &m, ium, ierr, &iter, &nsc, nda,
			ndja, mda, mdja, &msel, &fac, &resi[1], levels);
		cg_(&c__1, &icgr, &iter, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1],
			&imin[1], &imax[1], &iminw[1], &m, &time[1], ierr, ium);
		

		rms = 0.0;
		for (integer i23 = 0; i23 < n_a; i23++) {
			if (u[i23 + iadd1] > tmax) {
				tmax = u[i23 + iadd1];				
			}
			rms += (u[i23 + iadd1] - u_old[i23])*(u[i23 + iadd1] - u_old[i23]);
		}
		rms /= n_a;
		if (brms1) {
			rms1 = rms;
			brms1 = false;
		}
		std::cout << "tmax=" << tmax << " rms=" << rms << " " << rms1 / rms << "\n";

		if (iter > 0) {
			// установить 0 в случае отката на предыдущую стабильную локально-линейную версию алгоритма.
			// главная причина установки значения 1 является сокращение числа проходов для устранения
			// нелинейности в системе с 26 до 4. При установке 1 в данном месте кода надо в модуле
			// mysolver_v0_03 установить fHORF=1.0;
			//(iVar == TEMP) && всегда выполняется, данная нелинейная функция вызывается всегда только для температуры.
			if (/*(my_amg_manager.istabilization == 3)*/
				(AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5 == stabilization_amg1r5_algorithm))
			{
				if (bonly_solid_calculation  ) {
					if (bvacuumPrism) {
						// предполагается неизменный порядок следования позиций в x
						// и rthdsd.					
						
#pragma omp parallel for
						for (integer i23 = 0; i23 < n_a; i23++) {
							if (u[i23 + iadd1] < -272.15) u[i23 + iadd1] = -272.15;
							//x_temper[i23] = x[i23 + 1];
							// 0.01 параметр нижней релаксации.
							// 0.25; 0.2; 0.01.
							// 0.005
							// etalon 0.01 (1250it; 891it; 2it; 7s 770ms)
							// experimental 0.02 (1250it; 28it; 2it; 5s 120ms)
							// experimental 0.04 (740it; 2it; 3s 360ms)
							// experimental 0.05 (618it; 2it; 3s 50ms)
							// experimental 0.06 (533it; 2it; 2s 900ms)
							// experimental 0.07 (469it; 2it; 2s 800ms)
							// experimental 0.08 (420it; 2it; 2s 570ms)
							// experimental 0.09 (381 it; 2it; 2s 390ms)
							// experimental 0.1 (1250it и не сходится, переборщил).
							if (fabs(u[i23 + iadd1] - u_old[i23]) > 5.0) {
								// Порог 5С оптимум.
								// 0.04 14s 440ms 358it
								// 0.09 13s 290ms 314it optimum
								x_temper[i23] = u_old[i23] + 0.084 * (u[i23 + iadd1] - u_old[i23]);
							}
							else if (fabs(u[i23 + iadd1] - u_old[i23]) > 1.0) {
								// Порог 1С оптимум.
								// 0.09 13s 290ms 314it
								// 0.095 17s 10ms 318; 102; 22;
								// 0.091 316it; 10; 10; 13s 350ms;
								// 0.092 optimum
								x_temper[i23] = u_old[i23] + 0.092 * (u[i23 + iadd1] - u_old[i23]);
							}
							else {
								// 0.09 0.25 перебор.
								// 0,09; 0.11; 314it; 9it; 9it; 2s 310ms optimum
								// 0.09; 0.12; 300it; 33it; 20it; 2s 530ms; 
								// 0.09; 0.1; 308it; 35it; 2s 400ms;
								// 0.11 optimum
								x_temper[i23] = u_old[i23] + 0.11 * (u[i23 + iadd1] - u_old[i23]);
							}
							if (x_temper[i23] < -272.15) x_temper[i23] = -272.15;
							u[i23 + iadd1] = x_temper[i23];
						}

						// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
						// Вычисление осреднённых температур в К на границах вакуумных промежутков:
						for (integer i23 = 0; i23 < lb; i23++) {
							update_avg_temperatures(x_temper, my_body[i23]);
						}
						// Вычисление плотностей радиационных тепловых потоков:
						for (integer i23 = 0; i23 < lb; i23++) {
							calculation_density_radiation_heat_flux(my_body[i23]);
						}						

						for (integer i23 = 0; i23 < n_a; i23++) {
							rthdsd_loc123[i23] = rthdsd_no_radiosity_patch_1[i23];
							if ((i23 >= iadd_qnbc_maxelm) &&
								(qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
								(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on)) {
								// 0.25 13s 60ms 314; 10; 13; old value
								// 0.26 13s 0ms; 310; 21; 13;
								// 0.27 318; 13; 3;
								// 0.28 325; 42; 13;
								// 0.35 13s 560ms 347; //new optimum
								// 0.09 1250; 1250; время запредельно большое.
								doublereal alpha_relax142 = 0.35;// d_my_optimetric1_6_12_2019;// 0.35;
								//if (x_temper[i23] > qnbc[i23 - iadd_qnbc_maxelm].Tamb) {
									rthdsd_loc123[i23] = alpha_relax142 *
										(-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST *
										((273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) *
											(273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) -
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb))) +
											(1.0 - alpha_relax142) *
										(-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST *
										((273.15 + u_old[i23]) * (273.15 + u_old[i23]) *
											(273.15 + u_old[i23]) * (273.15 + u_old[i23]) -
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
											(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb)));
									rthdsd_loc123[i23] *= qnbc[i23 - iadd_qnbc_maxelm].dS;
								//}
								//else {
									// Модификация 17.06.2020. Мы не можем остыть меньше чем температура окружающей среды.
									// Тепловая мощность строго больше нуля - нагрев.
									// Не помогает. Сомнительная модификация.
									//rthdsd_loc123[i23] = 0.0;
								//}
							}
						}

						radiosity_patch_for_vacuum_Prism_Object_(rthdsd_loc123, my_body, lb, maxelm_out, whot_is_block);
#pragma omp parallel for
						for (integer i23 = 0; i23 < n_a; i23++) {
							u_old[i23] = x_temper[i23];
							//x_old[i23 + 1] = x[i23 + 1];
						}
#pragma omp parallel for
						for (integer i23 = 0; i23 < n_a; i23++) {
							f[i23 + iadd1] = rthdsd_loc123[i23];
						}

						
					}
					else if (b_sign_on_nonlinear_bc) {
						//  25 декабря 2015. Ускорение сходимости при использовании 
						// нелинейных граничных условий.
						
						bool bNewtonRichman = false;
						bool bStefanBolcman = false;

						for (integer i23 = 0; i23 < n_a; i23++) {
							if (i23 >= iadd_qnbc_maxelm) {
								if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on == false)) {
									// Стефан-Больцман.
									bStefanBolcman = true;
								}
								if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on == false)) {

									// Ньютон-Рихман.
									bNewtonRichman = true;
								}
								if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on)) {
									// Условие смешанного типа.
									bStefanBolcman = true;
									bNewtonRichman = true;
								}
							}

						}

#pragma omp parallel for
						for (integer i23 = 0; i23 < n_a; i23++) {
							if (u[i23 + iadd1] < -272.15) u[i23 + iadd1] = -272.15;
							//x_temper[i23] = x[i23 + 1];
							// 0.01 параметр нижней релаксации.
							// 0.25
							// 0.2
							// 10 июня 2018 года заменил на коэффициент нижней релаксации равный 0.9.
							// 0.01 <--
							//0.005
							//if (i23<maxelm_out)
							{//0.02
								doublereal alpha_relax_nonlinear_boundary_condition = 0.01;
								if (bNewtonRichman && (bStefanBolcman == false)) {
									alpha_relax_nonlinear_boundary_condition = 0.9;
								}
								x_temper[i23] = u_old[i23] + alpha_relax_nonlinear_boundary_condition * (u[i23 + iadd1] - u_old[i23]);
							}
							if (x_temper[i23] < -272.15) x_temper[i23] = -272.15;
							u[i23 + iadd1] = x_temper[i23];
						}


						
						if (iadd_qnbc_maxelm != maxelm_out) {
							std::cout << "error!!! iadd_qnbc_maxelm != maxelm_out \n"<< iadd_qnbc_maxelm<<" "<< maxelm_out;
							getchar();
						}
						for (integer i23 = 0; i23 < n_a; i23++) {
							rthdsd_loc123[i23] = rthdsd_no_radiosity_patch_1[i23];
							if (i23 >= iadd_qnbc_maxelm) {

								if (qnbc[i23 - iadd_qnbc_maxelm].dS < 0.0) {
									std::cout << "error!!! dS negative \n"<< qnbc[i23 - iadd_qnbc_maxelm].dS;
									getchar();
								}

								if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on == false)) {



									// Стефан Больцман.
									if (qnbc[i23 - iadd_qnbc_maxelm].dS <= 0.0) {
										std::cout << "error!!! dS negative \n" << qnbc[i23 - iadd_qnbc_maxelm].dS;
										getchar();
									}

									doublereal alpha_relax142 = 0.25;
									rthdsd_loc123[i23] = alpha_relax142 * (-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST * ((273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) - (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb))) +
										(1.0 - alpha_relax142) * (-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST * ((273.15 + u_old[i23 ]) * (273.15 + u_old[i23 ]) * (273.15 + u_old[i23]) * (273.15 + u_old[i23]) - (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb)));
									rthdsd_loc123[i23] *= qnbc[i23 - iadd_qnbc_maxelm].dS;
								}
								if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on == false)) {


									// Ньютон-Рихман.
									//doublereal alpha_relax142 = 0.25;
									//printf("Tamb=%e alpha=%e temp=%e dS=%e\n", qnbc[i23 - iadd_qnbc_maxelm].Tamb, qnbc[i23 - iadd_qnbc_maxelm].film_coefficient, x_temper[i23], qnbc[i23 - iadd_qnbc_maxelm].dS);

									if (qnbc[i23 - iadd_qnbc_maxelm].dS <= 0.0) {
										std::cout << "error!!! dS negative \n" << qnbc[i23 - iadd_qnbc_maxelm].dS;
										getchar();
									}

									// По физическому смыслу температура не может опуститься ниже температуры среды.
									doublereal dtmp= x_temper[i23] - qnbc[i23 - iadd_qnbc_maxelm].Tamb;
									if (dtmp<0.0) dtmp=0.0;
									rthdsd_loc123[i23] = -qnbc[i23 - iadd_qnbc_maxelm].film_coefficient * (dtmp);
									rthdsd_loc123[i23] *= qnbc[i23 - iadd_qnbc_maxelm].dS;
								}
								if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) && 
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on)) {
									// Условие смешанного типа.


									//doublerealT alpha_relax142 = 0.25;

									if (qnbc[i23 - iadd_qnbc_maxelm].dS <= 0.0) {
										std::cout << "error!!! dS negative \n" << qnbc[i23 - iadd_qnbc_maxelm].dS;
										getchar();
									}

									rthdsd_loc123[i23] = (-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST * ((273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) - (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb)));
									rthdsd_loc123[i23] += -qnbc[i23 - iadd_qnbc_maxelm].film_coefficient * (x_temper[i23] - qnbc[i23 - iadd_qnbc_maxelm].Tamb);
									rthdsd_loc123[i23] *= qnbc[i23 - iadd_qnbc_maxelm].dS;
								}
							}
						}
#pragma omp parallel for
						for (integer i23 = 0; i23 < n_a; i23++) {
							u_old[i23] = x_temper[i23];
							
						}
#pragma omp parallel for
						for (integer i23 = 0; i23 < n_a; i23++) {
							f[i23 + iadd1] = rthdsd_loc123[i23];
						}

						

					}
				}
			}
		}
		
		

		if (*ierr > 0) {
			printf("ierr>0\n");
			//return 0;
		}
		if (iter == 1) {
			mfirst = m;
			if (*iout == 3) {

				if (yes_print_amg) {
#if doubleintprecision == 1
					printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%lld,)\n", m);
#else
					printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%d,)\n", m);
#endif

				}
			}
		}
		else if (*iout == 3 && m != mfirst) {
			mfirst = m;

			if (yes_print_amg) {
#if doubleintprecision == 1
				printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%lld,)\n", m);
#else
				printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%d,)\n", m);
#endif

			}


		}
		if (*iout == 3 || iconv != 1) {
			resid_(&c__1, res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
				imin[1], &imax[1], &iminw[1]);
		}
		*ncyc0 = iter;
		if (*iout != 3) {
			goto L110;
		}
		cfac = *res / (resold + 1e-40);
		resold = *res;
		if (1 == m) {
			goto L150;
		}
		resid_(&c__2, &rescg, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
			imin[1], &imax[1], &iminw[1]);		

		if (yes_print_amg) {		
			std::cout << "( CYCLE ," << iter << ",:   RESCG=" << rescg << std::endl;
			std::cout << ",   RES=" << *res << ",   CFAC=" << cfac << ")" << std::endl;
		}
		goto L110;
	L150:
		if (yes_print_amg) {			
			std::cout << "( CYCLE ," << iter << ",:   RES=," << *res << "," << std::endl;			
			std::cout << "   CFAC=," << cfac << ")" << std::endl;
		}

	L110:
		if (iconv == 1) {
			goto L120;
		}
		epsil = epsi;
		if (iconv != 4) {
			goto L115;
		}
		umax = 0.;
		i__2 = imax[1];
		for (i__ = imin[1]; i__ <= i__2; ++i__) {
			/* Computing MAX */
			d__2 = umax, d__3 = (d__1 = u[i__], fabs(d__1));
			umax = myr_max(d__2, d__3);
			/* L160: */
		}
		epsil = epsi * umax;
	L115:
		if (*res < epsil) {
			goto L170;
		}
	L120:
		;
		
		
	}

	if (u_old != nullptr) {
		free(u_old);
	}
	u_old = nullptr;
	if (rthdsd_no_radiosity_patch_1!=nullptr) {
		free(rthdsd_no_radiosity_patch_1);
	}
	rthdsd_no_radiosity_patch_1 = nullptr;

	if (rthdsd_loc123 != nullptr) {
		free(rthdsd_loc123);
	}
	rthdsd_loc123 = nullptr;

	if (x_temper != nullptr) {
		free(x_temper);
	}
	x_temper = nullptr;

	

L170:
	if (*iout != 3 && *iout != 0) {
		resid_(&c__1, res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[
			1], &imax[1], &iminw[1]);
	}

	if (u_old != nullptr) {
		free(u_old);
	}
	u_old = nullptr;
	if (rthdsd_no_radiosity_patch_1!=nullptr) {
		free(rthdsd_no_radiosity_patch_1);
	}
	rthdsd_no_radiosity_patch_1 = nullptr;

	if (rthdsd_loc123 != nullptr) {
		free(rthdsd_loc123);
	}
	rthdsd_loc123 = nullptr;

	if (x_temper != nullptr) {
		free(x_temper);
	}
	x_temper = nullptr;

	return 0;

	/* L9030: */
} /* solve_non_linear_ */

  //#if AMG1R6_LABEL==1

  /*     Last change:  K    25 Jul 2002   12:55 pm */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R6                                        MAIN SUBROUTINE */

  /*     RELEASE 1.6, July 2002 */
  /* 1.  changed: value of ntrim in pcol */
  /* 2.  dimensioning (1) changed to (*) in some subroutines to avoid subscript */
  /*     range checks in sparse solvers */

  //#endif


  /*     Last change:  ERB  22 Aug 2000   10:31 am */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R5                                        MAIN SUBROUTINE */

  /*     RELEASE 1.5, OCTOBER 1990 */

  /*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

  /* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
  /*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
  /*     FOLLOWING LINE IN SUBROUTINE RESC: */
  /*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
  /*     HAS BEEN CHANGED TO: */
  /*     IW(IMAXW(KC-1)+1) = IAUX */

  /*
  1.  БЫЛА ОБНАРУЖЕНА ОШИБКА, КОТОРАЯ ПРИ ОПРЕДЕЛЕННЫХ ОБСТОЯТЕЛЬСТВАХ
  НЕМНОГО ПОВЛИЯЛА НА СКОРОСТЬ КОНВЕРГЕНЦИИ AMG1R1. ПО ЭТОЙ ПРИЧИНЕ
  СЛЕДУЮЩАЯ СТРОКА В ПОДПРОГРАММЕ RESC:
  IW(IMAXW(KC-1)+1) = IA(IMIN(KC))
  БЫЛ ИЗМЕНЕН НА:
  IW(IMAXW(KC-1)+1) = IAUX
  */

  /* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
  /*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
  /*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
  /*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
  /*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
  /*     IBLCK1 = IMINW(K). */

  /* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
  /*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
  /*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
  /*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
  /*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
  /*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
  /*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

  /*
  3.  ВВЕДЕН ПАРАМЕТР LRATIO, ОБОЗНАЧАЮЩИЙ ОТНОШЕНИЕ
  ПРОСТРАНСТВА ЗАНЯТОГО ПЕРЕМЕННОЙ ДВОЙНОЙ ТОЧНОСТИ РЕАЛЬНОЙ И
  ЦЕЛОГО ЧИСЛА. ДЛЯ IBM ЗАДАНО СООТНОШЕНИЕ LRATIO
  РАВНОЕ 2. ПРИ НЕОБХОДИМОСТИ ИЗМЕНИТЕ ЭТО ЗНАЧЕНИЕ. (ЕСЛИ, НАПРИМЕР, ВЫ
  ЧТОБЫ ИЗМЕНИТЬ ВЕКТОРЫ ДВОЙНОЙ ТОЧНОСТИ НА ОДИНАРНУЮ ТОЧНОСТЬ,
  НЕОБХОДИМО УСТАНОВИТЬ КОЭФФИЦИЕНТ РАВНЫМ 1. В ЙЕЛЬСКОМ SMP-РЕЖИМЕ
  NDRV СУЩЕСТВУЕТ ТАКЖЕ ПАРАМЕТР LRATIO.
  */

  /* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
  /*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
  /*     TIVELY. */

  /* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
  /*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
  /*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

  /* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

  /* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
  /*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
  /*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */
  /*
  1.  ОШИБКА В ПРОВЕРКЕ ПОДПРОГРАММЫ БЫЛА УДАЛЕНА. ЕСЛИ ИСХОДНАЯ МАТРИЦА
  БЫЛА СОХРАНЕНА НЕСИММЕТРИЧНЫМ СПОСОБОМ, СИММЕТРИЗАЦИЯ AMG1R3 МОЖЕТ
  ПОТЕРПЕТЬ НЕУДАЧУ ПРИ ОПРЕДЕЛЕННЫХ ОБСТОЯТЕЛЬСТВАХ. ДЛЯ ИСПРАВЛЕНИЯ
  БЫЛИ ИЗМЕНЕНЫ СЛЕДУЮЩИЕ ИНСТРУКЦИИ В ФУНКЦИИ CHECK:
  */

  /*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 450 J=IA(I)+1,ICG(I)-1 */

  /*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

  /*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 550 J=IA(I)+1,ICG(I)-1 */

  /*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

  /* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
  /*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */
  /*
  2.  ПОЯСНИТЕЛЬНАЯ ЧАСТЬ ПОДПРОГРАММЫ AMG1R5 БЫЛА РАСШИРЕНА
  ВО ИЗБЕЖАНИЕ НЕДОРАЗУМЕНИЙ В ОПРЕДЕЛЕНИИ СПИСКА АРГУМЕНТОВ.
  */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
  /*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
  /*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
  /*     ENTRIES. */
  /*
  1.  ОШИБКА В ПРОВЕРКЕ ПОДПРОГРАММЫ БЫЛА УДАЛЕНА. ЕСЛИ ИСХОДНАЯ МАТРИЦА
  БЫЛА СОХРАНЕНА НЕСИММЕТРИЧНЫМ СПОСОБОМ, СИММЕТРИЗАЦИЯ ПО AMG1R3 ВСЕ ЕЩЕ
  МОГЛА ПОТЕРПЕТЬ НЕУДАЧУ ПРИ ОПРЕДЕЛЕННЫХ ОБСТОЯТЕЛЬСТВАХ И НЕ БЫЛА
  ЗАФИКСИРОВАНА В ПРЕДЫДУЩЕЙ ВЕРСИИ. КРОМЕ ТОГО, ПРОЦЕДУРА БЫЛА ИЗМЕНЕНА,
  ЧТОБЫ ИЗБЕЖАТЬ НЕКОТОРЫХ НЕНУЖНЫХ ПОИСКОВ СТРОК ДЛЯ ТРАНСПОНИРОВАНИЯ ЗАПИСЕЙ.
  */

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* Subroutine */ integer amg1r5_non_linear_(doublereal *a, integer *ia, integer *ja,
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig,
	integer *nnu, integer *matrix, integer *iswtch, integer *iout,
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc,
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco,
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2,
	integer *nwt, integer *ntr, integer *ierr,
	BLOCK*& my_body, int& lb, integer maxelm_out, integer maxelm_plus_maxbound,
	int *& whot_is_block)
{
	/* Format strings */

	doublereal tmax = -1.0e30;
	for (integer i23 = 0; i23 < maxelm_plus_maxbound; i23++) {

		if (u[i23] > tmax) {
			tmax = u[i23];
		}
	}
	printf("amg1r5_non_linear_ tmax=%e\n", tmax);

	/* Builtin functions */


	/* Local variables */
	integer mda = 0, mdf = 0, mdu = 0;
	doublereal res = 0.0;
	integer ium = 0, iup = 0;
	doublereal res0 = 0.0;
	extern /* Subroutine */ integer idec_(integer *, integer *, integer *,
		integer *);
	integer mdia = 0, mdja = 0, mdig = 0, iarr[100] = { 0 };
	//static real time[20];
	unsigned int time[20] = { 0 };
	integer imin[100] = { 0 }, imax[100] = { 0 };
	doublereal resi[100] = { 0.0 };
	integer kout = 0, ncyc0 = 0, irow0 = 0, ndicg = 0, icgst = 0, iminw[100] = { 0 }, imaxw[100] = { 0 };
	extern /* Subroutine */ integer first_(integer *, doublereal *, integer *,
		integer *, integer *, integer *), solve_(integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *, integer *, integer *, integer *,
			doublereal *, integer *, integer *, integer *, integer *, integer
			*, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *), setup_(integer *, integer *,
				integer *, doublereal *, doublereal *, doublereal *, integer *,
				integer *, integer *, doublereal *, doublereal *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *);
	integer ndigit = 0, levels = 0, kevelx = 0, nstcol[100] = { 0 }, kswtch = 0;
	extern /* Subroutine */ integer wrkcnt_(integer *, integer *, integer *,
		integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
		integer *, integer *, integer *, integer *, integer *, integer *,
		integer *, doublereal *, doublereal *);

	/* Fortran I/O blocks */




	/*         ----------------------------------------------- */
	/*         | AMG-MODULE FOR SOLVING LINEAR SYSTEMS L*U=F | */
	/*         ----------------------------------------------- */

	/*         -------------------------------------------------------------- */

	/*     ASSUMPTIONS ON L: */

	/*         THE PROGRAM REQUIRES: */ // Требования к входным данным:

										/*             - DIAGONAL ENTRIES ARE ALWAYS POSITIVE (ON ALL GRIDS); */
										/*             - L IS A SQUARE MATRIX WHICH IS EITHER REGULAR OR SINGULAR */
										/*               WITH ROWSUMS=0. */
										// Матрица L это квадратная матрица с всегда положительными диагональными элементами с нулевой суммой коэффициентов в каждой строке (положительная определённость).

										/*         FOR THEORETICAL REASONS THE FOLLOWING SHOULD HOLD: */

										/*             - L POSITIVE DEFINITE (OR SEMI-DEFINITE WITH ROWSUM=0) */
										/*             - L "ESSENTIALLY" POSITIVE TYPE, I.E., */

										/*                  -- DIAGONAL ENTRIES MUST BE > 0 ; */
										/*                  -- MOST OF THE OFF-DIAGONAL ENTRIES <= 0 ; */
										/*                  -- ROWSUMS SHOULD BE >= 0 . */

										// Матрица L положительно определённая: диагональные элементы  строго больше нуля, большинство внедиагональных элементов <= 0.
										// Сумма коэффициентов в строке больше либо равна нулю - диагональное преобладание.


										/*     THE USER HAS TO PROVIDE THE MATRIX L, THE RIGHT HAND SIDE F AND */
										/*     CERTAIN POINTER VECTORS IA AND JA. */

										// Пользователь задаёт матрицу коэффициентов L, правую часть F, а также информацию о связях между элементами матрицы в специальных столбцах IA и JA.
										// IA - позиция первого элемента в строке. JA - номер столбца для каждого элемента, первым записывается диагональный элемент.
										// В фортране нумерация начинается с единицы.

										/*         -------------------------------------------------------------- */

										/*     STORAGE OF L: */ // Требования к хранению матрицы L.

										/*         THE NON-ZERO ENTRIES OF THE MATRIX L ARE STORED IN */
										/*         "COMPRESSED" SKY-LINE FASHION IN A 1-D VECTOR A, I.E., ROW */
										/*         AFTER ROW, EACH ROW STARTING WITH ITS DIAGONAL ELEMENT. THE */
										/*         OTHER NON-ZERO ROW ENTRIES FOLLOW THEIR DIAGONAL ENTRY IN ANY */
										/*         ORDER. */


										/*         IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO */
										/*         PROVIDE TWO POINTER ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL */
										/*         NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L */
										/*         (1.LE.I.LE.NNU) ARE STORED IN A(J) WHERE THE RANGE OF J */
										/*         IS GIVEN BY */

										/*                     IA(I) .LE. J .LE. IA(I+1)-1. */

										/*         THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF */
										/*         ROW I WITHIN THE VECTOR A. IN PARTICULAR, */

										/*                     IA(1) = 1 ,  IA(NNU+1) = 1 + NNA */

										/*         WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED. */
										/*         THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT */
										/*         ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E., */
										/*         JA(J) POINTS TO THE COLUMN INDEX OF A(J). */
										/*         IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I */
										/*         AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I. */

										/*         IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS */
										/*         (FOR ANY I WITH  1.LE.I.LE.NNU): */

										/*                  F(I) =        SUM      A(J) * U(JA(J)) */
										/*                           J1.LE.J.LE.J2 */

										/*         WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND */
										/*         SIDE AND */

										/*                     J1 = IA(I) ,  J2 = IA(I+1)-1. */

										/*         NOTES: THE ENTRY IA(NNU+1) HAS TO TOCHKA TO THE FIRST FREE */
										/*                ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, */
										/*                AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW. */

										/*                THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5. */
										/*                SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT */
										/*                BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED */
										/*                THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE */
										/*                SETUP PHASE WILL FAIL. */
										/*                  ON THE OTHER HAND, RUNNING AMG A SECOND TIME ON THE */
										/*                SAME INPUT DATA WITH ISWTCH=4 HAS NO SENSE, BECAUSE */
										/*                THE RESULTS OF THE FIRST SETUP PHASE ARE STILL STORED */
										/*                AND THUS THIS PHASE CAN BE SKIPPED IN A SECOND CALL. */
										/*                IN ORDER TO DO THIS, SET ISWTCH TO 1, 2 OR 3. */

										/* ----------------------------------------------------------------------- */

										/*         THE FORM OF THE CALLING PROGRAM HAS TO BE AS FOLLOWS: */

										/*               PROGRAM DRIVER */
										/*         C */
										/*               DOUBLE PRECISION A(#NDA),U(#NDU),F(#NDF) */
										/*               INTEGER IA(#NDIA),JA(#NDJA),IG(#NDIG) */
										/*         C */
										/*               NDA  = #NDA */
										/*               NDU  = #NDU */
										/*               NDF  = #NDF */
										/*               NDIA = #NDIA */
										/*               NDJA = #NDJA */
										/*               NDIG = #NDIG */
										/*         C */
										/*         C     SET UP A, F, IA, JA AND SPECIFY NECESSARY PARAMETERS */
										/*         C */
										/*               .... */
										/*               .... */
										/*         C */
										/*               CALL AMG1R5(A,IA,JA,U,F,IG, */
										/*        +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX, */
										/*        +                  ISWTCH,IOUT,IPRINT, */
										/*        +                LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU, */
										/*        +                  ECG1,ECG2,EWT2,NWT,NTR, */
										/*        +                  IERR) */
										/*         C */
										/*               .... */
										/*               .... */
										/*         C */
										//#if AMG1R6_LABEL==1
										/*               CALL USTOP(' ') */
										//#else
										/*               STOP */
										//#endif
										/*               END */

										/* ----------------------------------------------------------------------- */

										/*     INPUT VIA ARRAYS (SEE ABOVE): */

										/*     A        -   MATRIX L */

										/*     IA       -   POINTER VECTOR */

										/*     JA       -   POINTER VECTOR */

										/*     U        -   FIRST APPROXIMATION TO SOLUTION */

										/*     F        -   RIGHT HAND SIDE */


										/* ----------------------------------------------------------------------- */


										/*     SCALAR INPUT PARAMETERS OF AMG1R5: */

										/*     THE INPUT PARAMETERS OF AMG1R5 IN THE LIST BELOW ARE ARRANGED */
										/*     ACCORDING TO THEIR IMPORTANCE TO THE GENERAL USER. THE PARAMETERS */
										/*     PRECEEDED BY A * MUST BE SPECIFIED EXPLICITELY. ALL THE OTHER */
										/*     PARAMETERS ARE SET TO STANDARD VALUES IF ZERO ON INPUT. */

										/*     THERE ARE FOUR CLASSES OF INPUT PARAMETERS WITH DECREASING PRI- */
										/*     ORITY: */

										/*     1. PARAMETERS DESCRIBING THE USER-DEFINED PROBLEM AND DIMENSIONING */
										/*        OF VECTORS IN THE CALLING PROGRAM */

										/*     2. PARAMETERS SPECIFYING SOME GENERAL ALGORITHMIC ALTERNATIVES AND */
										/*        THE AMOUNT OF OUTPUT DURING SOLUTION */

										/*     3. PARAMETERS CONTROLLING THE MULTIGRID CYCLING DURING THE SOLU- */
										/*        TION PHASE */

										/*     4. PARAMETERS CONTROLLING THE CREATION OF COARSER GRIDS AND INTER- */
										/*        POLATION FORMULAS. */

										/*     ONLY THE CLASS 1 - PARAMETERS MUST BE SPECIFIED EXPLICITELY BY */
										/*     THE USER. CLASS 2 - PARAMETERS CONTROL THE GENERAL PERFORMANCE OF */
										/*     AMG1R5. CHANGING THEM DOESN'T REQUIRE UNDERSTANDING THE AMG - */
										/*     ALGORITHM. SPECIFYING NON-STANDARD-VALUES FOR CLASS 3 - PARAMETERS */
										/*     PRESUPPOSES A GENERAL KNOWLEDGE OF MULTIGRID METHODS, WHEREAS THE */
										/*     FUNCTION OF CLASS 4 - PARAMETERS IS ONLY UNDERSTANDABLE AFTER */
										/*     STUDYING THE AMG-ALGORITHM IN DETAIL. FORTUNATELY IN MOST CASES */
										/*     THE CHOICE OF CLASS 3 AND 4 - PARAMETERS ISN'T CRITICAL AND USING */
										/*     THE AMG1R5 - SUPPLIED STANDARD VALUES SHOULD GIVE SATISFACTORY */
										/*     RESULTS. */

										/*         -------------------------------------------------------------- */

										/*     CLASS 1 - PARAMETERS: */

										/*  *  NDA      -   DIMENSIONING OF VECTOR A IN CALLING PROGRAM */

										/*  *  NDIA     -   DIMENSIONING OF VECTOR IA IN CALLING PROGRAM */

										/*  *  NDJA     -   DIMENSIONING OF VECTOR JA IN CALLING PROGRAM */

										/*  *  NDU      -   DIMENSIONING OF VECTOR U IN CALLING PROGRAM */

										/*  *  NDF      -   DIMENSIONING OF VECTOR F IN CALLING PROGRAM */

										/*  *  NDIG     -   DIMENSIONING OF VECTOR IG IN CALLING PROGRAM */

										/*  *  NNU      -   NUMBER OF UNKNOWNS */

										/*  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L. */

										/*                  1ST DIGIT OF MATRIX  --  ISYM: */
										/*                    =1: L IS SYMMETRIC; */
										/*                    =2: L IS NOT SYMMETRIC. */

										/*                  2ND DIGIT OF MATRIX  --  IROW0: */
										/*                    =1: L HAS ROWSUM ZERO; */
										/*                    =2: L DOES NOT HAVE ROWSUM ZERO. */

										/*         -------------------------------------------------------------- */

										/*     CLASS 2 - PARAMETERS: */

										/*     ISWTCH   -   PARAMETER CONTROLLING WHICH MODULES OF AMG1R5 ARE TO */
										/*                  BE USED. */
										/*                    =1:   CALL FOR -----, -----, -----, WRKCNT. */
										/*                    =2:   CALL FOR -----, -----, SOLVE, WRKCNT. */
										/*                    =3:   CALL FOR -----, FIRST, SOLVE, WRKCNT. */
										/*                    =4:   CALL FOR SETUP, FIRST, SOLVE, WRKCNT. */
										/*                  SETUP DEFINES THE OPERATORS NEEDED IN THE SOLUTION */
										/*                         PHASE. */
										/*                  FIRST INITIALIZES THE SOLUTION VECTOR (SEE PARAMETER */
										/*                         IFIRST). */
										/*                  SOLVE COMPUTES THE SOLUTION BY AMG CYCLING (SEE */
										/*                         PARAMETER NCYC). */
										/*                  WRKCNT PROVIDES THE USER WITH INFORMATION ABOUT */
										/*                         RESIDUALS, STORAGE REQUIREMENTS AND CP-TIMES */
										/*                         (SEE PARAMETER IOUT). */
										/*                  IF AMG1R5 IS CALLED THE FIRST TIME, ISWTCH HAS TO */
										/*                  BE =4. INDEPENDENT OF ISWTCH, SINGLE MODULES CAN BE */
										/*                  BYPASSED BY A PROPER CHOICE OF THE CORRESPONDING */
										/*                  PARAMETER. */

										/*     IOUT     -   PARAMETER CONTROLLING THE AMOUNT OF OUTPUT DURING */
										/*                  SOLUTION PHASE: */

										/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO. */

										/*                  2ND DIGIT: */
										/*                    =0: NO OUTPUT (EXCEPT FOR MESSAGES) */
										/*                    =1: RESIDUAL BEFORE AND AFTER SOLUTION PROCESS */
										/*                    =2: ADD.: STATISTICS ON CP-TIMES AND STORAGE REQUI- */
										/*                        REMENTS */
										/*                    =3: ADD.: RESIDUAL AFTER EACH AMG-CYCLE */

										/*     IPRINT   -   PARAMETER SPECIFYING THE FORTRAN UNIT NUMBERS FOR */
										/*                  OUTPUT: */

										/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO */

										/*                  2ND AND 3RD DIGIT  --  IUP: UNIT NUMBER FOR RESULTS */

										/*                  4TH AND 5TH DIGIT  --  IUM: UNIT NUMBER FOR MESSAGES */

										/*         -------------------------------------------------------------- */

										/*     CLASS 3 - PARAMETERS: */

										/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

										/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

										/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

										/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
										/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
										/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
										/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
										/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
										/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
										/*                        BY THE FOLLWING DIGITS. */

										/*                  REST OF IFIRST  --  RNDU: */
										/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
										/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
										/*                    IFIRST=1372815) */

										/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
										/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

										/*                  1ST DIGIT OF NCYC  --  IGAM: */
										/*                    =1: V -CYCLE, */
										/*                    =2: V*-CYCLE, */
										/*                    =3: F -CYCLE, */
										/*                    =4: W -CYCLE. */
										/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
										/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
										/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

										/*                  2ND DIGIT OF NCYC  --  ICGR: */
										/*                    =0: NO CONJUGATE GRADIENT, */
										/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
										/*                    =2: CONJUGATE GRADIENT (FULL CG). */

										/*                  3RD DIGIT OF NCYC  --  ICONV: */
										/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
										/*                    (FINEST GRID): */
										/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
										/*                        NCYCLE (SEE BELOW) */
										/*                    =2: STOP, IF  ||RES|| < EPS */
										/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
										/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
										/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
										/*                           EPS     (SEE INPUT PARAMETER EPS) */
										/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
										/*                           |U|   = SUPREMUM NORM OF SOLUTION */
										/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
										/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
										/*                    AFTER AT MOST NCYCLE CYCLES. */

										/*                  REST OF NCYC  --  NCYCLE: */
										/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
										/*                    NCYCLE=0: NO CYCLING. */

										/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
										/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
										/*                  ARE PERFORMED, REGARDLESS OF EPS. */

										/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
										/*                  GRID IN CYCLING: */

										/*                  1ST DIGIT OF MADAPT  --  MSEL: */
										/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
										/*                        PHASE ARE USED WITHOUT CHECK. */
										/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
										/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
										/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
										/*                        (SEE BELOW). */

										/*                  REST OF MADAPT  --  FAC */
										/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
										/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
										/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
										/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
										/*                        BY DEFAULT. */


										/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

										/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

										/*                  2ND DIGIT OF NRD  --  NRDX: */
										/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
										/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

										/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
										/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
										/*                    =2: FULL GS SWEEP */
										/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
										/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

										/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

										/*                  1ST DIGIT  --  NSC: */
										/*                    =1: GAUSS-SEIDEL METHOD */
										/*                    =2: DIRECT SOLVER (YALE SMP) */

										/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
										/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
										/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
										/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
										/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

										/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

										/*         -------------------------------------------------------------- */

										/*     CLASS 4 - PARAMETERS: */

										/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
										/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
										/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
										/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

										/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
										/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
										/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
										/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

										/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
										/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
										/*                        OPERATORS */
										/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */


										/* ----------------------------------------------------------------------- */

										/*     OUTPUT: */

										/*     U        -   CONTAINS THE COMPUTED SOLUTION */


										/*     IERR     -   ERROR PARAMETER: */

										/*                    >0: FATAL ERROR (ABNORMAL TERMINATION OF AMG1R5) */
										/*                    <0: NON-FATAL ERROR (EXECUTION OF AMG1R5 CONTINUES) */

										/*                  ERROR CODES IN DETAIL: */

										/*                  1. DIMENSIONING TOO SMALL FOR VECTOR */
										/*                        A      (IERR = 1) */
										/*                        IA     (IERR = 2) */
										/*                        JA     (IERR = 3) */
										/*                        U      (IERR = 4) */
										/*                        F      (IERR = 5) */
										/*                        IG     (IERR = 6) */

										/*                     NO YALE-SMP BECAUSE OF STORAGE (NDA TOO SMALL): */
										/*                               (IERR = -1) */
										/*                     NO YALE-SMP BECAUSE OF STORAGE (NDJA TOO SMALL): */
										/*                               (IERR = -3) */
										/*                     NO CG BECAUSE OF STORAGE (NDU TOO SMALL): */
										/*                               (IERR = -4) */
										/*                     NO SPACE FOR TRANSPOSE OF INTERPOLATION (NDA OR */
										/*                                                     NDJA TOO SMALL): */
										/*                               (IERR = -1) */

										/*                  2. INPUT DATA ERRONEOUS: */

										/*                     A-ENTRY MISSING, ISYM = 1:           (IERR = -11) */
										/*                     PARAMETER MATRIX MAY BE ERRONEOUS:   (IERR = -12) */
										/*                     DIAGONAL ELEMENT NOT STORED FIRST:   (IERR =  13) */
										/*                     DIAGONAL ELEMENT NOT POSITIV:        (IERR =  14) */
										/*                     POINTER IA ERRONEOUS:                (IERR =  15) */
										/*                     POINTER JA ERRONEOUS:                (IERR =  16) */
										/*                     PARAMETER ISWTCH ERRONEOUS:          (IERR =  17) */
										/*                     PARAMETER LEVELX ERRONEOUS:          (IERR =  18) */

										/*                  3. ERRORS OF THE AMG1R5-SYSTEM (SHOULD NOT OCCUR): */

										/*                     TRANSPOSE A-ENTRY MISSING:           (IERR =  21) */
										/*                     INTERPOLATION ENTRY MISSING:         (IERR =  22) */

										/*                  4. ALGORITHMIC ERRORS: */

										/*                     CG-CORRECTION NOT DEFINED:           (IERR =  31) */
										/*                     NO YALE-SMP BECAUSE OF ERROR IN */
										/*                     FACTORIZATION:                       (IERR = -32) */

										/* ----------------------------------------------------------------------- */

										/*     WORK SPACE: */

										/*     THE INTEGER VECTOR IG HAS TO BE PASSED TO AMG1R5 AS WORK SPACE. */

										/* ----------------------------------------------------------------------- */

										/*     DIMENSIONING OF INPUT VECTORS AND WORK SPACE: */

										/*     IT'S IMPOSSIBLE TO TELL IN ADVANCE THE EXACT STORAGE REQUIREMENTS */
										/*     OF AMG. THUS, THE FOLLOWING FORMULAS GIVE ONLY REASONABLE GUESSES */
										/*     FOR THE VECTOR LENGTHS WHICH HAVE TO BE DECLARED IN THE CALLING */
										/*     PROGRAM. IN THESE FORMULAS NNA DENOTES THE NUMBER OF NON-ZERO */
										/*     ENTRIES IN THE INPUT-MATRIX L AND NNU IS THE NUMBER OF UNKNOWNS. */

										/*     VECTOR         NEEDED LENGTH (GUESS) */
										/*       A               3*NNA + 5*NNU */
										/*       JA              3*NNA + 5*NNU */
										/*       IA              2.2*NNU */
										/*       U               2.2*NNU */
										/*       F               2.2*NNU */
										/*       IG              5.4*NNU */

										/* ----------------------------------------------------------------------- */


										/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

										/*          ISWTCH = 4 */
										/*          IOUT   = 12 */
										/*          IPRINT = 10606 */

										/*          LEVELX = 100 */
										/*          IFIRST = 13 */
										/*          NCYC   = 10110 */
										/*          EPS    = 1.D-12 */
										/*          MADAPT = 27 */
										/*          NRD    = 1131 */
										/*          NSOLCO = 110 */
										/*          NRU    = 1131 */

										/*          ECG1   = 0. */
										/*          ECG2   = 0.25 */
										/*          EWT2   = 0.35 */
										/*          NWT    = 2 */
										/*          NTR    = 0 */

										/*     IF ANY ONE OF THESE PARAMETERS IS 0 ON INPUT, ITS CORRESPONDING */
										/*     STANDARD VALUE IS USED BY AMG1R5. */

										/* ----------------------------------------------------------------------- */

										/*     PORTABILITY RESTRICTIONS: */

										/*     1. ROUTINE CTIME IS MACHINE DEPENDENT AND HAS TO BE ADAPTED TO */
										/*        YOUR COMPUTER INSTALLATION OR REPLACED BY A DUMMY ROUTINE. */

										/*     2. MOST INPUT PARAMETERS ARE COMPOSED OF SEVERAL DIGITS, THEIR */
										/*        SIGNIFICANCE HAVING BEEN DESCRIBED ABOVE. BE SURE NOT TO ENTER */
										/*        MORE DIGITS THAN YOUR COMPUTER CAN STORE ON AN INTEGER VARI- */
										/*        ABLE. */

										/*     3. APART FROM FORTRAN INTRINSIC FUNCTIONS AND SERVICE ROUTINES, */
										/*        THERE IS ONLY ONE EXTERNAL REFERENCE TO A PROGRAM NOT CONTAINED */
										/*        IN THE AMG1R5 - SYSTEM, I.E. THE LINEAR SYSTEM SOLVER NDRV OF */
										/*        THE YALE SPARSE MATRIX PACKAGE. IF YOU HAVN'T ACCESS TO THIS */
										/*        PACKAGE, ENTER A DUMMY ROUTINE NDRV AND AVOID CHOOSING NSC=2 */
										/*        (SUBPARAMETER OF NSOLCO). THEN NDRV ISN'T CALLED BY AMG1R5. */
										/*        IN THIS CASE, HOWEVER, INDEFINITE PROBLEMS WILL NOT BE SOLV- */
										/*        ABLE. */
										/*          THE YALE SPARSE MATRIX PACKAGE IS FREELY AVAILABLE FOR NON- */
										/*        PROFIT PURPOSES. CONTACT THE DEPARTMENT OF COMPUTER SCIENCE, */
										/*        YALE UNITVERSITY. */

										/*     4. IN AMG1R5 THERE IS THE PARAMETER LRATIO, DENOTING THE RATIO */
										/*        OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
										/*        THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
										/*        TO 2. CHANGE THIS VALUE IF NECESSARY. (THE SAME HAS TO BE */
										/*        DONE WITH THE YALE SMP-ROUTINE NDRV.) */


										/* ----------------------------------------------------------------------- */

										/*     AUTHORS: */

										/*          JOHN RUGE, FORT COLLINS (USA), */
										/*              INSTITUTE FOR COMPUTATIONAL STUDIES AT CSU; */

										/*          KLAUS STUEBEN, D-5205 ST. AUGUSTIN (W.-GERMANY), */
										/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

										/*          ROLF HEMPEL, D-5205 ST. AUGUSTIN (W.-GERMANY), */
										/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

										/* ----------------------------------------------------------------------- */


										/* ===> LRATIO HAS TO BE SET TO THE NUMBER OF INTEGERS OCCUPYING THE SAME */
										/*     AMOUNT OF STORAGE AS ONE DOUBLE PRECISION REAL. */


										/* ===> MAXGR IS THE MAXIMAL NUMBER OF GRIDS. CHANGING THIS UPPER LIMIT */
										/*     JUST REQUIRES CHANGING THE PARAMETER STATEMENT. */


										/* Parameter adjustments */
	--ig;
	--f;
	--u;
	--ja;
	--ia;
	--a;

	/* Function Body */
	*ierr = 0;

	/* ===> SET PARAMETERS TO STANDARD VALUES, IF NECCESSARY */


	if (*iout != 0) {
		idec_(iout, &c__2, &ndigit, iarr);
		kout = iarr[1];
	}
	else {
		kout = 2;
	}

	if (*iswtch != 0) {
		kswtch = *iswtch;
	}
	else {
		kswtch = 4;
	}

	if (*levelx > 0) {
		kevelx = my_imin(*levelx, 100);
	}
	else if (*levelx < 0) {
		goto L70;
	}
	else {
		kevelx = 100;
	}

	if (*iprint != 0) {
		idec_(iprint, &c__4, &ndigit, iarr);
		iup = iarr[1] * 10 + iarr[2];
		ium = iarr[3];
	}
	else {
		iup = 6;
		ium = 6;
	}
	icgst = *nnu + 3;
	ndicg = (*ndig - icgst + 1) / 2;
	if (ndicg <= 0) {
		goto L60;
	}

	switch (kswtch) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	}
	//io___10.ciunit = ium;
	//s_wsfe(&io___10);
	//e_wsfe();


	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER I.  SWTCH ***\n");

	*ierr = 17;
	return 0;

L40:
	setup_(nnu, matrix, &kevelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[
		1], &ia[1], &ja[1], &ig[1], imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &levels, &irow0, nda,
			ndia, ndja, ndu, ndf, &ndicg, &ium, &mda, &mdia, &mdja, &mdu, &
			mdf, &mdig, &c__25, &c__2);
	if (*ierr > 0) {
		return 0;
	}
L30:
	//first_(ifirst, &u[1], imin, imax, iarr, &irow0);
L20:

	/*
	// Для корректности работы надо передавать информацию о модифицированном размере в вызывающие функции вверх по коду.
	if ((a != nullptr) && (ja != nullptr)) {

	// На данном этапе доступен истинный размер матрицы СЛАУ - хранимое число ненулевых элементов.
	// Вычислим реальное число ненулевых элементов матрицы a. Если номер столбца равен нулю то число ненулевых элементов
	// в матрице заведомо меньше чем позиция этого нуля.
	// С помощью низкоуровневой операции realloc ужимаем размер матрицы a.
	// Последующие выделения оперативной памяти для внешнего Крыловского итерационного процесса BiCGStab не приведут к ещё большему расходу
	// оперативной памяти, т.к. мы её освободили.

	printf("nda=%lld\n", nda[0]);
	integer isize97 = 0;
	for (integer i_95 = 1; i_95 < ndja[0]; i_95++) {
	if (ja[i_95] == 0) {
	isize97 = i_95;
	if (i_95 + 2 < ndja[0] && ja[i_95 + 1] == 0 && ja[i_95 + 2] == 0) {
	break;
	}
	}
	}
	printf("nda_new=%lld\n", isize97);
	++a;
	++ja;
	a = (doublereal*)realloc(a, ((integer)(isize97)+2) * sizeof(doublereal));
	ja = (integer*)realloc(ja, ((integer)(isize97)+2) * sizeof(integer));
	--a;
	--ja;
	//getchar();
	}
	*/
	//getchar();
	solve_non_linear_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
		ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
		&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
		nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
		res0, &res, my_body, lb, maxelm_out, maxelm_plus_maxbound, whot_is_block);
	if (*ierr > 0) {
		return 0;
	}
L10:
	wrkcnt_(&kout, &ia[1], &ig[1], imin, imax, iminw, &levels, time, &ncyc0, &
		iup, &mda, &mdia, &mdja, &mdu, &mdf, &mdig, &res0, &res);
	return 0;
L60:
	//io___29.ciunit = ium;
	//s_wsfe(&io___29);
	//e_wsfe();
	printf("*** ERROR IN AMG1R5: NDIG TOO SMALL ***\n");

	*ierr = 6;
	return 0;
L70:
	//io___30.ciunit = ium;
	//s_wsfe(&io___30);
	//e_wsfe();

	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER LEVELX ***\n");

	*ierr = 18;
	return 0;
}// non_linear_
//#if AMG1R6_LABEL==1
/* amg1r6_ */
//#else
/* amg1r5_ */
//#endif

#endif
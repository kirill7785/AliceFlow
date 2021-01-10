// Файл pseudo_time.cpp содержит функцию которая отвечает за вычисление пседовремени tau
// которое используется при коррекции скорости и участвует в качестве коэффициента диффузии 
// в уравнении на поправку давления.
// begin 20 июня 2012 года.

#ifndef MY_PSEUDO_TIME_CPP
#define MY_PSEUDO_TIME_CPP 1


// Сборка одного скалярного псевдовремени.
// На пространственно неоднородных областях параметр tau сильно неоднороден это приводит к ухудшению 
// обусловленности системы и в конечном счёте к расходимости. Здесь предпринята попытка регулизовать параметр tau.
void tau_calc(doublereal* &tau, integer maxelm, integer maxbound,
	          float** prop, float** prop_b, doublereal* alpha, 
			  int** nvtx, TOCHKA* pa, equation3D** sl,
			  int*** neighbors_for_the_internal_node, equation3D_bon** slb,
			  bool btimedep, doublereal dtime, doublereal CFL1,
			  integer inumberSIMPLEiteration, bool &bVERY_STABILITY_ON, bool boldscheme) {

	// В случае с сильно неоднородной по пространству сеткой, а также сильно меняющимися
	// коэффициентами диффузии в уравнениях на скорость распределение шага по псевдо-времени
	// будет также сильно неоднородно в пространстве, что приводит к расходимости вычислительного
	// процесса в рамках SIMPLE процедуры. Данная функция по вычислению псевдо временного шага 
	// специально написана для того чтобы сгладить сильно неоднородное распределение tau особенно
	// на первых итерациях SIMPLE алгоритма. Особой устойчивостью и стабильностью обладает просто постоянное 
	// значение шага по псевдовремени.
	

	// Судя по данным CFX если с солвером всё нормально то порядка за 200 итераций устанавливается ламинарное течение.
	// Эта информация может помочь для правильного определения отсечки по номеру итерации. Предполагается что на первых итерациях
	// мы будем использовать bVERY_STABILITY_ON==true, а дальше false.

	// inumberSIMPLEiteration - номер итерации SIMPLE алгоритма. Нужна для идентификации первых итераций алгоритма,
	// которые могут быть особенно неустойчивыми.
	//bool bVERY_STABILITY_ON=false; // если равно true то будет использоваться просто постоянный шаг по псевдовремени во всём пространстве,
	// т.к. это очень устойчиво.

	// dtime - шаг по времени,
	// если btimedep==true значит задача нестационарная.
	// default parameters CFL1=1.0, CFL2=alpha/(1-alpha) приблизительно 4.
	//doublereal CFL1=1.0; 

	bVERY_STABILITY_ON=false;

	// инициализация.
    for (integer i=0; i<maxelm+maxbound; i++) {
		tau[i]=0.0; // псевдовремя.
	}

	for (integer iP=0; iP<maxelm; iP++) {

        // iP - номер центрального контрольного объёма
	    integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	    iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
	


	    // если с одной из сторон граница расчётной области 
	    // то переменная равна true
	    bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

        if (iE>=maxelm) bE=true;
	    if (iN>=maxelm) bN=true;
	    if (iT>=maxelm) bT=true;
        if (iW>=maxelm) bW=true;
	    if (iS>=maxelm) bS=true;
	    if (iB>=maxelm) bB=true;

	    // вычисление размеров текущего контрольного объёма:
    	doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
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

	    //printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	    //getchar();
		
		// Вообще говоря tau принято собирать лишь по конвективно -диффузионной составляющей, так что диагональный
		// член матрицы лучше не использовать а использовать Sum(af,f==E,N,T,W,S,B).
	    doublereal apue=1.0, apuw=1.0, apvn=1.0, apvs=1.0, apwt=1.0, apwb=1.0;
	    if (!bE) apue=sl[VELOCITY_X_COMPONENT][iE].ap*sl[VELOCITY_X_COMPONENT][iP].ap/(feplus*sl[VELOCITY_X_COMPONENT][iE].ap+(1-feplus)*sl[VELOCITY_X_COMPONENT][iP].ap); else apue=slb[VELOCITY_X_COMPONENT][iE-maxelm].aw;
	    if (!bW) apuw=sl[VELOCITY_X_COMPONENT][iW].ap*sl[VELOCITY_X_COMPONENT][iP].ap/(fwplus*sl[VELOCITY_X_COMPONENT][iW].ap+(1-fwplus)*sl[VELOCITY_X_COMPONENT][iP].ap); else apuw=slb[VELOCITY_X_COMPONENT][iW-maxelm].aw;
	    if (!bN) apvn=sl[VELOCITY_Y_COMPONENT][iN].ap*sl[VELOCITY_Y_COMPONENT][iP].ap/(fnplus*sl[VELOCITY_Y_COMPONENT][iN].ap+(1-fnplus)*sl[VELOCITY_Y_COMPONENT][iP].ap); else apvn=slb[VELOCITY_Y_COMPONENT][iN-maxelm].aw;
	    if (!bS) apvs=sl[VELOCITY_Y_COMPONENT][iS].ap*sl[VELOCITY_Y_COMPONENT][iP].ap/(fsplus*sl[VELOCITY_Y_COMPONENT][iS].ap+(1-fsplus)*sl[VELOCITY_Y_COMPONENT][iP].ap); else apvs=slb[VELOCITY_Y_COMPONENT][iS-maxelm].aw;
	    if (!bT) apwt=sl[VELOCITY_Z_COMPONENT][iT].ap*sl[VELOCITY_Z_COMPONENT][iP].ap/(ftplus*sl[VELOCITY_Z_COMPONENT][iT].ap+(1-ftplus)*sl[VELOCITY_Z_COMPONENT][iP].ap); else apwt=slb[VELOCITY_Z_COMPONENT][iT-maxelm].aw;
	    if (!bB) apwb=sl[VELOCITY_Z_COMPONENT][iB].ap*sl[VELOCITY_Z_COMPONENT][iP].ap/(fbplus*sl[VELOCITY_Z_COMPONENT][iB].ap+(1-fbplus)*sl[VELOCITY_Z_COMPONENT][iP].ap); else apwb=slb[VELOCITY_Z_COMPONENT][iB-maxelm].aw;

		doublereal alpha_avg=(alpha[VELOCITY_X_COMPONENT]+alpha[VELOCITY_Y_COMPONENT]+alpha[VELOCITY_Z_COMPONENT])/3.0;
		doublereal CFL2=4.0; // соответствует SIMPLEC alpha=0.8; 
		if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
		    // SIMPLEC 
		    CFL2=alpha_avg/(1.0-alpha_avg);
	    }
	    if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
		    // SIMPLE
		    CFL2=alpha_avg;
	    }
		doublereal tau_loc=CFL2*prop[RHO][iP]*dx*dy*dz/(apue+apuw+apvn+apvs+apwt+apwb);
		
		if (boldscheme) {
			tau[iP]=tau_loc;
		}
		else {
		// Внимание для повышения устойчивости рекомендуется увеличить параметр CFL1.
		// Увеличивая CFL2 можно уменьшать влияние локальных параметров течения на поправку давления,
		// чтобы увеличить CFL2 рекомендуется сделать параметр нижней релаксации для скорости сделать более близким к 1.0.
		// Например:
		// alpha CFL2
		// 0.7 2.33333
		// 0.5 1.0 // Рекомендовано С. Патанкаром по умолчанию.
		// 0.8 4.0 // Рекомендовано Гавриловым Андреем по умолчанию.
		// 0.85 5.666
		if (bVERY_STABILITY_ON==true) {
			// пока просто запоминаем, основное вычисление псевдовремени будет позже.
			tau[iP]=tau_loc;
		}
		else {
		   if (btimedep) {
			   // нестационарная задача.
			   tau[iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc));
		   }
		   else {
			   // стационарный солвер.
			   // в этом случае также можно применить формулу для нестационарного случая,
			   // вопрос в том что шаг по времени dtime теперь надо подбирать.
			   if (dtime<=0.0) {
				   // возьмём среднее арифметическое от всех существующих tau.
			       tau[iP]=tau_loc; // пока просто запомним чтобы потом можно было вычислить среднее арифметическое.
			   }
			   else {
				   // если параметр dtime положителен то мы считаем что он задан пользователем из соображений пользователя.
				   tau[iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc));
			   }
		   }
		}
		}

	} // next iP

	if (!boldscheme) {
	if (bVERY_STABILITY_ON==true) {
		// во всей области просто постоянное значение шага по псевдовремени.
		doublereal tau_avg=0.0;
		for (integer iP=0; iP<maxelm; iP++) {
		   tau_avg+=tau[iP];
		}
		tau_avg=tau_avg/((doublereal)(maxelm));
		if (btimedep) {
			// нестационарная задача.
			for (integer iP=0; iP<maxelm; iP++) {
			    tau[iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg)); // постоянный шаг по псевдовремени во всей расчётной области.
			}
		}
		else {
			if (dtime<=0.0) {
				// пользователь не определил дополнительный временной параметр в стационарном случае.
				for (integer iP=0; iP<maxelm; iP++) {
					tau[iP]=CFL1*tau_avg;
				}
			}
			else {
				// пользователь определил дополнительный опорный стабилизирующий временной параметр.
				for (integer iP=0; iP<maxelm; iP++) {
			        tau[iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg)); // постоянный шаг по псевдовремени во всей расчётной области.
			    }
			}
		}
	}
	else {
		// bVERY_STABILITY_ON==false.
		// Обычный рабочий режим, всё-таки более стабильный чем простое вычисление tau ~ rho*Vol*alpha/ap или tau ~ rho*Vol*alpha/(ap*(1.0-alpha)).

		if (dtime<=0.0) {
		   doublereal tau_avg=0.0;
		   for (integer iP=0; iP<maxelm; iP++) {
			   tau_avg+=tau[iP];
		   }
		   tau_avg=tau_avg/((doublereal)(maxelm));

		   for (integer iP=0; iP<maxelm; iP++) {
			  tau[iP]=CFL1*(1.0/(1.0/tau_avg+1.0/tau[iP]));
		   }
	    }
	}
	}

	// продолжение параметра tau на границу расчётной области.
	// Просто сносим значение tau из ближайшего внутреннего узла.
	// iP - номер центрального контрольного объёма
	for (integer iP=0; iP<maxelm; iP++) {
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	    iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

        if (iE>=maxelm) {
			// граничный узел
			tau[iE]=tau[iP];
		}

		if (iW>=maxelm) {
			// boundary node
			tau[iW]=tau[iP];
		}

		if (iN>=maxelm) {
			// граничный узел
			tau[iN]=tau[iP];
		}

		if (iS>=maxelm) {
			// boundary node
			tau[iS]=tau[iP];
		}
		if (iT>=maxelm) {
			// граничный узел
			tau[iT]=tau[iP];
		}

		if (iB>=maxelm) {
			// boundary node
			tau[iB]=tau[iP];
		}
	}

} // tau_calc

// На пространственно неоднородных областях параметр tau сильно неоднороден это приводит к ухудшению 
// обусловленности системы и в конечном счёте к расходимости. Здесь предпринята попытка регулизовать параметр tau.
// Реализовано 23 июня 2012 года.
// Модифицировано для АЛИС сеток 5 декабря 2018.
void tau_calc3(doublereal** &tau, integer maxelm, integer maxbound,
	          float** prop, float** prop_b, doublereal* alpha, 
			  int** nvtx, TOCHKA* pa, 
			  int*** neighbors_for_the_internal_node, doublereal** sumanb,
			  bool btimedep, doublereal dtime, doublereal CFL1,
			  integer inumberSIMPLEiteration, 
			  bool &bVERY_STABILITY_ON, bool boldscheme) {

	// В случае с сильно неоднородной по пространству сетки, а также сильно меняющимися
	// коэффициентами диффузии в уравнениях на скорость распределение шага по псевдо-времени
	// будет также сильно неоднородно в пространстве, что приводит к расходимости вычислительного
	// процесса в рамках SIMPLE процедуры. Данная функция по вычислению псевдо временного шага 
	// специально написана для того чтобы сгладить сильно неоднородное распределение tau особенно
	// на первых итерациях SIMPLE алгоритма. Особой устойчивостью и стабильностью обладает просто постоянное 
	// значение шага по псевдовремени.
	// tau[VX or VY or VZ][iP] - для каждой компоненты скорости своё псевдовремя.

	// Судя по данным CFX если с солвером всё нормально то порядка за 200 итераций устанавливается ламинарное течение.
	// Эта информация может помочь для правильного определения отсечки по номеру итерации. Предполагается что на первых итерациях
	// мы будем использовать bVERY_STABILITY_ON==true, а дальше false.

	// inumberSIMPLEiteration - номер итерации SIMPLE алгоритма. Нужна для идентификации первых итераций алгоритма,
	// которые могут быть особенно неустойчивыми.
	//bool bVERY_STABILITY_ON=false; // если равно true то будет использоваться просто постоянный шаг по псевдовремени во всём пространстве,
	// т.к. это очень устойчиво.

	// dtime - шаг по времени,
	// если btimedep==true значит задача нестационарная.
	// default parameters CFL1=1.0, CFL2=alpha/(1-alpha) приблизительно 4.
	//doublereal CFL1=1.0; 

	bVERY_STABILITY_ON=false;

	// инициализация.
#pragma omp parallel for
    for (integer i=0; i<maxelm+maxbound; i++) {
		tau[VELOCITY_X_COMPONENT][i]=0.0; // псевдовремя.
		tau[VELOCITY_Y_COMPONENT][i]=0.0;
		tau[VELOCITY_Z_COMPONENT][i]=0.0;
	}

#pragma omp parallel for
	for (integer iP=0; iP<maxelm; iP++) {

        // iP - номер центрального контрольного объёма
	    integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	    iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
	
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
	    // если с одной из сторон граница расчётной области 
	    // то переменная равна true
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
    	doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
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

	    //printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	    //getchar();
		
		// ВНИМАНИЕ ! БУДЕМ СОБИРАТЬ tau ТОЛЬКО ПО КОНВЕКТИВНО ДИФФУЗИОННОЙ СОСТАВЛЯЮЩЕЙ.
		// Если конвективная составляющая полностью равна нулю то только по полностью незанулённой  диффузионной.
		doublereal apX=0.0, apY=0.0, apZ=0.0;
		apX=sumanb[VELOCITY_X_COMPONENT][iP];
		apY=sumanb[VELOCITY_Y_COMPONENT][iP];
		apZ=sumanb[VELOCITY_Z_COMPONENT][iP];

		/*
		if ((inumberSIMPLEiteration==2)&&(iP==3477)) {
			printf("%e, %e, %e\n",sumanb[VX][iP],sumanb[VY][iP],sumanb[VZ][iP]);
			getchar();
		}
		*/

		//doublereal CFL2[3]={4.0,4.0,4.0}; // соответствует SIMPLEC alpha=0.8; 
		doublereal CFL2[3];
		for (integer i9=0; i9<=2; i9++) {
			CFL2[i9]=4.0;
		}
		if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
		    // SIMPLEC 
		    CFL2[VELOCITY_X_COMPONENT]=alpha[VELOCITY_X_COMPONENT]/(1.0-alpha[VELOCITY_X_COMPONENT]);
			CFL2[VELOCITY_Y_COMPONENT]=alpha[VELOCITY_Y_COMPONENT]/(1.0-alpha[VELOCITY_Y_COMPONENT]);
			CFL2[VELOCITY_Z_COMPONENT]=alpha[VELOCITY_Z_COMPONENT]/(1.0-alpha[VELOCITY_Z_COMPONENT]);
	    }
	    if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
		    // SIMPLE
		    CFL2[VELOCITY_X_COMPONENT]=alpha[VELOCITY_X_COMPONENT];
			CFL2[VELOCITY_Y_COMPONENT]=alpha[VELOCITY_Y_COMPONENT];
			CFL2[VELOCITY_Z_COMPONENT]=alpha[VELOCITY_Z_COMPONENT];
	    }
		//doublereal tau_loc[3]={CFL2[VX]*prop[RHO][iP]*dx*dy*dz/apX,CFL2[VY]*prop[RHO][iP]*dx*dy*dz/apY,CFL2[VZ]*prop[RHO][iP]*dx*dy*dz/apZ};
		doublereal tau_loc[3];
		for (integer i9=0; i9<=2; i9++) {
			tau_loc[i9]=CFL2[i9]*prop[RHO][iP]*dx*dy*dz/sumanb[i9][iP];
		}
		
		if (boldscheme) {
			
			if (!btimedep) {
				// Стационарный решатель.
			   // работающий в данный момент вариант. 10 апреля 2012.
			   tau[VELOCITY_X_COMPONENT][iP]=tau_loc[VELOCITY_X_COMPONENT];
			   tau[VELOCITY_Y_COMPONENT][iP]=tau_loc[VELOCITY_Y_COMPONENT];
               tau[VELOCITY_Z_COMPONENT][iP]=tau_loc[VELOCITY_Z_COMPONENT];
			}
			else {
				// нестационарный решатель.
				// 7 мая 2013 года.
				// Так рекомендует делать Гаврилов Андрей в нестационарной постановке.
				// dtime - реальный шаг по времени.
				// tau_loc - псевдовремя.
				if (dtime<=0.0) {
					printf("negativ or zero dtime in unsteady solver tau calculation...\n");
					printf("please, press any key to exit...\n");
					//getchar();
					system("pause");
					exit(0);
				}
				// Так делать не следует так как в файле pamendment всё завязано на стандартное tau равное tau_loc.
				// Внимание !!! нужно быть очень осторожным. Нельзя допускать внутренних противоречий.
               //tau[VX][iP]=CFL1/(1.0/dtime + 1.0/tau_loc[VX]); // нельзя.
			   //tau[VY][iP]=CFL1/(1.0/dtime + 1.0/tau_loc[VY]); // нельзя.
               //tau[VZ][iP]=CFL1/(1.0/dtime + 1.0/tau_loc[VZ]); // нельзя.
			   // Согласованный с файлом pamendment вариант. 7 мая 2013.
			   tau[VELOCITY_X_COMPONENT][iP]=tau_loc[VELOCITY_X_COMPONENT];
			   tau[VELOCITY_Y_COMPONENT][iP]=tau_loc[VELOCITY_Y_COMPONENT];
               tau[VELOCITY_Z_COMPONENT][iP]=tau_loc[VELOCITY_Z_COMPONENT];
			   
			}
		}
		else {
		    // Внимание для повышения устойчивости рекомендуется увеличить параметр CFL1.
		    // Увеличивая CFL2 можно уменьшать влияние локальных параметров течения на поправку давления,
		    // чтобы увеличить CFL2 рекомендуется сделать параметр нижней релаксации для скорости сделать более близким к 1.0.
		    // Например:
		    // alpha CFL2
		    // 0.7 2.33333
		    // 0.5 1.0 // Рекомендовано С. Патанкаром по умолчанию.
		    // 0.8 4.0 // Рекомендовано Гавриловым Андреем по умолчанию.
		    // 0.85 5.666
		    if (bVERY_STABILITY_ON==true) {
			   // пока просто запоминаем, оновное вычисление псевдовремени будет позже.
			   tau[VELOCITY_X_COMPONENT][iP]=tau_loc[VELOCITY_X_COMPONENT];
			   tau[VELOCITY_Y_COMPONENT][iP]=tau_loc[VELOCITY_Y_COMPONENT];
			   tau[VELOCITY_Z_COMPONENT][iP]=tau_loc[VELOCITY_Z_COMPONENT];
	 	    }
	    	else {
		       if (btimedep) {
			       // нестационарная задача.
				   tau[VELOCITY_X_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc[VELOCITY_X_COMPONENT]));
				   tau[VELOCITY_Y_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc[VELOCITY_Y_COMPONENT]));
				   tau[VELOCITY_Z_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc[VELOCITY_Z_COMPONENT]));
		       }
		   else {
			   // стационарный солвер.
			   // в этом случае также можно применить формулу для нестационарного случая,
			   // вопрос в том что шаг по времени dtime теперь надо подбирать.
			   if (dtime<=0.0) {
				   // возьмём среднее арифметическое от всех существующих tau.
			       tau[VELOCITY_X_COMPONENT][iP]=tau_loc[VELOCITY_X_COMPONENT]; // пока просто запомним чтобы потом можно было вычислить среднее арифметическое.
                   tau[VELOCITY_Y_COMPONENT][iP]=tau_loc[VELOCITY_Y_COMPONENT];
				   tau[VELOCITY_Z_COMPONENT][iP]=tau_loc[VELOCITY_Z_COMPONENT];
			   }
			   else {
				   // если параметр dtime положителен то мы считаем что он задан пользователем из соображений пользователя.
				   tau[VELOCITY_X_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc[VELOCITY_X_COMPONENT]));
                   tau[VELOCITY_Y_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc[VELOCITY_Y_COMPONENT]));
                   tau[VELOCITY_Z_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_loc[VELOCITY_Z_COMPONENT]));
			   }
		   }
		}
		}

	} // next iP

	if (!boldscheme) {
	   if (bVERY_STABILITY_ON==true) {
		   // во всей области просто постоянное значение шага по псевдовремени.
		   doublereal tau_avg[3]={0.0,0.0,0.0};
		   for (integer iP=0; iP<maxelm; iP++) {
		       tau_avg[VELOCITY_X_COMPONENT]+=tau[VELOCITY_X_COMPONENT][iP];
               tau_avg[VELOCITY_Y_COMPONENT]+=tau[VELOCITY_Y_COMPONENT][iP];
               tau_avg[VELOCITY_Z_COMPONENT]+=tau[VELOCITY_Z_COMPONENT][iP];
		   }
		   tau_avg[VELOCITY_X_COMPONENT]=tau_avg[VELOCITY_X_COMPONENT]/((doublereal)(maxelm));
		   tau_avg[VELOCITY_Y_COMPONENT]=tau_avg[VELOCITY_Y_COMPONENT]/((doublereal)(maxelm));
		   tau_avg[VELOCITY_Z_COMPONENT]=tau_avg[VELOCITY_Z_COMPONENT]/((doublereal)(maxelm));

		   if (btimedep) {
			   // нестационарная задача.
			   for (integer iP=0; iP<maxelm; iP++) {
			       tau[VELOCITY_X_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg[VELOCITY_X_COMPONENT])); // постоянный шаг по псевдовремени во всей расчётной области.
				   tau[VELOCITY_Y_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg[VELOCITY_Y_COMPONENT]));
				   tau[VELOCITY_Z_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg[VELOCITY_Z_COMPONENT]));
			   }
		   }
		   else {
			   if (dtime<=0.0) {
				  // пользователь не определил дополнительный временной параметр в стационарном случае.
				  for (integer iP=0; iP<maxelm; iP++) {
					  tau[VELOCITY_X_COMPONENT][iP]=CFL1*tau_avg[VELOCITY_X_COMPONENT];
					  tau[VELOCITY_Y_COMPONENT][iP]=CFL1*tau_avg[VELOCITY_Y_COMPONENT];
					  tau[VELOCITY_Z_COMPONENT][iP]=CFL1*tau_avg[VELOCITY_Z_COMPONENT];

				  }
			   }
			   else {
				    // пользователь определил дополнительный опорный стабилизирующий временной параметр.
				    for (integer iP=0; iP<maxelm; iP++) {
			            tau[VELOCITY_X_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg[VELOCITY_X_COMPONENT])); // постоянный шаг по псевдовремени во всей расчётной области.
					    tau[VELOCITY_Y_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg[VELOCITY_Y_COMPONENT]));
					    tau[VELOCITY_Z_COMPONENT][iP]=CFL1*(1.0/(1.0/dtime+1.0/tau_avg[VELOCITY_Z_COMPONENT]));
			        }
			   }
		   }
	   }
	   else {
		   // bVERY_STABILITY_ON==false.
		   // Обычный рабочий режим, всё-таки более стабильный чем простое вычисление tau ~ rho*Vol*alpha/ap или tau ~ rho*Vol*alpha/(ap*(1.0-alpha)).

		   if (dtime<=0.0) {

			   //printf("relax tau incomming...");
			   //getchar();


			  doublereal tau_avg[3]={0.0,0.0,0.0};
		      for (integer iP=0; iP<maxelm; iP++) {
			      tau_avg[VELOCITY_X_COMPONENT]+=tau[VELOCITY_X_COMPONENT][iP];
			      tau_avg[VELOCITY_Y_COMPONENT]+=tau[VELOCITY_Y_COMPONENT][iP];
			      tau_avg[VELOCITY_Z_COMPONENT]+=tau[VELOCITY_Z_COMPONENT][iP];
		      }
		      tau_avg[VELOCITY_X_COMPONENT]=tau_avg[VELOCITY_X_COMPONENT]/((doublereal)(maxelm));
              tau_avg[VELOCITY_Y_COMPONENT]=tau_avg[VELOCITY_Y_COMPONENT]/((doublereal)(maxelm));
		      tau_avg[VELOCITY_Z_COMPONENT]=tau_avg[VELOCITY_Z_COMPONENT]/((doublereal)(maxelm));

		      for (integer iP=0; iP<maxelm; iP++) {
			      tau[VELOCITY_X_COMPONENT][iP]=CFL1*(1.0/(1.0/tau_avg[VELOCITY_X_COMPONENT]+1.0/tau[VELOCITY_X_COMPONENT][iP]));
			      tau[VELOCITY_Y_COMPONENT][iP]=CFL1*(1.0/(1.0/tau_avg[VELOCITY_Y_COMPONENT]+1.0/tau[VELOCITY_Y_COMPONENT][iP]));
			      tau[VELOCITY_Z_COMPONENT][iP]=CFL1*(1.0/(1.0/tau_avg[VELOCITY_Z_COMPONENT]+1.0/tau[VELOCITY_Z_COMPONENT][iP]));
		      }
	       }
	   }
	}

	if (1) {
		// продолжение параметра tau на границу расчётной области.
		// Просто сносим значение tau из ближайшего внутреннего узла.
		// iP - номер центрального контрольного объёма
#pragma omp parallel for
		for (integer iP = 0; iP < maxelm; iP++) {
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

			if (iE > -1) {
				if (iE >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iE] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iE] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iW > -1) {
				if (iW >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iW] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iW] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iN > -1) {
				if (iN >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iN] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iN] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iS > -1) {
				if (iS >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iS] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iS] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}
			if (iT > -1) {
				if (iT >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iT] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iT] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iB > -1) {
				if (iB >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iB] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iB] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iE2 > -1) {
				if (iE2 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE2] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iE2] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iE2] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iW2 > -1) {
				if (iW2 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW2] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iW2] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iW2] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iN2 > -1) {
				if (iN2 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN2] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iN2] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iN2] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iS2 > -1) {
				if (iS2 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS2] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iS2] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iS2] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}
			if (iT2 > -1) {
				if (iT2 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT2] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iT2] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iT2] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iB2 > -1) {
				if (iB2 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB2] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iB2] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iB2] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iE3 > -1) {
				if (iE3 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE3] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iE3] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iE3] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iW3 > -1) {
				if (iW3 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW3] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iW3] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iW3] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iN3 > -1) {
				if (iN3 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN3] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iN3] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iN3] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iS3 > -1) {
				if (iS3 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS3] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iS3] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iS3] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}
			if (iT3 > -1) {
				if (iT3 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT3] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iT3] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iT3] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iB3 > -1) {
				if (iB3 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB3] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iB3] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iB3] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iE4 > -1) {
				if (iE4 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE4] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iE4] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iE4] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iW4 > -1) {
				if (iW4 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW4] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iW4] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iW4] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iN4 > -1) {
				if (iN4 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN4] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iN4] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iN4] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iS4 > -1) {
				if (iS4 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS4] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iS4] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iS4] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}
			if (iT4 > -1) {
				if (iT4 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT4] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iT4] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iT4] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

			if (iB4 > -1) {
				if (iB4 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB4] = tau[VELOCITY_X_COMPONENT][iP];
					tau[VELOCITY_Y_COMPONENT][iB4] = tau[VELOCITY_Y_COMPONENT][iP];
					tau[VELOCITY_Z_COMPONENT][iB4] = tau[VELOCITY_Z_COMPONENT][iP];
				}
			}

		}
	}
	else {
		// на основе диагональных коэффициентов в граничном условии.
#pragma omp parallel for
		for (integer iP = 0; iP < maxelm; iP++) {
			integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
			iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
			iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

			integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
			iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
			iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

			integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
			iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
			iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

			integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
			iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
			iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0; // размеры контрольного объёма
			volume3D(iP, nvtx, pa, dx, dy, dz);

			doublereal Vol = dx*dy*dz;
			doublereal CFL2[3] = {4.0, 4.0, 4.0};
			for (integer i9 = 0; i9 <= 2; i9++) {
				CFL2[i9] = 4.0;
			}

			if (iSIMPLE_alg == SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
				// SIMPLEC 
				CFL2[VELOCITY_X_COMPONENT] = alpha[VELOCITY_X_COMPONENT] / (1.0 - alpha[VELOCITY_X_COMPONENT]);
				CFL2[VELOCITY_Y_COMPONENT] = alpha[VELOCITY_Y_COMPONENT] / (1.0 - alpha[VELOCITY_Y_COMPONENT]);
				CFL2[VELOCITY_Z_COMPONENT] = alpha[VELOCITY_Z_COMPONENT] / (1.0 - alpha[VELOCITY_Z_COMPONENT]);
			}
			if (iSIMPLE_alg == SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
				// SIMPLE
				CFL2[VELOCITY_X_COMPONENT] = alpha[VELOCITY_X_COMPONENT];
				CFL2[VELOCITY_Y_COMPONENT] = alpha[VELOCITY_Y_COMPONENT];
				CFL2[VELOCITY_Z_COMPONENT] = alpha[VELOCITY_Z_COMPONENT];
			}

			if (iE > -1) {
				if (iE >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iE - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iE];
					tau[VELOCITY_Y_COMPONENT][iE] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iE - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iE];
					tau[VELOCITY_Z_COMPONENT][iE] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iE - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iE];
				}
			}

			if (iW > -1) {
				if (iW >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iW - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iW];
					tau[VELOCITY_Y_COMPONENT][iW] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iW - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iW];
					tau[VELOCITY_Z_COMPONENT][iW] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iW - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iW];
				}
			}

			if (iN > -1) {
				if (iN >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iN - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iN];
					tau[VELOCITY_Y_COMPONENT][iN] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iN - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iN];
					tau[VELOCITY_Z_COMPONENT][iN] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iN - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iN];
				}
			}

			if (iS > -1) {
				if (iS >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iS - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iS];
					tau[VELOCITY_Y_COMPONENT][iS] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iS - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iS];
					tau[VELOCITY_Z_COMPONENT][iS] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iS - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iS];
				}
			}

			if (iT > -1) {
				if (iT >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iT - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iT];
					tau[VELOCITY_Y_COMPONENT][iT] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iT - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iT];
					tau[VELOCITY_Z_COMPONENT][iT] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iT - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iT];
				}
			}

			if (iB > -1) {
				if (iB >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iB - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iB];
					tau[VELOCITY_Y_COMPONENT][iB] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iB - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iB];
					tau[VELOCITY_Z_COMPONENT][iB] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iB - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iB];
				}
			}

			if (iE2 > -1) {
				if (iE2 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE2] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iE2 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iE2];
					tau[VELOCITY_Y_COMPONENT][iE2] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iE2 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iE2];
					tau[VELOCITY_Z_COMPONENT][iE2] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iE2 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iE2];
				}
			}

			if (iW2 > -1) {
				if (iW2 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW2] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iW2 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iW2];
					tau[VELOCITY_Y_COMPONENT][iW2] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iW2 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iW2];
					tau[VELOCITY_Z_COMPONENT][iW2] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iW2 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iW2];
				}
			}

			if (iN2 > -1) {
				if (iN2 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN2] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iN2 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iN2];
					tau[VELOCITY_Y_COMPONENT][iN2] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iN2 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iN2];
					tau[VELOCITY_Z_COMPONENT][iN2] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iN2 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iN2];
				}
			}

			if (iS2 > -1) {
				if (iS2 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS2] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iS2 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iS2];
					tau[VELOCITY_Y_COMPONENT][iS2] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iS2 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iS2];
					tau[VELOCITY_Z_COMPONENT][iS2] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iS2 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iS2];
				}
			}

			if (iT2 > -1) {
				if (iT2 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT2] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iT2 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iT2];
					tau[VELOCITY_Y_COMPONENT][iT2] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iT2 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iT2];
					tau[VELOCITY_Z_COMPONENT][iT2] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iT2 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iT2];
				}
			}

			if (iB2 > -1) {
				if (iB2 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB2] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iB2 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iB2];
					tau[VELOCITY_Y_COMPONENT][iB2] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iB2 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iB2];
					tau[VELOCITY_Z_COMPONENT][iB2] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iB2 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iB2];
				}
			}

			if (iE3 > -1) {
				if (iE3 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE3] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iE3 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iE3];
					tau[VELOCITY_Y_COMPONENT][iE3] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iE3 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iE3];
					tau[VELOCITY_Z_COMPONENT][iE3] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iE3 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iE3];
				}
			}

			if (iW3 > -1) {
				if (iW3 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW3] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iW3 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iW3];
					tau[VELOCITY_Y_COMPONENT][iW3] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iW3 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iW3];
					tau[VELOCITY_Z_COMPONENT][iW3] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iW3 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iW3];
				}
			}

			if (iN3 > -1) {
				if (iN3 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN3] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iN3 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iN3];
					tau[VELOCITY_Y_COMPONENT][iN3] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iN3 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iN3];
					tau[VELOCITY_Z_COMPONENT][iN3] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iN3 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iN3];
				}
			}

			if (iS3 > -1) {
				if (iS3 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS3] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iS3 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iS3];
					tau[VELOCITY_Y_COMPONENT][iS3] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iS3 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iS3];
					tau[VELOCITY_Z_COMPONENT][iS3] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iS3 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iS3];
				}
			}

			if (iT3 > -1) {
				if (iT3 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT3] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iT3 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iT3];
					tau[VELOCITY_Y_COMPONENT][iT3] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iT3 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iT3];
					tau[VELOCITY_Z_COMPONENT][iT3] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iT3 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iT3];
				}
			}

			if (iB3 > -1) {
				if (iB3 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB3] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iB3 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iB3];
					tau[VELOCITY_Y_COMPONENT][iB3] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iB3 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iB3];
					tau[VELOCITY_Z_COMPONENT][iB3] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iB3 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iB3];
				}
			}

			if (iE4 > -1) {
				if (iE4 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iE4] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iE4 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iE4];
					tau[VELOCITY_Y_COMPONENT][iE4] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iE4 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iE4];
					tau[VELOCITY_Z_COMPONENT][iE4] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iE4 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iE4];
				}
			}

			if (iW4 > -1) {
				if (iW4 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iW4] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iW4 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iW4];
					tau[VELOCITY_Y_COMPONENT][iW4] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iW4 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iW4];
					tau[VELOCITY_Z_COMPONENT][iW4] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iW4 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iW4];
				}
			}

			if (iN4 > -1) {
				if (iN4 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iN4] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iN4 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iN4];
					tau[VELOCITY_Y_COMPONENT][iN4] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iN4 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iN4];
					tau[VELOCITY_Z_COMPONENT][iN4] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iN4 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iN4];
				}
			}

			if (iS4 > -1) {
				if (iS4 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iS4] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iS4 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iS4];
					tau[VELOCITY_Y_COMPONENT][iS4] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iS4 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iS4];
					tau[VELOCITY_Z_COMPONENT][iS4] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iS4 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iS4];
				}
			}

			if (iT4 > -1) {
				if (iT4 >= maxelm) {
					// граничный узел
					tau[VELOCITY_X_COMPONENT][iT4] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iT4 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iT4];
					tau[VELOCITY_Y_COMPONENT][iT4] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iT4 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iT4];
					tau[VELOCITY_Z_COMPONENT][iT4] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iT4 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iT4];
				}
			}

			if (iB4 > -1) {
				if (iB4 >= maxelm) {
					// boundary node
					tau[VELOCITY_X_COMPONENT][iB4] = CFL2[VELOCITY_X_COMPONENT] * prop_b[RHO][iB4 - maxelm] * Vol / sumanb[VELOCITY_X_COMPONENT][iB4];
					tau[VELOCITY_Y_COMPONENT][iB4] = CFL2[VELOCITY_Y_COMPONENT] * prop_b[RHO][iB4 - maxelm] * Vol / sumanb[VELOCITY_Y_COMPONENT][iB4];
					tau[VELOCITY_Z_COMPONENT][iB4] = CFL2[VELOCITY_Z_COMPONENT] * prop_b[RHO][iB4 - maxelm] * Vol / sumanb[VELOCITY_Z_COMPONENT][iB4];
				}
			}
		
		}
	}

} // tau_calc3

#endif
// В файле correct_velocity.cpp содержится код который применяется
// для коррекции скорости так чтобы они стали удовлетворять уравнению 
// неразрывности.

// Данные функции выделены в отдельный файл 31 марта 2012 года.

// 
// Коррекцию скорости надо разделить на две процедуры.
// Первая процедура для внутренних КО. При этом используется интерполяция для 
// значений поправки давления на границе расчётной области. 
// 
// Вторая процедура коррекции для граничных КО. С интерполяцией скорости изнутри области на границу.
// Коррекцию граничных КО можно и отключить, т.к. она может не понадобиться.
// Возможно расчёт будет правилен даже в случае отсутствия коррекции граничных узлов.
// Это нужно проверить непосредственно путём численного расчёта.
// На данный момент решено включить коррекцию скорости на границе путём интерполяции изнутри области.

// Коррекция граничных КО включена.
// По поводу коррекции граничных КО смотри mysolver SIMPLE Algorithm.
// В этом файле исправлены все интерполяционные формулы с 
// учётом предельного перехода к предельным случаям.
// Квадратичная интерполяция, увеличивающая скорость сходимости, введена в этом модуле 13,14,15 мая 2012 года.

#ifndef CORRECT_VELOCITYv_0_07_CPP
#define CORRECT_VELOCITYv_0_07_CPP 1



// коррекция скорости для внутренних КО.
// скорректированные скорости удовлетворяют уравнению 
// несжимаемости.
void correct_internal_volume(integer iP, integer iVar, equation3D** sl,   
			 int** nvtx, doublereal** &potent, integer maxelm, doublereal* alpha,
			 TOCHKA* pa, int*** neighbors_for_the_internal_node, integer iternumber,
	ORDER_DERIVATIVE iderivative_pressure) {

    // квадратичная интерполяция добавлена 14 мая 2012 года.
	// Ещё как вариант можно применить МНК интерполяцию (построить прямую по трём точкам).
	integer interpol=0; // 0 - без всякой интерполяции, 1- линейная интерполяция, 2- квадратичная интерполяция.

	// Внимание можно использовать и второй порядок тоже. формула из книги Г.З.Гарбера  
	// работает. Формула из книги Г.З. Гарбера обеспечивает точность производной от 
	// поправки давления второго порядка на неравномерной сетке (она найдена как производная 
	// от параболы построенная по трём точкам).
	// порядок точности нахождения первой производной от давления.
	//integer iderivative_pressure=1; // 1 - первый порядок, 2 - второй порядок.


    // Внутренний узел и его соседи:
    // iP принадлежит интервалу 0..maxelm-1

    // iP - номер центрального контрольного объёма
	integer iE=-1, iN=-1, iT=-1, iW=-1, iS=-1, iB=-1; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
	// индексировать соседей для поправки давления 
	// ненужно, т.к. они уже проиндексированы в 
	// my_elmatr_quad_PAm...
	//sl[PAM][iP].iP=iP;
	//sl[PAM][iP].iE=iE; sl[PAM][iP].iN=iN; 
	//sl[PAM][iP].iS=iS; sl[PAM][iP].iW=iW;
    //sl[PAM][iP].iT=iT; sl[PAM][iP].iB=iB;

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


	doublereal feplus=0.0, fwplus=0.0, fnplus=0.0, fsplus=0.0, ftplus=0.0, fbplus=0.0;
	// x-direction
	feplus=0.5*dx/dxe; // if(bE) then feplus=1.0;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;

	//printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	//getchar();


	doublereal ds=0.0, dl=0.0, dv=0.0; // площадь грани, длина интервала и объём контрольного объёма.
	switch (iVar) {
       case VELOCITY_X_COMPONENT: ds=dy*dz; dv=ds*dx; dl=dx; break;
       case VELOCITY_Y_COMPONENT: ds=dx*dz; dv=ds*dy; dl=dy;  break;
       case VELOCITY_Z_COMPONENT: ds=dx*dy; dv=ds*dz; dl=dz;  break;
	}

	
	doublereal PAmP=0.0, PAmW=0.0, PAmE=0.0, PAmS=0.0, PAmN=0.0, PAmB=0.0, PAmT=0.0;
	PAmP=potent[PAM][iP];

	// поправка давления известна
	// как во внутренних так и в 
	// граничных узлах.
	// поправка давления и само давление должны быть интерполированы на границу
	// с помощью линейной или квадратичной интерполяции (квадратичная наверно лучше).
	// Однородное условие Неймана на границе для поправки давления не оправдывает себя,
	// т.к. в приграничном узле возникает нефизичный максимум скорости.
	// Проблему нужно решать по-видимому не в блоке коррекции скорости, а сразу после решения
	// уравнения для поправки давления. Тогда линейная интерполяция здесь и везде будет выполнена автоматически.
	/*
	PAmW=potent[PAM][iW];
	PAmE=potent[PAM][iE];
	PAmS=potent[PAM][iS];
	PAmN=potent[PAM][iN];
	PAmB=potent[PAM][iB];
	PAmT=potent[PAM][iT];//*/

	doublereal hxminus=0.0, hxplus=0.0, hyminus=0.0, hyplus=0.0, hzminus=0.0, hzplus=0.0;

	//*
    if (!bW) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iW, nvtx, pa, pb,W_SIDE);
		hxminus=fabs(pp.x-pb.x);

		PAmW=potent[PAM][iW];
	} else {
		// узел W граничный
		TOCHKA pp, pb,pbb;
		center_cord3D(iP, nvtx, pa, pp,100);
		hxminus=0.5*dx;

		switch (interpol) {
		case 0: // значение которое получено после решения СЛАУ.
			    PAmW=potent[PAM][iW]; 
			    break;
		case 1: // линейная интерполяция
		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iE, nvtx, pa, pb,E_SIDE);
		         PAmW=my_linear_interpolation('-', potent[PAM][iP], potent[PAM][iE], pp.x, pb.x, pp.x-0.5*dx); 
			    break;
		case 2: // квадратичная интерполяция.

		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iE, nvtx, pa, pb,E_SIDE);
			     center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb,EE_SIDE);
					
			     PAmW=my_quadratic_interpolation('-', potent[PAM][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[PAM][iE], potent[PAM][iP], pbb.x , pb.x, pp.x, pp.x-0.5*dx); 
			    break;
		default: // значение которое получено после решения СЛАУ.
			      PAmW=potent[PAM][iW];
			    break;
		} // end switch
		
	}
	if (!bE) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iE, nvtx, pa, pb,E_SIDE);
		hxplus=fabs(pb.x-pp.x);

		PAmE=potent[PAM][iE];
	} else {
        // узел E граничный
		hxplus=0.5*dx;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0: // значение которое получено после решения СЛАУ. 
			    PAmE=potent[PAM][iE];
			    break;
		case 1: // линейная интерполяция

		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iW, nvtx, pa, pb,W_SIDE);
		         PAmE=my_linear_interpolation('+', potent[PAM][iP], potent[PAM][iW], pp.x, pb.x, pp.x+0.5*dx);

			    break;
		case 2: // квадратичная интерполяция.
						     
		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iW, nvtx, pa, pb,W_SIDE);
				 center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb,WW_SIDE);
					
				 PAmE = my_quadratic_interpolation('+', potent[PAM][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[PAM][iW], potent[PAM][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			    break;
		default: // значение которое получено после решения СЛАУ.
			     PAmE=potent[PAM][iE]; 
			     break;
		}

	}
	if (!bS) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iS, nvtx, pa, pb,S_SIDE);
		hyminus=fabs(pp.y-pb.y);

		PAmS=potent[PAM][iS];
	} else {
		// узел S граничный
		hyminus=0.5*dy;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0: // значение которое получено после решения СЛАУ.  
			    PAmS=potent[PAM][iS];
			    break;
		case 1:
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iN, nvtx, pa, pb,N_SIDE);
		        PAmS=my_linear_interpolation('-', potent[PAM][iP], potent[PAM][iN], pp.y, pb.y, pp.y-0.5*dy);

			    break;
		case 2:  
			    // квадратичная интерполяция.

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iN, nvtx, pa, pb,N_SIDE);
				center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb,NN_SIDE);
					
				PAmS = my_quadratic_interpolation('-', potent[PAM][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[PAM][iN], potent[PAM][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);

			    break;
		default: // значение которое получено после решения СЛАУ. 
			    PAmS=potent[PAM][iS];
			    break;
		}

	}
	if (!bN) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iN, nvtx, pa, pb,N_SIDE);
		hyplus=fabs(pb.y-pp.y);

		PAmN=potent[PAM][iN];
	} else {
		// узел N граничный
		hyplus=0.5*dy;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0: // значение которое получено после решения СЛАУ. 
			    PAmN=potent[PAM][iN];
			    break;
		case 1:
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iS, nvtx, pa, pb,S_SIDE);
		        PAmN=my_linear_interpolation('+', potent[PAM][iP], potent[PAM][iS], pp.y, pb.y, pp.y+0.5*dy);
			    break;
		case 2:
			    // квадратичная интерполяция.

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iS, nvtx, pa, pb,S_SIDE);
				center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb,SS_SIDE);
					
				PAmN = my_quadratic_interpolation('+', potent[PAM][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[PAM][iS], potent[PAM][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			    break;
		default: 
			    // значение которое получено после решения СЛАУ. 
			    PAmN=potent[PAM][iN];
			    break;
		}

		 
	} 

	if (!bB) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iB, nvtx, pa, pb,B_SIDE);
		hzminus=fabs(pp.z-pb.z);

		PAmB=potent[PAM][iB];
	} else {
		// узел B граничный
		hzminus=0.5*dz;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0:
			    // значение которое получено после решения СЛАУ. 
			    PAmB=potent[PAM][iB];
			    break;
		case 1: 
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iT, nvtx, pa, pb,T_SIDE);
		        PAmB=my_linear_interpolation('-', potent[PAM][iP], potent[PAM][iT], pp.z, pb.z, pp.z-0.5*dz);
			    break;
		case 2: // квадратичная интерполяция.

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iT, nvtx, pa, pb,T_SIDE);
				center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb,TT_SIDE);
					
				PAmB = my_quadratic_interpolation('-', potent[PAM][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[PAM][iT], potent[PAM][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			    break;
		default:
			    // значение которое получено после решения СЛАУ. 
			    PAmB=potent[PAM][iB];
			    break;
		}

	}
	if (!bT) { 
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iT, nvtx, pa, pb,T_SIDE);
		hzplus=fabs(pb.z-pp.z);

		PAmT=potent[PAM][iT];
	} else {
		// узел T граничный
		hzplus=0.5*dz;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0: // значение которое получено после решения СЛАУ. 
			    PAmT=potent[PAM][iT];
			    break;
		case 1: 
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iB, nvtx, pa, pb,B_SIDE);
		        PAmT=my_linear_interpolation('+', potent[PAM][iP], potent[PAM][iB], pp.z, pb.z, pp.z+0.5*dz);
			    break;
		case 2:// квадратичная интерполяция.
                
		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iB, nvtx, pa, pb,B_SIDE);
				center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb,BB_SIDE);
					
				PAmT = my_quadratic_interpolation('+', potent[PAM][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[PAM][iB], potent[PAM][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			    break;
		default: // значение которое получено после решения СЛАУ. 
			     PAmT=potent[PAM][iT];
			    break;
		} 

	}
	//*/

	// Случай граничного узла G учитывается правильно,
	// т.к. в этом случае fgplus==1.0; 

	// Линейная интерполяция давления на грань КО.
	doublereal deltaP=0.0, gradP=0.0;
	switch (iVar) {
	case VELOCITY_X_COMPONENT: if (iderivative_pressure== ORDER_DERIVATIVE::FIRST_ORDER) {
			          // естественная аппроксимация первого порядка.
				      deltaP=(fwplus*PAmW+(1-fwplus)*PAmP);
			          deltaP-=(feplus*PAmE+(1-feplus)*PAmP); 
				      gradP=deltaP/dl; // первая производная от давления.
				  }
				  else if (iderivative_pressure== ORDER_DERIVATIVE::SECOND_ORDER) {
					  // аппроксимация второго порядка точности.
					  // обязательно нужен знак минус иначе скорость
					  // будет направлена в другую сторону.
					  gradP=-rgradF(PAmW, PAmP, PAmE, hxminus, hxplus);
				  }
				  break;
		case VELOCITY_Y_COMPONENT: if (iderivative_pressure== ORDER_DERIVATIVE::FIRST_ORDER) {
			          // естественная аппроксимация первого порядка.
				      deltaP=(fsplus*PAmS+(1-fsplus)*PAmP);
			          deltaP-=(fnplus*PAmN+(1-fnplus)*PAmP);
				      gradP=deltaP/dl; // первая производная от давления.
				  }
				  else if (iderivative_pressure== ORDER_DERIVATIVE::SECOND_ORDER) {
					  // аппроксимация второго порядка точности.
					  // обязательно нужен знак минус иначе скорость
					  // будет направлена в другую сторону.
					  gradP=-rgradF(PAmS, PAmP, PAmN, hyminus, hyplus);
				  }
			      break;
        case VELOCITY_Z_COMPONENT: if (iderivative_pressure== ORDER_DERIVATIVE::FIRST_ORDER) {
			          // естественная аппроксимация первого порядка.
			          deltaP=(fbplus*PAmB+(1-fbplus)*PAmP);
			          deltaP-=(ftplus*PAmT+(1-ftplus)*PAmP);
				      gradP=deltaP/dl; // первая производная от давления.
				  }
				  else if (iderivative_pressure== ORDER_DERIVATIVE::SECOND_ORDER) {
					  // аппроксимация второго порядка точности.
					  // обязательно нужен знак минус иначе скорость
					  // будет направлена в другую сторону.
				      gradP=-rgradF(PAmB, PAmP, PAmT, hzminus, hzplus);
				      //printf("%e\n",PAmB); getchar();
				  }
				  break;
	}

    // коррекция компоненты скорости во внутреннем узле iP. 

	// коррекция скорости не должна подвергаться нижней релаксации.
	// Так предлагает делать С. Патанкар. Это должно быть согласовано с составлением уравнения
	// для поправки давления.
	//potent[iVar][iP]+=ds*(deltaP)/sl[iVar][iP].ap;//alpha[iVar]*
	// Так предлагает делать Гаврилов Андрей в расчётном комплексе Sigma-flow.
	// Это должно быть согласовано с составлением уравнения для поправки давления.
	// potent[iVar][iP]+=(tau/rho)*deltaP;
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
		// SIMPLE алгоритм: Carretto et al., 1973.
		// tau ~ (alpha[iVar]*rho*ds)/(sl[iVar][iP].ap);
		//potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
	}
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC алгоритм: Van Doormal and Raithby., 1984
		// tau ~ (alpha[iVar]*rho*ds)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
        //potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
	}

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скорректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume

// коррекция скорости для внутренних КО.
// скорректированные скорости удовлетворяют уравнению 
// несжимаемости.
void correct_internal_volume2(integer iP, integer iVar, equation3D** sl,   
			 int** nvtx, doublereal** &potent, integer maxelm, doublereal* alpha,
			 TOCHKA* pa, int*** neighbors_for_the_internal_node, integer iternumber) {

	// Процедура коррекции скорости во внутренних контрольных объёмах на 
	// основе осреднённого градиента давления.
    
    // Внутренний узел и его соседи:
    // iP принадлежит интервалу 0..maxelm-1

	
	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
    volume3D(iP, nvtx, pa, dx, dy, dz);
	    

	doublereal ds=0.0, dv=0.0; // площадь грани, длина интервала и объём контрольного объёма.
	switch (iVar) {
       case VELOCITY_X_COMPONENT: ds=dy*dz; dv=ds*dx; break;
       case VELOCITY_Y_COMPONENT: ds=dx*dz; dv=ds*dy; break;
       case VELOCITY_Z_COMPONENT: ds=dx*dy; dv=ds*dz; break;
	}
	
	// Случай граничного узла G учитывается правильно,
	// т.к. в этом случае fgplus==1.0; 

	// Линейная интерполяция давления на грань КО.
	doublereal deltaP=0.0, gradP=0.0;
	switch (iVar) {
		case VELOCITY_X_COMPONENT: gradP=-potent[GRADXPAM][iP];
				  break;
		case VELOCITY_Y_COMPONENT: gradP=-potent[GRADYPAM][iP];
				  break;
        case VELOCITY_Z_COMPONENT: gradP=-potent[GRADZPAM][iP];
				  break;
	}


    // коррекция компоненты скорости во внутреннем узле iP. 

	// коррекция скорости не должна подвергаться нижней релаксации.
	// Так предлагает делать С. Патанкар. Это должно быть согласовано с составлением уравнения
	// для поправки давления.
	//potent[iVar][iP]+=ds*(deltaP)/sl[iVar][iP].ap;//alpha[iVar]*
	// Так предлагает делать Гаврилов Андрей в расчётном комплексе Sigma-flow.
	// Это должно быть согласовано с составлением уравнения для поправки давления.
	// potent[iVar][iP]+=(tau/rho)*deltaP;
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
		// SIMPLE алгоритм: Carretto et al., 1973.
		// tau ~ (alpha[iVar]*rho*ds)/(sl[iVar][iP].ap);
		//potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
	}
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC алгоритм: Van Doormal and Raithby., 1984
		// tau ~ (alpha[iVar]*rho*ds)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
        //potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
	}

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скорректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume2

// коррекция скорости для внутренних КО.
// скорректированные скорости удовлетворяют уравнению 
// несжимаемости (неразрывности).
// begin 20 июня 2012 года.
void correct_internal_volume3(integer iP, integer iVar, float** prop,   
			  doublereal** &potent,  doublereal* tau) {

	// Процедура коррекции скорости во внутренних контрольных объёмах на 
	// основе предварительно вычисленного по формуле Грена-Гаусса градиента поправки давления.
    
    // Внутренний узел 
    // iP принадлежит интервалу 0..maxelm-1

	
	// Запоминание градиента от поправки давления
	// в зависимости от выбранной компоненты скорости или
	// координатного направления.
	doublereal gradP=0.0;
	switch (iVar) {
		case VELOCITY_X_COMPONENT: gradP=-potent[GRADXPAM][iP];
				  break;
		case VELOCITY_Y_COMPONENT: gradP=-potent[GRADYPAM][iP];
				  break;
        case VELOCITY_Z_COMPONENT: gradP=-potent[GRADZPAM][iP];
				  break;
	}


    // коррекция компоненты скорости во внутреннем узле iP. 

	// коррекция скорости не должна подвергаться нижней релаксации.
	// Так предлагает делать С. Патанкар. Это должно быть согласовано с составлением уравнения
	// для поправки давления.
	//potent[iVar][iP]+=ds*(deltaP)/sl[iVar][iP].ap;//alpha[iVar]*
	// Так предлагает делать Гаврилов Андрей в расчётном комплексе Sigma-flow.
	// Это должно быть согласовано с составлением уравнения для поправки давления.
	// potent[iVar][iP]+=(tau/rho)*deltaP;
	/*
	if (iSIMPLE_alg==SIMPLE_Carretto) {
		// SIMPLE алгоритм: Carretto et al., 1973.
		// tau ~ (alpha[iVar]*rho*ds)/(sl[iVar][iP].ap);
		//potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
	}
	if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC алгоритм: Van Doormal and Raithby., 1984
		// tau ~ (alpha[iVar]*rho*ds)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
        //potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
	}
	*/

	// Предполагается что шаг по псевдовремени tau предварительно сглажен,
	// применением нестационарной формулы расчёта где шаг по времени в стационарном случае
	// играет роль гладкого ориентира к которому можно стремиться. См. Гаврилов Андрей.
	// Если tau не сглаживать то возможна и весьма вероятна расходимость вычислительного
	// процесса или утрата сходимости.
	// В формуле для tau уже учитывается два варианта рабочего алгоритма: SIMPLE и SIMPLEC.
	potent[iVar][iP]+=tau[iP]*gradP/prop[RHO][iP];

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скорректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume3

// коррекция скорости для внутренних КО.
// скорректированные скорости удовлетворяют уравнению 
// несжимаемости (неразрывности).
// реализовано 23 июня 2012 года.
// На основе сглаженного псевдовремени, а точнее
// трёх скалярных полей псевдовремени
void correct_internal_volume4(integer iP, integer iVar, float** prop,   
			  doublereal** &potent,  doublereal** tau) {

	// Процедура коррекции скорости во внутренних контрольных объёмах на 
	// основе предварительно вычисленного по формуле Грена-Гаусса градиента поправки давления.
    
    // Внутренний узел 
    // iP принадлежит интервалу 0..maxelm-1

	
	// Запоминание градиента от поправки давления
	// в зависимости от выбранной компоненты скорости или
	// координатного направления.
	doublereal gradPAM=0.0, tauP=0.0;
	switch (iVar) {
		// всё правильно перед градиентом именно знак минус.
		case VELOCITY_X_COMPONENT: gradPAM=-potent[GRADXPAM][iP];
			      tauP=tau[VELOCITY_X_COMPONENT][iP];
				  break;
		case VELOCITY_Y_COMPONENT: gradPAM=-potent[GRADYPAM][iP];
			      tauP=tau[VELOCITY_Y_COMPONENT][iP];
				  break;
        case VELOCITY_Z_COMPONENT: gradPAM=-potent[GRADZPAM][iP];
			      tauP=tau[VELOCITY_Z_COMPONENT][iP];
				  break;
	}


    // коррекция компоненты скорости во внутреннем узле iP. 

	// коррекция скорости не должна подвергаться нижней релаксации.
	// Так предлагает делать С. Патанкар. Это должно быть согласовано с составлением уравнения
	// для поправки давления.
	//potent[iVar][iP]+=ds*(deltaP)/sl[iVar][iP].ap;//alpha[iVar]*
	// Так предлагает делать Гаврилов Андрей в расчётном комплексе Sigma-flow.
	// Это должно быть согласовано с составлением уравнения для поправки давления.
	// potent[iVar][iP]+=(tau/rho)*deltaP;
	/*
	if (iSIMPLE_alg==SIMPLE_Carretto) {
		// SIMPLE алгоритм: Carretto et al., 1973.
		// tau ~ (alpha[iVar]*rho*ds)/(sl[iVar][iP].ap);
		//potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/sl[iVar][iP].ap;// см. статьи Гаврилова Андрея.
	}
	if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC алгоритм: Van Doormal and Raithby., 1984
		// tau ~ (alpha[iVar]*rho*ds)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
        //potent[iVar][iP]+=alpha[iVar]*ds*(deltaP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);// см. статьи Гаврилова Андрея.
		potent[iVar][iP]+=alpha[iVar]*dv*(gradP)/((1.0-alpha[iVar])*sl[iVar][iP].ap);
	}
	*/

	// Предполагается что шаг по псевдовремени tau предварительно сглажен,
	// применением нестационарной формулы расчёта где шаг по времени в стационарном случае
	// играет роль гладкого ориентира к которому можно стремиться. См. Гаврилов Андрей.
	// Если tau не сглаживать то возможна и весьма вероятна расходимость вычислительного
	// процесса или утрата сходимости.
	// В формуле для tau уже учитывается два варианта рабочего алгоритма: SIMPLE и SIMPLEC.
	potent[iVar][iP]+=tauP*gradPAM/prop[RHO][iP];

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скорректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume4

// коррекция массового потока на грани КО.
// begin 25 июня 2012 года.
// 8.12.2018 адаптация кода к АЛИС.
void correct_mf(doublereal** &mfcurrentretune, doublereal** potent,  doublereal** tau,
	TOCHKA* pa, int*** neighbors_for_the_internal_node, int** nvtx, integer maxelm,
				BOUND* &border_neighbor, integer ls, integer lw, WALL* w, float** prop_b,
	integer *ilevel_alice, int* ptr) {

					

	doublereal** mfloc = NULL;
	mfloc=new doublereal*[maxelm];
	for (integer i=0; i<maxelm; i++) {
		mfloc[i]=new doublereal[6];
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
	    doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
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



	    doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	    // x-direction
	    feplus=0.5*dx/dxe; // if(bE) then feplus=1.0;
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE2]]) {
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

			if (bW) {
				// граничный узел.
				dSqw2 = border_neighbor[iW - maxelm].dS;
			}
			else {
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW]]) {
					dSqw2 = dy * dz;
				}
				else {
					// вычисление размеров соседнего контрольного объёма:
					doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
					volume3D(iW, nvtx, pa, dx_loc, dy_loc, dz_loc);

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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS2]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT2]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB2]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE3]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW3]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN3]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS3]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT3]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB3]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iE4]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iW4]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iN4]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iS4]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iT4]]) {
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
				if (ilevel_alice[ptr[iP]] >= ilevel_alice[ptr[iB4]]) {
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


		// значение псевдовремени на грани контрольного объёма.
		doublereal taue = 0.0, tauw = 0.0, taun = 0.0, taus = 0.0, taut = 0.0, taub = 0.0;
		doublereal taue2 = 0.0, tauw2 = 0.0, taun2 = 0.0, taus2 = 0.0, taut2 = 0.0, taub2 = 0.0;
		doublereal taue3 = 0.0, tauw3 = 0.0, taun3 = 0.0, taus3 = 0.0, taut3 = 0.0, taub3 = 0.0;
		doublereal taue4 = 0.0, tauw4 = 0.0, taun4 = 0.0, taus4 = 0.0, taut4 = 0.0, taub4 = 0.0;
        // интерполяция псевдовремени сделана так, чтобы выполнялись 
	    // предельные соотношения.
		if (iE > -1) {
			if (!bE) taue = tau[VELOCITY_X_COMPONENT][iE] * tau[VELOCITY_X_COMPONENT][iP] / (feplus*tau[VELOCITY_X_COMPONENT][iE] + (1.0 - feplus)*tau[VELOCITY_X_COMPONENT][iP]); else taue = tau[VELOCITY_X_COMPONENT][iE];
		}
		if (iW > -1) {
			if (!bW) tauw = tau[VELOCITY_X_COMPONENT][iW] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus*tau[VELOCITY_X_COMPONENT][iW] + (1.0 - fwplus)*tau[VELOCITY_X_COMPONENT][iP]); else tauw = tau[VELOCITY_X_COMPONENT][iW];
		}
		if (iN > -1) {
			if (!bN) taun = tau[VELOCITY_Y_COMPONENT][iN] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus*tau[VELOCITY_Y_COMPONENT][iN] + (1.0 - fnplus)*tau[VELOCITY_Y_COMPONENT][iP]); else taun = tau[VELOCITY_Y_COMPONENT][iN];
		}
		if (iS > -1) {
			if (!bS) taus = tau[VELOCITY_Y_COMPONENT][iS] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus*tau[VELOCITY_Y_COMPONENT][iS] + (1.0 - fsplus)*tau[VELOCITY_Y_COMPONENT][iP]); else taus = tau[VELOCITY_Y_COMPONENT][iS];
		}
		if (iT > -1) {
			if (!bT) taut = tau[VELOCITY_Z_COMPONENT][iT] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus*tau[VELOCITY_Z_COMPONENT][iT] + (1.0 - ftplus)*tau[VELOCITY_Z_COMPONENT][iP]); else taut = tau[VELOCITY_Z_COMPONENT][iT];
		}
		if (iB > -1) {
			if (!bB) taub = tau[VELOCITY_Z_COMPONENT][iB] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus*tau[VELOCITY_Z_COMPONENT][iB] + (1.0 - fbplus)*tau[VELOCITY_Z_COMPONENT][iP]); else taub = tau[VELOCITY_Z_COMPONENT][iB];
		}

		if (iE2 > -1) {
			if (!bE2) taue2 = tau[VELOCITY_X_COMPONENT][iE2] * tau[VELOCITY_X_COMPONENT][iP] / (feplus2*tau[VELOCITY_X_COMPONENT][iE2] + (1.0 - feplus2)*tau[VELOCITY_X_COMPONENT][iP]); else taue2 = tau[VELOCITY_X_COMPONENT][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) tauw2 = tau[VELOCITY_X_COMPONENT][iW2] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus2*tau[VELOCITY_X_COMPONENT][iW2] + (1.0 - fwplus2)*tau[VELOCITY_X_COMPONENT][iP]); else tauw2 = tau[VELOCITY_X_COMPONENT][iW2];
		}
		if (iN2 > -1) {
			if (!bN2) taun2 = tau[VELOCITY_Y_COMPONENT][iN2] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus2*tau[VELOCITY_Y_COMPONENT][iN2] + (1.0 - fnplus2)*tau[VELOCITY_Y_COMPONENT][iP]); else taun2 = tau[VELOCITY_Y_COMPONENT][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) taus2 = tau[VELOCITY_Y_COMPONENT][iS2] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus2*tau[VELOCITY_Y_COMPONENT][iS2] + (1.0 - fsplus2)*tau[VELOCITY_Y_COMPONENT][iP]); else taus2 = tau[VELOCITY_Y_COMPONENT][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) taut2 = tau[VELOCITY_Z_COMPONENT][iT2] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus2*tau[VELOCITY_Z_COMPONENT][iT2] + (1.0 - ftplus2)*tau[VELOCITY_Z_COMPONENT][iP]); else taut2 = tau[VELOCITY_Z_COMPONENT][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) taub2 = tau[VELOCITY_Z_COMPONENT][iB2] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus2*tau[VELOCITY_Z_COMPONENT][iB2] + (1.0 - fbplus2)*tau[VELOCITY_Z_COMPONENT][iP]); else taub2 = tau[VELOCITY_Z_COMPONENT][iB2];
		}

		if (iE3 > -1) {
			if (!bE3) taue3 = tau[VELOCITY_X_COMPONENT][iE3] * tau[VELOCITY_X_COMPONENT][iP] / (feplus3*tau[VELOCITY_X_COMPONENT][iE3] + (1.0 - feplus3)*tau[VELOCITY_X_COMPONENT][iP]); else taue3 = tau[VELOCITY_X_COMPONENT][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) tauw3 = tau[VELOCITY_X_COMPONENT][iW3] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus3*tau[VELOCITY_X_COMPONENT][iW3] + (1.0 - fwplus3)*tau[VELOCITY_X_COMPONENT][iP]); else tauw3 = tau[VELOCITY_X_COMPONENT][iW3];
		}
		if (iN3 > -1) {
			if (!bN3) taun3 = tau[VELOCITY_Y_COMPONENT][iN3] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus3*tau[VELOCITY_Y_COMPONENT][iN3] + (1.0 - fnplus3)*tau[VELOCITY_Y_COMPONENT][iP]); else taun3 = tau[VELOCITY_Y_COMPONENT][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) taus3 = tau[VELOCITY_Y_COMPONENT][iS3] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus3*tau[VELOCITY_Y_COMPONENT][iS3] + (1.0 - fsplus3)*tau[VELOCITY_Y_COMPONENT][iP]); else taus3 = tau[VELOCITY_Y_COMPONENT][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) taut3 = tau[VELOCITY_Z_COMPONENT][iT3] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus3*tau[VELOCITY_Z_COMPONENT][iT3] + (1.0 - ftplus3)*tau[VELOCITY_Z_COMPONENT][iP]); else taut3 = tau[VELOCITY_Z_COMPONENT][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) taub3 = tau[VELOCITY_Z_COMPONENT][iB3] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus3*tau[VELOCITY_Z_COMPONENT][iB3] + (1.0 - fbplus3)*tau[VELOCITY_Z_COMPONENT][iP]); else taub3 = tau[VELOCITY_Z_COMPONENT][iB3];
		}

		if (iE4 > -1) {
			if (!bE4) taue4 = tau[VELOCITY_X_COMPONENT][iE4] * tau[VELOCITY_X_COMPONENT][iP] / (feplus4*tau[VELOCITY_X_COMPONENT][iE4] + (1.0 - feplus4)*tau[VELOCITY_X_COMPONENT][iP]); else taue4 = tau[VELOCITY_X_COMPONENT][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) tauw4 = tau[VELOCITY_X_COMPONENT][iW4] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus4*tau[VELOCITY_X_COMPONENT][iW4] + (1.0 - fwplus4)*tau[VELOCITY_X_COMPONENT][iP]); else tauw4 = tau[VELOCITY_X_COMPONENT][iW4];
		}
		if (iN4 > -1) {
			if (!bN4) taun4 = tau[VELOCITY_Y_COMPONENT][iN4] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus4*tau[VELOCITY_Y_COMPONENT][iN4] + (1.0 - fnplus4)*tau[VELOCITY_Y_COMPONENT][iP]); else taun4 = tau[VELOCITY_Y_COMPONENT][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) taus4 = tau[VELOCITY_Y_COMPONENT][iS4] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus4*tau[VELOCITY_Y_COMPONENT][iS4] + (1.0 - fsplus4)*tau[VELOCITY_Y_COMPONENT][iP]); else taus4 = tau[VELOCITY_Y_COMPONENT][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) taut4 = tau[VELOCITY_Z_COMPONENT][iT4] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus4*tau[VELOCITY_Z_COMPONENT][iT4] + (1.0 - ftplus4)*tau[VELOCITY_Z_COMPONENT][iP]); else taut4 = tau[VELOCITY_Z_COMPONENT][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) taub4 = tau[VELOCITY_Z_COMPONENT][iB4] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus4*tau[VELOCITY_Z_COMPONENT][iB4] + (1.0 - fbplus4)*tau[VELOCITY_Z_COMPONENT][iP]); else taub4 = tau[VELOCITY_Z_COMPONENT][iB4];
		}

		// Градиент поправки давления на грани контрольного объёма.
		doublereal gradpame = 0.0, gradpamw = 0.0, gradpamn = 0.0, gradpams = 0.0, gradpamt = 0.0, gradpamb = 0.0;
		doublereal gradpame2 = 0.0, gradpamw2 = 0.0, gradpamn2 = 0.0, gradpams2 = 0.0, gradpamt2 = 0.0, gradpamb2 = 0.0;
		doublereal gradpame3 = 0.0, gradpamw3 = 0.0, gradpamn3 = 0.0, gradpams3 = 0.0, gradpamt3 = 0.0, gradpamb3 = 0.0;
		doublereal gradpame4 = 0.0, gradpamw4 = 0.0, gradpamn4 = 0.0, gradpams4 = 0.0, gradpamt4 = 0.0, gradpamb4 = 0.0;

		if (iE > -1) {
			if (!bE) gradpame = feplus * potent[GRADXPAM][iE] + (1.0 - feplus)*potent[GRADXPAM][iP]; else gradpame = potent[GRADXPAM][iE];
		}
		if (iW > -1) {
			if (!bW) gradpamw = fwplus * potent[GRADXPAM][iW] + (1.0 - fwplus)*potent[GRADXPAM][iP]; else gradpamw = potent[GRADXPAM][iW];
		}
		if (iN > -1) {
			if (!bN) gradpamn = fnplus * potent[GRADYPAM][iN] + (1.0 - fnplus)*potent[GRADYPAM][iP]; else gradpamn = potent[GRADYPAM][iN];
		}
		if (iS > -1) {
			if (!bS) gradpams = fsplus * potent[GRADYPAM][iS] + (1.0 - fsplus)*potent[GRADYPAM][iP]; else gradpams = potent[GRADYPAM][iS];
		}
		if (iT > -1) {
			if (!bT) gradpamt = ftplus * potent[GRADZPAM][iT] + (1.0 - ftplus)*potent[GRADZPAM][iP]; else gradpamt = potent[GRADZPAM][iT];
		}
		if (iB > -1) {
			if (!bB) gradpamb = fbplus * potent[GRADZPAM][iB] + (1.0 - fbplus)*potent[GRADZPAM][iP]; else gradpamb = potent[GRADZPAM][iB];
		}

		if (iE2 > -1) {
			if (!bE2) gradpame2 = feplus2 * potent[GRADXPAM][iE2] + (1.0 - feplus2)*potent[GRADXPAM][iP]; else gradpame2 = potent[GRADXPAM][iE2];
		}
		if (iW2 > -1) {
			if (!bW2) gradpamw2 = fwplus2 * potent[GRADXPAM][iW2] + (1.0 - fwplus2)*potent[GRADXPAM][iP]; else gradpamw2 = potent[GRADXPAM][iW2];
		}
		if (iN2 > -1) {
			if (!bN2) gradpamn2 = fnplus2 * potent[GRADYPAM][iN2] + (1.0 - fnplus2)*potent[GRADYPAM][iP]; else gradpamn2 = potent[GRADYPAM][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) gradpams2 = fsplus2 * potent[GRADYPAM][iS2] + (1.0 - fsplus2)*potent[GRADYPAM][iP]; else gradpams2 = potent[GRADYPAM][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) gradpamt2 = ftplus2 * potent[GRADZPAM][iT2] + (1.0 - ftplus2)*potent[GRADZPAM][iP]; else gradpamt2 = potent[GRADZPAM][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) gradpamb2 = fbplus2 * potent[GRADZPAM][iB2] + (1.0 - fbplus2)*potent[GRADZPAM][iP]; else gradpamb2 = potent[GRADZPAM][iB2];
		}

		if (iE3 > -1) {
			if (!bE3) gradpame3 = feplus3 * potent[GRADXPAM][iE3] + (1.0 - feplus3)*potent[GRADXPAM][iP]; else gradpame3 = potent[GRADXPAM][iE3];
		}
		if (iW3 > -1) {
			if (!bW3) gradpamw3 = fwplus3 * potent[GRADXPAM][iW3] + (1.0 - fwplus3)*potent[GRADXPAM][iP]; else gradpamw3 = potent[GRADXPAM][iW3];
		}
		if (iN3 > -1) {
			if (!bN3) gradpamn3 = fnplus3 * potent[GRADYPAM][iN3] + (1.0 - fnplus3)*potent[GRADYPAM][iP]; else gradpamn3 = potent[GRADYPAM][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) gradpams3 = fsplus3 * potent[GRADYPAM][iS3] + (1.0 - fsplus3)*potent[GRADYPAM][iP]; else gradpams3 = potent[GRADYPAM][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) gradpamt3 = ftplus3 * potent[GRADZPAM][iT3] + (1.0 - ftplus3)*potent[GRADZPAM][iP]; else gradpamt3 = potent[GRADZPAM][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) gradpamb3 = fbplus3 * potent[GRADZPAM][iB3] + (1.0 - fbplus3)*potent[GRADZPAM][iP]; else gradpamb3 = potent[GRADZPAM][iB3];
		}

		if (iE4 > -1) {
			if (!bE4) gradpame4 = feplus4 * potent[GRADXPAM][iE4] + (1.0 - feplus4)*potent[GRADXPAM][iP]; else gradpame4 = potent[GRADXPAM][iE4];
		}
		if (iW4 > -1) {
			if (!bW4) gradpamw4 = fwplus4 * potent[GRADXPAM][iW4] + (1.0 - fwplus4)*potent[GRADXPAM][iP]; else gradpamw4 = potent[GRADXPAM][iW4];
		}
		if (iN4 > -1) {
			if (!bN4) gradpamn4 = fnplus4 * potent[GRADYPAM][iN4] + (1.0 - fnplus4)*potent[GRADYPAM][iP]; else gradpamn4 = potent[GRADYPAM][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) gradpams4 = fsplus4 * potent[GRADYPAM][iS4] + (1.0 - fsplus4)*potent[GRADYPAM][iP]; else gradpams4 = potent[GRADYPAM][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) gradpamt4 = ftplus4 * potent[GRADZPAM][iT4] + (1.0 - ftplus4)*potent[GRADZPAM][iP]; else gradpamt4 = potent[GRADZPAM][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) gradpamb4 = fbplus4 * potent[GRADZPAM][iB4] + (1.0 - fbplus4)*potent[GRADZPAM][iP]; else gradpamb4 = potent[GRADZPAM][iB4];
		}


		// Наконец вычисление скорректированного массового 
		// потока на грани КО.
		mfloc[iP][E_SIDE]=mfcurrentretune[iP][E_SIDE]-taue*gradpame*dSqe - taue2 * gradpame2*dSqe2 - taue3 * gradpame3*dSqe3 - taue4 * gradpame4*dSqe4;
		mfloc[iP][W_SIDE]=mfcurrentretune[iP][W_SIDE]-tauw*gradpamw*dSqw - tauw2 * gradpamw2*dSqw2 - tauw3 * gradpamw3*dSqw3 - tauw4 * gradpamw4*dSqw4;
		mfloc[iP][N_SIDE]=mfcurrentretune[iP][N_SIDE]-taun*gradpamn*dSqn - taun2 * gradpamn2*dSqn2 - taun3 * gradpamn3*dSqn3 - taun4 * gradpamn4*dSqn4;
		mfloc[iP][S_SIDE]=mfcurrentretune[iP][S_SIDE]-taus*gradpams*dSqs - taus2 * gradpams2*dSqs2 - taus3 * gradpams3*dSqs3 - taus4 * gradpams4*dSqs4;
		mfloc[iP][T_SIDE]=mfcurrentretune[iP][T_SIDE]-taut*gradpamt*dSqt - taut2 * gradpamt2*dSqt2 - taut3 * gradpamt3*dSqt3 - taut4 * gradpamt4*dSqt4;
		mfloc[iP][B_SIDE]=mfcurrentretune[iP][B_SIDE]-taub*gradpamb*dSqb - taub2 * gradpamb2*dSqb2 - taub3 * gradpamb3*dSqb3 - taub4 * gradpamb4*dSqb4;

		if (mfloc[iP][E_SIDE]!= mfloc[iP][E_SIDE]) {
			printf("mfcurrentretune[%lld][ESIDE]=%e  taue=%e, gradpame=%e, dSqe=%e\n",iP, mfcurrentretune[iP][E_SIDE], taue, gradpame, dSqe);
			printf("taue2=%e, gradpame2=%e, dSqe2=%e\n", taue2, gradpame2, dSqe2);
			printf("taue3=%e, gradpame3=%e, dSqe3=%e\n", taue3, gradpame3, dSqe3);
			printf("taue4=%e, gradpame4=%e, dSqe4=%e\n", taue4, gradpame4, dSqe4);
			system("pause");
		}

		if (mfloc[iP][W_SIDE] != mfloc[iP][W_SIDE]) {
			printf("mfcurrentretune[%lld][WSIDE]=%e  tauw=%e, gradpamw=%e, dSqw=%e\n", iP, mfcurrentretune[iP][W_SIDE], tauw, gradpamw, dSqw);
			printf("tauw2=%e, gradpamw2=%e, dSqw2=%e\n", tauw2, gradpamw2, dSqw2);
			printf("tauw3=%e, gradpamw3=%e, dSqw3=%e\n", tauw3, gradpamw3, dSqw3);
			printf("tauw4=%e, gradpamw4=%e, dSqw4=%e\n", tauw4, gradpamw4, dSqw4);
			system("pause");
		}

		if (mfloc[iP][N_SIDE] != mfloc[iP][N_SIDE]) {
			printf("mfcurrentretune[%lld][NSIDE]=%e  taun=%e, gradpamn=%e, dSqn=%e\n", iP, mfcurrentretune[iP][N_SIDE], taun, gradpamn, dSqn);
			printf("taun2=%e, gradpamn2=%e, dSqn2=%e\n", taun2, gradpamn2, dSqn2);
			printf("taun3=%e, gradpamn3=%e, dSqn3=%e\n", taun3, gradpamn3, dSqn3);
			printf("taun4=%e, gradpamn4=%e, dSqn4=%e\n", taun4, gradpamn4, dSqn4);
			system("pause");
		}

		if (mfloc[iP][S_SIDE] != mfloc[iP][S_SIDE]) {
			printf("mfcurrentretune[%lld][SSIDE]=%e  taus=%e, gradpams=%e, dSqs=%e\n", iP, mfcurrentretune[iP][S_SIDE], taus, gradpams, dSqs);
			printf("taus2=%e, gradpams2=%e, dSqs2=%e\n", taus2, gradpams2, dSqs2);
			printf("taus3=%e, gradpams3=%e, dSqs3=%e\n", taus3, gradpams3, dSqs3);
			printf("taus4=%e, gradpams4=%e, dSqs4=%e\n", taus4, gradpams4, dSqs4);
			system("pause");
		}

		if (mfloc[iP][T_SIDE] != mfloc[iP][T_SIDE]) {
			printf("mfcurrentretune[%lld][TSIDE]=%e  taut=%e, gradpamt=%e, dSqt=%e\n", iP, mfcurrentretune[iP][T_SIDE], taut, gradpamt, dSqt);
			printf("taut2=%e, gradpamt2=%e, dSqt2=%e\n", taut2, gradpamt2, dSqt2);
			printf("taut3=%e, gradpamt3=%e, dSqt3=%e\n", taut3, gradpamt3, dSqt3);
			printf("taut4=%e, gradpamt4=%e, dSqt4=%e\n", taut4, gradpamt4, dSqt4);
			system("pause");
		}

		if (mfloc[iP][B_SIDE] != mfloc[iP][B_SIDE]) {
			printf("mfcurrentretune[%lld][BSIDE]=%e  taub=%e, gradpamb=%e, dSqb=%e\n", iP, mfcurrentretune[iP][B_SIDE], taub, gradpamb, dSqb);
			printf("taub2=%e, gradpamb2=%e, dSqb2=%e\n", taub2, gradpamb2, dSqb2);
			printf("taub3=%e, gradpamb3=%e, dSqb3=%e\n", taub3, gradpamb3, dSqb3);
			printf("taub4=%e, gradpamb4=%e, dSqb4=%e\n", taub4, gradpamb4, dSqb4);
			system("pause");
		}

}

	// Однако есть границы где массовый поток задан пользователем,
	// очевидно это надо как-то учесть.
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
		const doublereal relax_bound = 1.0;

		if (bE||bE2||bE3||bE4) {
			integer inumber = -1;
			if (bE) {
				inumber = iE - maxelm;
			}
			else if (bE2) {
				inumber = iE2 - maxelm;
			}
			else if (bE3) {
				inumber = iE3 - maxelm;
			}
			else if (bE4) {
				inumber = iE4 - maxelm;
			}
			if (inumber > -1) {
				if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
					// Выходная граница оставляем всё как есть
					mfloc[iP][E_SIDE] = relax_bound * (mfloc[iP][E_SIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][E_SIDE];
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry) {
					mfloc[iP][E_SIDE] = 0.0;
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					// заданная скорость на входной границе.
					mfloc[iP][E_SIDE] = prop_b[RHO][inumber] * w[border_neighbor[inumber].MCB - ls].Vx*dy*dz; // заданный массовый поток.
				}
				else {
					// твёрдая неподвижная стенка по умолчанию
					mfloc[iP][E_SIDE] = 0.0;
				}
			}
		}

		if (bW||bW2||bW3||bW4) {
			integer inumber = -1;
			if (bW) {
				inumber = iW - maxelm;
			}
			else if (bW2) {
				inumber = iW2 - maxelm;
			}
			else if (bW3) {
				inumber = iW3 - maxelm;
			}
			else if (bW4) {
				inumber = iW4 - maxelm;
			}
			if (inumber > -1) {

				if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
					// Выходная граница оставляем всё как есть
					mfloc[iP][W_SIDE] = relax_bound * (mfloc[iP][W_SIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][W_SIDE];
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry) {
					mfloc[iP][W_SIDE] = 0.0;
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					// заданная скорость на входной границе.
					mfloc[iP][W_SIDE] = prop_b[RHO][inumber] * w[border_neighbor[inumber].MCB - ls].Vx*dy*dz; // заданный массовый поток.
				}
				else {
					// твёрдая неподвижная стенка по умолчанию
					mfloc[iP][W_SIDE] = 0.0;
				}
			}
		}

		if (bN||bN2||bN3||bN4) {
			integer inumber=-1;
			if (bN) {
				inumber = iN - maxelm;
			}
			else if (bN2) {
				inumber = iN2 - maxelm;
			}
			else if (bN3) {
				inumber = iN3 - maxelm;
			}
			else if (bN4) {
				inumber = iN4 - maxelm;
			}
			if (inumber > -1) {

				if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
					// Выходная граница оставляем всё как есть
					mfloc[iP][N_SIDE] = relax_bound * (mfloc[iP][N_SIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][N_SIDE];
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry) {
					mfloc[iP][N_SIDE] = 0.0;
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					// заданная скорость на входной границе.
					mfloc[iP][N_SIDE] = prop_b[RHO][inumber] * w[border_neighbor[inumber].MCB - ls].Vy*dx*dz; // заданный массовый поток.
				}
				else {
					// твёрдая неподвижная стенка по умолчанию
					mfloc[iP][N_SIDE] = 0.0;
				}
			}
		}

		if (bS||bS2||bS3||bS4) {
			integer inumber = -1;
			if (bS) {
				inumber = iS - maxelm;
			}
			else if (bS2) {
				inumber = iS2 - maxelm;
			}
			else if (bS3) {
				inumber = iS3 - maxelm;
			}
			else if (bS4) {
				inumber = iS4 - maxelm;
			}
			if (inumber > -1) {

				if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
					// Выходная граница оставляем всё как есть
					mfloc[iP][S_SIDE] = relax_bound * (mfloc[iP][S_SIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][S_SIDE];
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry) {
					mfloc[iP][S_SIDE] = 0.0;
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					// заданная скорость на входной границе.
					mfloc[iP][S_SIDE] = prop_b[RHO][inumber] * w[border_neighbor[inumber].MCB - ls].Vy*dx*dz; // заданный массовый поток.
				}
				else {
					// твёрдая неподвижная стенка по умолчанию
					mfloc[iP][S_SIDE] = 0.0;
				}
			}
		}

		if (bT||bT2||bT3||bT4) {
			integer inumber = -1;
			if (bT) {
				inumber = iT - maxelm;
			}
			else if (bT2) {
				inumber = iT2 - maxelm;
			}
			else if (bT3) {
				inumber = iT3 - maxelm;
			}
			else if (bT4) {
				inumber = iT4 - maxelm;
			}
			if (inumber > -1) {

				if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
					// Выходная граница оставляем всё как есть
					mfloc[iP][T_SIDE] = relax_bound * (mfloc[iP][T_SIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][T_SIDE];
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry) {
					mfloc[iP][T_SIDE] = 0.0;
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					// заданная скорость на входной границе.
					mfloc[iP][T_SIDE] = prop_b[RHO][inumber] * w[border_neighbor[inumber].MCB - ls].Vz*dx*dy; // заданный массовый поток.
				}
				else {
					// твёрдая неподвижная стенка по умолчанию
					mfloc[iP][T_SIDE] = 0.0;
				}
			}
		}

		if (bB||bB2||bB3||bB4) {

			integer inumber = -1;
			if (bB) {
				inumber = iB - maxelm;
			}
			else if (bB2) {
				inumber = iB2 - maxelm;
			}
			else if (bB3) {
				inumber = iB3 - maxelm;
			}
			else if (bB4) {
				inumber = iB4 - maxelm;
			}
			if (inumber > -1) {

				if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
					// Выходная граница оставляем всё как есть
					mfloc[iP][B_SIDE] = relax_bound * (mfloc[iP][B_SIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][B_SIDE];
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry) {
					mfloc[iP][B_SIDE] = 0.0;
				}
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					// заданная скорость на входной границе.
					mfloc[iP][B_SIDE] = prop_b[RHO][inumber] * w[border_neighbor[inumber].MCB - ls].Vz*dx*dy; // заданный массовый поток.
				}
				else {
					// твёрдая неподвижная стенка по умолчанию
					mfloc[iP][B_SIDE] = 0.0;
				}
			}
		}

	}

	// Обратное копирование.
	for (integer iG=0; iG<6; iG++) {
#pragma omp parallel for
		for (integer iP=0; iP<maxelm; iP++) {
			mfcurrentretune[iP][iG]=mfloc[iP][iG];
		}
	}

	// Освобождение памяти.
	
	if (mfloc != NULL) {
		for (integer i = 0; i < maxelm; i++) {
			if (mfloc[i] != NULL) {
				delete[] mfloc[i];
			}
		}
		delete[] mfloc;
	}

} // correct_mf

// 02.04.2019 Теперь данная функция не просто 
// сообщает о том что на твердой стенке ошибочно  ненулевой поток массы вещества,
// но и автоматом корректирует его в ноль. Это проявилось при расчёте теплопередачи на АЛИС
// сетки со скоростью, считанной из load.txt файла.
void iscorrectmf(doublereal** &mf,
	integer maxelm,
	int*** neighbors_for_the_internal_node, BOUND* &border_neighbor,
	integer ls, integer lw, WALL* w) {

	integer iP = 0;
	integer inumber;
	bool biscorrectmf = false;
	// 02.04.2019
	bool bdiagnostic_message = (!(bonly_solid_calculation&&b_on_adaptive_local_refinement_mesh));
	// iP - номер центрального контрольного объёма
	for (iP = 0; iP < maxelm; iP++) {
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
				inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][E_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iE=%lld\n", iE);
#else
							printf("wall mf flux velocity non zero iE=%d\n", iE);
#endif
						}


						mf[iP][E_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}


				}
			}
		}
		if (iW > -1) {
			if (iW >= maxelm) {
				// граничный узел
				inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][W_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iW=%lld\n", iW);
#else
							printf("wall mf flux velocity non zero iW=%d\n", iW);
#endif
						}
						mf[iP][W_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iN > -1) {
			if (iN >= maxelm) {
				// граничный узел
				inumber = iN - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][N_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iN=%lld\n", iN);
#else
							printf("wall mf flux velocity non zero iN=%d\n", iN);
#endif
						}

						mf[iP][N_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
        			}
				}
			}
		}

		if (iS >-1) {
		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				if (fabs(mf[iP][S_SIDE]) > admission) {
					if (bdiagnostic_message) {
#if doubleintprecision == 1
						printf("wall mf flux velocity non zero iS=%lld\n", iS);
#else
						printf("wall mf flux velocity non zero iS=%d\n", iS);
#endif
			}
					
					mf[iP][S_SIDE] = 0.0;
					biscorrectmf = true;
					//system("pause");
				}
				
			}
		}
		}

		if (iT > -1) {
			if (iT >= maxelm) {
				// граничный узел
				inumber = iT - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][T_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iT=%lld\n", iT);
#else
							printf("wall mf flux velocity non zero iT=%d\n", iT);
#endif
						}
						mf[iP][T_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iB > -1) {
			if (iB >= maxelm) {
				// граничный узел
				inumber = iB - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][B_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iB=%lld\n", iB);
#else
							printf("wall mf flux velocity non zero iB=%d\n", iB);
#endif
						}

						mf[iP][B_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iE2 > -1) {
			if (iE2 >= maxelm) {
				// граничный узел
				inumber = iE2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][E_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iE2=%lld\n", iE2);
#else
							printf("wall mf flux velocity non zero iE2=%d\n", iE2);
#endif
						}

						mf[iP][E_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iW2 > -1) {
			if (iW2 >= maxelm) {
				// граничный узел
				inumber = iW2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][W_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iW2=%lld\n", iW2);
#else
							printf("wall mf flux velocity non zero iW2=%d\n", iW2);
#endif
						}

						mf[iP][W_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iN2 > -1) {
			if (iN2 >= maxelm) {
				// граничный узел
				inumber = iN2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][N_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iN2=%lld\n", iN2);
#else
							printf("wall mf flux velocity non zero iN2=%d\n", iN2);
#endif
						}

						mf[iP][N_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iS2 >-1) {
			if (iS2 >= maxelm) {
				// граничный узел
				inumber = iS2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][S_SIDE])>admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iS2=%lld\n", iS2);
#else
							printf("wall mf flux velocity non zero iS2=%d\n", iS2);
#endif
						}

						mf[iP][S_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iT2 > -1) {
			if (iT2 >= maxelm) {
				// граничный узел
				inumber = iT2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][T_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iT2=%lld\n", iT2);
#else
							printf("wall mf flux velocity non zero iT2=%d\n", iT2);
#endif
						}

						mf[iP][T_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iB2 > -1) {
			if (iB2 >= maxelm) {
				// граничный узел
				inumber = iB2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][B_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iB2=%lld\n", iB2);
#else
							printf("wall mf flux velocity non zero iB2=%d\n", iB2);
#endif
						}

						mf[iP][B_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iE3 > -1) {
			if (iE3 >= maxelm) {
				// граничный узел
				inumber = iE3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][E_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iE3=%lld\n", iE3);
#else
							printf("wall mf flux velocity non zero iE3=%d\n", iE3);
#endif
						}


						mf[iP][E_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iW3 > -1) {
			if (iW3 >= maxelm) {
				// граничный узел
				inumber = iW3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][W_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iW=%lld\n", iW3);
#else
							printf("wall mf flux velocity non zero iW=%d\n", iW3);
#endif
				}

						mf[iP][W_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iN3 > -1) {
			if (iN3 >= maxelm) {
				// граничный узел
				inumber = iN3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][N_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iN3=%lld\n", iN3);
#else
							printf("wall mf flux velocity non zero iN3=%d\n", iN3);
#endif
						}

						mf[iP][N_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iS3 >-1) {
			if (iS3 >= maxelm) {
				// граничный узел
				inumber = iS3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][S_SIDE])>admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iS=%lld\n", iS3);
#else
							printf("wall mf flux velocity non zero iS=%d\n", iS3);
#endif
						}

						mf[iP][S_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iT3 > -1) {
			if (iT3 >= maxelm) {
				// граничный узел
				inumber = iT3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][T_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iT3=%lld\n", iT3);
#else
							printf("wall mf flux velocity non zero iT3=%d\n", iT3);
#endif
						}

						mf[iP][T_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iB3 > -1) {
			if (iB3 >= maxelm) {
				// граничный узел
				inumber = iB3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][B_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iB3=%lld\n", iB3);
#else
							printf("wall mf flux velocity non zero iB3=%d\n", iB3);
#endif
						}

						mf[iP][B_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iE4 > -1) {
			if (iE4 >= maxelm) {
				// граничный узел
				inumber = iE4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][E_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iE4=%lld\n", iE4);
#else
							printf("wall mf flux velocity non zero iE4=%d\n", iE4);
#endif
				}


						mf[iP][E_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iW4 > -1) {
			if (iW4 >= maxelm) {
				// граничный узел
				inumber = iW4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][W_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iW4=%lld\n", iW4);
#else
							printf("wall mf flux velocity non zero iW4=%d\n", iW4);
#endif
						}

						mf[iP][W_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iN4 > -1) {
			if (iN4 >= maxelm) {
				// граничный узел
				inumber = iN4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][N_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iN4=%lld\n", iN4);
#else
							printf("wall mf flux velocity non zero iN4=%d\n", iN4);
#endif
						}

						mf[iP][N_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iS4 >-1) {
			if (iS4 >= maxelm) {
				// граничный узел
				inumber = iS4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][S_SIDE])>admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iS4=%lld\n", iS4);
#else
							printf("wall mf flux velocity non zero iS4=%d\n", iS4);
#endif
						}

						mf[iP][S_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}

		if (iT4 > -1) {
			if (iT4 >= maxelm) {
				// граничный узел
				inumber = iT4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][T_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iT4=%lld\n", iT4);
#else
							printf("wall mf flux velocity non zero iT4=%d\n", iT4);
#endif
						}

						mf[iP][T_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}
		if (iB4 > -1) {
			if (iB4 >= maxelm) {
				// граничный узел
				inumber = iB4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(mf[iP][B_SIDE]) > admission) {
						if (bdiagnostic_message) {
#if doubleintprecision == 1
							printf("wall mf flux velocity non zero iB4=%lld\n", iB4);
#else
							printf("wall mf flux velocity non zero iB4=%d\n", iB4);
#endif
						}

						mf[iP][B_SIDE] = 0.0;
						biscorrectmf = true;
						//system("pause");
					}

				}
			}
		}


	}

	if (bdiagnostic_message) {
		if (biscorrectmf  ) {
			system("PAUSE");
		}
	}

}

void iscorrectOk(doublereal** &potent,
	integer maxelm,
	int*** neighbors_for_the_internal_node, BOUND* &border_neighbor,
	integer ls, integer lw, WALL* w)
{
	
	
	// iP - номер центрального контрольного объёма
#pragma omp parallel for
	for (integer iP = 0; iP < maxelm; iP++) {

		integer inumber;

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
				inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iE]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iE=%lld", iE);
#else
						printf("wall VX velocity non zero iE=%d", iE);
#endif


						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iE]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iE=%lld", iE);
#else
						printf("wall VY velocity non zero iE=%d", iE);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iE]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iE=%lld", iE);
#else
						printf("wall VZ velocity non zero iE=%d", iE);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW > -1) {
			if (iW >= maxelm) {
				// граничный узел
				inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iW]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iW=%lld", iW);
#else
						printf("wall VX velocity non zero iW=%d", iW);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iW]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iW=%lld", iW);
#else
						printf("wall VY velocity non zero iW=%d", iW);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iW]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iW=%lld", iW);
#else
						printf("wall VZ velocity non zero iW=%d", iW);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN > -1) {
			if (iN >= maxelm) {
				// граничный узел
				inumber = iN - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iN]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iN=%lld", iN);
#else
						printf("wall VX velocity non zero iN=%d", iN);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iN]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iN=%lld", iN);
#else
						printf("wall VY velocity non zero iN=%d", iN);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iN]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iN=%lld", iN);
#else
						printf("wall VZ velocity non zero iN=%d", iN);
#endif

						//getchar();
						system("pause");
			}
		}
	}
}

		if (iS>-1) {
		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (border_neighbor[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VELOCITY_X_COMPONENT][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iS=%lld", iS);
#else
					printf("wall VX velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VELOCITY_Y_COMPONENT][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iS=%lld", iS);
#else
					printf("wall VY velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VELOCITY_Z_COMPONENT][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VZ velocity non zero iS=%lld", iS);
#else
					printf("wall VZ velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
			}
		}
		}


		if (iT > -1) {
			if (iT >= maxelm) {
				// граничный узел
				inumber = iT - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iT]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iT=%lld", iT);
#else
						printf("wall VX velocity non zero iT=%d", iT);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iT]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iT=%lld", iT);
#else
						printf("wall VY velocity non zero iT=%d", iT);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iT]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iT=%lld", iT);
#else
						printf("wall VZ velocity non zero iT=%d", iT);
#endif

						//getchar();
						system("pause");
					}
				}
					}
				}

		if (iB > -1) {
			if (iB >= maxelm) {
				// граничный узел
				inumber = iB - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iB]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iB=%lld", iB);
#else
						printf("wall VX velocity non zero iB=%d", iB);
#endif

						//getchar();
						system("pause");
				}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iB]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iB=%lld", iB);
#else
						printf("wall VY velocity non zero iB=%d", iB);
#endif

						//getchar();
						system("pause");
			}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iB]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iB=%lld", iB);
#else
						printf("wall VZ velocity non zero iB=%d", iB);
#endif

						//getchar();
						system("pause");
					}
		}
			}
		}

		if (iE2 > -1) {
			if (iE2 >= maxelm) {
				// граничный узел
				inumber = iE2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iE2]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iE2=%lld", iE2);
#else
						printf("wall VX velocity non zero iE2=%d", iE2);
#endif


						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iE2]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iE2=%lld", iE2);
#else
						printf("wall VY velocity non zero iE2=%d", iE2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iE2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iE2=%lld", iE2);
#else
						printf("wall VZ velocity non zero iE2=%d", iE2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW2 > -1) {
			if (iW2 >= maxelm) {
				// граничный узел
				inumber = iW2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iW2]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iW2=%lld", iW2);
#else
						printf("wall VX velocity non zero iW2=%d", iW2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iW2]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iW2=%lld", iW2);
#else
						printf("wall VY velocity non zero iW2=%d", iW2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iW2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iW2=%lld", iW2);
#else
						printf("wall VZ velocity non zero iW2=%d", iW2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN2 > -1) {
			if (iN2 >= maxelm) {
				// граничный узел
				inumber = iN2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iN2]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iN=%lld", iN2);
#else
						printf("wall VX velocity non zero iN=%d", iN2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iN2]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iN=%lld", iN2);
#else
						printf("wall VY velocity non zero iN=%d", iN2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iN2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iN=%lld", iN2);
#else
						printf("wall VZ velocity non zero iN=%d", iN2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iS2>-1) {
			if (iS2 >= maxelm) {
				// граничный узел
				inumber = iS2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iS2])>admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iS2=%lld", iS2);
#else
						printf("wall VX velocity non zero iS2=%d", iS2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iS2])>admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iS2=%lld", iS2);
#else
						printf("wall VY velocity non zero iS2=%d", iS2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iS2])>admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iS2=%lld", iS2);
#else
						printf("wall VZ velocity non zero iS2=%d", iS2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}


		if (iT2 > -1) {
			if (iT2 >= maxelm) {
				// граничный узел
				inumber = iT2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iT2]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iT2=%lld", iT2);
#else
						printf("wall VX velocity non zero iT2=%d", iT2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iT2]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iT2=%lld", iT2);
#else
						printf("wall VY velocity non zero iT2=%d", iT2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iT2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iT2=%lld", iT2);
#else
						printf("wall VZ velocity non zero iT2=%d", iT2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iB2 > -1) {
			if (iB2 >= maxelm) {
				// граничный узел
				inumber = iB2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iB2]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iB2=%lld", iB2);
#else
						printf("wall VX velocity non zero iB2=%d", iB2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iB2]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iB2=%lld", iB2);
#else
						printf("wall VY velocity non zero iB2=%d", iB2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iB2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iB2=%lld", iB2);
#else
						printf("wall VZ velocity non zero iB2=%d", iB2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iE3 > -1) {
			if (iE3 >= maxelm) {
				// граничный узел
				inumber = iE3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iE3]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iE3=%lld", iE3);
#else
						printf("wall VX velocity non zero iE3=%d", iE3);
#endif


						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iE3]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iE3=%lld", iE3);
#else
						printf("wall VY velocity non zero iE3=%d", iE3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iE3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iE3=%lld", iE3);
#else
						printf("wall VZ velocity non zero iE3=%d", iE3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW3 > -1) {
			if (iW3 >= maxelm) {
				// граничный узел
				inumber = iW3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iW3]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iW3=%lld", iW3);
#else
						printf("wall VX velocity non zero iW3=%d", iW3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iW3]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iW3=%lld", iW3);
#else
						printf("wall VY velocity non zero iW3=%d", iW3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iW3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iW3=%lld", iW3);
#else
						printf("wall VZ velocity non zero iW3=%d", iW3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN3 > -1) {
			if (iN3 >= maxelm) {
				// граничный узел
				inumber = iN3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iN3]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iN3=%lld", iN3);
#else
						printf("wall VX velocity non zero iN3=%d", iN3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iN3]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iN3=%lld", iN3);
#else
						printf("wall VY velocity non zero iN3=%d", iN3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iN3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iN3=%lld", iN3);
#else
						printf("wall VZ velocity non zero iN3=%d", iN3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iS3>-1) {
			if (iS3 >= maxelm) {
				// граничный узел
				inumber = iS3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iS3])>admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iS3=%lld", iS3);
#else
						printf("wall VX velocity non zero iS3=%d", iS3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iS3])>admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iS3=%lld", iS3);
#else
						printf("wall VY velocity non zero iS3=%d", iS3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iS3])>admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iS3=%lld", iS3);
#else
						printf("wall VZ velocity non zero iS3=%d", iS3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}


		if (iT3 > -1) {
			if (iT3 >= maxelm) {
				// граничный узел
				inumber = iT3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iT3]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iT3=%lld", iT3);
#else
						printf("wall VX velocity non zero iT3=%d", iT3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iT3]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iT3=%lld", iT3);
#else
						printf("wall VY velocity non zero iT3=%d", iT3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iT3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iT3=%lld", iT3);
#else
						printf("wall VZ velocity non zero iT3=%d", iT3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iB3 > -1) {
			if (iB3 >= maxelm) {
				// граничный узел
				inumber = iB3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iB3]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iB3=%lld", iB3);
#else
						printf("wall VX velocity non zero iB3=%d", iB3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iB3]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iB3=%lld", iB3);
#else
						printf("wall VY velocity non zero iB3=%d", iB3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iB3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iB3=%lld", iB3);
#else
						printf("wall VZ velocity non zero iB3=%d", iB3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iE4 > -1) {
			if (iE4 >= maxelm) {
				// граничный узел
				inumber = iE4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iE4]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iE4=%lld", iE4);
#else
						printf("wall VX velocity non zero iE4=%d", iE4);
#endif


						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iE4]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iE4=%lld", iE4);
#else
						printf("wall VY velocity non zero iE4=%d", iE4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iE4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iE4=%lld", iE4);
#else
						printf("wall VZ velocity non zero iE4=%d", iE4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW4 > -1) {
			if (iW4 >= maxelm) {
				// граничный узел
				inumber = iW4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iW4]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iW4=%lld", iW4);
#else
						printf("wall VX velocity non zero iW4=%d", iW4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iW4]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iW4=%lld", iW4);
#else
						printf("wall VY velocity non zero iW4=%d", iW4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iW4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iW4=%lld", iW4);
#else
						printf("wall VZ velocity non zero iW4=%d", iW4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN4 > -1) {
			if (iN4 >= maxelm) {
				// граничный узел
				inumber = iN4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iN4]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iN4=%lld", iN4);
#else
						printf("wall VX velocity non zero iN4=%d", iN4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iN4]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iN4=%lld", iN4);
#else
						printf("wall VY velocity non zero iN4=%d", iN4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iN4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iN4=%lld", iN4);
#else
						printf("wall VZ velocity non zero iN4=%d", iN4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iS4>-1) {
			if (iS4 >= maxelm) {
				// граничный узел
				inumber = iS4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iS4])>admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iS4=%lld", iS4);
#else
						printf("wall VX velocity non zero iS4=%d", iS4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iS4])>admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iS4=%lld", iS4);
#else
						printf("wall VY velocity non zero iS4=%d", iS4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iS4])>admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iS4=%lld", iS4);
#else
						printf("wall VZ velocity non zero iS4=%d", iS4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}


		if (iT4 > -1) {
			if (iT4 >= maxelm) {
				// граничный узел
				inumber = iT4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iT4]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iT4=%lld", iT4);
#else
						printf("wall VX velocity non zero iT4=%d", iT4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iT4]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iT4=%lld", iT4);
#else
						printf("wall VY velocity non zero iT4=%d", iT4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iT4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iT4=%lld", iT4);
#else
						printf("wall VZ velocity non zero iT4=%d", iT4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iB4 > -1) {
			if (iB4 >= maxelm) {
				// граничный узел
				inumber = iB4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VELOCITY_X_COMPONENT][iB4]) > admission) {
#if doubleintprecision == 1
						printf("wall VX velocity non zero iB4=%lld", iB4);
#else
						printf("wall VX velocity non zero iB4=%d", iB4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Y_COMPONENT][iB4]) > admission) {
#if doubleintprecision == 1
						printf("wall VY velocity non zero iB4=%lld", iB4);
#else
						printf("wall VY velocity non zero iB4=%d", iB4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VELOCITY_Z_COMPONENT][iB4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZ velocity non zero iB4=%lld", iB4);
#else
						printf("wall VZ velocity non zero iB4=%d", iB4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iE > -1) {
			if (iE >= maxelm) {
				// граничный узел
				inumber = iE - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iE]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iE=%lld", iE);
#else
						printf("wall VXCOR velocity non zero iE=%d", iE);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iE]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iE=%lld", iE);
#else
						printf("wall VYCOR velocity non zero iE=%d", iE);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iE]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iE=%lld", iE);
#else
						printf("wall VZCOR velocity non zero iE=%d", iE);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW > -1) {
			if (iW >= maxelm) {
				// граничный узел
				inumber = iW - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iW]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iW=%lld", iW);
#else
						printf("wall VXCOR velocity non zero iW=%d", iW);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iW]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iW=%lld", iW);
#else
						printf("wall VYCOR velocity non zero iW=%d", iW);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iW]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iW=%lld", iW);
#else
						printf("wall VZCOR velocity non zero iW=%d", iW);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN>-1) {
		if (iN>=maxelm) {
			// граничный узел
			inumber=iN-maxelm;
			if (border_neighbor[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VXCOR][iN])>admission) {
#if doubleintprecision == 1
					printf("wall VXCOR velocity non zero iN=%lld", iN);
#else
					printf("wall VXCOR velocity non zero iN=%d", iN);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VYCOR][iN])>admission) {
#if doubleintprecision == 1
					printf("wall VYCOR velocity non zero iN=%lld", iN);
#else
					printf("wall VYCOR velocity non zero iN=%d", iN);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZCOR][iN])>admission) {
#if doubleintprecision == 1
					printf("wall VZCOR velocity non zero iN=%lld", iN);
#else
					printf("wall VZCOR velocity non zero iN=%d", iN);
#endif
					
					//getchar();
					system("pause");
				}
			}
		}
		}

		if (iS>-1) {

		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (border_neighbor[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VXCOR][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VXCOR velocity non zero iS=%lld", iS);
#else
					printf("wall VXCOR velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VYCOR][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VYCOR velocity non zero iS=%lld", iS);
#else
					printf("wall VYCOR velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZCOR][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VZCOR velocity non zero iS=%lld", iS);
#else
					printf("wall VZCOR velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
			}
		}
		}

		if (iT>-1) {
		if (iT>=maxelm) {
			// граничный узел
			inumber=iT-maxelm;
			if (border_neighbor[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VXCOR][iT])>admission) {
#if doubleintprecision == 1
					printf("wall VXCOR velocity non zero iT=%lld", iT);
#else
					printf("wall VXCOR velocity non zero iT=%d", iT);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VYCOR][iT])>admission) {
#if doubleintprecision == 1
					printf("wall VYCOR velocity non zero iT=%lld", iT);
#else
					printf("wall VYCOR velocity non zero iT=%d", iT);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZCOR][iT])>admission) {
#if doubleintprecision == 1
					printf("wall VZCOR velocity non zero iT=%lld", iT);
#else
					printf("wall VZCOR velocity non zero iT=%d", iT);
#endif
					
					//getchar();
					system("pause");
				}
			}
		}
		}


		if (iB>-1) {
		if (iB>=maxelm) {
			// граничный узел
			inumber=iB-maxelm;
			if (border_neighbor[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VXCOR][iB])>admission) {
#if doubleintprecision == 1
					printf("wall VXCOR velocity non zero iB=%lld", iB);
#else
					printf("wall VXCOR velocity non zero iB=%d", iB);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VYCOR][iB])>admission) {
#if doubleintprecision == 1
					printf("wall VYCOR velocity non zero iB=%lld", iB);
#else
					printf("wall VYCOR velocity non zero iB=%d", iB);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZCOR][iB])>admission) {
#if doubleintprecision == 1
					printf("wall VZCOR velocity non zero iB=%lld", iB);
#else
					printf("wall VZCOR velocity non zero iB=%d", iB);
#endif
					
					//getchar();
					system("pause");
				}
			}
		}

	}

		if (iE2 > -1) {
			if (iE2 >= maxelm) {
				// граничный узел
				inumber = iE2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iE2]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iE2=%lld", iE2);
#else
						printf("wall VXCOR velocity non zero iE2=%d", iE2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iE2]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iE2=%lld", iE2);
#else
						printf("wall VYCOR velocity non zero iE2=%d", iE2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iE2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iE2=%lld", iE2);
#else
						printf("wall VZCOR velocity non zero iE2=%d", iE2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW2 > -1) {
			if (iW2 >= maxelm) {
				// граничный узел
				inumber = iW2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iW2]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iW2=%lld", iW2);
#else
						printf("wall VXCOR velocity non zero iW2=%d", iW2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iW2]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iW2=%lld", iW2);
#else
						printf("wall VYCOR velocity non zero iW2=%d", iW2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iW2]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iW2=%lld", iW2);
#else
						printf("wall VZCOR velocity non zero iW2=%d", iW2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN2>-1) {
			if (iN2 >= maxelm) {
				// граничный узел
				inumber = iN2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iN2])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iN2=%lld", iN2);
#else
						printf("wall VXCOR velocity non zero iN2=%d", iN2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iN2])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iN2=%lld", iN2);
#else
						printf("wall VYCOR velocity non zero iN2=%d", iN2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iN2])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iN2=%lld", iN2);
#else
						printf("wall VZCOR velocity non zero iN2=%d", iN2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iS2>-1) {

			if (iS2 >= maxelm) {
				// граничный узел
				inumber = iS2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iS2])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iS2=%lld", iS2);
#else
						printf("wall VXCOR velocity non zero iS2=%d", iS2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iS2])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iS2=%lld", iS2);
#else
						printf("wall VYCOR velocity non zero iS2=%d", iS2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iS2])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iS2=%lld", iS2);
#else
						printf("wall VZCOR velocity non zero iS2=%d", iS2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iT2>-1) {
			if (iT2 >= maxelm) {
				// граничный узел
				inumber = iT2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iT2])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iT2=%lld", iT2);
#else
						printf("wall VXCOR velocity non zero iT2=%d", iT2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iT2])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iT2=%lld", iT2);
#else
						printf("wall VYCOR velocity non zero iT2=%d", iT2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iT2])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iT2=%lld", iT2);
#else
						printf("wall VZCOR velocity non zero iT2=%d", iT2);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}


		if (iB2>-1) {
			if (iB2 >= maxelm) {
				// граничный узел
				inumber = iB2 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iB2])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iB2=%lld", iB2);
#else
						printf("wall VXCOR velocity non zero iB2=%d", iB2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iB2])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iB2=%lld", iB2);
#else
						printf("wall VYCOR velocity non zero iB2=%d", iB2);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iB2])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iB2=%lld", iB2);
#else
						printf("wall VZCOR velocity non zero iB2=%d", iB2);
#endif

						//getchar();
						system("pause");
					}
				}
			}

		}
		

		if (iE3 > -1) {
			if (iE3 >= maxelm) {
				// граничный узел
				inumber = iE3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iE3]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iE3=%lld", iE3);
#else
						printf("wall VXCOR velocity non zero iE3=%d", iE3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iE3]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iE3=%lld", iE3);
#else
						printf("wall VYCOR velocity non zero iE3=%d", iE3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iE3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iE3=%lld", iE3);
#else
						printf("wall VZCOR velocity non zero iE3=%d", iE3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW3 > -1) {
			if (iW3 >= maxelm) {
				// граничный узел
				inumber = iW3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iW3]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iW3=%lld", iW3);
#else
						printf("wall VXCOR velocity non zero iW3=%d", iW3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iW3]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iW3=%lld", iW3);
#else
						printf("wall VYCOR velocity non zero iW3=%d", iW3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iW3]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iW3=%lld", iW3);
#else
						printf("wall VZCOR velocity non zero iW3=%d", iW3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN3>-1) {
			if (iN3 >= maxelm) {
				// граничный узел
				inumber = iN3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iN3])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iN3=%lld", iN3);
#else
						printf("wall VXCOR velocity non zero iN3=%d", iN3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iN3])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iN3=%lld", iN3);
#else
						printf("wall VYCOR velocity non zero iN3=%d", iN3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iN3])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iN3=%lld", iN3);
#else
						printf("wall VZCOR velocity non zero iN3=%d", iN3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iS3>-1) {

			if (iS3 >= maxelm) {
				// граничный узел
				inumber = iS3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iS3])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iS=%lld", iS3);
#else
						printf("wall VXCOR velocity non zero iS=%d", iS3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iS3])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iS3=%lld", iS3);
#else
						printf("wall VYCOR velocity non zero iS3=%d", iS3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iS3])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iS3=%lld", iS3);
#else
						printf("wall VZCOR velocity non zero iS3=%d", iS3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iT3>-1) {
			if (iT3 >= maxelm) {
				// граничный узел
				inumber = iT3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iT3])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iT3=%lld", iT3);
#else
						printf("wall VXCOR velocity non zero iT3=%d", iT3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iT3])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iT3=%lld", iT3);
#else
						printf("wall VYCOR velocity non zero iT3=%d", iT3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iT3])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iT3=%lld", iT3);
#else
						printf("wall VZCOR velocity non zero iT3=%d", iT3);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}


		if (iB3>-1) {
			if (iB3 >= maxelm) {
				// граничный узел
				inumber = iB3 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iB3])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iB3=%lld", iB3);
#else
						printf("wall VXCOR velocity non zero iB3=%d", iB3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iB3])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iB3=%lld", iB3);
#else
						printf("wall VYCOR velocity non zero iB3=%d", iB3);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iB3])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iB3=%lld", iB3);
#else
						printf("wall VZCOR velocity non zero iB3=%d", iB3);
#endif

						//getchar();
						system("pause");
					}
				}
			}

		}
		

		if (iE4 > -1) {
			if (iE4 >= maxelm) {
				// граничный узел
				inumber = iE4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iE4]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iE4=%lld", iE4);
#else
						printf("wall VXCOR velocity non zero iE4=%d", iE4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iE4]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iE4=%lld", iE4);
#else
						printf("wall VYCOR velocity non zero iE4=%d", iE4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iE4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iE4=%lld", iE4);
#else
						printf("wall VZCOR velocity non zero iE4=%d", iE4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iW4 > -1) {
			if (iW4 >= maxelm) {
				// граничный узел
				inumber = iW4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iW4]) > admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iW4=%lld", iW4);
#else
						printf("wall VXCOR velocity non zero iW4=%d", iW4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iW4]) > admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iW4=%lld", iW4);
#else
						printf("wall VYCOR velocity non zero iW4=%d", iW4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iW4]) > admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iW4=%lld", iW4);
#else
						printf("wall VZCOR velocity non zero iW4=%d", iW4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iN4>-1) {
			if (iN4 >= maxelm) {
				// граничный узел
				inumber = iN4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iN4])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iN4=%lld", iN4);
#else
						printf("wall VXCOR velocity non zero iN4=%d", iN4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iN4])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iN4=%lld", iN4);
#else
						printf("wall VYCOR velocity non zero iN4=%d", iN4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iN4])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iN4=%lld", iN4);
#else
						printf("wall VZCOR velocity non zero iN4=%d", iN4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iS4>-1) {

			if (iS4 >= maxelm) {
				// граничный узел
				inumber = iS4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iS4])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iS4=%lld", iS4);
#else
						printf("wall VXCOR velocity non zero iS4=%d", iS4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iS4])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iS4=%lld", iS4);
#else
						printf("wall VYCOR velocity non zero iS4=%d", iS4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iS4])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iS4=%lld", iS4);
#else
						printf("wall VZCOR velocity non zero iS4=%d", iS4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}

		if (iT4>-1) {
			if (iT4 >= maxelm) {
				// граничный узел
				inumber = iT4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iT4])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iT4=%lld", iT4);
#else
						printf("wall VXCOR velocity non zero iT4=%d", iT4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iT4])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iT4=%lld", iT4);
#else
						printf("wall VYCOR velocity non zero iT4=%d", iT4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iT4])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iT4=%lld", iT4);
#else
						printf("wall VZCOR velocity non zero iT4=%d", iT4);
#endif

						//getchar();
						system("pause");
					}
				}
			}
		}


		if (iB4>-1) {
			if (iB4 >= maxelm) {
				// граничный узел
				inumber = iB4 - maxelm;
				if (border_neighbor[inumber].MCB == (ls + lw)) {
					if (fabs(potent[VXCOR][iB4])>admission) {
#if doubleintprecision == 1
						printf("wall VXCOR velocity non zero iB4=%lld", iB4);
#else
						printf("wall VXCOR velocity non zero iB4=%d", iB4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VYCOR][iB4])>admission) {
#if doubleintprecision == 1
						printf("wall VYCOR velocity non zero iB4=%lld", iB4);
#else
						printf("wall VYCOR velocity non zero iB4=%d", iB4);
#endif

						//getchar();
						system("pause");
					}
					if (fabs(potent[VZCOR][iB4])>admission) {
#if doubleintprecision == 1
						printf("wall VZCOR velocity non zero iB4=%lld", iB4);
#else
						printf("wall VZCOR velocity non zero iB4=%d", iB4);
#endif

						//getchar();
						system("pause");
					}
				}
			}

		}
		}
	
}

// коррекция граничных узлов если на них стоит условие Неймана.
void correct_boundary_volume(integer iVar, doublereal** &potent,
							 integer maxelm, int** nvtx, TOCHKA* pa, 
							 int*** neighbors_for_the_internal_node, BOUND* &border_neighbor,
							 integer ls, integer lw, WALL* w, doublereal* &relax_value, integer binterpol) {

   // Алгоритм:
   // проход по всем внутренним КО. 
   // для каждого внутреннего КО проход по всем его соседям:
   // если соседний узел окажется граничным, принадлежащим стенке wall
   // на которой стоит условие равенства давления нулю, то для всех
   // компонент скорости на этой границе стоит однородное условие Неймана,
   // и поэтому нужно вновь соблюсти граничное условие для скорректированной скорости
   // приравняв компоненты скорости на границе к значениям соответствующих компонент скорости
   // в ближайшем внутреннем узле.

    // iVar - определитель (идентефикатор) коррекируемой компоненты скорости.

	// 14 мая в код добавлена квадратичная интерполяция на границе.
	// она положительно сказывается на сходимости алгоритма SIMPLEC (немного ускоряя её).
	// Вообще говоря, главная цель граничных значений удовлетворить граничным условиям.
	// А на границе стоит однородное условие Неймана, что означает равенство значений величины F
	// в граничном и ближайшем приграничном узле. Поэтому, исходя из этих рассуждений, точное выполнение граничных
	// условий обеспечивается величиной binterpol==0. Однако значение 0 замедляет сходимость. Наибольшая скорость 
	// сходимости обеспечивается значением 2 (квадратичная интерполяция).
	// 1 использовать линейную интерполяцию на границе. (0 не использовать).
	//const integer binterpol=0; // 2 - использовать квадратичную интерполяцию.
	const bool brelax_bound = false;
	const bool brelax_val2 = true;
	const doublereal relaxboundconstvel = 1.0;

    
	
    // iP - номер центрального контрольного объёма
#pragma omp parallel for 
	for (integer iP=0; iP<maxelm; iP++) {

		integer inumber;

		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	    iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

		integer iE2=-1, iN2=-1, iT2=-1, iW2=-1, iS2=-1, iB2=-1; // номера соседних контрольных объёмов
		integer iE3=-1, iN3=-1, iT3=-1, iW3=-1, iS3=-1, iB3=-1; // номера соседних контрольных объёмов
		integer iE4=-1, iN4=-1, iT4=-1, iW4=-1, iS4=-1, iB4=-1; // номера соседних контрольных объёмов

		if (b_on_adaptive_local_refinement_mesh) {

			iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
			iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

			iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
			iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

			iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
			iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
		}
		// вычисление размеров текущего контрольного объёма:
	    doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
        volume3D(iP, nvtx, pa, dx, dy, dz);
		
		if (iE > -1) {
			if (iE >= maxelm) {
				// граничный узел
				inumber = iE - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE];
							}
							else {
								potent[iVar][iE] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE];
							}
						}
						else {
							potent[iVar][iE] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// линейная интерполяция:
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

						potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE];
							}
							else {
								potent[iVar][iE] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE];
							}
						}
						else {
							potent[iVar][iE] = potent[iVar][iP];
						}
					}
							 else if (binterpol == 1) {
								 // линейная интерполяция:
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

								 potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE] = 0.0; break;
					}
				}


			} // iE
		}

		if (iW > -1) {
			if (iW >= maxelm) {
				// граничный узел
				inumber = iW - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW];
							}
							else {
								potent[iVar][iW] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW];
							}
						}
						else {
							potent[iVar][iW] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

						potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW];
							}
							else {
								potent[iVar][iW] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW];
							}
						}
						else {
							potent[iVar][iW] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

								 potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW] = 0.0; break;
					}
				}

			} // iW
		}

		if (iN > -1) {
			if (iN >= maxelm) {
				// граничный узел
				inumber = iN - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN];
							}
							else {
								potent[iVar][iN] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN];
							}
						}
						else {
							potent[iVar][iN] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

						potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN];
							}
							else {
								potent[iVar][iN] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN];
							}
						}
						else {
							potent[iVar][iN] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }
								 
								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.

								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

								 potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iN] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN] = 0.0; break;
					}
				}

			} // iN
		}

		if (iS > -1) {
			if (iS >= maxelm) {
				// граничный узел
				inumber = iS - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS];
							}
							else {
								potent[iVar][iS] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS];
							}
						}
						else {
							potent[iVar][iS] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

						potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
					}
					//if (iVar==VY) { printf("Vs==%e, Vp==%e\n",potent[iVar][iS],potent[iVar][iP]); getchar(); } // debug
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS];
							}
							else {
								potent[iVar][iS] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS];
							}
						}
						else {
							potent[iVar][iS] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

								 potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iS] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS] = 0.0; break;
					}
				}

			} // iS
		}

		if (iT > -1) {
			if (iT >= maxelm) {
				// граничный узел
				inumber = iT - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT];
							}
							else {
								potent[iVar][iT] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT];
							}
						}
						else {
							potent[iVar][iT] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

						potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
					}
				} // pressure outlet
				else  if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT];
							}
							else {
								potent[iVar][iT] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT];
							}
						}
						else {
							potent[iVar][iT] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }


								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

								 potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iT] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT] = 0.0; break;
					}
				}

			} // iT
		}

		if (iB > -1) {
			if (iB >= maxelm) {
				// граничный узел
				inumber = iB - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB];
							}
							else {
								potent[iVar][iB] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB];
							}
						}
						else {
							potent[iVar][iB] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

						potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB];
							}
							else {
								potent[iVar][iB] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB];
							}
						}
						else {
							potent[iVar][iB] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

								 potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iB] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB] = 0.0; break;
					}
				}

			} // iB
		}


		if (iE2 > -1) {
			if (iE2 >= maxelm) {
				// граничный узел
				inumber = iE2 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE2];
							}
							else {
								potent[iVar][iE2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE2];
							}
						}
						else {
							potent[iVar][iE2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// линейная интерполяция:
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

						potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE2] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE2];
							}
							else {
								potent[iVar][iE2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE2];
							}
						}
						else {
							potent[iVar][iE2] = potent[iVar][iP];
						}
					}
							 else if (binterpol == 1) {
								 // линейная интерполяция:
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

								 potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE2] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE2] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE2] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE2] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE2] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE2] = 0.0; break;
					}
				}


			} // iE2
		}

		if (iW2 > -1) {
			if (iW2 >= maxelm) {
				// граничный узел
				inumber = iW2 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW2];
							}
							else {
								potent[iVar][iW2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW2];
							}
						}
						else {
							potent[iVar][iW2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

						potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW2] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW2];
							}
							else {
								potent[iVar][iW2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW2];
							}
						}
						else {
							potent[iVar][iW2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

								 potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW2] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW2] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW2] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW2] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW2] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW2] = 0.0; break;
					}
				}

			} // iW2
		}

		if (iN2 > -1) {
			if (iN2 >= maxelm) {
				// граничный узел
				inumber = iN2 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN2];
							}
							else {
								potent[iVar][iN2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN2];
							}
						}
						else {
							potent[iVar][iN2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

						potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN2];
							}
							else {
								potent[iVar][iN2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN2];
							}
						}
						else {
							potent[iVar][iN2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.

								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

								 potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iN2] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN2] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN2] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN2] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN2] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN2] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN2] = 0.0; break;
					}
				}

			} // iN2
		}

		if (iS2 > -1) {
			if (iS2 >= maxelm) {
				// граничный узел
				inumber = iS2 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS2];
							}
							else {
								potent[iVar][iS2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS2];
							}
						}
						else {
							potent[iVar][iS2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

						potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
					}
					//if (iVar==VY) { printf("Vs==%e, Vp==%e\n",potent[iVar][iS],potent[iVar][iP]); getchar(); } // debug
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS2];
							}
							else {
								potent[iVar][iS2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS2];
							}
						}
						else {
							potent[iVar][iS2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

								 potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iS2] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS2] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS2] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS2] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS2] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS2] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS2] = 0.0; break;
					}
				}

			} // iS2
		}

		if (iT2 > -1) {
			if (iT2 >= maxelm) {
				// граничный узел
				inumber = iT2 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT2];
							}
							else {
								potent[iVar][iT2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT2];
							}
						}
						else {
							potent[iVar][iT2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

						potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
					}
				} // pressure outlet
				else  if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT2];
							}
							else {
								potent[iVar][iT2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT2];
							}
						}
						else {
							potent[iVar][iT2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }


								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

								 potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iT2] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT2] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT2] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT2] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT2] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT2] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT2] = 0.0; break;
					}
				}

			} // iT2
		}

		if (iB2 > -1) {
			if (iB2 >= maxelm) {
				// граничный узел
				inumber = iB2 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB2];
							}
							else {
								potent[iVar][iB2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB2];
							}
						}
						else {
							potent[iVar][iB2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

						potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB2];
							}
							else {
								potent[iVar][iB2] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB2];
							}
						}
						else {
							potent[iVar][iB2] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

								 potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iB2] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB2] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB2] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB2] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB2] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB2] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB2] = 0.0; break;
					}
				}

			} // iB
		}

		if (iE3 > -1) {
			if (iE3 >= maxelm) {
				// граничный узел
				inumber = iE3 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE3];
							}
							else {
								potent[iVar][iE3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE3];
							}
						}
						else {
							potent[iVar][iE3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// линейная интерполяция:
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

						potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE3] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE3];
							}
							else {
								potent[iVar][iE3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE3];
							}
						}
						else {
							potent[iVar][iE3] = potent[iVar][iP];
						}
					}
							 else if (binterpol == 1) {
								 // линейная интерполяция:
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

								 potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE3] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE3] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE3] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE3] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE3] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE3] = 0.0; break;
					}
				}


			} // iE3
		}

		if (iW3 > -1) {
			if (iW3 >= maxelm) {
				// граничный узел
				inumber = iW3 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW3];
							}
							else {
								potent[iVar][iW3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW3];
							}
						}
						else {
							potent[iVar][iW3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

						potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW3] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW3];
							}
							else {
								potent[iVar][iW3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW3];
							}
						}
						else {
							potent[iVar][iW3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

								 potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW3] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW3] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW3] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW3] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW3] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW3] = 0.0; break;
					}
				}

			} // iW3
		}

		if (iN3 > -1) {
			if (iN3 >= maxelm) {
				// граничный узел
				inumber = iN3 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN3];
							}
							else {
								potent[iVar][iN3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN3];
							}
						}
						else {
							potent[iVar][iN3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

						potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN3];
							}
							else {
								potent[iVar][iN3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN3];
							}
						}
						else {
							potent[iVar][iN3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.

								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

								 potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iN3] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN3] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN3] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN3] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN3] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN3] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN3] = 0.0; break;
					}
				}

			} // iN3
		}

		if (iS3 > -1) {
			if (iS3 >= maxelm) {
				// граничный узел
				inumber = iS3 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS3];
							}
							else {
								potent[iVar][iS3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS3];
							}
						}
						else {
							potent[iVar][iS3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

						potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
					}
					//if (iVar==VY) { printf("Vs==%e, Vp==%e\n",potent[iVar][iS],potent[iVar][iP]); getchar(); } // debug
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS3];
							}
							else {
								potent[iVar][iS3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS3];
							}
						}
						else {
							potent[iVar][iS3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

								 potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iS3] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS3] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS3] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS3] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS3] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS3] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS3] = 0.0; break;
					}
				}

			} // iS3
		}

		if (iT3 > -1) {
			if (iT3 >= maxelm) {
				// граничный узел
				inumber = iT3 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT3];
							}
							else {
								potent[iVar][iT3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT3];
							}
						}
						else {
							potent[iVar][iT3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

						potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
					}
				} // pressure outlet
				else  if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT3];
							}
							else {
								potent[iVar][iT3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT3];
							}
						}
						else {
							potent[iVar][iT3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }


								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

								 potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iT3] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT3] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT3] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT3] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT3] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT3] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT3] = 0.0; break;
					}
				}

			} // iT3
		}

		if (iB3 > -1) {
			if (iB3 >= maxelm) {
				// граничный узел
				inumber = iB3 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB3];
							}
							else {
								potent[iVar][iB3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB3];
							}
						}
						else {
							potent[iVar][iB3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

						potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB3];
							}
							else {
								potent[iVar][iB3] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB3];
							}
						}
						else {
							potent[iVar][iB3] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

								 potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iB3] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB3] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB3] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB3] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB3] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB3] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB3] = 0.0; break;
					}
				}

			} // iB3
		}

		if (iE4 > -1) {
			if (iE4 >= maxelm) {
				// граничный узел
				inumber = iE4 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE4];
							}
							else {
								potent[iVar][iE4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE4];
							}
						}
						else {
							potent[iVar][iE4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// линейная интерполяция:
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iW, nvtx, pa, pb, W_SIDE);
						center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

						potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE4] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iE4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE4];
							}
							else {
								potent[iVar][iE4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE4];
							}
						}
						else {
							potent[iVar][iE4] = potent[iVar][iP];
						}
					}
							 else if (binterpol == 1) {
								 // линейная интерполяция:
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 potent[iVar][iE] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x + 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iE correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iW, nvtx, pa, pb, W_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][iW], nvtx, pa, pbb, WW_SIDE);

								 potent[iVar][iE] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[W_SIDE][0][iW]], potent[iVar][iW], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE4] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE4] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE4] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iE4] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iE4] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iE4] = 0.0; break;
					}
				}


			} // iE
		}

		if (iW4 > -1) {
			if (iW4 >= maxelm) {
				// граничный узел
				inumber = iW4 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW4];
							}
							else {
								potent[iVar][iW4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW4];
							}
						}
						else {
							potent[iVar][iW4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iE, nvtx, pa, pb, E_SIDE);
						center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

						potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
					// Значит скорость VY и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW4] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iW4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW4];
							}
							else {
								potent[iVar][iW4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW4];
							}
						}
						else {
							potent[iVar][iW4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 potent[iVar][iW] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x - 0.5*dx);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iW correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iE, nvtx, pa, pb, E_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][iE], nvtx, pa, pbb, EE_SIDE);

								 potent[iVar][iW] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[E_SIDE][0][iE]], potent[iVar][iE], potent[iVar][iP], pbb.x, pb.x, pp.x, pp.x - 0.5*dx);
							 }
							 break; // корректируем скорость.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW4] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW4] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW4] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iW4] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iW4] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iW4] = 0.0; break;
					}
				}

			} // iW4
		}

		if (iN4 > -1) {
			if (iN4 >= maxelm) {
				// граничный узел
				inumber = iN4 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN4];
							}
							else {
								potent[iVar][iN4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN4];
							}
						}
						else {
							potent[iVar][iN4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iS, nvtx, pa, pb, S_SIDE);
						center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

						potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iN4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN4];
							}
							else {
								potent[iVar][iN4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN4];
							}
						}
						else {
							potent[iVar][iN4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 potent[iVar][iN] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y + 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.

								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iN correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iS, nvtx, pa, pb, S_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][iS], nvtx, pa, pbb, SS_SIDE);

								 potent[iVar][iN] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[S_SIDE][0][iS]], potent[iVar][iS], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iN4] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN4] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN4] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN4] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iN4] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iN4] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iN4] = 0.0; break;
					}
				}

			} // iN4
		}

		if (iS4 > -1) {
			if (iS4 >= maxelm) {
				// граничный узел
				inumber = iS4 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS4];
							}
							else {
								potent[iVar][iS4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS4];
							}
						}
						else {
							potent[iVar][iS4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iN, nvtx, pa, pb, N_SIDE);
						center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

						potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
					}
					//if (iVar==VY) { printf("Vs==%e, Vp==%e\n",potent[iVar][iS],potent[iVar][iP]); getchar(); } // debug
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
					// Значит скорость VX и VZ в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Z_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iS4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS4];
							}
							else {
								potent[iVar][iS4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS4];
							}
						}
						else {
							potent[iVar][iS4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 potent[iVar][iS] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y - 0.5*dy);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iS correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iN, nvtx, pa, pb, N_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][iN], nvtx, pa, pbb, NN_SIDE);

								 potent[iVar][iS] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[N_SIDE][0][iN]], potent[iVar][iN], potent[iVar][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Y_COMPONENT: potent[iVar][iS4] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS4] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS4] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS4] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iS4] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iS4] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iS4] = 0.0; break;
					}
				}

			} // iS4
		}

		if (iT4 > -1) {
			if (iT4 >= maxelm) {
				// граничный узел
				inumber = iT4 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT4];
							}
							else {
								potent[iVar][iT4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT4];
							}
						}
						else {
							potent[iVar][iT4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.

						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iB, nvtx, pa, pb, B_SIDE);
						center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

						potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
					}
				} // pressure outlet
				else  if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iT4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT4];
							}
							else {
								potent[iVar][iT4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT4];
							}
						}
						else {
							potent[iVar][iT4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 potent[iVar][iT] = my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z + 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iT correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }


								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iB, nvtx, pa, pb, B_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][iB], nvtx, pa, pbb, BB_SIDE);

								 potent[iVar][iT] = my_quadratic_interpolation('+', potent[iVar][neighbors_for_the_internal_node[B_SIDE][0][iB]], potent[iVar][iB], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iT4] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT4] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT4] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT4] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iT4] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iT4] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iT4] = 0.0; break;
					}
				}

			} // iT4
		}

		if (iB4 > -1) {
			if (iB4 >= maxelm) {
				// граничный узел
				inumber = iB4 - maxelm;
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening))) {
					// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
					// Значит скорость в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла.
					if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB4];
							}
							else {
								potent[iVar][iB4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB4];
							}
						}
						else {
							potent[iVar][iB4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
					else if (binterpol == 1) {
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}

						TOCHKA pp, pb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
					}
					else if (binterpol == 2) {
						// квадратичная интерполяция.
						// не работает на АЛИС.
						if (b_on_adaptive_local_refinement_mesh) {
							printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
							system("pause");
							exit(1);
						}


						TOCHKA pp, pb, pbb;
						center_cord3D(iP, nvtx, pa, pp, 100);
						center_cord3D(iT, nvtx, pa, pb, T_SIDE);
						center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

						potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
					}
				} // pressure outlet
				else if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bsymmetry)) {
					// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
					// Значит скорость VX и VY в граничном узле нужно скорректировать записав в неё значение из ближайшего внутреннего узла,
					// так чтобы выполнялось граничное условие для скорректированной скорости.
					switch (iVar) {
					case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: if (binterpol == 0) {
						if (brelax_bound) {
							// Здесь возможно надо релаксировать к скорректированной скорости удовлетворяющей уравнению неразрывности.
							if (brelax_val2) {
								potent[iVar][iB4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB4];
							}
							else {
								potent[iVar][iB4] = relaxboundconstvel * potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB4];
							}
						}
						else {
							potent[iVar][iB4] = potent[iVar][iP]; // корректируем скорость.
						}
					}
							 else if (binterpol == 1) {
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 1) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 potent[iVar][iB] = my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z - 0.5*dz);
							 }
							 else if (binterpol == 2) {
								 // квадратичная интерполяция.
								 // не работает на АЛИС.
								 if (b_on_adaptive_local_refinement_mesh) {
									 printf("function iB correct_boundary_volume in module correct_velocity.cpp if (binterpol == 2) not worked in ALICE mesh...\n ");
									 system("pause");
									 exit(1);
								 }

								 TOCHKA pp, pb, pbb;
								 center_cord3D(iP, nvtx, pa, pp, 100);
								 center_cord3D(iT, nvtx, pa, pb, T_SIDE);
								 center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][iT], nvtx, pa, pbb, TT_SIDE);

								 potent[iVar][iB] = my_quadratic_interpolation('-', potent[iVar][neighbors_for_the_internal_node[T_SIDE][0][iT]], potent[iVar][iT], potent[iVar][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
							 }
							 break; // корректируем скорость.
					case VELOCITY_Z_COMPONENT: potent[iVar][iB4] = 0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
					}

				} // symmetry
				else if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw))) {
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB4] = w[border_neighbor[inumber].MCB - ls].Vx; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB4] = w[border_neighbor[inumber].MCB - ls].Vy; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB4] = w[border_neighbor[inumber].MCB - ls].Vz; break;
					}
				}
				else {
					// Твёрдая неподвижная стенка Stacionary WALL
					switch (iVar) {
					case VELOCITY_X_COMPONENT: potent[iVar][iB4] = 0.0; break;
					case VELOCITY_Y_COMPONENT: potent[iVar][iB4] = 0.0; break;
					case VELOCITY_Z_COMPONENT: potent[iVar][iB4] = 0.0; break;
					}
				}

			} // iB4
		}

	}

} // correct_boundary_volume

#endif
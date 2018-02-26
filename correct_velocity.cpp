// В файле correct_velocity.cpp содержится код который применяется
// для коррекции скорости так чтобы они стали удовлетворять уравнению 
// неразрывности.

// Данные функции выделены в отдельный файл 31 марта 2012 года.

// 
// Коррекцию скорости надо разделить на две процедуры.
// Первая процедура для внутренних КО. При этом используется интерполляция для 
// значений поправки давления на границе расчётной области. 
// 
// Вторая процедура коррекции для граничных КО. С интерполяцией скорости изнутри области на границу.
// Коррекцию граничных КО можно и отключить, т.к. она может не понадобиться.
// Возможно расчёт будет правилен даже в случае отсутствия коррекции граничных узлов.
// Это нужно проверить непосредственно путём численного расчёта.
// На данный момент решено включить коррекцию скорости на границе путём интерполляции изнутри области.

// Коррекция граничных КО включена.
// По поводу корекции граничных КО смотри mysolver SIMPLE Algorithm.
// В этом файле исправлены все интерполяционные формулы с 
// учётом предельного перехода к предельным случаям.
// Квадратичная интерполляция, увеличивающая скорость сходимости, введена в этом модуле 13,14,15 мая 2012 года.

#ifndef CORRECT_VELOCITYv_0_07_CPP
#define CORRECT_VELOCITYv_0_07_CPP 1



// коррекция скорости для внутренних КО.
// Скоректированные скорости удовлетворяют уравнению 
// несжимаемости.
void correct_internal_volume(integer iP, integer iVar, equation3D** sl,   
			 integer** nvtx, doublereal** &potent, integer maxelm, doublereal* alpha,
			 TOCHKA* pa, ALICE_PARTITION** sosedi, integer iternumber) {

    // квадратичная интерполляция добавлена 14 мая 2012 года.
	// Ещё как вариант можно применить МНК интерполляцию (построить прямую по трём точкам).
	integer interpol=0; // 0 - без всякой интерполляции, 1- линейная интерполляция, 2- квадратичная интерполляция.

	// Внимание можно использовать и второй порядок тоже. формула из книги Г.З.Гарбера  
	// работает. Формула из книги Г.З. Гарбера обеспечивает точность производной от 
	// поправки давления второго порядка на неравномерной сетке (она найдена как производная 
	// от параболы построенная по трём точкам).
	// порядок точности нахождения первой производной от давления.
	integer iderivative_pressure=1; // 1 - первый порядок, 2 - второй порядок.


    // Внутренний узел и его соседи:
    // iP принадлежит интервалу 0..maxelm-1

    // iP - номер центрального контрольного объёма
	integer iE=-1, iN=-1, iT=-1, iW=-1, iS=-1, iB=-1; // номера соседних контрольных объёмов
	iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;
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
       case VX : ds=dy*dz; dv=ds*dx; dl=dx; break;
       case VY : ds=dx*dz; dv=ds*dy; dl=dy;  break;
       case VZ : ds=dx*dy; dv=ds*dz; dl=dz;  break;
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
		center_cord3D(iW, nvtx, pa, pb,WSIDE);
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
		case 1 : // линейная интерполляция
		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iE, nvtx, pa, pb,ESIDE);
		         PAmW=my_linear_interpolation('-', potent[PAM][iP], potent[PAM][iE], pp.x, pb.x, pp.x-0.5*dx); 
			    break;
		case 2 : // квадратичная интерполляция.

		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iE, nvtx, pa, pb,ESIDE);
			     center_cord3D(sosedi[ESIDE][iE].iNODE1, nvtx, pa, pbb,EE);
					
			     PAmW=my_quadratic_interpolation('-', potent[PAM][sosedi[ESIDE][iE].iNODE1], potent[PAM][iE], potent[PAM][iP], pbb.x , pb.x, pp.x, pp.x-0.5*dx); 
			    break;
		default : // значение которое получено после решения СЛАУ.
			      PAmW=potent[PAM][iW];
			    break;
		} // end switch
		
	}
	if (!bE) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iE, nvtx, pa, pb,ESIDE);
		hxplus=fabs(pb.x-pp.x);

		PAmE=potent[PAM][iE];
	} else {
        // узел E граничный
		hxplus=0.5*dx;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0 : // значение которое получено после решения СЛАУ. 
			    PAmE=potent[PAM][iE];
			    break;
		case 1 : // линейная интерполяция

		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iW, nvtx, pa, pb,WSIDE);
		         PAmE=my_linear_interpolation('+', potent[PAM][iP], potent[PAM][iW], pp.x, pb.x, pp.x+0.5*dx);

			    break;
		case 2 : // квадратичная интерполляция.
						     
		         center_cord3D(iP, nvtx, pa, pp,100);
		         center_cord3D(iW, nvtx, pa, pb,WSIDE);
				 center_cord3D(sosedi[WSIDE][iW].iNODE1, nvtx, pa, pbb,WW);
					
				 PAmE = my_quadratic_interpolation('+', potent[PAM][sosedi[WSIDE][iW].iNODE1], potent[PAM][iW], potent[PAM][iP], pbb.x, pb.x, pp.x, pp.x + 0.5*dx);
			    break;
		default : // значение которое получено после решения СЛАУ.
			     PAmE=potent[PAM][iE]; 
			     break;
		}

	}
	if (!bS) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iS, nvtx, pa, pb,SSIDE);
		hyminus=fabs(pp.y-pb.y);

		PAmS=potent[PAM][iS];
	} else {
		// узел S граничный
		hyminus=0.5*dy;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0 : // значение которое получено после решения СЛАУ.  
			    PAmS=potent[PAM][iS];
			    break;
		case 1 :
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iN, nvtx, pa, pb,NSIDE);
		        PAmS=my_linear_interpolation('-', potent[PAM][iP], potent[PAM][iN], pp.y, pb.y, pp.y-0.5*dy);

			    break;
		case 2 :  
			    // квадратичная интерполляция.

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iN, nvtx, pa, pb,NSIDE);
				center_cord3D(sosedi[NSIDE][iN].iNODE1, nvtx, pa, pbb,NN);
					
				PAmS = my_quadratic_interpolation('-', potent[PAM][sosedi[NSIDE][iN].iNODE1], potent[PAM][iN], potent[PAM][iP], pbb.y, pb.y, pp.y, pp.y - 0.5*dy);

			    break;
		default : // значение которое получено после решения СЛАУ. 
			    PAmS=potent[PAM][iS];
			    break;
		}

	}
	if (!bN) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iN, nvtx, pa, pb,NSIDE);
		hyplus=fabs(pb.y-pp.y);

		PAmN=potent[PAM][iN];
	} else {
		// узел N граничный
		hyplus=0.5*dy;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0 : // значение которое получено после решения СЛАУ. 
			    PAmN=potent[PAM][iN];
			    break;
		case 1 :
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iS, nvtx, pa, pb,SSIDE);
		        PAmN=my_linear_interpolation('+', potent[PAM][iP], potent[PAM][iS], pp.y, pb.y, pp.y+0.5*dy);
			    break;
		case 2 :
			    // квадратичная интерполляция.

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iS, nvtx, pa, pb,SSIDE);
				center_cord3D(sosedi[SSIDE][iS].iNODE1, nvtx, pa, pbb,SS);
					
				PAmN = my_quadratic_interpolation('+', potent[PAM][sosedi[SSIDE][iS].iNODE1], potent[PAM][iS], potent[PAM][iP], pbb.y, pb.y, pp.y, pp.y + 0.5*dy);
			    break;
		default : 
			    // значение которое получено после решения СЛАУ. 
			    PAmN=potent[PAM][iN];
			    break;
		}

		 
	} 

	if (!bB) {
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iB, nvtx, pa, pb,BSIDE);
		hzminus=fabs(pp.z-pb.z);

		PAmB=potent[PAM][iB];
	} else {
		// узел B граничный
		hzminus=0.5*dz;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0 :
			    // значение которое получено после решения СЛАУ. 
			    PAmB=potent[PAM][iB];
			    break;
		case 1 : 
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iT, nvtx, pa, pb,TSIDE);
		        PAmB=my_linear_interpolation('-', potent[PAM][iP], potent[PAM][iT], pp.z, pb.z, pp.z-0.5*dz);
			    break;
		case 2 : // квадратичная интерполляция.

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iT, nvtx, pa, pb,TSIDE);
				center_cord3D(sosedi[TSIDE][iT].iNODE1, nvtx, pa, pbb,TTSIDE);
					
				PAmB = my_quadratic_interpolation('-', potent[PAM][sosedi[TSIDE][iT].iNODE1], potent[PAM][iT], potent[PAM][iP], pbb.z, pb.z, pp.z, pp.z - 0.5*dz);
			    break;
		default :
			    // значение которое получено после решения СЛАУ. 
			    PAmB=potent[PAM][iB];
			    break;
		}

	}
	if (!bT) { 
		TOCHKA pp, pb;
		center_cord3D(iP, nvtx, pa, pp,100);
		center_cord3D(iT, nvtx, pa, pb,TSIDE);
		hzplus=fabs(pb.z-pp.z);

		PAmT=potent[PAM][iT];
	} else {
		// узел T граничный
		hzplus=0.5*dz;
		TOCHKA pp,pb,pbb;

		switch (interpol) {
		case 0 : // значение которое получено после решения СЛАУ. 
			    PAmT=potent[PAM][iT];
			    break;
		case 1 : 
			    // линейная интерполяция

		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iB, nvtx, pa, pb,BSIDE);
		        PAmT=my_linear_interpolation('+', potent[PAM][iP], potent[PAM][iB], pp.z, pb.z, pp.z+0.5*dz);
			    break;
		case 2 :// квадратичная интерполляция.
                
		        center_cord3D(iP, nvtx, pa, pp,100);
		        center_cord3D(iB, nvtx, pa, pb,BSIDE);
				center_cord3D(sosedi[BSIDE][iB].iNODE1, nvtx, pa, pbb,BB);
					
				PAmT = my_quadratic_interpolation('+', potent[PAM][sosedi[BSIDE][iB].iNODE1], potent[PAM][iB], potent[PAM][iP], pbb.z, pb.z, pp.z, pp.z + 0.5*dz);
			    break;
		default : // значение которое получено после решения СЛАУ. 
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
		case VX : if (iderivative_pressure==1) {
			          // естественная аппроксимация первого порядка.
				      deltaP=(fwplus*PAmW+(1-fwplus)*PAmP);
			          deltaP-=(feplus*PAmE+(1-feplus)*PAmP); 
				      gradP=deltaP/dl; // первая производная от давления.
				  }
				  else if (iderivative_pressure==2) {
					  // аппроксимация второго порядка точности.
					  // обязательно нужен знак минус иначе скорость
					  // будет направлена в другую сторону.
					  gradP=-rgradF(PAmW, PAmP, PAmE, hxminus, hxplus);
				  }
				  break;
		case VY : if (iderivative_pressure==1) {
			          // естественная аппроксимация первого порядка.
				      deltaP=(fsplus*PAmS+(1-fsplus)*PAmP);
			          deltaP-=(fnplus*PAmN+(1-fnplus)*PAmP);
				      gradP=deltaP/dl; // первая производная от давления.
				  }
				  else if (iderivative_pressure==2) {
					  // аппроксимация второго порядка точности.
					  // обязательно нужен знак минус иначе скорость
					  // будет направлена в другую сторону.
					  gradP=-rgradF(PAmS, PAmP, PAmN, hyminus, hyplus);
				  }
			      break;
        case VZ : if (iderivative_pressure==1) {
			          // естественная аппроксимация первого порядка.
			          deltaP=(fbplus*PAmB+(1-fbplus)*PAmP);
			          deltaP-=(ftplus*PAmT+(1-ftplus)*PAmP);
				      gradP=deltaP/dl; // первая производная от давления.
				  }
				  else if (iderivative_pressure==2) {
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

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скоректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume

// коррекция скорости для внутренних КО.
// Скоректированные скорости удовлетворяют уравнению 
// несжимаемости.
void correct_internal_volume2(integer iP, integer iVar, equation3D** sl,   
			 integer** nvtx, doublereal** &potent, integer maxelm, doublereal* alpha,
			 TOCHKA* pa, ALICE_PARTITION** sosedi, integer iternumber) {

	// Процедура коррекции скорости во внутренних контрольных объёмах на 
	// основе осреднённого градиента давления.
    
    // Внутренний узел и его соседи:
    // iP принадлежит интервалу 0..maxelm-1

	
	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
    volume3D(iP, nvtx, pa, dx, dy, dz);
	    

	doublereal ds=0.0, dv=0.0; // площадь грани, длина интервала и объём контрольного объёма.
	switch (iVar) {
       case VX : ds=dy*dz; dv=ds*dx; break;
       case VY : ds=dx*dz; dv=ds*dy; break;
       case VZ : ds=dx*dy; dv=ds*dz; break;
	}
	
	// Случай граничного узла G учитывается правильно,
	// т.к. в этом случае fgplus==1.0; 

	// Линейная интерполяция давления на грань КО.
	doublereal deltaP=0.0, gradP=0.0;
	switch (iVar) {
		case VX : gradP=-potent[GRADXPAM][iP];
				  break;
		case VY : gradP=-potent[GRADYPAM][iP];
				  break;
        case VZ : gradP=-potent[GRADZPAM][iP];
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

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скоректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume2

// коррекция скорости для внутренних КО.
// Скоректированные скорости удовлетворяют уравнению 
// несжимаемости (неразрывности).
// begin 20 июня 2012 года.
void correct_internal_volume3(integer iP, integer iVar, doublereal** prop,   
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
		case VX : gradP=-potent[GRADXPAM][iP];
				  break;
		case VY : gradP=-potent[GRADYPAM][iP];
				  break;
        case VZ : gradP=-potent[GRADZPAM][iP];
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
	// В формуле для tau уже учитывается два варианта рабочего алгоритма : SIMPLE и SIMPLEC.
	potent[iVar][iP]+=tau[iP]*gradP/prop[RHO][iP];

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скоректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume3

// коррекция скорости для внутренних КО.
// Скоректированные скорости удовлетворяют уравнению 
// несжимаемости (неразрывности).
// реализовано 23 июня 2012 года.
// На основе сглаженного псевдовремени, а точнее
// трёх скалярных полей псевдовремени
void correct_internal_volume4(integer iP, integer iVar, doublereal** prop,   
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
		case VX : gradPAM=-potent[GRADXPAM][iP];
			      tauP=tau[VX][iP];
				  break;
		case VY : gradPAM=-potent[GRADYPAM][iP];
			      tauP=tau[VY][iP];
				  break;
        case VZ : gradPAM=-potent[GRADZPAM][iP];
			      tauP=tau[VZ][iP];
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
	// В формуле для tau уже учитывается два варианта рабочего алгоритма : SIMPLE и SIMPLEC.
	potent[iVar][iP]+=tauP*gradPAM/prop[RHO][iP];

	// Коррекция скорости в граничных узлах должна быть всецело подчинена
	// выполнению граничных условий. А именно 1. если в граничном узле
	// по скорости стоит условие Дирихле то никакой коррекции скорости 
	// в граничном узле не требуется. 2. если в граничном узле по скорости
	// стоит однородное условие Неймана то необходимо снести скоректированную 
	// скорость из ближайшего внутреннего узла в граничный.
	
    
} // correct_internal_volume4

// коррекция массового потока на грани КО.
// begin 25 июня 2012 года.
void correct_mf(doublereal** &mfcurrentretune, doublereal** potent,  doublereal** tau,
	TOCHKA* pa, ALICE_PARTITION** sosedi, integer** nvtx, integer maxelm,
				BOUND* &sosedb, integer ls, integer lw, WALL* w, doublereal** prop_b) {

					

	doublereal** mfloc = NULL;
	mfloc=new doublereal*[maxelm];
	for (integer i=0; i<maxelm; i++) {
		mfloc[i]=new doublereal[6];
	}

	

	for (integer iP=0; iP<maxelm; iP++) {

		// iP - номер центрального контрольного объёма
	    integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	    iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;
	    

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

		// значение псевдовремени на грани контрольного объёма.
		doublereal taue, tauw, taun, taus, taut, taub;
        // интерполяция псевдовремени сделана так, чтобы выполнялись 
	    // предельные соотношения.
	    if (!bE) taue=tau[VX][iE]*tau[VX][iP]/(feplus*tau[VX][iE]+(1.0-feplus)*tau[VX][iP]); else taue=tau[VX][iE];
	    if (!bW) tauw=tau[VX][iW]*tau[VX][iP]/(fwplus*tau[VX][iW]+(1.0-fwplus)*tau[VX][iP]); else tauw=tau[VX][iW];
	    if (!bN) taun=tau[VY][iN]*tau[VY][iP]/(fnplus*tau[VY][iN]+(1.0-fnplus)*tau[VY][iP]); else taun=tau[VY][iN];
	    if (!bS) taus=tau[VY][iS]*tau[VY][iP]/(fsplus*tau[VY][iS]+(1.0-fsplus)*tau[VY][iP]); else taus=tau[VY][iS];
        if (!bT) taut=tau[VZ][iT]*tau[VZ][iP]/(ftplus*tau[VZ][iT]+(1.0-ftplus)*tau[VZ][iP]); else taut=tau[VZ][iT];
	    if (!bB) taub=tau[VZ][iB]*tau[VZ][iP]/(fbplus*tau[VZ][iB]+(1.0-fbplus)*tau[VZ][iP]); else taub=tau[VZ][iB];

		// Градиент поправки давления на грани контрольного объёма.
		doublereal gradpame, gradpamw, gradpamn, gradpams, gradpamt, gradpamb;
		if (!bE) gradpame=feplus*potent[GRADXPAM][iE]+(1.0-feplus)*potent[GRADXPAM][iP]; else gradpame=potent[GRADXPAM][iE];
        if (!bW) gradpamw=fwplus*potent[GRADXPAM][iW]+(1.0-fwplus)*potent[GRADXPAM][iP]; else gradpamw=potent[GRADXPAM][iW];
	    if (!bN) gradpamn=fnplus*potent[GRADYPAM][iN]+(1.0-fnplus)*potent[GRADYPAM][iP]; else gradpamn=potent[GRADYPAM][iN];
        if (!bS) gradpams=fsplus*potent[GRADYPAM][iS]+(1.0-fsplus)*potent[GRADYPAM][iP]; else gradpams=potent[GRADYPAM][iS];
        if (!bT) gradpamt=ftplus*potent[GRADZPAM][iT]+(1.0-ftplus)*potent[GRADZPAM][iP]; else gradpamt=potent[GRADZPAM][iT];
        if (!bB) gradpamb=fbplus*potent[GRADZPAM][iB]+(1.0-fbplus)*potent[GRADZPAM][iP]; else gradpamb=potent[GRADZPAM][iB];

		// Наконец вычисление скорректированного массового 
		// потока на грани КО.
		mfloc[iP][ESIDE]=mfcurrentretune[iP][ESIDE]-taue*gradpame*dy*dz;
		mfloc[iP][WSIDE]=mfcurrentretune[iP][WSIDE]-tauw*gradpamw*dy*dz;
		mfloc[iP][NSIDE]=mfcurrentretune[iP][NSIDE]-taun*gradpamn*dx*dz;
		mfloc[iP][SSIDE]=mfcurrentretune[iP][SSIDE]-taus*gradpams*dx*dz;
		mfloc[iP][TSIDE]=mfcurrentretune[iP][TSIDE]-taut*gradpamt*dx*dy;
		mfloc[iP][BSIDE]=mfcurrentretune[iP][BSIDE]-taub*gradpamb*dx*dy;

	}

	// Однако есть границы где массовый поток задан пользователем,
	// очевидно это надо как-то учесть.
	for (integer iP=0; iP<maxelm; iP++) {
		// iP - номер центрального контрольного объёма
	    integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	    iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;
	    

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
		const doublereal relax_bound = 1.0;

		if (bE) {
			integer inumber=iE-maxelm;
			if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening)) {
				// Выходная граница оставляем всё как есть
				mfloc[iP][ESIDE] = relax_bound*(mfloc[iP][ESIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][ESIDE];
			}
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry) {
				mfloc[iP][ESIDE]=0.0;
			} else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				// заданная скорость на входной границе.
				mfloc[iP][ESIDE]=prop_b[RHO][inumber]*w[sosedb[inumber].MCB-ls].Vx*dy*dz; // заданный массовый поток.
			}
			else {
				// твёрдая неподвижная стенка по умолчанию
				mfloc[iP][ESIDE]=0.0;
			}
		}

		if (bW) {
			integer inumber=iW-maxelm;
			if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure ||  w[sosedb[inumber].MCB - ls].bopening)) {
				// Выходная граница оставляем всё как есть
				mfloc[iP][WSIDE] = relax_bound*(mfloc[iP][WSIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][WSIDE];
			}
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry) {
				mfloc[iP][WSIDE]=0.0;
			} else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				// заданная скорость на входной границе.
				mfloc[iP][WSIDE]=prop_b[RHO][inumber]*w[sosedb[inumber].MCB-ls].Vx*dy*dz; // заданный массовый поток.
			}
			else {
				// твёрдая неподвижная стенка по умолчанию
				mfloc[iP][WSIDE]=0.0;
			}
		}

		if (bN) {
			integer inumber=iN-maxelm;
			if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure ||  w[sosedb[inumber].MCB - ls].bopening)) {
				// Выходная граница оставляем всё как есть
				mfloc[iP][NSIDE] = relax_bound*(mfloc[iP][NSIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][NSIDE];
			}
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry) {
				mfloc[iP][NSIDE]=0.0;
			} else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				// заданная скорость на входной границе.
				mfloc[iP][NSIDE]=prop_b[RHO][inumber]*w[sosedb[inumber].MCB-ls].Vy*dx*dz; // заданный массовый поток.
			}
			else {
				// твёрдая неподвижная стенка по умолчанию
				mfloc[iP][NSIDE]=0.0;
			}
		}

		if (bS) {
			integer inumber=iS-maxelm;
			if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure ||  w[sosedb[inumber].MCB - ls].bopening)) {
				// Выходная граница оставляем всё как есть
				mfloc[iP][SSIDE] = relax_bound*(mfloc[iP][SSIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][SSIDE];
			}
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry) {
				mfloc[iP][SSIDE]=0.0;
			} else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				// заданная скорость на входной границе.
				mfloc[iP][SSIDE]=prop_b[RHO][inumber]*w[sosedb[inumber].MCB-ls].Vy*dx*dz; // заданный массовый поток.
			}
			else {
				// твёрдая неподвижная стенка по умолчанию
				mfloc[iP][SSIDE]=0.0;
			}
		}

		if (bT) {
			integer inumber=iT-maxelm;
			if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure ||  w[sosedb[inumber].MCB - ls].bopening)) {
				// Выходная граница оставляем всё как есть
				mfloc[iP][TSIDE] = relax_bound*(mfloc[iP][TSIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][TSIDE];
			}
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry) {
				mfloc[iP][TSIDE]=0.0;
			} else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				// заданная скорость на входной границе.
				mfloc[iP][TSIDE]=prop_b[RHO][inumber]*w[sosedb[inumber].MCB-ls].Vz*dx*dy; // заданный массовый поток.
			}
			else {
				// твёрдая неподвижная стенка по умолчанию
				mfloc[iP][TSIDE]=0.0;
			}
		}

		if (bB) {
			integer inumber=iB-maxelm;
			if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure ||  w[sosedb[inumber].MCB - ls].bopening)) {
				// Выходная граница оставляем всё как есть
				mfloc[iP][BSIDE] = relax_bound*(mfloc[iP][BSIDE]) + (1.0 - relax_bound)*mfcurrentretune[iP][BSIDE];
			}
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry) {
				mfloc[iP][BSIDE]=0.0;
			} else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				// заданная скорость на входной границе.
				mfloc[iP][BSIDE]=prop_b[RHO][inumber]*w[sosedb[inumber].MCB-ls].Vz*dx*dy; // заданный массовый поток.
			}
			else {
				// твёрдая неподвижная стенка по умолчанию
				mfloc[iP][BSIDE]=0.0;
			}
		}

	}

	// Обратное копирование.
	for (integer iG=0; iG<6; iG++) {
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

void iscorrectmf(doublereal** &mf,
							 integer maxelm, 
							 ALICE_PARTITION** sosedi, BOUND* &sosedb,
							 integer ls, integer lw, WALL* w) {
    integer iP=0;
	integer inumber;
    // iP - номер центрального контрольного объёма
	for (iP=0; iP<maxelm; iP++) {
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	    iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;

		if (iE>=maxelm) {
			// граничный узел
			inumber=iE-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(mf[iP][ESIDE])>admission) {
#if doubleintprecision == 1
					printf("wall mf flux velocity non zero iE=%lld", iE);
#else
					printf("wall mf flux velocity non zero iE=%d", iE);
#endif

					
					//getchar();
					system("pause");
				}
				
			}
		}
		if (iW>=maxelm) {
			// граничный узел
			inumber=iW-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(mf[iP][WSIDE])>admission) {
#if doubleintprecision == 1
					printf("wall mf flux velocity non zero iW=%lld", iW);
#else
					printf("wall mf flux velocity non zero iW=%d", iW);
#endif
					
					//getchar();
					system("pause");
				}
				
			}
		}

		if (iN>=maxelm) {
			// граничный узел
			inumber=iN-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(mf[iP][NSIDE])>admission) {
#if doubleintprecision == 1
					printf("wall mf flux velocity non zero iN=%lld", iN);
#else
					printf("wall mf flux velocity non zero iN=%d", iN);
#endif
					
					//getchar();
					system("pause");
				}
				
			}
		}
		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(mf[iP][SSIDE])>admission) {
#if doubleintprecision == 1
					printf("wall mf flux velocity non zero iS=%lld", iS);
#else
					printf("wall mf flux velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				
			}
		}
		if (iT>=maxelm) {
			// граничный узел
			inumber=iT-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(mf[iP][TSIDE])>admission) {
#if doubleintprecision == 1
					printf("wall mf flux velocity non zero iT=%lld", iT);
#else
					printf("wall mf flux velocity non zero iT=%d", iT);
#endif
					
					//getchar();
					system("pause");
				}
				
			}
		}
		if (iB>=maxelm) {
			// граничный узел
			inumber=iB-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(mf[iP][BSIDE])>admission) {
#if doubleintprecision == 1
					printf("wall mf flux velocity non zero iB=%lld", iB);
#else
					printf("wall mf flux velocity non zero iB=%d", iB);
#endif
					
					//getchar();
					system("pause");
				}
				
			}
		}


	}

}

void iscorrectOk(doublereal** &potent,
							 integer maxelm, 
							 ALICE_PARTITION** sosedi, BOUND* &sosedb,
							 integer ls, integer lw, WALL* w)
{
	integer iP=0;
	integer inumber;
    // iP - номер центрального контрольного объёма
	for (iP=0; iP<maxelm; iP++) {
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	    iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;

		if (iE>=maxelm) {
			// граничный узел
			inumber=iE-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VX][iE])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iE=%lld", iE);
#else
					printf("wall VX velocity non zero iE=%d", iE);
#endif

					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VY][iE])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iE=%lld", iE);
#else
					printf("wall VY velocity non zero iE=%d", iE);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZ][iE])>admission) {
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
		if (iW>=maxelm) {
			// граничный узел
			inumber=iW-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VX][iW])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iW=%lld", iW);
#else
					printf("wall VX velocity non zero iW=%d", iW);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VY][iW])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iW=%lld", iW);
#else
					printf("wall VY velocity non zero iW=%d", iW);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZ][iW])>admission) {
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

		if (iN>=maxelm) {
			// граничный узел
			inumber=iN-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VX][iN])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iN=%lld", iN);
#else
					printf("wall VX velocity non zero iN=%d", iN);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VY][iN])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iN=%lld", iN);
#else
					printf("wall VY velocity non zero iN=%d", iN);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZ][iN])>admission) {
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
		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VX][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iS=%lld", iS);
#else
					printf("wall VX velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VY][iS])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iS=%lld", iS);
#else
					printf("wall VY velocity non zero iS=%d", iS);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZ][iS])>admission) {
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

		if (iT>=maxelm) {
			// граничный узел
			inumber=iT-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VX][iT])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iT=%lld", iT);
#else
					printf("wall VX velocity non zero iT=%d", iT);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VY][iT])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iT=%lld", iT);
#else
					printf("wall VY velocity non zero iT=%d", iT);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZ][iT])>admission) {
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
		if (iB>=maxelm) {
			// граничный узел
			inumber=iB-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VX][iB])>admission) {
#if doubleintprecision == 1
					printf("wall VX velocity non zero iB=%lld", iB);
#else
					printf("wall VX velocity non zero iB=%d", iB);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VY][iB])>admission) {
#if doubleintprecision == 1
					printf("wall VY velocity non zero iB=%lld", iB);
#else
					printf("wall VY velocity non zero iB=%d", iB);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZ][iB])>admission) {
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

		if (iE>=maxelm) {
			// граничный узел
			inumber=iE-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VXCOR][iE])>admission) {
#if doubleintprecision == 1
					printf("wall VXCOR velocity non zero iE=%lld", iE);
#else
					printf("wall VXCOR velocity non zero iE=%d", iE);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VYCOR][iE])>admission) {
#if doubleintprecision == 1
					printf("wall VYCOR velocity non zero iE=%lld", iE);
#else
					printf("wall VYCOR velocity non zero iE=%d", iE);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZCOR][iE])>admission) {
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
		if (iW>=maxelm) {
			// граничный узел
			inumber=iW-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
				if (fabs(potent[VXCOR][iW])>admission) {
#if doubleintprecision == 1
					printf("wall VXCOR velocity non zero iW=%lld", iW);
#else
					printf("wall VXCOR velocity non zero iW=%d", iW);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VYCOR][iW])>admission) {
#if doubleintprecision == 1
					printf("wall VYCOR velocity non zero iW=%lld", iW);
#else
					printf("wall VYCOR velocity non zero iW=%d", iW);
#endif
					
					//getchar();
					system("pause");
				}
				if (fabs(potent[VZCOR][iW])>admission) {
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

		if (iN>=maxelm) {
			// граничный узел
			inumber=iN-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
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
		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
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

		if (iT>=maxelm) {
			// граничный узел
			inumber=iT-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
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
		if (iB>=maxelm) {
			// граничный узел
			inumber=iB-maxelm;
			if (sosedb[inumber].MCB==(ls+lw)) {
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
	
}

// коррекция граничных узлов если на них стоит условие Неймана.
void correct_boundary_volume(integer iVar, doublereal** &potent,
							 integer maxelm, integer** nvtx, TOCHKA* pa, 
							 ALICE_PARTITION** sosedi, BOUND* &sosedb,
							 integer ls, integer lw, WALL* w, doublereal* &relax_value) {

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
	// сходимости обеспечивается значением 2 (квадратичная интерполляция).
	// 1 использовать линейную интерполяцию на границе. (0 не использовать).
	integer binterpol=0; // 2 - использовать квадратичную интерполляцию.
	bool brelax_bound = false;
	bool brelax_val2 = true;
	const doublereal relaxboundconstvel = 1.0;

    integer iP=0;
	integer inumber;
    // iP - номер центрального контрольного объёма
	for (iP=0; iP<maxelm; iP++) {
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	    iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	    iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;

		// вычисление размеров текущего контрольного объёма:
	    doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
        volume3D(iP, nvtx, pa, dx, dy, dz);
		
		if (iE>=maxelm) {
			// граничный узел
			inumber=iE-maxelm;
			if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening))) {
				// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
				// Значит скорость в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла.
				if (binterpol==0) {
					if (brelax_bound) {
						// Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						if (brelax_val2) {
							potent[iVar][iE] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE];
						}
						else {
							potent[iVar][iE] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE];
						}
					}
					else {
						potent[iVar][iE] = potent[iVar][iP]; // корректируем скорость.
					}
				}
				else if (binterpol==1) {
					// линейная интерполяция:
					TOCHKA pp,pb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iW, nvtx, pa, pb,WSIDE);
		            potent[iVar][iE]=my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x+0.5*dx);
				}
				else if (binterpol==2) {
					// квадратичная интерполляция.

					TOCHKA pp,pb,pbb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iW, nvtx, pa, pb,WSIDE);
					center_cord3D(sosedi[WSIDE][iW].iNODE1, nvtx, pa, pbb,WW);
					
					potent[iVar][iE]=my_quadratic_interpolation('+', potent[iVar][sosedi[WSIDE][iW].iNODE1], potent[iVar][iW], potent[iVar][iP], pbb.x , pb.x, pp.x, pp.x+0.5*dx);
				}
			} // pressure outlet
			else if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry)) {
				// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
				// Значит скорость VY и VZ в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла,
				// так чтобы выполнялось граничное условие для скорректированной скорости.
				switch (iVar) {
				   case VX : potent[iVar][iE]=0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
				   case VY : case VZ : if (binterpol==0) {
					   if (brelax_bound) {
						   // Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						   if (brelax_val2) {
							   potent[iVar][iE] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iE];
						   }
						   else {
							   potent[iVar][iE] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iE];
						   }
					   }
					   else {
						   potent[iVar][iE] = potent[iVar][iP];
					   }
							 }
							 else if (binterpol==1) {
								 // линейная интерполяция:
					             TOCHKA pp,pb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iW, nvtx, pa, pb,WSIDE);
		                         potent[iVar][iE]=my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iW], pp.x, pb.x, pp.x+0.5*dx);
							 }
							 else if (binterpol==2) {
								 // квадратичная интерполляция.

					             TOCHKA pp,pb,pbb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iW, nvtx, pa, pb,WSIDE);
					             center_cord3D(sosedi[WSIDE][iW].iNODE1, nvtx, pa, pbb,WW);
					
					             potent[iVar][iE]=my_quadratic_interpolation('+', potent[iVar][sosedi[WSIDE][iW].iNODE1], potent[iVar][iW], potent[iVar][iP], pbb.x , pb.x, pp.x, pp.x+0.5*dx);
				             }
							 break; // корректируем скорость.
				}
				 
			} // symmetry
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				switch (iVar) {
				  case VX : potent[iVar][iE]=w[sosedb[inumber].MCB-ls].Vx; break;
				  case VY : potent[iVar][iE]=w[sosedb[inumber].MCB-ls].Vy; break;
				  case VZ : potent[iVar][iE]=w[sosedb[inumber].MCB-ls].Vz; break;
				}				
			}
			else {
				// Твёрдая неподвижная стенка Stacionary WALL
                switch (iVar) {
				  case VX : potent[iVar][iE]=0.0; break;
				  case VY : potent[iVar][iE]=0.0; break;
				  case VZ : potent[iVar][iE]=0.0; break;
				}				
			}


		} // iE

		if (iW>=maxelm) {
			// граничный узел
			inumber=iW-maxelm;
			if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening))) {
				// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
				// Значит скорость в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла.
				if (binterpol==0) {
					if (brelax_bound) {
						// Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						if (brelax_val2) {
							potent[iVar][iW] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW];
						}
						else {
							potent[iVar][iW] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW];
						}
					}
					else {
						potent[iVar][iW] = potent[iVar][iP]; // корректируем скорость.
					}
				}
				else if (binterpol==1) {
                    TOCHKA pp,pb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iE, nvtx, pa, pb,ESIDE);
		            potent[iVar][iW]=my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x-0.5*dx);
				}
				else if (binterpol==2) {
					// квадратичная интерполляция.

					TOCHKA pp,pb,pbb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iE, nvtx, pa, pb,ESIDE);
					center_cord3D(sosedi[ESIDE][iE].iNODE1, nvtx, pa, pbb,EE);
					
					potent[iVar][iW]=my_quadratic_interpolation('-', potent[iVar][sosedi[ESIDE][iE].iNODE1], potent[iVar][iE], potent[iVar][iP], pbb.x , pb.x, pp.x, pp.x-0.5*dx);
				}
			} // pressure outlet
			else if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry)) {
				// граница симметрии: по VY и VZ стоит однородное условие Неймана, а для VX==0.0;
				// Значит скорость VY и VZ в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла,
				// так чтобы выполнялось граничное условие для скорректированной скорости.
				switch (iVar) {
				   case VX : potent[iVar][iW]=0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
				   case VY : case VZ : if (binterpol==0) {
					          if (brelax_bound) {
						          // Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
								  if (brelax_val2) {
									  potent[iVar][iW] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iW];
								  }
								  else {
									  potent[iVar][iW] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iW];
								  }
							}
					          else {
						          potent[iVar][iW] = potent[iVar][iP]; // корректируем скорость.
					          }
					      }
				             else if (binterpol==1) {
                                TOCHKA pp,pb;
		                        center_cord3D(iP, nvtx, pa, pp,100);
		                        center_cord3D(iE, nvtx, pa, pb,ESIDE);
		                        potent[iVar][iW]=my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iE], pp.x, pb.x, pp.x-0.5*dx);
				             }
							 else if (binterpol==2) {
								 // квадратичная интерполляция.

					             TOCHKA pp,pb,pbb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iE, nvtx, pa, pb,ESIDE);
					             center_cord3D(sosedi[ESIDE][iE].iNODE1, nvtx, pa, pbb,EE);
					
					             potent[iVar][iW]=my_quadratic_interpolation('-', potent[iVar][sosedi[ESIDE][iE].iNODE1], potent[iVar][iE], potent[iVar][iP], pbb.x , pb.x, pp.x, pp.x-0.5*dx);
				             }
					         break; // корректируем скорость.
				}
				 
			} // symmetry
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				switch (iVar) {
				  case VX : potent[iVar][iW]=w[sosedb[inumber].MCB-ls].Vx; break;
				  case VY : potent[iVar][iW]=w[sosedb[inumber].MCB-ls].Vy; break;
				  case VZ : potent[iVar][iW]=w[sosedb[inumber].MCB-ls].Vz; break;
				}
			}
			else {
				// Твёрдая неподвижная стенка Stacionary WALL
                switch (iVar) {
				  case VX : potent[iVar][iW]=0.0; break;
				  case VY : potent[iVar][iW]=0.0; break;
				  case VZ : potent[iVar][iW]=0.0; break;
				}
			}

		} // iW

		if (iN>=maxelm) {
			// граничный узел
			inumber=iN-maxelm;
			if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening))) {
				// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
				// Значит скорость в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла.
				if (binterpol==0) {
					if (brelax_bound) {
						// Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						if (brelax_val2) {
							potent[iVar][iN] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN];
						}
						else {
							potent[iVar][iN] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN];
						}
					}
					else {
						potent[iVar][iN] = potent[iVar][iP]; // корректируем скорость.
					}
				}
				else if (binterpol==1) {
					 TOCHKA pp,pb;
		             center_cord3D(iP, nvtx, pa, pp,100);
		             center_cord3D(iS, nvtx, pa, pb,SSIDE);
		             potent[iVar][iN]=my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y+0.5*dy);
				}
				else if (binterpol==2) {
					// квадратичная интерполляция.

					TOCHKA pp,pb,pbb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iS, nvtx, pa, pb,SSIDE);
					center_cord3D(sosedi[SSIDE][iS].iNODE1, nvtx, pa, pbb,SS);
					
					potent[iVar][iN]=my_quadratic_interpolation('+', potent[iVar][sosedi[SSIDE][iS].iNODE1], potent[iVar][iS], potent[iVar][iP], pbb.y , pb.y, pp.y, pp.y+0.5*dy);
				}
			} // pressure outlet
			else if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry)) {
				// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
				// Значит скорость VX и VZ в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла,
				// так чтобы выполнялось граничное условие для скорректированной скорости.
				switch (iVar) {
				   case VX : case VZ : if (binterpol==0) {
					   if (brelax_bound) {
						   // Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						   if (brelax_val2) {
							   potent[iVar][iN] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iN];
						   }
						   else {
							   potent[iVar][iN] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iN];
						   }
					   }
					   else {
						   potent[iVar][iN] = potent[iVar][iP]; // корректируем скорость.
					   }
				             }
				             else if (binterpol==1) {
					             TOCHKA pp,pb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iS, nvtx, pa, pb,SSIDE);
		                         potent[iVar][iN]=my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iS], pp.y, pb.y, pp.y+0.5*dy);
				             }
							 else if (binterpol==2) {
								 // квадратичная интерполляция.

					             TOCHKA pp,pb,pbb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iS, nvtx, pa, pb,SSIDE);
					             center_cord3D(sosedi[SSIDE][iS].iNODE1, nvtx, pa, pbb,SS);
					
					             potent[iVar][iN]=my_quadratic_interpolation('+', potent[iVar][sosedi[SSIDE][iS].iNODE1], potent[iVar][iS], potent[iVar][iP], pbb.y , pb.y, pp.y, pp.y+0.5*dy);
				             }
					         break; // корректируем скорость.
				   case VY : potent[iVar][iN]=0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
				}
				 
			} // symmetry
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				switch (iVar) {
				  case VX : potent[iVar][iN]=w[sosedb[inumber].MCB-ls].Vx; break;
				  case VY : potent[iVar][iN]=w[sosedb[inumber].MCB-ls].Vy; break;
				  case VZ : potent[iVar][iN]=w[sosedb[inumber].MCB-ls].Vz; break;
				}
			}
			else {
				// Твёрдая неподвижная стенка Stacionary WALL
                switch (iVar) {
				  case VX : potent[iVar][iN]=0.0; break;
				  case VY : potent[iVar][iN]=0.0; break;
				  case VZ : potent[iVar][iN]=0.0; break;
				}
			}

		} // iN

		if (iS>=maxelm) {
			// граничный узел
			inumber=iS-maxelm;
			if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening))) {
				// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
				// Значит скорость в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла.
				if (binterpol==0) {
					if (brelax_bound) {
						// Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						if (brelax_val2) {
							potent[iVar][iS] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS];
						}
						else {
							potent[iVar][iS] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS];
						}
					}
					else {
						potent[iVar][iS] = potent[iVar][iP]; // корректируем скорость.
					}
				}
				else if (binterpol==1) {
					TOCHKA pp,pb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iN, nvtx, pa, pb,NSIDE);
		            potent[iVar][iS]=my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y-0.5*dy);
				}
				else if (binterpol==2) {
					// квадратичная интерполляция.

					TOCHKA pp,pb,pbb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iN, nvtx, pa, pb,NSIDE);
					center_cord3D(sosedi[NSIDE][iN].iNODE1, nvtx, pa, pbb,NN);
					
					potent[iVar][iS]=my_quadratic_interpolation('-', potent[iVar][sosedi[NSIDE][iN].iNODE1], potent[iVar][iN], potent[iVar][iP], pbb.y , pb.y, pp.y, pp.y-0.5*dy);
				}
				//if (iVar==VY) { printf("Vs==%e, Vp==%e\n",potent[iVar][iS],potent[iVar][iP]); getchar(); } // debug
			} // pressure outlet
			else if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry)) {
				// граница симметрии: по VX и VZ стоит однородное условие Неймана, а для VY==0.0;
				// Значит скорость VX и VZ в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла,
				// так чтобы выполнялось граничное условие для скорректированной скорости.
				switch (iVar) {
				   case VX : case VZ : if (binterpol==0) {
					   if (brelax_bound) {
						      // Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						   if (brelax_val2) {
							   potent[iVar][iS] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iS];
						   }
						   else {
							   potent[iVar][iS] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iS];
						   }
					        }
					        else {
						        potent[iVar][iS] = potent[iVar][iP]; // корректируем скорость.
					          }
				             }
				             else if (binterpol==1) {
					             TOCHKA pp,pb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iN, nvtx, pa, pb,NSIDE);
		                         potent[iVar][iS]=my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iN], pp.y, pb.y, pp.y-0.5*dy);
				             }
							 else if (binterpol==2) {
								 // квадратичная интерполляция.

					             TOCHKA pp,pb,pbb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iN, nvtx, pa, pb,NSIDE);
					             center_cord3D(sosedi[NSIDE][iN].iNODE1, nvtx, pa, pbb,NN);
					
					             potent[iVar][iS]=my_quadratic_interpolation('-', potent[iVar][sosedi[NSIDE][iN].iNODE1], potent[iVar][iN], potent[iVar][iP], pbb.y , pb.y, pp.y, pp.y-0.5*dy);
				             }
					         break; // корректируем скорость.
				   case VY : potent[iVar][iS]=0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
				}
				 
			} // symmetry
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				switch (iVar) {
				  case VX : potent[iVar][iS]=w[sosedb[inumber].MCB-ls].Vx; break;
				  case VY : potent[iVar][iS]=w[sosedb[inumber].MCB-ls].Vy; break;
				  case VZ : potent[iVar][iS]=w[sosedb[inumber].MCB-ls].Vz; break;
				}
			}
			else {
				// Твёрдая неподвижная стенка Stacionary WALL
                switch (iVar) {
				  case VX : potent[iVar][iS]=0.0; break;
				  case VY : potent[iVar][iS]=0.0; break;
				  case VZ : potent[iVar][iS]=0.0; break;
				}
			}

		} // iS

		if (iT>=maxelm) {
			// граничный узел
			inumber=iT-maxelm;
			if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening))) {
				// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
				// Значит скорость в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла.
				if (binterpol==0) {
					if (brelax_bound) {
						// Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						if (brelax_val2) {
							potent[iVar][iT] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT];
						}
						else {
							potent[iVar][iT] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT];
						}
					}
					else {
						potent[iVar][iT] = potent[iVar][iP]; // корректируем скорость.
					}
				}
				else if (binterpol==1) {
					TOCHKA pp,pb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iB, nvtx, pa, pb,BSIDE);
		            potent[iVar][iT]=my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z+0.5*dz);
				}
				else if (binterpol==2) {
					// квадратичная интерполляция.

					TOCHKA pp,pb,pbb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iB, nvtx, pa, pb,BSIDE);
					center_cord3D(sosedi[BSIDE][iB].iNODE1, nvtx, pa, pbb,BB);
					
					potent[iVar][iT]=my_quadratic_interpolation('+', potent[iVar][sosedi[BSIDE][iB].iNODE1], potent[iVar][iB], potent[iVar][iP], pbb.z , pb.z, pp.z, pp.z+0.5*dz);
				}
			} // pressure outlet
			else  if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry)) {
				// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
				// Значит скорость VX и VY в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла,
				// так чтобы выполнялось граничное условие для скорректированной скорости.
				switch (iVar) {
				   case VX : case VY : if (binterpol==0) {
					   if (brelax_bound) {
						   // Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						   if (brelax_val2) {
							   potent[iVar][iT] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iT];
						   }
						   else {
							   potent[iVar][iT] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iT];
						   }
					   }
					   else {
						   potent[iVar][iT] = potent[iVar][iP]; // корректируем скорость.
					   }
				             }
				             else if (binterpol==1) {
					            TOCHKA pp,pb;
		                        center_cord3D(iP, nvtx, pa, pp,100);
		                        center_cord3D(iB, nvtx, pa, pb,BSIDE);
		                        potent[iVar][iT]=my_linear_interpolation('+', potent[iVar][iP], potent[iVar][iB], pp.z, pb.z, pp.z+0.5*dz);
				             }
							 else if (binterpol==2) {
								 // квадратичная интерполляция.

								 TOCHKA pp,pb,pbb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iB, nvtx, pa, pb,BSIDE);
								 center_cord3D(sosedi[BSIDE][iB].iNODE1, nvtx, pa, pbb,BB);

								 potent[iVar][iT]=my_quadratic_interpolation('+', potent[iVar][sosedi[BSIDE][iB].iNODE1], potent[iVar][iB], potent[iVar][iP], pbb.z , pb.z, pp.z, pp.z+0.5*dz);
							 }
					         break; // корректируем скорость.
				   case VZ : potent[iVar][iT]=0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
				}
				 
			} // symmetry
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				switch (iVar) {
				  case VX : potent[iVar][iT]=w[sosedb[inumber].MCB-ls].Vx; break;
				  case VY : potent[iVar][iT]=w[sosedb[inumber].MCB-ls].Vy; break;
				  case VZ : potent[iVar][iT]=w[sosedb[inumber].MCB-ls].Vz; break;
				}
			}
			else {
				// Твёрдая неподвижная стенка Stacionary WALL
                switch (iVar) {
				  case VX : potent[iVar][iT]=0.0; break;
				  case VY : potent[iVar][iT]=0.0; break;
				  case VZ : potent[iVar][iT]=0.0; break;
				}
			}

		} // iT

		if (iB>=maxelm) {
			// граничный узел
			inumber=iB-maxelm;
			if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure ||  w[sosedb[inumber].MCB - ls].bopening))) {
				// на этой границе фиксировано давление значит по всем скоростям стоят условия Неймана.
				// Значит скорость в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла.
				if (binterpol==0) {
					if (brelax_bound) {
						// Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						if (brelax_val2) {
							potent[iVar][iB] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB];
						}
						else {
							potent[iVar][iB] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB];
						}
					}
					else {
						potent[iVar][iB] = potent[iVar][iP]; // корректируем скорость.
					}
				}
				else if (binterpol==1) {
					TOCHKA pp,pb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iT, nvtx, pa, pb,TSIDE);
		            potent[iVar][iB]=my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z-0.5*dz);
				}
				else if (binterpol==2) {
					// квадратичная интерполляция.

					TOCHKA pp,pb,pbb;
		            center_cord3D(iP, nvtx, pa, pp,100);
		            center_cord3D(iT, nvtx, pa, pb,TSIDE);
					center_cord3D(sosedi[TSIDE][iT].iNODE1, nvtx, pa, pbb,TTSIDE);
					
					potent[iVar][iB]=my_quadratic_interpolation('-', potent[iVar][sosedi[TSIDE][iT].iNODE1], potent[iVar][iT], potent[iVar][iP], pbb.z , pb.z, pp.z, pp.z-0.5*dz);
				}
			} // pressure outlet
			else if (((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bsymmetry)) {
				// граница симметрии: по VX и VY стоит однородное условие Неймана, а для VZ==0.0;
				// Значит скорость VX и VY в граничном узле нужно скоректировать записав в неё значение из ближайшего внутреннего узла,
				// так чтобы выполнялось граничное условие для скорректированной скорости.
				switch (iVar) {
				   case VX : case VY : if (binterpol==0) {
					   if (brelax_bound) {
						   // Здесь возможно надо релаксировать к скоректированной скорости удовлетворяющей уравнению неразрывности.
						   if (brelax_val2) {
							   potent[iVar][iB] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*relax_value[iB];
						   }
						   else {
							   potent[iVar][iB] = relaxboundconstvel*potent[iVar][iP] + (1.0 - relaxboundconstvel)*potent[iVar][iB];
						   }
					   }
					   else {
						   potent[iVar][iB] = potent[iVar][iP]; // корректируем скорость.
					   }
				             }
				             else if (binterpol==1) {
					             TOCHKA pp,pb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iT, nvtx, pa, pb,TSIDE);
		                         potent[iVar][iB]=my_linear_interpolation('-', potent[iVar][iP], potent[iVar][iT], pp.z, pb.z, pp.z-0.5*dz);
				             }
							 else if (binterpol==2) {
								 // квадратичная интерполляция.

					             TOCHKA pp,pb,pbb;
		                         center_cord3D(iP, nvtx, pa, pp,100);
		                         center_cord3D(iT, nvtx, pa, pb,TSIDE);
					             center_cord3D(sosedi[TSIDE][iT].iNODE1, nvtx, pa, pbb,TTSIDE);
					
					             potent[iVar][iB]=my_quadratic_interpolation('-', potent[iVar][sosedi[TSIDE][iT].iNODE1], potent[iVar][iT], potent[iVar][iP], pbb.z , pb.z, pp.z, pp.z-0.5*dz);
				             }
					         break; // корректируем скорость.
				   case VZ : potent[iVar][iB]=0.0; break; // по физическому смыслу эта компонента скорости равна нулю.
				}
				 
			} // symmetry
			else if ((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw))) {
				switch (iVar) {
				  case VX : potent[iVar][iB]=w[sosedb[inumber].MCB-ls].Vx; break;
				  case VY : potent[iVar][iB]=w[sosedb[inumber].MCB-ls].Vy; break;
				  case VZ : potent[iVar][iB]=w[sosedb[inumber].MCB-ls].Vz; break;
				}
			}
			else {
				// Твёрдая неподвижная стенка Stacionary WALL
                switch (iVar) {
				  case VX : potent[iVar][iB]=0.0; break;
				  case VY : potent[iVar][iB]=0.0; break;
				  case VZ : potent[iVar][iB]=0.0; break;
				}
			}

		} // iB

	}

} // correct_boundary_volume

#endif
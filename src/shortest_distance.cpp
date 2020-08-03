// Файл shortest_distance.cpp содержит код для вычисления
// кратчайшего расстояния до ближайшей стенки.

// Суть модуля: в приложении 2 книги А.Ю.Снегирёва
// Высокопроизводительные вычисления в технической физике.
// Численное моделирование турбулентных течений. 
// Санкт-Петербург, Издательство Политехнического университета 2009.
// говорится, что в ANSYS CFX для вычисления кратчайшего расстояния до
// ближайшей стенки решается эллиптическое уравнение для функции FI в правой
// части уравнения стоит -1. На твердой стенке расстояние до которой требуется
// вычислить стоит нулевое условие Дирихле, а на остальных границах которые не
// являются стенкой стоят однородные условия Неймана. Само кратчайшее расстояние 
// вычисляется как некая функция зависящая от градиентов найденной функции FI.
// Данное уравнение эвляется эллиптическим с положительно определённой симметричной
// матрицей для обращения которой может быть применён ICCG solver.
// В данном модуле реализуется программный код находящий кратчайшее расстояние до стенки
// путём решения эллиптического уравнения в частных производных.
// begin 10 апреля 2012 года.

// Кратчайшее расстояние до стенки используется в частности в Zero Equation Turbulence Model.

#ifndef SHORTEST_DISTANCE_TO_THE_WALL_CPP
#define SHORTEST_DISTANCE_TO_THE_WALL_CPP 1

#define FIDISTW 0 // потенциал FI
#define GRADXFI 1 // и его градиенты
#define GRADYFI 2
#define GRADZFI 3

// Вычисление градиентов потенциала FI в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
void green_gauss_FI(integer iP, doublereal** &potent, integer** nvtx, TOCHKA* pa,
	ALICE_PARTITION** neighbors_for_the_internal_node, integer maxelm, bool bbond) {
	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond == true то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][iP].iNODE1; iN=neighbors_for_the_internal_node[N_SIDE][iP].iNODE1; iT=neighbors_for_the_internal_node[T_SIDE][iP].iNODE1; iW=neighbors_for_the_internal_node[W_SIDE][iP].iNODE1; iS=neighbors_for_the_internal_node[S_SIDE][iP].iNODE1; iB=neighbors_for_the_internal_node[B_SIDE][iP].iNODE1;

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

	// линейная интерполяция скорости VX на грань КО.
    doublereal fe, fw, fn, fs, ft, fb;
	if (!bbond) {
		// внутренние КО.

	    if (!bE) fe=feplus*potent[FIDISTW][iE]+(1.0-feplus)*potent[FIDISTW][iP]; else fe=potent[FIDISTW][iE];
        if (!bW) fw=fwplus*potent[FIDISTW][iW]+(1.0-fwplus)*potent[FIDISTW][iP]; else fw=potent[FIDISTW][iW];
	    if (!bN) fn=fnplus*potent[FIDISTW][iN]+(1.0-fnplus)*potent[FIDISTW][iP]; else fn=potent[FIDISTW][iN];
        if (!bS) fs=fsplus*potent[FIDISTW][iS]+(1.0-fsplus)*potent[FIDISTW][iP]; else fs=potent[FIDISTW][iS];
        if (!bT) ft=ftplus*potent[FIDISTW][iT]+(1.0-ftplus)*potent[FIDISTW][iP]; else ft=potent[FIDISTW][iT];
        if (!bB) fb=fbplus*potent[FIDISTW][iB]+(1.0-fbplus)*potent[FIDISTW][iP]; else fb=potent[FIDISTW][iB];

	    potent[GRADXFI][iP]=(fe-fw)/dx;
	    potent[GRADYFI][iP]=(fn-fs)/dy;
	    potent[GRADZFI][iP]=(ft-fb)/dz;
	}
	else {
		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

		if (bE) {
			potent[GRADXFI][iE]=potent[GRADXFI][iP]+(dxe/dxw)*(potent[GRADXFI][iP]-potent[GRADXFI][iW]);
			potent[GRADYFI][iE]=potent[GRADYFI][iP]+(dxe/dxw)*(potent[GRADYFI][iP]-potent[GRADYFI][iW]);
			potent[GRADZFI][iE]=potent[GRADZFI][iP]+(dxe/dxw)*(potent[GRADZFI][iP]-potent[GRADZFI][iW]);
		}

		if (bW) {
			potent[GRADXFI][iW]=potent[GRADXFI][iP]+(dxw/dxe)*(potent[GRADXFI][iP]-potent[GRADXFI][iE]);
			potent[GRADYFI][iW]=potent[GRADYFI][iP]+(dxw/dxe)*(potent[GRADYFI][iP]-potent[GRADYFI][iE]);
			potent[GRADZFI][iW]=potent[GRADZFI][iP]+(dxw/dxe)*(potent[GRADZFI][iP]-potent[GRADZFI][iE]);
		}

		if (bN) {
			potent[GRADXFI][iN]=potent[GRADXFI][iP]+(dyn/dys)*(potent[GRADXFI][iP]-potent[GRADXFI][iS]);
			potent[GRADYFI][iN]=potent[GRADYFI][iP]+(dyn/dys)*(potent[GRADYFI][iP]-potent[GRADYFI][iS]);
			potent[GRADZFI][iN]=potent[GRADZFI][iP]+(dyn/dys)*(potent[GRADZFI][iP]-potent[GRADZFI][iS]);
		}

		if (bS) {
			potent[GRADXFI][iS]=potent[GRADXFI][iP]+(dys/dyn)*(potent[GRADXFI][iP]-potent[GRADXFI][iN]);
			potent[GRADYFI][iS]=potent[GRADYFI][iP]+(dys/dyn)*(potent[GRADYFI][iP]-potent[GRADYFI][iN]);
			potent[GRADZFI][iS]=potent[GRADZFI][iP]+(dys/dyn)*(potent[GRADZFI][iP]-potent[GRADZFI][iN]);
		}

		if (bT) {
			potent[GRADXFI][iT]=potent[GRADXFI][iP]+(dzt/dzb)*(potent[GRADXFI][iP]-potent[GRADXFI][iB]);
			potent[GRADYFI][iT]=potent[GRADYFI][iP]+(dzt/dzb)*(potent[GRADYFI][iP]-potent[GRADYFI][iB]);
			potent[GRADZFI][iT]=potent[GRADZFI][iP]+(dzt/dzb)*(potent[GRADZFI][iP]-potent[GRADZFI][iB]);
		}

		if (bB) {
			potent[GRADXFI][iB]=potent[GRADXFI][iP]+(dzb/dzt)*(potent[GRADXFI][iP]-potent[GRADXFI][iT]);
			potent[GRADYFI][iB]=potent[GRADYFI][iP]+(dzb/dzt)*(potent[GRADYFI][iP]-potent[GRADYFI][iT]);
			potent[GRADZFI][iB]=potent[GRADZFI][iP]+(dzb/dzt)*(potent[GRADZFI][iP]-potent[GRADZFI][iT]);
		}
	}

} // green_gauss_FI

// учёт граничных условий для вычисления расстояния до ближайшей
// твёрдой неподвижной поверхности.
void my_elmatr_quad_short_dist_bound(integer inumber, integer maxelm, 
							  bool bDirichlet, BOUND* border_neighbor, integer ls, integer lw,
							  WALL* w, equation3D_bon* &slb, doublereal dbeta,
							  TOCHKA* pa, integer** nvtx, doublereal* potent
							  ) 
{

	 // potent - потенциал FI на основе градиентов которого 
	 // вычисляется кратчайшее расстояние до ближайшей стенки.
	 // В данном методе значения potent могут использоваться для граничного условия
	 // Неймана повышенного порядка точности при соответствующем значении dbeta.

     // bDirichlet == true осуществляется сборка только граничных условий Дирихле.
     // bDirichlet == false осуществляется сборка только однородных условий Неймана.

     // inumber - номер граничного КО.
	 // inumber изменяется от 0..maxbound-1

     bool bDirichleti=false; // bDirichleti,  i - internal


     // Сначала запишем граничные условия Дирихле
	 if ((border_neighbor[inumber].MCB<(ls+lw)) && (border_neighbor[inumber].MCB>=ls) && (!w[border_neighbor[inumber].MCB-ls].bsymmetry) && (!w[border_neighbor[inumber].MCB-ls].bpressure)&&(!w[border_neighbor[inumber].MCB-ls].bopening)) {
		// граничное условие Дирихле
		// Задана скорость на границе
        // Это не граница симметрии и не выходная граница.

		 doublereal vel_mag=sqrt(w[border_neighbor[inumber].MCB-ls].Vx*w[border_neighbor[inumber].MCB-ls].Vx+
			          w[border_neighbor[inumber].MCB-ls].Vy*w[border_neighbor[inumber].MCB-ls].Vy+
					  w[border_neighbor[inumber].MCB-ls].Vz*w[border_neighbor[inumber].MCB-ls].Vz);

		 doublereal epsilon0=1e-32;

		 

		 if (fabs(vel_mag)<epsilon0) {
			 // скорость на стенке равна нулю.
			 bDirichleti=true;
		 } else bDirichleti=false;
        
	}
	else if (( (border_neighbor[inumber].MCB==(ls+lw)) ||(border_neighbor[inumber].MCB<ls)) ) { // 
		// источник тоже является твёрдой неподвижной стенкой.

        // граничное условие Дирихле
		bDirichleti=true;

	}
	else  {
		// Все остальные границы: выходная, с условием симметрии, входная.
		// На них должно стоять однородное условие Неймана.
        bDirichleti=false;
	}

	 if (bDirichlet&&bDirichleti) {
			 // условие Дирихле функция равна нулю.
			 slb[inumber].aw=1.0;
		     slb[inumber].ai=0.0;
			 slb[inumber].b=0.0; // на твёрдой неподвижной стенке функция принимает значение 0.
			 slb[inumber].iI=-1; // не присутствует в матрице
		     slb[inumber].iW=border_neighbor[inumber].iB;
	}
	if ((!bDirichlet)&&(!bDirichleti)) {
		// однородное условие Неймана.
		doublereal dl, dS, deltal, fiplus;

		switch (border_neighbor[inumber].Norm) {
		case E_SIDE: dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
                 dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
                 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
                 deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
				 deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                 fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
			    break;
		case N_SIDE: 
			     dl = pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
				 deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                 fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				 break;
       case T_SIDE:  
			     dl = pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
				 deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                 fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
                 break;
		case W_SIDE: 
			     dl = pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
                 dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
    			 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
				 deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                 fiplus=0.5*dl/deltal;
     			 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
    			 break;
         case S_SIDE:
			     dl = pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
				 deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                 fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				 break;
		 case B_SIDE:  
			     dl = pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
                 deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
    			 deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                 fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				break;
		} // end switch
	}

}

// Составляет матрицу для уравнения 
// поправки давления
void my_elmatr_quad_short_dist(integer iP, equation3D* &sl, equation3D_bon* &slb,  
	TOCHKA* pa, integer** nvtx, ALICE_PARTITION** neighbors_for_the_internal_node, integer maxelm, doublereal dbeta) {

   
	doublereal eps=1e-37; // для отделения вещественного нуля.
	

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][iP].iNODE1; iN=neighbors_for_the_internal_node[N_SIDE][iP].iNODE1; iT=neighbors_for_the_internal_node[T_SIDE][iP].iNODE1;
	iW=neighbors_for_the_internal_node[W_SIDE][iP].iNODE1; iS=neighbors_for_the_internal_node[S_SIDE][iP].iNODE1; iB=neighbors_for_the_internal_node[B_SIDE][iP].iNODE1;
	sl[iP].iP=iP;
	sl[iP].iE=iE; sl[iP].iN=iN; 
	sl[iP].iS=iS; sl[iP].iW=iW;
    sl[iP].iT=iT; sl[iP].iB=iB;


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


    doublereal De, Dw, Ds, Dn, Dt, Db; // диффузионный поток через грани КО.
	

	if (!bE) {
		if (bW) De=dbeta*dy*dz/dxe;
		else De=dy*dz/dxe;
	}
	else De=dbeta*dy*dz/dxe;

	if (!bW) {
		if (bE) Dw=dbeta*dy*dz/dxw;
		else Dw=dy*dz/dxw; 
	}
	else Dw=dbeta*dy*dz/dxw;

	if (!bN) {
		if (bS) Dn=dbeta*dx*dz/dyn;
		else Dn=dx*dz/dyn; 
	}
	else Dn=dbeta*dx*dz/dyn;


	if (!bS) {
		if (bN) Ds=dbeta*dx*dz/dys;
		else Ds=dx*dz/dys;
	}
	else Ds=dbeta*dx*dz/dys;


	if (!bT) {
		if (bB) Dt=dbeta*dx*dy/dzt;
		else Dt=dx*dy/dzt;
	} 
	else Dt=dbeta*dx*dy/dzt;


	if (!bB) {
		if (bT) Db=dbeta*dx*dy/dzb;
		else Db=dx*dy/dzb;
	}
	else Db=dbeta*dx*dy/dzb;

	// Число Пекле равно нулю т.к. это чисто эллиптическое уравнение.
	sl[iP].ae=De; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	sl[iP].ap=sl[iP].ae+sl[iP].aw+sl[iP].an+sl[iP].as+sl[iP].at+sl[iP].ab;

	doublereal dSc=-1.0;

	sl[iP].b=dSc*dx*dy*dz;

	// Симметризация СЛАУ:

	// для потенциала FI получается эллиптическое уравнение с SPD матрицей.
	
	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

	    // 1. Учёт условия Дирихле:
		if ((iE>=maxelm) && (fabs(slb[iE-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].ae*slb[iE-maxelm].b/slb[iE-maxelm].aw;
			sl[iP].ae=0.0;
			sl[iP].iE=-1; // не входит в матрицу СЛАУ.
		}
		if ((iW>=maxelm) && (fabs(slb[iW-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].aw*slb[iW-maxelm].b/slb[iW-maxelm].aw;
			sl[iP].aw=0.0;
			sl[iP].iW=-1; // не входит в матрицу СЛАУ.
		}
		if ((iN>=maxelm) && (fabs(slb[iN-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].an*slb[iN-maxelm].b/slb[iN-maxelm].aw;
			sl[iP].an=0.0;
			sl[iP].iN=-1; // не входит в матрицу СЛАУ.
		}
		if ((iS>=maxelm) && (fabs(slb[iS-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].as*slb[iS-maxelm].b/slb[iS-maxelm].aw;
			sl[iP].as=0.0;
			sl[iP].iS=-1; // не входит в матрицу СЛАУ.
		}
		if ((iT>=maxelm) && (fabs(slb[iT-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].at*slb[iT-maxelm].b/slb[iT-maxelm].aw;
			sl[iP].at=0.0;
			sl[iP].iT=-1; // не входит в матрицу СЛАУ.
		}
		if ((iB>=maxelm) && (fabs(slb[iB-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].ab*slb[iB-maxelm].b/slb[iB-maxelm].aw;
			sl[iP].ab=0.0;
			sl[iP].iB=-1; // не входит в матрицу СЛАУ.
		}

} // внутренность матрицы для вычисления кратчайшего расстояния до стенки.

// Вычисление расстояния до ближайшей стенки.
void calcdistwallCFX(FLOW &f, integer ls, integer lw, WALL* w) {

	bool bprintmessage=true; // печатать ли сообщения на консоль.

	// Выделение памяти
	doublereal** potent=NULL;
	potent=new doublereal*[4];
	for (integer i=0; i<4; i++) potent[i]=new doublereal[f.maxelm+f.maxbound];
	// матрица
	equation3D_bon* slb = NULL;
		slb = new equation3D_bon[f.maxbound];
	equation3D* sl = NULL;
		sl=new equation3D[f.maxelm];
	doublereal* rthdsd = NULL;
		rthdsd=new doublereal[f.maxelm + f.maxbound];
	// инициализация
	for (integer i=0; i<4; i++) for (integer j=0; j<f.maxelm+f.maxbound; j++) potent[i][j]=0.0;
	// сборка матрицы

	// Граничные условия Дирихле обязательно 
	// должны собираться в первую очередь
	doublereal dbeta=1.0; // первый порядок точности.
    for (integer i=0; i<f.maxbound; i++) {
		my_elmatr_quad_short_dist_bound(i, f.maxelm, 
			              true, f.border_neighbor, ls, lw, w,
						  slb, dbeta, f.pa, f.nvtx,
						  potent[FIDISTW]);
	}

	// Собираем однородные условия Неймана.
    // последний параметр bool bDirichlet равен false.
    for (integer i=0; i<f.maxbound; i++) {
		my_elmatr_quad_short_dist_bound(i, f.maxelm, 
			              false, f.border_neighbor, ls, lw, w,
						  slb, dbeta, f.pa, f.nvtx,
						  potent[FIDISTW]);
	}

	// Собираем матрицу для внутренних контрольных объёмов.
    for (integer iP=0; iP<f.maxelm; iP++) {
		my_elmatr_quad_short_dist(iP, sl, slb,  
						f.pa, f.nvtx, f.neighbors_for_the_internal_node,
						f.maxelm, dbeta);
	}

	SIMPLESPARSE sparseM; // разреженная матрица
	IMatrix sparseS; // разреженная матрица в формате IMatrix

	// выделение памяти и инициализация для 
	// простейшей разреженной матрицы.
	initsimplesparse(sparseM, f.maxelm + f.maxbound );
	initIMatrix(&sparseS, f.maxelm + f.maxbound);

	// Для внутренних узлов расчётной сетки:
	for (integer i=0; i<f.maxelm; i++) {
		addelmsimplesparse(sparseM, sl[i].ap, sl[i].iP, sl[i].iP, true);
		setValueIMatrix(&sparseS, sl[i].iP, sl[i].iP, sl[i].ap);
		rthdsd[sl[i].iP]=sl[i].b;

		const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

		     if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) {
				 addelmsimplesparse(sparseM, -sl[i].ae, sl[i].iP, sl[i].iE, true);
				 setValueIMatrix(&sparseS, sl[i].iP, sl[i].iE, -sl[i].ae);
			 }
			 if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) {
				 addelmsimplesparse(sparseM, -sl[i].an, sl[i].iP, sl[i].iN, true);
				 setValueIMatrix(&sparseS, sl[i].iP, sl[i].iN, -sl[i].an);
			 }
			 if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) {
				 addelmsimplesparse(sparseM, -sl[i].at, sl[i].iP, sl[i].iT, true);
                 setValueIMatrix(&sparseS, sl[i].iP, sl[i].iT, -sl[i].at);
			 }
			 if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) {
				 addelmsimplesparse(sparseM, -sl[i].as, sl[i].iP, sl[i].iS, true);
                 setValueIMatrix(&sparseS, sl[i].iP, sl[i].iS, -sl[i].as);
			 }
			 if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) {
				 addelmsimplesparse(sparseM, -sl[i].aw, sl[i].iP, sl[i].iW, true);
                 setValueIMatrix(&sparseS, sl[i].iP, sl[i].iW, -sl[i].aw);
			 }
			 if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) {
				 addelmsimplesparse(sparseM, -sl[i].ab, sl[i].iP, sl[i].iB, true);
				 setValueIMatrix(&sparseS, sl[i].iP, sl[i].iB, -sl[i].ab);
			 }
	}

	// Для граничных узлов расчётной сетки:
	for (integer i=0; i<f.maxbound; i++) {
		addelmsimplesparse(sparseM, slb[i].aw, slb[i].iW, slb[i].iW, true);
		setValueIMatrix(&sparseS, slb[i].iW, slb[i].iW, slb[i].aw);
		rthdsd[slb[i].iW]=slb[i].b;

		const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			 addelmsimplesparse(sparseM, -slb[i].ai, slb[i].iW, slb[i].iI, true);
             setValueIMatrix(&sparseS, slb[i].iW, slb[i].iI, -slb[i].ai);
		}
	}


	// решение СЛАУ
	for (integer i=0; i<(f.maxelm+f.maxbound); i++) potent[FIDISTW][i]=0.0;
    integer maxiter=2000; //120
	//doublereal *val;
    //integer *col_ind, *row_ptr;

	freeIMatrix(&sparseS); // освобождение оперативной памяти.
	//simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (f.maxelm+f.maxbound)); // преобразование матрицы из одного формата хранения в другой.

	//delete val; delete col_ind; delete row_ptr;
	//SOR3D(sl, slb, potent[FIDISTW], f.maxelm, f.maxbound, PAM); // проверка кода
	// getchar(); отладка
	// MICCG решатель.
	ICCG(FIDISTW,sparseM, rthdsd, potent[FIDISTW], f.maxelm + f.maxbound ,bprintmessage,true,2000); //->//
	// при использовании ICCG память из под sparseM освобождается автоматом

	// Вычисление градиентов.
	// на основе поля скорости удовлетворяющего уравнению неразрывности.
	for (integer i=0; i<f.maxelm; i++) {
		// градиенты потенциала FI для внутренних КО.
	    green_gauss_FI(i, potent, f.nvtx, f.pa,
	                   f.neighbors_for_the_internal_node, f.maxelm, false);
	}
	for (integer i=0; i<f.maxelm; i++) {
		// градиенты потенциала FI для граничных КО.
	    green_gauss_FI(i, potent, f.nvtx, f.pa,
	                   f.neighbors_for_the_internal_node, f.maxelm, true);
    }

	// Вычисление расстояния до ближайшей стенки.
	f.rdistWall=new doublereal[f.maxelm+f.maxbound];
	for (integer iP=0; iP<f.maxelm+f.maxbound; iP++) {
		// формула (223) из приложения 2 в книге Снегирёва.
		// А также смотри: Transport equation based wall distance calculation
		// www.cfd-online.com/Wiki/
		// В обоих источниках этой формулы содержится неточность.
		// Так как мы берём отрицательный источниковый член -1.0*dx*dy*dz
		// То это эквивалентно "холодильнику" и поэтому внутри расчётной области
		// значения потенциала FI будут отрицателдьными. Но подкоренное выражение
		// не может быть отрицательным,  поэтому нужно перед 2.0 поставить минус,
		// тем самым -2.0*FI будет положительной величиной и с подкоренным выражением 
		// будет всё впорядке.
		f.rdistWall[iP]=-sqrt(potent[GRADXFI][iP]*potent[GRADXFI][iP]+
			                  potent[GRADYFI][iP]*potent[GRADYFI][iP]+
							  potent[GRADZFI][iP]*potent[GRADZFI][iP])
					    +sqrt(potent[GRADXFI][iP]*potent[GRADXFI][iP]+
							  potent[GRADYFI][iP]*potent[GRADYFI][iP]+
							  potent[GRADZFI][iP]*potent[GRADZFI][iP]-
							  2.0*potent[FIDISTW][iP]);
	}
	// Вычисление толщины пограничного слоя в формуле Эскудиера
	f.rdistWallmax=-1e20; // очень маленькое значение
	for (integer i=0; i<(f.maxelm+f.maxbound); i++) {
		if (f.rdistWall[i]>f.rdistWallmax) f.rdistWallmax=f.rdistWall[i];
	}
	// Освобождение памяти.
	if (potent != NULL) {
		for (integer i = 0; i < 4; i++) {
			if (potent[i] != NULL) {
				delete[] potent[i];
			}
		}
		delete[] potent;
	}
	// освобождение памяти из под матрицы.
	if (sl != NULL) {
		delete[] sl;
	}
	if (slb != NULL) {
		delete[] slb;
	}
	if (rthdsd != NULL) {
		delete[] rthdsd;
	}
}

#endif
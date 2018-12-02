// pamendment3.c 
// текущий экспериментальный вариант 
// файла в котором составляется дискретный
// аналог для уравнения поправки давления.


// В случае если на всей границе FLUID стоит условие Неймана
// давление фиксируется в последней граничной точке.
// Вопрос о симметричности портрета СЛАУ, а также о самой
// симметричности СЛАУ требует рассмотрения в специальном модуле.
// Возможно с графической визуализацией матрицы. 
// Здесь сделано всё по масимуму для того чтобы матрица
// уравнения для поправки давления была симметричной и 
// положительно определённой. Возможно для целей SPD 
// в матрицу потребуется ввести несколько нулевых элементов.
// Вопрос остаётся открытым и требующим изучения свойств реальных матриц.



#ifndef PAMENDMENTv_0_07_C
#define PAMENDMENTv_0_07_C 1


// аппроксимация обобщённого уравнения конвекции-диффузии
// на совмещённой сетке
#include "my_elmatr_quad_f3D.c"
// вычисление шага по псевдовремени так как рекомендовал Гаврилов Андрей.
#include "pseudo_time.cpp"
// коррекция скорости на совмещённой сетке
// для внутренних и граничных контрольных объёмов.
#include "correct_velocity.cpp"

#define PRESSUREOUTLET 0 // давление заданное на выходной границе
#define BERNULLI 1 // если жидкость вытекает из расчётной области то давление полагается равным нулю, а если втекает то работает закон Бернулли.


// Заполнение граничных условий.
// В уравнении для поправки давления.
// Метод отложенной коррекции.
void my_elmatr_quad_PAm_bon(equation3D_bon** &slb, equation3D** sl, 
							integer inumber, integer maxelm, integer maxbound, 
							BOUND* &sosedb, integer** nvtx, bool bPfix,
							doublereal dbeta, TOCHKA* pa, doublereal** potent,
							doublereal** prop, doublereal** prop_b, doublereal* alpha,
							integer ls, integer lw, WALL* w, bool bDirichlet, 
							ALICE_PARTITION** sosedi, doublereal **diag_coef, doublereal RCh,
							bool &breversedflow) {
	
    // inumber - номер граничного узла, в порядке нумерации.
	bool bFReeStyle=false; // вообще отказ от фиксации поправки давления, поправка давления если ==true теперь везде плавающая.
	// Если bFReeStyle==false то поправка давления фиксируется либо на выходной границе либо в точке.

    
	

	// bDirichlet==true значит собираем только краевые условия Дирихле.

	if ( (inumber==(maxbound-1)) && bPfix && (!bFReeStyle)) {
        // Для корректного решения СЛАУ 
		// иногда возникает необходимость 
		// фиксировать давление в одной точке расчётной области.

		if (bDirichlet) {
           // поправка давления фиксированно равна нулю:
		   slb[PAM][inumber].aw=1.0;
		   slb[PAM][inumber].ai=0.0;
		   slb[PAM][inumber].b=0.0;
		   slb[PAM][inumber].iI=-1; // не присутствует в матрице
		   slb[PAM][inumber].iW=sosedb[inumber].iB; 

		   // Это условие Дирихле:
		   // только диагональный элемент
		   // не равен нулю.
		   slb[PAM][inumber].iW1=-1;
		   slb[PAM][inumber].iW2=-1;
		   slb[PAM][inumber].iW3=-1;
		   slb[PAM][inumber].iW4=-1;
		   //printf("pressure fix\n"); getchar(); // debug
		}

	}
	 else 
	{

		if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure)) {
            
			if (!bFReeStyle) {
				doublereal rCe=0.0; // значение давления на выходной границе
			    integer ioutflowcondition=BERNULLI; // PRESSUREOUTLET // BERNULLI И.Ю.Чумаков
			    doublereal kineticenergy=0.0; // кинетическая энергия потока: 0.5*rho*Vel!2.
			    doublereal rsign=1.0; // >0 скорость направлена из расчётной области наружу, если <0 то внутрь и используется соотношение Бернулли.
			
			    if (bDirichlet) {

			    	switch (ioutflowcondition) {
				     case PRESSUREOUTLET : // поправка давления фиксированна и равна нулю:
				                      // простейший способ - традиционно используется во многих расчётах.
				                      // Но только не при наличии рециркуляционных зон на границе.
		                              slb[PAM][inumber].aw=1.0;
		                              slb[PAM][inumber].ai=0.0;
		                              slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                              slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                              slb[PAM][inumber].iW=sosedb[inumber].iB;
					                  break;
			      	 case BERNULLI : // Рассматривается выходная граница потока и условие для поправки давления для неё.
					      // Если жидкость вытекает через выходную границу то используется постоянное давление rCe на выходной границе.
					      // Т.е. при наличии вытекания жидкости из расчётной области через выходную границу используется стандартное 
					      // условие PRESSUREOUTLET. Но если через границу проходят зоны возвратно циркуляционного течения и жидкость
                          // втекает в расчётную область через выходную границу то используется соотношение Бернулли для определения 
					      // давления на выходной границе. Соотношение Бернулли должно обеспечивать наилучшую сходимость решения в
					      // в ущерб его точности. Данное условие предложено в статью И.Ю.Чумакова Использование различных условий
					      // для давления на выходной границе при расчёте сложных внутренних течений несжимаемой жидкости на совмещённых
					      // сетках. 1997 год. с.48-54. Вестник молодых учёных. серия прикладная математика и механика.
					      
					      // определение знака величины rsign:
					      // внутренняя нормаль на выходной границе расчётной области.
	                      switch (sosedb[inumber].Norm) {
						  case ESIDE : if (potent[VX][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case WSIDE : if (potent[VX][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case NSIDE : if (potent[VY][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case SSIDE :if (potent[VY][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case TSIDE : if (potent[VZ][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case BSIDE : if (potent[VZ][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
			              } // end switch
					      // вычисление поправки давления на границе расчётной области:
					      if (rsign>=0.0) {
                             slb[PAM][inumber].aw=1.0;
		                     slb[PAM][inumber].ai=0.0;
		                     slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                     slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                     slb[PAM][inumber].iW=sosedb[inumber].iB; 
						  }
						  else {
							  kineticenergy=0.5*prop_b[RHO][sosedb[inumber].iB-maxelm];
							  kineticenergy*=(potent[VX][sosedb[inumber].iB]*potent[VX][sosedb[inumber].iB]+
								              potent[VY][sosedb[inumber].iB]*potent[VY][sosedb[inumber].iB]+
											  potent[VZ][sosedb[inumber].iB]*potent[VZ][sosedb[inumber].iB]);
							  slb[PAM][inumber].aw=1.0;
		                      slb[PAM][inumber].ai=0.0;
		                      slb[PAM][inumber].b=rCe-kineticenergy-potent[PRESS][sosedb[inumber].iB]; // Соотношение Бернулли
		                      slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                      slb[PAM][inumber].iW=sosedb[inumber].iB; 
						  }
					      break;
				     default : // поправка давления фиксированна и равна нулю:
				          // простейший способ - традиционно используется во многих расчётах.
				          // Но только не при наличии рециркуляционных зон на границе.
		                  slb[PAM][inumber].aw=1.0;
		                  slb[PAM][inumber].ai=0.0;
		                  slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                  slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                  slb[PAM][inumber].iW=sosedb[inumber].iB;
					      break;
				     } // end switch

					// Это условие Дирихле:
		            // только диагональный элемент
		            // не равен нулю.
		            slb[PAM][inumber].iW1=-1;
		            slb[PAM][inumber].iW2=-1;
		            slb[PAM][inumber].iW3=-1;
		            slb[PAM][inumber].iW4=-1;
			        //printf("pressure outlet\n"); // debug
			    
			     }	
						 

                 
			}
			else { // bFReeStyle==true
				if (!bDirichlet) {
					// однородное условие Неймана
		    doublereal dl, deltal, dS;
	        doublereal rhoi; // плотность на грани КО.
		    doublereal aUPi; // интерполированный диагональный коэффициент по скорости
	        doublereal fiplus; // учёт неравномерности сетки
			//doublereal FgRhieChow; // поправка Рхи-Чоу

			// внутренняя нормаль
	        switch (sosedb[inumber].Norm) {
		        case ESIDE :
			    
				    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[1][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VX]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VX][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VX]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					deltal-=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI] + fiplus*prop[RHO][sosedb[inumber].iII]);// проверено !
					aUPi=(sl[VX][sosedb[inumber].iI].ap*sl[VX][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VX][sosedb[inumber].iI].ap + fiplus*sl[VX][sosedb[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VX]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VX]);

					/*
					*  Итак лапласиан от поправки давления равен дивергенции скорости.
					*  Дивергенция скорости преобразует векторное поле скоростей
					*  в скалярное поле. Смысл дивергенции скорости для данного 
					*  контророльного объёма (КО) она показывает насколько расходятся
					*  входящий и исходящий поток в данном КО.
					*  Здесь идёт речь о граничном КО. Он является плоскостью (имеет нулевой объём). 
					*  Такой граничный  КО характеризуется нулевой расходимостью т.к. в нём по определению
					*  что вошло то и вышло, т.к. геометрически положение точки входа совпадает с геометрическим 
					*  положением точки выхода.
					*  Вывод: никакой дополнительной добавки в источниковый член не требуется!!! И не должно быть.
					*  Последний вывод проверен рядом тестов.
					*/

					//slb[PAM][inumber].b+=rhoi*dS*(potent[VX][sosedb[inumber].iB]-potent[VX][sosedb[inumber].iI]);
					//printf("iII=%d, iI=%d, iB=%d, maxelm=%d",sosedb[inumber].iII,sosedb[inumber].iI,sosedb[inumber].iB,maxelm);
					//getchar(); // debug внутренняя нормаль
					  //FgRhieChow=rFgRhieChow_internal_border(sosedb[inumber].iI, WSIDE, prop_b[RHO][sosedb[inumber].iB-maxelm], alpha[VX], nvtx, sosedi, maxelm, potent[PRESS], pa, diag_coef);
					//slb[PAM][inumber].b+=prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*0.5*(potent[VX][sosedb[inumber].iB]+potent[VX][sosedb[inumber].iI]);
					  //slb[PAM][inumber].b+=prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*potent[VX][sosedb[inumber].iB];
					  //slb[PAM][inumber].b+=RCh*FgRhieChow;
					  // данная добавка по-видимому приводит к расходимости.
					   //printf("RCh=%e, Fe=%e\n",RCh*FgRhieChow,slb[PAM][inumber].b+=prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*potent[VX][sosedb[inumber].iB]);
					   // getchar();
				    break;			
			
		        case NSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VY]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VY][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VY]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					deltal-=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI] + fiplus*prop[RHO][sosedb[inumber].iII]);// проверено !
					aUPi=(sl[VY][sosedb[inumber].iI].ap*sl[VY][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VY][sosedb[inumber].iI].ap + fiplus*sl[VY][sosedb[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VY]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VY]);

				    break;

			    case TSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VZ]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VZ][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VZ]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					deltal-=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI] + fiplus*prop[RHO][sosedb[inumber].iII]);// проверено !
					aUPi=(sl[VZ][sosedb[inumber].iI].ap*sl[VZ][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VZ][sosedb[inumber].iI].ap + fiplus*sl[VZ][sosedb[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VZ]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VZ]);
					
				    break;

			    case WSIDE :

                    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VX]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VX][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VX]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					deltal+=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]);// проверено !
					aUPi=(sl[VX][sosedb[inumber].iI].ap*sl[VX][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VX][sosedb[inumber].iI].ap + fiplus*sl[VX][sosedb[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VX]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VX]);
					
					break;

		        case SSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VY]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VY][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VY]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					deltal+=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI] + fiplus*prop[RHO][sosedb[inumber].iII]);// проверено !
					aUPi=(sl[VY][sosedb[inumber].iI].ap*sl[VY][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VY][sosedb[inumber].iI].ap + fiplus*sl[VY][sosedb[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VY]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VY]);
					
				   break;

		        case BSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VZ]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VZ][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VZ]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					deltal+=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VZ][sosedb[inumber].iI].ap*sl[VZ][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VZ][sosedb[inumber].iI].ap + fiplus*sl[VZ][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VZ]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VZ]);
					
				    break;
	        } // switch

			//*/
            integer j,l,xitem,k;
			// сортировка по возрастанию
			for (j=0; j<5; j++) {
				k=j; xitem=sosedb[inumber].iW[j];
				for (l=j+1; l<6; l++) {
					if (sosedb[inumber].iW[l] < xitem) {
						k=l; xitem=sosedb[inumber].iW[k];
					}
				}
                sosedb[inumber].iW[k]=sosedb[inumber].iW[j];
				sosedb[inumber].iW[j]=xitem;
			}

            j=0; l=0;
			while (sosedb[inumber].iW[j]==(-1)) j++;

			if (j<6) { slb[PAM][inumber].iW1=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW2=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW3=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW4=sosedb[inumber].iW[j++]; l++; } 

			switch (l) {
				case 0 : slb[PAM][inumber].iW1=-1;
		                 slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
		                 break;
				case 1 : slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 2 : slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 3 : slb[PAM][inumber].iW4=-1;
						 break;
			}

			
				}
			}

	    }
		 else if ((!bDirichlet) && !((inumber == (maxbound - 1)) && bPfix) && !((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure )))
	    {  
			// Не условие Дирихле, Не фиксация давления в точке, Не выходная граница
			//printf("neiman!\n"); getchar(); // debug

			/*
            // поправка давления фиксированно равна нулю:
		    slb[PAM][inumber].aw=1.0;
		    slb[PAM][inumber].ai=0.0;
		    slb[PAM][inumber].b=0.0;
			slb[PAM][inumber].iI=-1;//sosedb[inumber].iI; //  присутствует в матрице
		    slb[PAM][inumber].iW=sosedb[inumber].iB;
			*/

            // Здесь предполагается, что узлы iI и iII внутренние, иначе поведение
		    // программы будет неправильным. НЕ МЕНЕЕ 2-ух ненулевых КО в зазоре между рёбрами.
			//printf("neiman pressure...\n"); // debug
			//getchar();

			///*
		    // однородное условие Неймана
		    doublereal dl, deltal, dS;
	        doublereal rhoi; // плотность на грани КО.
		    doublereal aUPi; // интерполированный диагональный коэффициент по скорости
	        doublereal fiplus; // учёт неравномерности сетки
			//doublereal FgRhieChow; // поправка Рхи-Чоу

			// внутренняя нормаль
	        switch (sosedb[inumber].Norm) {
		        case ESIDE :
			    
				    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[1][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VX]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VX][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VX]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					deltal-=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VX][sosedb[inumber].iI].ap*sl[VX][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VX][sosedb[inumber].iI].ap + fiplus*sl[VX][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VX]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VX]);

					/*
					*  Итак лапласиан от поправки давления равен дивергенции скорости.
					*  Дивергенция скорости преобразует векторное поле скоростей
					*  в скалярное поле. Смысл дивергенции скорости для данного 
					*  контророльного объёма (КО) она показывает насколько расходятся
					*  входящий и исходящий поток в данном КО.
					*  Здесь идёт речь о граничном КО. Он является плоскостью (имеет нулевой объём). 
					*  Такой граничный  КО характеризуется нулевой расходимостью т.к. в нём по определению
					*  что вошло то и вышло, т.к. геометрически положение точки входа совпадает с геометрическим 
					*  положением точки выхода.
					*  Вывод: никакой дополнительной добавки в источниковый член не требуется!!! И не должно быть.
					*  Последний вывод проверен рядом тестов.
					*/

					//slb[PAM][inumber].b+=rhoi*dS*(potent[VX][sosedb[inumber].iB]-potent[VX][sosedb[inumber].iI]);
					//printf("iII=%d, iI=%d, iB=%d, maxelm=%d",sosedb[inumber].iII,sosedb[inumber].iI,sosedb[inumber].iB,maxelm);
					//getchar(); // debug внутренняя нормаль
					  //FgRhieChow=rFgRhieChow_internal_border(sosedb[inumber].iI, WSIDE, prop_b[RHO][sosedb[inumber].iB-maxelm], alpha[VX], nvtx, sosedi, maxelm, potent[PRESS], pa, diag_coef);
					//slb[PAM][inumber].b+=prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*0.5*(potent[VX][sosedb[inumber].iB]+potent[VX][sosedb[inumber].iI]);
					  //slb[PAM][inumber].b+=prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*potent[VX][sosedb[inumber].iB];
					  //slb[PAM][inumber].b+=RCh*FgRhieChow;
					  // данная добавка по-видимому приводит к расходимости.
					   //printf("RCh=%e, Fe=%e\n",RCh*FgRhieChow,slb[PAM][inumber].b+=prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*potent[VX][sosedb[inumber].iB]);
					   // getchar();
				    break;			
			
		        case NSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VY]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VY][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VY]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					deltal-=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VY][sosedb[inumber].iI].ap*sl[VY][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VY][sosedb[inumber].iI].ap + fiplus*sl[VY][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VY]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VY]);

				    break;

			    case TSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VZ]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VZ][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VZ]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					deltal-=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VZ][sosedb[inumber].iI].ap*sl[VZ][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VZ][sosedb[inumber].iI].ap + fiplus*sl[VZ][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VZ]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VZ]);
					
				    break;

			    case WSIDE :

                    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VX]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VX][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VX]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					deltal+=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI] + fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VX][sosedb[inumber].iI].ap*sl[VX][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VX][sosedb[inumber].iI].ap + fiplus*sl[VX][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VX]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VX]);
					
					break;

		        case SSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VY]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VY][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VY]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					deltal+=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VY][sosedb[inumber].iI].ap*sl[VY][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VY][sosedb[inumber].iI].ap + fiplus*sl[VY][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VY]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VY]);
					
				   break;

		        case BSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VZ]*prop_b[RHO][sosedb[inumber].iB-maxelm]*dS*dS/slb[VZ][sosedb[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VZ]);
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					deltal+=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][sosedb[inumber].iI]*prop[RHO][sosedb[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][sosedb[inumber].iI]+ fiplus*prop[RHO][sosedb[inumber].iII]); // проверено !
					aUPi=(sl[VZ][sosedb[inumber].iI].ap*sl[VZ][sosedb[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VZ][sosedb[inumber].iI].ap + fiplus*sl[VZ][sosedb[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VZ]*rhoi*dS*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/aUPi;
					if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VZ]);
					
				    break;
	        } // switch

			//*/
            integer j,l,xitem,k;
			// сортировка по возрастанию
			for (j=0; j<5; j++) {
				k=j; xitem=sosedb[inumber].iW[j];
				for (l=j+1; l<6; l++) {
					if (sosedb[inumber].iW[l] < xitem) {
						k=l; xitem=sosedb[inumber].iW[k];
					}
				}
                sosedb[inumber].iW[k]=sosedb[inumber].iW[j];
				sosedb[inumber].iW[j]=xitem;
			}

            j=0; l=0;
			while (sosedb[inumber].iW[j]==(-1)) j++;

			if (j<6) { slb[PAM][inumber].iW1=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW2=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW3=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW4=sosedb[inumber].iW[j++]; l++; } 

			switch (l) {
				case 0 : slb[PAM][inumber].iW1=-1;
		                 slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
		                 break;
				case 1 : slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 2 : slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 3 : slb[PAM][inumber].iW4=-1;
						 break;
			}

			// Наложения исключения в случае совпадения узла iI с диагональным 
			// элементом узла для котторого стоит условие Дирихле не требуется. 
			// т.к. узлы iI строго внутренние для которых iI < maxelm, а диагональный 
			// элемент с условием Дирихле стоит в позициях >= maxelm поэтому они не
			// пересекаются.
	    }
	}
	
	
} // my_elmatr_quad_PAm_bon

// Заполнение граничных условий.
// В уравнении для поправки давления
// На основе сглаженного псевдовремени.
// Метод отложенной коррекции.
// реализовано 22 июня 2012 года.
void my_elmatr_quad_PAm_bon2(equation3D_bon** &slb, equation3D** sl, 
							integer inumber, integer maxelm, integer maxbound, 
							BOUND* &sosedb, integer** nvtx, bool bPfix,
							doublereal dbeta, TOCHKA* pa, doublereal** potent,
							doublereal** prop, doublereal** prop_b, doublereal* alpha,
							integer ls, integer lw, WALL* w, bool bDirichlet, 
							ALICE_PARTITION** sosedi, doublereal RCh,
							bool &breversedflow, doublereal* tau) {
	
    // inumber - номер граничного узла, в порядке нумерации.
	bool bFReeStyle=false; // вообще отказ от фиксации поправки давления, поправка давления если ==true теперь везде плавающая.
	// Если bFReeStyle==false то поправка давления фиксируется либо на выходной границе либо в точке.

   
	// bDirichlet==true значит собираем только краевые условия Дирихле.

	if ( (inumber==(maxbound-1)) && bPfix && (!bFReeStyle)) {
        // Для корректного решения СЛАУ 
		// иногда возникает необходимость 
		// фиксировать давление в одной точке расчётной области.

		if (bDirichlet) {
           // поправка давления фиксированно равна нулю:
		   slb[PAM][inumber].aw=1.0;
		   slb[PAM][inumber].ai=0.0;
		   slb[PAM][inumber].b=0.0;
		   slb[PAM][inumber].iI=-1; // не присутствует в матрице
		   slb[PAM][inumber].iW=sosedb[inumber].iB; 

		   // Это условие Дирихле:
		   // только диагональный элемент
		   // не равен нулю.
		   slb[PAM][inumber].iW1=-1;
		   slb[PAM][inumber].iW2=-1;
		   slb[PAM][inumber].iW3=-1;
		   slb[PAM][inumber].iW4=-1;
		   //printf("pressure fix\n"); getchar(); // debug
		}

	}
	 else 
	{

		if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure )) {
            
			if (!bFReeStyle) {
				doublereal rCe=0.0; // значение давления на выходной границе
			    integer ioutflowcondition=PRESSUREOUTLET; // PRESSUREOUTLET // BERNULLI И.Ю.Чумаков
			    doublereal kineticenergy=0.0; // кинетическая энергия потока: 0.5*rho*Vel!2.
			    doublereal rsign=1.0; // >0 скорость направлена из расчётной области наружу, если <0 то внутрь и используется соотношение Бернулли.
			
			    if (bDirichlet) {

			    	switch (ioutflowcondition) {
				     case PRESSUREOUTLET : // поправка давления фиксированна и равна нулю:
				                      // простейший способ - традиционно используется во многих расчётах.
				                      // Но только не при наличии рециркуляционных зон на границе.
		                              slb[PAM][inumber].aw=1.0;
		                              slb[PAM][inumber].ai=0.0;
		                              slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                              slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                              slb[PAM][inumber].iW=sosedb[inumber].iB;
					                  break;
			      	 case BERNULLI : // Рассматривается выходная граница потока и условие для поправки давления для неё.
					      // Если жидкость вытекает через выходную границу то используется постоянное давление rCe на выходной границе.
					      // Т.е. при наличии вытекания жидкости из расчётной области через выходную границу используется стандартное 
					      // условие PRESSUREOUTLET. Но если через границу проходят зоны возвратно циркуляционного течения и жидкость
                          // втекает в расчётную область через выходную границу то используется соотношение Бернулли для определения 
					      // давления на выходной границе. Соотношение Бернулли должно обеспечивать наилучшую сходимость решения в
					      // в ущерб его точности. Данное условие предложено в статью И.Ю.Чумакова Использование различных условий
					      // для давления на выходной границе при расчёте сложных внутренних течений несжимаемой жидкости на совмещённых
					      // сетках. 1997 год. с.48-54. Вестник молодых учёных. серия прикладная математика и механика.
					      
					      // определение знака величины rsign:
					      // внутренняя нормаль на выходной границе расчётной области.
	                      switch (sosedb[inumber].Norm) {
						  case ESIDE : if (potent[VX][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case WSIDE : if (potent[VX][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case NSIDE : if (potent[VY][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case SSIDE :if (potent[VY][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case TSIDE : if (potent[VZ][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case BSIDE : if (potent[VZ][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
			              } // end switch
					      // вычисление поправки давления на границе расчётной области:
					      if (rsign>=0.0) {
                             slb[PAM][inumber].aw=1.0;
		                     slb[PAM][inumber].ai=0.0;
		                     slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                     slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                     slb[PAM][inumber].iW=sosedb[inumber].iB; 
						  }
						  else {
							  kineticenergy=0.5*prop_b[RHO][sosedb[inumber].iB-maxelm];
							  kineticenergy*=(potent[VX][sosedb[inumber].iB]*potent[VX][sosedb[inumber].iB]+
								              potent[VY][sosedb[inumber].iB]*potent[VY][sosedb[inumber].iB]+
											  potent[VZ][sosedb[inumber].iB]*potent[VZ][sosedb[inumber].iB]);
							  slb[PAM][inumber].aw=1.0;
		                      slb[PAM][inumber].ai=0.0;
		                      slb[PAM][inumber].b=rCe-kineticenergy-potent[PRESS][sosedb[inumber].iB]; // Соотношение Бернулли
		                      slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                      slb[PAM][inumber].iW=sosedb[inumber].iB; 
						  }
					      break;
				     default : // поправка давления фиксированна и равна нулю:
				          // простейший способ - традиционно используется во многих расчётах.
				          // Но только не при наличии рециркуляционных зон на границе.
		                  slb[PAM][inumber].aw=1.0;
		                  slb[PAM][inumber].ai=0.0;
		                  slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                  slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                  slb[PAM][inumber].iW=sosedb[inumber].iB;
					      break;
				     } // end switch

					// Это условие Дирихле:
		            // только диагональный элемент
		            // не равен нулю.
		            slb[PAM][inumber].iW1=-1;
		            slb[PAM][inumber].iW2=-1;
		            slb[PAM][inumber].iW3=-1;
		            slb[PAM][inumber].iW4=-1;
			        //printf("pressure outlet\n"); // debug
			    
			     }	
						 

                 
			}
			else { // bFReeStyle==true
				if (!bDirichlet) {
					// однородное условие Неймана
		            doublereal dl, deltal, dS;
	                doublereal taui; // псевдовремя на грани КО.
	                doublereal fiplus; // учёт неравномерности сетки

        			// внутренняя нормаль
	                switch (sosedb[inumber].Norm) {
		              case ESIDE :
			    
				          dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					      dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[1][sosedb[inumber].iI]-1].y; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					      deltal=0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					      deltal-=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;

					      /*
					      *  Итак лапласиан от поправки давления равен дивергенции скорости.
					      *  Дивергенция скорости преобразует векторное поле скоростей
					      *  в скалярное поле. Смысл дивергенции скорости для данного 
					      *  контророльного объёма (КО) она показывает насколько расходятся
					      *  входящий и исходящий поток в данном КО.
					      *  Здесь идёт речь о граничном КО. Он является плоскостью (имеет нулевой объём). 
					      *  Такой граничный  КО характеризуется нулевой расходимостью т.к. в нём по определению
					      *  что вошло то и вышло, т.к. геометрически положение точки входа совпадает с геометрическим 
					      *  положением точки выхода.
				 	      *  Вывод: никакой дополнительной добавки в источниковый член не требуется!!! И не должно быть.
					      *  Последний вывод проверен рядом тестов.
					      */
				          break;			
			
		              case NSIDE :

                          dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					      deltal=0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					      deltal-=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+ fiplus*tau[sosedb[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					

				          break;

			          case TSIDE : 

                          dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					      deltal=0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					      deltal-=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
				          break;

			          case WSIDE :

                          dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					      dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					      deltal=-0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					      deltal+=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
					      break;

		              case SSIDE :

                          dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					      deltal=-0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					      deltal+=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
					
				          break;

		              case BSIDE : 

                          dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					      deltal=-0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					      deltal+=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					      
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					      
					
				          break;
	                } // switch

			        //*/
                    integer j,l,xitem,k;
			        // сортировка по возрастанию
			        for (j=0; j<5; j++) {
				        k=j; xitem=sosedb[inumber].iW[j];
				        for (l=j+1; l<6; l++) {
					        if (sosedb[inumber].iW[l] < xitem) {
						       k=l; xitem=sosedb[inumber].iW[k];
					        }
				        }
                        sosedb[inumber].iW[k]=sosedb[inumber].iW[j];
				        sosedb[inumber].iW[j]=xitem;
			        }

                    j=0; l=0;
			        while (sosedb[inumber].iW[j]==(-1)) j++;

			        if (j<6) { slb[PAM][inumber].iW1=sosedb[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW2=sosedb[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW3=sosedb[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW4=sosedb[inumber].iW[j++]; l++; } 

			        switch (l) {
				       case 0 : slb[PAM][inumber].iW1=-1;
		                        slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
		                        break;
				       case 1 : slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 2 : slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 3 : slb[PAM][inumber].iW4=-1;
						        break;
			        }

			
				}
			}

	    }
	     else if ((!bDirichlet) && !((inumber==(maxbound-1)) && bPfix) && !((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure) )
	    {  
			// Не условие Дирихле, Не фиксация давления в точке, Не выходная граница
			//printf("neiman!\n"); getchar(); // debug


			/*
            // поправка давления фиксированно равна нулю:
		    slb[PAM][inumber].aw=1.0;
		    slb[PAM][inumber].ai=0.0;
		    slb[PAM][inumber].b=0.0;
			slb[PAM][inumber].iI=-1;//sosedb[inumber].iI; //  присутствует в матрице
		    slb[PAM][inumber].iW=sosedb[inumber].iB;
			*/

            // Здесь предполагается, что узлы iI и iII внутренние, иначе поведение
		    // программы будет неправильным. НЕ МЕНЕЕ 2-ух ненулевых КО в зазоре между рёбрами.
			//printf("neiman pressure...\n"); // debug
			//getchar();

			///*
		    // однородное условие Неймана
		    doublereal dl, deltal, dS;
	        doublereal taui; // псевдовремя на грани КО.
	        doublereal fiplus; // учёт неравномерности сетки
			//doublereal FgRhieChow; // поправка Рхи-Чоу

			// внутренняя нормаль
	        switch (sosedb[inumber].Norm) {
		        case ESIDE :
			    
				    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[1][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					deltal-=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;

					/*
					*  Итак лапласиан от поправки давления равен дивергенции скорости.
					*  Дивергенция скорости преобразует векторное поле скоростей
					*  в скалярное поле. Смысл дивергенции скорости для данного 
					*  контророльного объёма (КО) она показывает насколько расходятся
					*  входящий и исходящий поток в данном КО.
					*  Здесь идёт речь о граничном КО. Он является плоскостью (имеет нулевой объём). 
					*  Такой граничный  КО характеризуется нулевой расходимостью т.к. в нём по определению
					*  что вошло то и вышло, т.к. геометрически положение точки входа совпадает с геометрическим 
					*  положением точки выхода.
					*  Вывод: никакой дополнительной добавки в источниковый член не требуется!!! И не должно быть.
					*  Последний вывод проверен рядом тестов.
					*/
				    break;			
			
		        case NSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					deltal-=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;

				    break;

			    case TSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					deltal-=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
								        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
				    break;

			    case WSIDE :

                    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[1][sosedb[inumber].iII]-1].x+pa[nvtx[0][sosedb[inumber].iII]-1].x);
					deltal+=0.5*(pa[nvtx[1][sosedb[inumber].iI]-1].x+pa[nvtx[0][sosedb[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
					break;

		        case SSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[2][sosedb[inumber].iII]-1].y+pa[nvtx[0][sosedb[inumber].iII]-1].y);
					deltal+=0.5*(pa[nvtx[2][sosedb[inumber].iI]-1].y+pa[nvtx[0][sosedb[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
				   break;

		        case BSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[4][sosedb[inumber].iII]-1].z+pa[nvtx[0][sosedb[inumber].iII]-1].z);
					deltal+=0.5*(pa[nvtx[4][sosedb[inumber].iI]-1].z+pa[nvtx[0][sosedb[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					taui=(tau[sosedb[inumber].iI]*tau[sosedb[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[sosedb[inumber].iI]+fiplus*tau[sosedb[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI]-potent[PAM][sosedb[inumber].iII])/deltal;
					
				    break;
	        } // switch

			//*/
            integer j,l,xitem,k;
			// сортировка по возрастанию
			for (j=0; j<5; j++) {
				k=j; xitem=sosedb[inumber].iW[j];
				for (l=j+1; l<6; l++) {
					if (sosedb[inumber].iW[l] < xitem) {
						k=l; xitem=sosedb[inumber].iW[k];
					}
				}
                sosedb[inumber].iW[k]=sosedb[inumber].iW[j];
				sosedb[inumber].iW[j]=xitem;
			}

            j=0; l=0;
			while (sosedb[inumber].iW[j]==(-1)) j++;

			if (j<6) { slb[PAM][inumber].iW1=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW2=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW3=sosedb[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW4=sosedb[inumber].iW[j++]; l++; } 

			switch (l) {
				case 0 : slb[PAM][inumber].iW1=-1;
		                 slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
		                 break;
				case 1 : slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 2 : slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 3 : slb[PAM][inumber].iW4=-1;
						 break;
			}

			// Наложения исключения в случае совпадения узла iI с диагональным 
			// элементом узла для котторого стоит условие Дирихле не требуется. 
			// т.к. узлы iI строго внутренние для которых iI < maxelm, а диагональный 
			// элемент с условием Дирихле стоит в позициях >= maxelm поэтому они не
			// пересекаются.
	    }
	}
	
	
} // my_elmatr_quad_PAm_bon2

// Заполнение граничных условий.
// В уравнении для поправки давления
// На основе сглаженного псевдовремени.
// Метод отложенной коррекции.
// реализовано 22 июня 2012 года.
// добавление условия opening не потребовало модификации. 27.07.2016.
void my_elmatr_quad_PAm_bon3(equation3D_bon** &slb, equation3D** sl, 
							integer inumber, integer maxelm, integer maxbound, 
							BOUND* &sosedb, integer** nvtx, bool bPfix,
							doublereal dbeta, TOCHKA* pa, doublereal** potent,
							doublereal** prop, doublereal** prop_b, doublereal* alpha,
							integer ls, integer lw, WALL* w, bool bDirichlet, 
							ALICE_PARTITION** sosedi, doublereal RCh,
							bool &breversedflow, doublereal** tau) {
	
	bool bcheck = false;

    // inumber - номер граничного узла, в порядке нумерации.
	bool bFReeStyle=false; // вообще отказ от фиксации поправки давления, поправка давления если ==true теперь везде плавающая.
	// Если bFReeStyle==false то поправка давления фиксируется либо на выходной границе либо в точке.
	///bFReeStyle = true; // 5.05.2017 Так не работает вообще.

	bool b_prosto = true;
	// bDirichlet==true значит собираем только краевые условия Дирихле.

	/*
	if ( (inumber==(maxbound-1)) && bPfix && (!bFReeStyle)) {
        // Для корректного решения СЛАУ 
		// иногда возникает необходимость 
		// фиксировать давление в одной точке расчётной области.
		/*
		printf("pressure fix Ok...");
		getchar();
		*//*

		if (bDirichlet) {
           // поправка давления фиксированно равна нулю:
		   slb[PAM][inumber].aw=1.0;
		   slb[PAM][inumber].ai=0.0;
		   slb[PAM][inumber].b=0.0;
		   slb[PAM][inumber].iI=-1; // не присутствует в матрице
		   slb[PAM][inumber].iW=sosedb[inumber].iB; 

		   // Это условие Дирихле:
		   // только диагональный элемент
		   // не равен нулю.
		   slb[PAM][inumber].iW1=-1;
		   slb[PAM][inumber].iW2=-1;
		   slb[PAM][inumber].iW3=-1;
		   slb[PAM][inumber].iW4=-1;
		   //printf("pressure fix\n"); getchar(); // debug
		}

	}
	 else 
	{*/

	if ((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB<(ls + lw)) && (w[sosedb[inumber].MCB - ls].bpressure || w[sosedb[inumber].MCB - ls].bopening)) {
            
			if (!bFReeStyle) {
				


				doublereal rCe=w[sosedb[inumber].MCB-ls].P; // значение давления на выходной границе
				//printf("rCe=%e %d \n",rCe,sosedb[inumber].MCB-ls);
				//getchar();

				integer ioutflowcondition =  BERNULLI; // PRESSUREOUTLET // BERNULLI И.Ю.Чумаков
			    doublereal kineticenergy=0.0; // кинетическая энергия потока: 0.5*rho*Vel!2.
			    doublereal rsign=1.0; // >0 скорость направлена из расчётной области наружу, если <0 то внутрь и используется соотношение Бернулли.
			
			    if (bDirichlet) {

					bcheck = true;

			    	switch (ioutflowcondition) {
				     case PRESSUREOUTLET : // поправка давления фиксированна и равна нулю:
				                      // простейший способ - традиционно используется во многих расчётах.
				                      // Но только не при наличии рециркуляционных зон на границе.
		                              slb[PAM][inumber].aw=1.0;
		                              slb[PAM][inumber].ai=0.0;
		                              slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                              slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                              slb[PAM][inumber].iW=sosedb[inumber].iB;
									 // printf("%d\n",slb[PAM][inumber].iW); getchar();// больше либо равно maxelm
					                  break;
			      	 case BERNULLI : // Рассматривается выходная граница потока и условие для поправки давления для неё.
					      // Если жидкость вытекает через выходную границу то используется постоянное давление rCe на выходной границе.
					      // Т.е. при наличии вытекания жидкости из расчётной области через выходную границу используется стандартное 
					      // условие PRESSUREOUTLET. Но если через границу проходят зоны возвратно циркуляционного течения и жидкость
                          // втекает в расчётную область через выходную границу то используется соотношение Бернулли для определения 
					      // давления на выходной границе. Соотношение Бернулли должно обеспечивать наилучшую сходимость решения в
					      // в ущерб его точности. Данное условие предложено в статью И.Ю.Чумакова Использование различных условий
					      // для давления на выходной границе при расчёте сложных внутренних течений несжимаемой жидкости на совмещённых
					      // сетках. 1997 год. с.48-54. Вестник молодых учёных. серия прикладная математика и механика.
					      
					      // определение знака величины rsign:
					      // внутренняя нормаль на выходной границе расчётной области.
	                      switch (sosedb[inumber].Norm) {
						  case ESIDE : if (potent[VX][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case WSIDE : if (potent[VX][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case NSIDE : if (potent[VY][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case SSIDE :if (potent[VY][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case TSIDE : if (potent[VZ][sosedb[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case BSIDE : if (potent[VZ][sosedb[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
			              } // end switch
					      // вычисление поправки давления на границе расчётной области:
					      if (rsign>=0.0) {
							  // Жидкость покидает расчётную область.

                             slb[PAM][inumber].aw=1.0;
		                     slb[PAM][inumber].ai=0.0;
		                     slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                     slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                     slb[PAM][inumber].iW=sosedb[inumber].iB; 
						  }
						  else {
							  // Жидкость проникает внутрь расчётной области.

							  kineticenergy=0.5*prop_b[RHO][sosedb[inumber].iB-maxelm];
							  // VComponentCOR вместо VComponent не повлияло.
							  kineticenergy*=(potent[VXCOR][sosedb[inumber].iB]*potent[VXCOR][sosedb[inumber].iB]+
								              potent[VYCOR][sosedb[inumber].iB]*potent[VYCOR][sosedb[inumber].iB]+
											  potent[VZCOR][sosedb[inumber].iB]*potent[VZCOR][sosedb[inumber].iB]);
							  slb[PAM][inumber].aw=1.0;
		                      slb[PAM][inumber].ai=0.0;
		                      slb[PAM][inumber].b=rCe-kineticenergy-potent[PRESS][sosedb[inumber].iB]; // Соотношение Бернулли
		                      slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                      slb[PAM][inumber].iW=sosedb[inumber].iB; 
						  }
					      break;
				     default : // поправка давления фиксированна и равна нулю:
				          // простейший способ - традиционно используется во многих расчётах.
				          // Но только не при наличии рециркуляционных зон на границе.
		                  slb[PAM][inumber].aw=1.0;
		                  slb[PAM][inumber].ai=0.0;
		                  slb[PAM][inumber].b=rCe-potent[PRESS][sosedb[inumber].iB];
		                  slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                  slb[PAM][inumber].iW=sosedb[inumber].iB;
					      break;
				     } // end switch

					// Это условие Дирихле:
		            // только диагональный элемент
		            // не равен нулю.
		            slb[PAM][inumber].iW1=-1;
		            slb[PAM][inumber].iW2=-1;
		            slb[PAM][inumber].iW3=-1;
		            slb[PAM][inumber].iW4=-1;
			        //printf("pressure outlet\n"); // debug
			    
			     }	
				 else if (1) {
					 // !bDirihlet (not bDirichlet)
					 if (w[sosedb[inumber].MCB - ls].bopening) {
						 // ВНИМАНИЕ!!! только в случае opening границы расчётной области.

						 // 5.05.2017. Обязательно оставить включённым ==1.
						 // Включение данного участка кода существенно улучшает сходимость на 
						 // поздних итерациях.
						 // Однородное условие Неймана для поправки давления в ячейках границы,
						 // через которые теплоноситель (воздух) покидает расчётную область.
						 // полезная модификация в случае Opening границы.
						 switch (ioutflowcondition) {
						 case BERNULLI:
							 // определение знака величины rsign:
							 // внутренняя нормаль на выходной границе расчётной области.
							 switch (sosedb[inumber].Norm) {
							 case ESIDE:  if (potent[VX][sosedb[inumber].iB] > 0.0) {
								        // жидкость втекает внутрь расчётной области через выходную границу.
								        rsign = -1.0;
										// Повторно не увеличиваем счетчик КО с рециркуляцией.
								        // breversedflow = true;
							         }
									 else rsign = 1.0;
									 break;
							 case WSIDE: if (potent[VX][sosedb[inumber].iB] < 0.0) {
								          // жидкость втекает внутрь расчётной области через выходную границу.
								          rsign = -1.0;
										  // Повторно не увеличиваем счетчик КО с рециркуляцией.
								          //breversedflow = true;
							         }
									 else rsign = 1.0;
									 break;
							 case NSIDE: if (potent[VY][sosedb[inumber].iB] > 0.0) {
								             // жидкость втекает внутрь расчётной области через выходную границу.
								             rsign = -1.0;
											 // Повторно не увеличиваем счетчик КО с рециркуляцией.
								             //breversedflow = true;
							             }
										 else rsign = 1.0;
										 break;
							 case SSIDE:if (potent[VY][sosedb[inumber].iB] < 0.0) {
								            // жидкость втекает внутрь расчётной области через выходную границу.
								            rsign = -1.0;
											// Повторно не увеличиваем счетчик КО с рециркуляцией.
								            //breversedflow = true;
							            }
										else rsign = 1.0;
										break;
							 case TSIDE: if (potent[VZ][sosedb[inumber].iB] > 0.0) {
								             // жидкость втекает внутрь расчётной области через выходную границу.
								             rsign = -1.0;
											 // Повторно не увеличиваем счетчик КО с рециркуляцией.				             
											 //breversedflow = true;
							             }
										 else rsign = 1.0;
										 break;
							 case BSIDE: if (potent[VZ][sosedb[inumber].iB] < 0.0) {
								             // жидкость втекает внутрь расчётной области через выходную границу.
								             rsign = -1.0;
											 // Повторно не увеличиваем счетчик КО с рециркуляцией.
											 //breversedflow = true;
							             }
										 else rsign = 1.0;
										 break;
							 } // end switch
							   // вычисление поправки давления на границе расчётной области:
							 if (rsign >= 0.0) {
								 // Жидкость покидает расчётную область.
								 // Однородные условия Неймана для поправки давления 
								 // при вытекании воздуха из расчётной области.
								 slb[PAM][inumber].aw = 1.0;
								 slb[PAM][inumber].ai = 1.0;
								 slb[PAM][inumber].b = 0.0;
								 slb[PAM][inumber].iI = sosedb[inumber].iI;
								 slb[PAM][inumber].iW = sosedb[inumber].iB;
							 }
							 break;
						 }
					 }
				 }
						 

                 
			}
			else { // bFReeStyle==true
				if (!bDirichlet) {
					// однородное условие Неймана
		            doublereal dl, deltal, dS;
	                doublereal taui; // псевдовремя на грани КО.
	                doublereal fiplus; // учёт неравномерности сетки

					bcheck = true;

        			// внутренняя нормаль
	                switch (sosedb[inumber].Norm) {
		              case ESIDE :
			    
				          dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					      dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[1][sosedb[inumber].iI]-1].y; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[VX][sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
						  if (!b_prosto) {


							  deltal = 0.5*(pa[nvtx[1][sosedb[inumber].iII] - 1].x + pa[nvtx[0][sosedb[inumber].iII] - 1].x);
							  deltal -= 0.5*(pa[nvtx[1][sosedb[inumber].iI] - 1].x + pa[nvtx[0][sosedb[inumber].iI] - 1].x);
							  fiplus = 0.5*dl / deltal;

							  taui = (tau[VX][sosedb[inumber].iI] * tau[VX][sosedb[inumber].iII]);
							  taui = taui / ((1.0 - fiplus)*tau[VX][sosedb[inumber].iI] + fiplus*tau[VX][sosedb[inumber].iII]); // проверено !


							  // правая часть:  
							  slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
						  }
						  else slb[PAM][inumber].b = 0.0;
					      /*
					      *  Итак лапласиан от поправки давления равен дивергенции скорости.
					      *  Дивергенция скорости преобразует векторное поле скоростей
					      *  в скалярное поле. Смысл дивергенции скорости для данного 
					      *  контророльного объёма (КО) она показывает насколько расходятся
					      *  входящий и исходящий поток в данном КО.
					      *  Здесь идёт речь о граничном КО. Он является плоскостью (имеет нулевой объём). 
					      *  Такой граничный  КО характеризуется нулевой расходимостью т.к. в нём по определению
					      *  что вошло то и вышло, т.к. геометрически положение точки входа совпадает с геометрическим 
					      *  положением точки выхода.
				 	      *  Вывод: никакой дополнительной добавки в источниковый член не требуется!!! И не должно быть.
					      *  Последний вывод проверен рядом тестов.
					      */
				          break;			
			
		              case NSIDE :

                          dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[VY][sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
						  if (!b_prosto) {

							  deltal = 0.5*(pa[nvtx[2][sosedb[inumber].iII] - 1].y + pa[nvtx[0][sosedb[inumber].iII] - 1].y);
							  deltal -= 0.5*(pa[nvtx[2][sosedb[inumber].iI] - 1].y + pa[nvtx[0][sosedb[inumber].iI] - 1].y);
							  fiplus = 0.5*dl / deltal;

							  taui = (tau[VY][sosedb[inumber].iI] * tau[VY][sosedb[inumber].iII]);
							  taui = taui / ((1.0 - fiplus)*tau[VY][sosedb[inumber].iI] + fiplus*tau[VY][sosedb[inumber].iII]); // проверено !


							  // правая часть:  
							  slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
						  }
						  else slb[PAM][inumber].b = 0.0;

				          break;

			          case TSIDE : 

                          dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[VZ][sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
						  if (!b_prosto) {

							  deltal = 0.5*(pa[nvtx[4][sosedb[inumber].iII] - 1].z + pa[nvtx[0][sosedb[inumber].iII] - 1].z);
							  deltal -= 0.5*(pa[nvtx[4][sosedb[inumber].iI] - 1].z + pa[nvtx[0][sosedb[inumber].iI] - 1].z);
							  fiplus = 0.5*dl / deltal;

							  taui = (tau[VZ][sosedb[inumber].iI] * tau[VZ][sosedb[inumber].iII]);
							  taui = taui / ((1.0 - fiplus)*tau[VZ][sosedb[inumber].iI] + fiplus*tau[VZ][sosedb[inumber].iII]);// проверено !


							  // правая часть:  
							  slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
						  }
						  else slb[PAM][inumber].b = 0.0;

				          break;

			          case WSIDE :

                          dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					      dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[VX][sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
						  if (!b_prosto) {

							  deltal = -0.5*(pa[nvtx[1][sosedb[inumber].iII] - 1].x + pa[nvtx[0][sosedb[inumber].iII] - 1].x);
							  deltal += 0.5*(pa[nvtx[1][sosedb[inumber].iI] - 1].x + pa[nvtx[0][sosedb[inumber].iI] - 1].x);
							  fiplus = 0.5*dl / deltal;

							  taui = (tau[VX][sosedb[inumber].iI] * tau[VX][sosedb[inumber].iII]);
							  taui = taui / ((1.0 - fiplus)*tau[VX][sosedb[inumber].iI] + fiplus*tau[VX][sosedb[inumber].iII]);// проверено !


							  // правая часть:  
							  slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
						  }
						  else slb[PAM][inumber].b = 0.0;

					      break;

		              case SSIDE :

                          dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[VY][sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
						  if (!b_prosto) {

							  deltal = -0.5*(pa[nvtx[2][sosedb[inumber].iII] - 1].y + pa[nvtx[0][sosedb[inumber].iII] - 1].y);
							  deltal += 0.5*(pa[nvtx[2][sosedb[inumber].iI] - 1].y + pa[nvtx[0][sosedb[inumber].iI] - 1].y);
							  fiplus = 0.5*dl / deltal;

							  taui = (tau[VY][sosedb[inumber].iI] * tau[VY][sosedb[inumber].iII]);
							  taui = taui / ((1.0 - fiplus)*tau[VY][sosedb[inumber].iI] + fiplus*tau[VY][sosedb[inumber].iII]);// проверено !


							  // правая часть:  
							  slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
						  }
						  else slb[PAM][inumber].b = 0.0;

				          break;

		              case BSIDE : 

                          dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					      dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[VZ][sosedb[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=sosedb[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=sosedb[inumber].iB;
					
						  if (!b_prosto) {

							  deltal = -0.5*(pa[nvtx[4][sosedb[inumber].iII] - 1].z + pa[nvtx[0][sosedb[inumber].iII] - 1].z);
							  deltal += 0.5*(pa[nvtx[4][sosedb[inumber].iI] - 1].z + pa[nvtx[0][sosedb[inumber].iI] - 1].z);
							  fiplus = 0.5*dl / deltal;

							  taui = (tau[VZ][sosedb[inumber].iI] * tau[VZ][sosedb[inumber].iII]);
							  taui = taui / ((1.0 - fiplus)*tau[VZ][sosedb[inumber].iI] + fiplus*tau[VZ][sosedb[inumber].iII]);// проверено !


							  // правая часть:  
							  slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;

						  }
						  else slb[PAM][inumber].b = 0.0;

				          break;
	                } // switch

			        //*/
                    integer j,l,xitem,k;
			        // сортировка по возрастанию
			        for (j=0; j<5; j++) {
				        k=j; xitem=sosedb[inumber].iW[j];
				        for (l=j+1; l<6; l++) {
					        if (sosedb[inumber].iW[l] < xitem) {
						       k=l; xitem=sosedb[inumber].iW[k];
					        }
				        }
                        sosedb[inumber].iW[k]=sosedb[inumber].iW[j];
				        sosedb[inumber].iW[j]=xitem;
			        }

                    j=0; l=0;
			        while (sosedb[inumber].iW[j]==(-1)) j++;

			        if (j<6) { slb[PAM][inumber].iW1=sosedb[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW2=sosedb[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW3=sosedb[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW4=sosedb[inumber].iW[j++]; l++; } 

			        switch (l) {
				       case 0 : slb[PAM][inumber].iW1=-1;
		                        slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
		                        break;
				       case 1 : slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 2 : slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 3 : slb[PAM][inumber].iW4=-1;
						        break;
			        }

			
				}
			}

	    }
	     else if ((!bDirichlet) /*&& !((inumber==(maxbound-1)) && bPfix && (!bFReeStyle))*/ && !((sosedb[inumber].MCB>=ls) && (sosedb[inumber].MCB<(ls+lw)) && w[sosedb[inumber].MCB-ls].bpressure) )
	    {  
			// Не условие Дирихле, Не фиксация давления в точке, Не выходная граница
			//printf("neiman!\n"); getchar(); // debug

		

			/*
            // поправка давления фиксированно равна нулю:
		    slb[PAM][inumber].aw=1.0;
		    slb[PAM][inumber].ai=0.0;
		    slb[PAM][inumber].b=0.0;
			slb[PAM][inumber].iI=-1;//sosedb[inumber].iI; //  присутствует в матрице
		    slb[PAM][inumber].iW=sosedb[inumber].iB;
			*/

            // Здесь предполагается, что узлы iI и iII внутренние, иначе поведение
		    // программы будет неправильным. НЕ МЕНЕЕ 2-ух ненулевых КО в зазоре между рёбрами.
			//printf("neiman pressure...\n"); // debug
			//getchar();

			///*
		    // однородное условие Неймана
		    doublereal dl, deltal, dS;
	        doublereal taui; // псевдовремя на грани КО.
	        doublereal fiplus; // учёт неравномерности сетки
			//doublereal FgRhieChow; // поправка Рхи-Чоу

			bcheck = true;

			// внутренняя нормаль
	        switch (sosedb[inumber].Norm) {
		        case ESIDE :
			    
				    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[1][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[VX][sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;

					if (!b_prosto) {

						deltal = 0.5*(pa[nvtx[1][sosedb[inumber].iII] - 1].x + pa[nvtx[0][sosedb[inumber].iII] - 1].x);
						deltal -= 0.5*(pa[nvtx[1][sosedb[inumber].iI] - 1].x + pa[nvtx[0][sosedb[inumber].iI] - 1].x);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VX][sosedb[inumber].iI] * tau[VX][sosedb[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VX][sosedb[inumber].iI] + fiplus*tau[VX][sosedb[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;

					}
					else slb[PAM][inumber].b = 0.0;
					/*
					*  Итак лапласиан от поправки давления равен дивергенции скорости.
					*  Дивергенция скорости преобразует векторное поле скоростей
					*  в скалярное поле. Смысл дивергенции скорости для данного 
					*  контророльного объёма (КО) она показывает насколько расходятся
					*  входящий и исходящий поток в данном КО.
					*  Здесь идёт речь о граничном КО. Он является плоскостью (имеет нулевой объём). 
					*  Такой граничный  КО характеризуется нулевой расходимостью т.к. в нём по определению
					*  что вошло то и вышло, т.к. геометрически положение точки входа совпадает с геометрическим 
					*  положением точки выхода.
					*  Вывод: никакой дополнительной добавки в источниковый член не требуется!!! И не должно быть.
					*  Последний вывод проверен рядом тестов.
					*/
				    break;			
			
		        case NSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[VY][sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					if (!b_prosto) {


						deltal = 0.5*(pa[nvtx[2][sosedb[inumber].iII] - 1].y + pa[nvtx[0][sosedb[inumber].iII] - 1].y);
						deltal -= 0.5*(pa[nvtx[2][sosedb[inumber].iI] - 1].y + pa[nvtx[0][sosedb[inumber].iI] - 1].y);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VY][sosedb[inumber].iI] * tau[VY][sosedb[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VY][sosedb[inumber].iI] + fiplus*tau[VY][sosedb[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

				    break;

			    case TSIDE : 

                    dl=pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[VZ][sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					if (!b_prosto) {


						deltal = 0.5*(pa[nvtx[4][sosedb[inumber].iII] - 1].z + pa[nvtx[0][sosedb[inumber].iII] - 1].z);
						deltal -= 0.5*(pa[nvtx[4][sosedb[inumber].iI] - 1].z + pa[nvtx[0][sosedb[inumber].iI] - 1].z);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VZ][sosedb[inumber].iI] * tau[VZ][sosedb[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VZ][sosedb[inumber].iI] + fiplus*tau[VZ][sosedb[inumber].iII]);// проверено !

						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

				    break;

			    case WSIDE :

                    dl=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x;
					dS=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[VX][sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;
					
					if (!b_prosto) {

						deltal = -0.5*(pa[nvtx[1][sosedb[inumber].iII] - 1].x + pa[nvtx[0][sosedb[inumber].iII] - 1].x);
						deltal += 0.5*(pa[nvtx[1][sosedb[inumber].iI] - 1].x + pa[nvtx[0][sosedb[inumber].iI] - 1].x);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VX][sosedb[inumber].iI] * tau[VX][sosedb[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VX][sosedb[inumber].iI] + fiplus*tau[VX][sosedb[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;

		        case SSIDE :

                    dl=pa[nvtx[2][sosedb[inumber].iI]-1].y-pa[nvtx[0][sosedb[inumber].iI]-1].y;
					dS=pa[nvtx[1][sosedb[inumber].iI]-1].x-pa[nvtx[0][sosedb[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][sosedb[inumber].iI]-1].z-pa[nvtx[0][sosedb[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[VY][sosedb[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=sosedb[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=sosedb[inumber].iB;

                    if (!b_prosto) {


	                    deltal = -0.5*(pa[nvtx[2][sosedb[inumber].iII] - 1].y + pa[nvtx[0][sosedb[inumber].iII] - 1].y);
	                    deltal += 0.5*(pa[nvtx[2][sosedb[inumber].iI] - 1].y + pa[nvtx[0][sosedb[inumber].iI] - 1].y);
	                    fiplus = 0.5*dl / deltal;

	                    taui = (tau[VY][sosedb[inumber].iI] * tau[VY][sosedb[inumber].iII]);
	                    taui = taui / ((1.0 - fiplus)*tau[VY][sosedb[inumber].iI] + fiplus*tau[VY][sosedb[inumber].iII]);// проверено !


	                    // правая часть:  
	                    slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
                    }
                    else slb[PAM][inumber].b = 0.0;

                    break;

		        case BSIDE:

					dl = pa[nvtx[4][sosedb[inumber].iI] - 1].z - pa[nvtx[0][sosedb[inumber].iI] - 1].z;
					dS = pa[nvtx[1][sosedb[inumber].iI] - 1].x - pa[nvtx[0][sosedb[inumber].iI] - 1].x;
					dS *= (pa[nvtx[2][sosedb[inumber].iI] - 1].y - pa[nvtx[0][sosedb[inumber].iI] - 1].y); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VZ][sosedb[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = sosedb[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = sosedb[inumber].iB;
					if (!b_prosto) {


						deltal = -0.5*(pa[nvtx[4][sosedb[inumber].iII] - 1].z + pa[nvtx[0][sosedb[inumber].iII] - 1].z);
						deltal += 0.5*(pa[nvtx[4][sosedb[inumber].iI] - 1].z + pa[nvtx[0][sosedb[inumber].iI] - 1].z);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VZ][sosedb[inumber].iI] * tau[VZ][sosedb[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VZ][sosedb[inumber].iI] + fiplus*tau[VZ][sosedb[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][sosedb[inumber].iI] - potent[PAM][sosedb[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;
			} // switch

			//*/
			integer j, l, xitem, k;
			// сортировка по возрастанию
			for (j = 0; j < 5; j++) {
				k = j; xitem = sosedb[inumber].iW[j];
				for (l = j + 1; l < 6; l++) {
					if (sosedb[inumber].iW[l] < xitem) {
						k = l; xitem = sosedb[inumber].iW[k];
					}
				}
				sosedb[inumber].iW[k] = sosedb[inumber].iW[j];
				sosedb[inumber].iW[j] = xitem;
			}

			j = 0; l = 0;
			while (sosedb[inumber].iW[j] == (-1)) j++;

			if (j < 6) { slb[PAM][inumber].iW1 = sosedb[inumber].iW[j++]; l++; }
			if (j < 6) { slb[PAM][inumber].iW2 = sosedb[inumber].iW[j++]; l++; }
			if (j < 6) { slb[PAM][inumber].iW3 = sosedb[inumber].iW[j++]; l++; }
			if (j < 6) { slb[PAM][inumber].iW4 = sosedb[inumber].iW[j++]; l++; }

			switch (l) {
			case 0: slb[PAM][inumber].iW1 = -1;
				slb[PAM][inumber].iW2 = -1;
				slb[PAM][inumber].iW3 = -1;
				slb[PAM][inumber].iW4 = -1;
				break;
			case 1: slb[PAM][inumber].iW2 = -1;
				slb[PAM][inumber].iW3 = -1;
				slb[PAM][inumber].iW4 = -1;
				break;
			case 2: slb[PAM][inumber].iW3 = -1;
				slb[PAM][inumber].iW4 = -1;
				break;
			case 3: slb[PAM][inumber].iW4 = -1;
				break;
			}

			// Наложения исключения в случае совпадения узла iI с диагональным 
			// элементом узла для котторого стоит условие Дирихле не требуется. 
			// т.к. узлы iI строго внутренние для которых iI < maxelm, а диагональный 
			// элемент с условием Дирихле стоит в позициях >= maxelm поэтому они не
			// пересекаются.
		}
		//}

		if (!bDirichlet) {
			if (!bcheck) {
				if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
					//printf("OPENING");
				}
				else {
					/*
					if (fabs(slb[PAM][inumber].ai - slb[PAM][inumber].aw) > 1.0e-9) {
						printf(" ai_ap=%e, b=%e\n", slb[PAM][inumber].ai / slb[PAM][inumber].aw, slb[PAM][inumber].b);
						printf("aw=%e ai=%e b=%e\n", slb[PAM][inumber].aw, slb[PAM][inumber].ai, slb[PAM][inumber].b);

						printf("my_elmatr_quad_PAm_bon3\n");
						system("pause");
					}
					*/
				}
		}
		//if (fabs(slb[PAM][inumber].b) > 1.0e-9) {
			//printf("non zero b: ai_ap=%e, b=%e\n", slb[PAM][inumber].ai / slb[PAM][inumber].aw, slb[PAM][inumber].b);
			//printf("my_elmatr_quad_PAm_bon3\n");
			//system("pause");
		//}

	}

		if (!bcheck) {
			if (!bDirichlet) {
				if (((sosedb[inumber].MCB >= ls) && (sosedb[inumber].MCB < (ls + lw)) && w[sosedb[inumber].MCB - ls].bopening)) {
					//printf("OPENING");
				}
				else {
					/*
					printf("undefined situation for pam found. 10.02.2017\n");
					printf("my_elmatr_quad_PAm_bon3 method.\n");
					//if (bDirichlet) {
						//printf("Dirichlet\n");
					//}
					//else {
						//printf("Neiman\n");
					//}

					system("pause");


					slb[PAM][inumber].ai = 1.0;
					slb[PAM][inumber].iI = sosedb[inumber].iI;
					slb[PAM][inumber].aw = 1.0;
					slb[PAM][inumber].iW = sosedb[inumber].iB;
					slb[PAM][inumber].b = 0.0;
					*/
				}
			}
	}


} // my_elmatr_quad_PAm_bon3


  // возвращает скорректированный массовый поток.
  // скоректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
  // скорректированной скорости или просто скорости на основе простейшей интерполляции.
// Это нужно для отдельного решения уравнения конвекции-диффузии где задана пользовательская скорость,
// никакой поправки Рхи-Чоу просто интерполляция.
// 26.03.2017
void return_calc_correct_mass_flux_only_interpolation(integer iP, doublereal** potent,
	TOCHKA* pa, doublereal** prop, doublereal** prop_b,
	integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm, 
	doublereal* &mfcurrentretune, BOUND* &sosedb)
{

	// SpeedCorOld - скоректированная скорость на предыдущей итерации.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	//bool bRhieChowi = true, bRhieChowb = false; // i - internal, b - border.

	// iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1;
	iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;


	// 26.09.2016 Добавок для АЛИС сетки.
	integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
	integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
	integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов

	// -1 если не используется и [0..maxelm+maxbound-1] если используется.

	iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
	iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;
	iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
	iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;
	iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
	iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


	// Присутствует в СЛАУ
	bool bslE = false, bslN = false, bslT = false, bslW = false, bslS = false, bslB = false;
	bool bslE2 = false, bslN2 = false, bslT2 = false, bslW2 = false, bslS2 = false, bslB2 = false;
	bool bslE3 = false, bslN3 = false, bslT3 = false, bslW3 = false, bslS3 = false, bslB3 = false;
	bool bslE4 = false, bslN4 = false, bslT4 = false, bslW4 = false, bslS4 = false, bslB4 = false;

	// если с одной из сторон граница расчётной области 
	// то переменная равна true
	bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;
	bool bE2 = false, bN2 = false, bT2 = false, bW2 = false, bS2 = false, bB2 = false;
	bool bE3 = false, bN3 = false, bT3 = false, bW3 = false, bS3 = false, bB3 = false;
	bool bE4 = false, bN4 = false, bT4 = false, bW4 = false, bS4 = false, bB4 = false;

	if (iE >= maxelm) bE = true;
	if (iN >= maxelm) bN = true;
	if (iT >= maxelm) bT = true;
	if (iW >= maxelm) bW = true;
	if (iS >= maxelm) bS = true;
	if (iB >= maxelm) bB = true;

	if (iE > -1) bslE = true;
	if (iW > -1) bslW = true;
	if (iN > -1) bslN = true;
	if (iS > -1) bslS = true;
	if (iT > -1) bslT = true;
	if (iB > -1) bslB = true;

	if (iE2 > -1) bslE2 = true;
	if (iW2 > -1) bslW2 = true;
	if (iN2 > -1) bslN2 = true;
	if (iS2 > -1) bslS2 = true;
	if (iT2 > -1) bslT2 = true;
	if (iB2 > -1) bslB2 = true;

	if (iE3 > -1) bslE3 = true;
	if (iW3 > -1) bslW3 = true;
	if (iN3 > -1) bslN3 = true;
	if (iS3 > -1) bslS3 = true;
	if (iT3 > -1) bslT3 = true;
	if (iB3 > -1) bslB3 = true;

	if (iE4 > -1) bslE4 = true;
	if (iW4 > -1) bslW4 = true;
	if (iN4 > -1) bslN4 = true;
	if (iS4 > -1) bslS4 = true;
	if (iT4 > -1) bslT4 = true;
	if (iB4 > -1) bslB4 = true;

	if (iE2 >= maxelm) bE2 = true;
	if (iW2 >= maxelm) bW2 = true;
	if (iN2 >= maxelm) bN2 = true;
	if (iS2 >= maxelm) bS2 = true;
	if (iT2 >= maxelm) bT2 = true;
	if (iB2 >= maxelm) bB2 = true;

	if (iE3 >= maxelm) bE3 = true;
	if (iW3 >= maxelm) bW3 = true;
	if (iN3 >= maxelm) bN3 = true;
	if (iS3 >= maxelm) bS3 = true;
	if (iT3 >= maxelm) bT3 = true;
	if (iB3 >= maxelm) bB3 = true;

	if (iE4 >= maxelm) bE4 = true;
	if (iW4 >= maxelm) bW4 = true;
	if (iN4 >= maxelm) bN4 = true;
	if (iS4 >= maxelm) bS4 = true;
	if (iT4 >= maxelm) bT4 = true;
	if (iB4 >= maxelm) bB4 = true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx = 0.0, dy = 0.0, dz = 0.0; // размеры контрольного объёма
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

	//printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	//getchar();

	// плотность аппроксимируется средним гармоническим
	doublereal rhoe=0.0, rhow = 0.0, rhon = 0.0, rhos = 0.0, rhot = 0.0, rhob = 0.0;
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


	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE > -1) {
		if (!bE) rhoe = rE * rP / (feplus*rE + (1.0 - feplus)*rP); else rhoe = rE;
	}
	if (iW > -1) {
		if (!bW) rhow = rW * rP / (fwplus*rW + (1.0 - fwplus)*rP); else rhow = rW;
	}
	if (iN > -1) {
		if (!bN) rhon = rN * rP / (fnplus*rN + (1.0 - fnplus)*rP); else rhon = rN;
	}
	if (iS > -1) {
		if (!bS) rhos = rS * rP / (fsplus*rS + (1.0 - fsplus)*rP); else rhos = rS;
	}
	if (iT > -1) {
		if (!bT) rhot = rT * rP / (ftplus*rT + (1.0 - ftplus)*rP); else rhot = rT;
	}
	if (iB > -1) {
		if (!bB) rhob = rB * rP / (fbplus*rB + (1.0 - fbplus)*rP); else rhob = rB;
	}


	doublereal rhoe2 = 0.0, rhow2 = 0.0, rhon2 = 0.0, rhos2 = 0.0, rhot2 = 0.0, rhob2 = 0.0;
	doublereal rhoe3 = 0.0, rhow3 = 0.0, rhon3 = 0.0, rhos3 = 0.0, rhot3 = 0.0, rhob3 = 0.0;
	doublereal rhoe4 = 0.0, rhow4 = 0.0, rhon4 = 0.0, rhos4 = 0.0, rhot4 = 0.0, rhob4 = 0.0;

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE2 > -1) {
		if (!bE2) rhoe2 = rE2 * rP / (feplus2*rE2 + (1.0 - feplus2)*rP); else rhoe2 = rE2;
	}
	if (iW2 > -1) {
		if (!bW2) rhow2 = rW2 * rP / (fwplus2*rW2 + (1.0 - fwplus2)*rP); else rhow2 = rW2;
	}
	if (iN2 > -1) {
		if (!bN2) rhon2 = rN2 * rP / (fnplus2*rN2 + (1.0 - fnplus2)*rP); else rhon2 = rN2;
	}
	if (iS2 > -1) {
		if (!bS2) rhos2 = rS2 * rP / (fsplus2*rS2 + (1.0 - fsplus2)*rP); else rhos2 = rS2;
	}
	if (iT2 > -1) {
		if (!bT2) rhot2 = rT2 * rP / (ftplus2*rT2 + (1.0 - ftplus2)*rP); else rhot2 = rT2;
	}
	if (iB2 > -1) {
		if (!bB2) rhob2 = rB2 * rP / (fbplus2*rB2 + (1.0 - fbplus2)*rP); else rhob2 = rB2;
	}

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE3 > -1) {
		if (!bE3) rhoe3 = rE3 * rP / (feplus3*rE3 + (1.0 - feplus3)*rP); else rhoe3 = rE3;
	}
	if (iW3 > -1) {
		if (!bW3) rhow3 = rW3 * rP / (fwplus3*rW3 + (1.0 - fwplus3)*rP); else rhow3 = rW3;
	}
	if (iN3 > -1) {
		if (!bN3) rhon3 = rN3 * rP / (fnplus3*rN3 + (1.0 - fnplus3)*rP); else rhon3 = rN3;
	}
	if (iS3 > -1) {
		if (!bS3) rhos3 = rS3 * rP / (fsplus3*rS3 + (1.0 - fsplus3)*rP); else rhos3 = rS3;
	}
	if (iT3 > -1) {
		if (!bT3) rhot3 = rT3 * rP / (ftplus3*rT3 + (1.0 - ftplus3)*rP); else rhot3 = rT3;
	}
	if (iB3 > -1) {
		if (!bB3) rhob3 = rB3 * rP / (fbplus3*rB3 + (1.0 - fbplus3)*rP); else rhob3 = rB3;
	}

	if (iE4 > -1) {
		if (!bE4) rhoe4 = rE4 * rP / (feplus4*rE4 + (1.0 - feplus4)*rP); else rhoe4 = rE4;
	}
	if (iW4 > -1) {
		if (!bW4) rhow4 = rW4 * rP / (fwplus4*rW4 + (1.0 - fwplus4)*rP); else rhow4 = rW4;
	}
	if (iN4 > -1) {
		if (!bN4) rhon4 = rN4 * rP / (fnplus4*rN4 + (1.0 - fnplus4)*rP); else rhon4 = rN4;
	}
	if (iS4 > -1) {
		if (!bS4) rhos4 = rS4 * rP / (fsplus4*rS4 + (1.0 - fsplus4)*rP); else rhos4 = rS4;
	}
	if (iT4 > -1) {
		if (!bT4) rhot4 = rT4 * rP / (ftplus4*rT4 + (1.0 - ftplus4)*rP); else rhot4 = rT4;
	}
	if (iB4 > -1) {
		if (!bB4) rhob4 = rB4 * rP / (fbplus4*rB4 + (1.0 - fbplus4)*rP); else rhob4 = rB4;
	}

	

	doublereal Fw = 0.0, Fe = 0.0, Fs = 0.0, Fn = 0.0, Ft = 0.0, Fb = 0.0;

	// Для АЛИС сетки.
	//doublereal Fe1 = 0.0, Fe2 = 0.0, Fe3 = 0.0, Fe4 = 0.0;
	//doublereal Fw1 = 0.0, Fw2 = 0.0, Fw3 = 0.0, Fw4 = 0.0;
	//doublereal Fn1 = 0.0, Fn2 = 0.0, Fn3 = 0.0, Fn4 = 0.0;
	//doublereal Fs1 = 0.0, Fs2 = 0.0, Fs3 = 0.0, Fs4 = 0.0;
	//doublereal Ft1 = 0.0, Ft2 = 0.0, Ft3 = 0.0, Ft4 = 0.0;
	//doublereal Fb1 = 0.0, Fb2 = 0.0, Fb3 = 0.0, Fb4 = 0.0;

	

		doublereal SpeedCorOlde=0.0, SpeedCorOldw = 0.0, SpeedCorOldn = 0.0, SpeedCorOlds = 0.0, SpeedCorOldt = 0.0, SpeedCorOldb = 0.0;
		if (iE > -1) {
			if (!bE) {
				SpeedCorOlde = feplus * potent[VX][iE] + (1.0 - feplus)*potent[VX][iP];
			}
			else {
				SpeedCorOlde = potent[VX][iE];
			}
		}
		if (iN > -1) {
			if (!bN) {
				SpeedCorOldn = fnplus * potent[VY][iN] + (1.0 - fnplus)*potent[VY][iP];
			}
			else {
				SpeedCorOldn = potent[VY][iN];
			}
		}
		if (iT > -1) {
			if (!bT) {
				SpeedCorOldt = ftplus * potent[VZ][iT] + (1.0 - ftplus)*potent[VZ][iP];
			}
			else {
				SpeedCorOldt = potent[VZ][iT];
			}
		}
		if (iW > -1) {
			if (!bW) {
				SpeedCorOldw = fwplus * potent[VX][iW] + (1.0 - fwplus)*potent[VX][iP];
			}
			else {
				SpeedCorOldw = potent[VX][iW];
			}
		}
		if (iS > -1) {
			if (!bS) {
				SpeedCorOlds = fsplus * potent[VY][iS] + (1.0 - fsplus)*potent[VY][iP];
			}
			else {
				SpeedCorOlds = potent[VY][iS];
			}
		}
		if (iB > -1) {
			if (!bB) {
				SpeedCorOldb = fbplus * potent[VZ][iB] + (1.0 - fbplus)*potent[VZ][iP];
			}
			else {
				SpeedCorOldb = potent[VZ][iB];
			}
		}
		/*
		if (fabs(potent[VZ][iP]) + fabs(potent[VZ][iB]) + fabs(potent[VZ][iT]) > 0.0) {
			printf("Vz non zero Ok: %e %e %e\n", potent[VZ][iP], potent[VZ][iB], potent[VZ][iT]);
			getchar();
		}
		if (fabs(SpeedCorOldt) + fabs(SpeedCorOldb) > 0.0) {
			printf("non zero mf. Ok.\n");
			getchar();
		}
		*/

		doublereal SpeedCorOlde2 = 0.0, SpeedCorOldw2 = 0.0, SpeedCorOldn2 = 0.0, SpeedCorOlds2 = 0.0, SpeedCorOldt2 = 0.0, SpeedCorOldb2 = 0.0;
		doublereal SpeedCorOlde3 = 0.0, SpeedCorOldw3 = 0.0, SpeedCorOldn3 = 0.0, SpeedCorOlds3 = 0.0, SpeedCorOldt3 = 0.0, SpeedCorOldb3 = 0.0;
		doublereal SpeedCorOlde4 = 0.0, SpeedCorOldw4 = 0.0, SpeedCorOldn4 = 0.0, SpeedCorOlds4 = 0.0, SpeedCorOldt4 = 0.0, SpeedCorOldb4 = 0.0;

		if (iE2 > -1) {
			if (!bE2) {
				SpeedCorOlde2 = feplus2 * potent[VX][iE2] + (1.0 - feplus2)*potent[VX][iP];
			}
			else {
				SpeedCorOlde2 = potent[VX][iE2];
			}
		}
		if (iN2 > -1) {
			if (!bN2) {
				SpeedCorOldn2 = fnplus2 * potent[VY][iN2] + (1.0 - fnplus2)*potent[VY][iP];
			}
			else {
				SpeedCorOldn2 = potent[VY][iN2];
			}
		}
		if (iT2 > -1) {
			if (!bT2) {
				SpeedCorOldt2 = ftplus2 * potent[VZ][iT2] + (1.0 - ftplus2)*potent[VZ][iP];
			}
			else {
				SpeedCorOldt2 = potent[VZ][iT2];
			}
		}
		if (iW2 > -1) {
			if (!bW2) {
				SpeedCorOldw2 = fwplus2 * potent[VX][iW2] + (1.0 - fwplus2)*potent[VX][iP];
			}
			else {
				SpeedCorOldw2 = potent[VX][iW2];
			}
		}
		if (iS2 > -1) {
			if (!bS2) {
				SpeedCorOlds2 = fsplus2 * potent[VY][iS2] + (1.0 - fsplus2)*potent[VY][iP];
			}
			else {
				SpeedCorOlds2 = potent[VY][iS2];
			}
		}
		if (iB2 > -1) {
			if (!bB2) {
				SpeedCorOldb2 = fbplus2 * potent[VZ][iB2] + (1.0 - fbplus2)*potent[VZ][iP];
			}
			else {
				SpeedCorOldb2 = potent[VZ][iB2];
			}
		}

		if (iE3 > -1) {
			if (!bE3) {
				SpeedCorOlde3 = feplus3 * potent[VX][iE3] + (1.0 - feplus3)*potent[VX][iP];
			}
			else {
				SpeedCorOlde3 = potent[VX][iE3];
			}
		}
		if (iN3 > -1) {
			if (!bN3) {
				SpeedCorOldn3 = fnplus3 * potent[VY][iN3] + (1.0 - fnplus3)*potent[VY][iP];
			}
			else {
				SpeedCorOldn3 = potent[VY][iN3];
			}
		}
		if (iT3 > -1) {
			if (!bT3) {
				SpeedCorOldt3 = ftplus3 * potent[VZ][iT3] + (1.0 - ftplus3)*potent[VZ][iP];
			}
			else {
				SpeedCorOldt3 = potent[VZ][iT3];
			}
		}
		if (iW3 > -1) {
			if (!bW3) {
				SpeedCorOldw3 = fwplus3 * potent[VX][iW3] + (1.0 - fwplus3)*potent[VX][iP];
			}
			else {
				SpeedCorOldw3 = potent[VX][iW3];
			}
		}
		if (iS3 > -1) {
			if (!bS3) {
				SpeedCorOlds3 = fsplus3 * potent[VY][iS3] + (1.0 - fsplus3)*potent[VY][iP];
			}
			else {
				SpeedCorOlds3 = potent[VY][iS3];
			}
		}
		if (iB3 > -1) {
			if (!bB3) {
				SpeedCorOldb3 = fbplus3 * potent[VZ][iB3] + (1.0 - fbplus3)*potent[VZ][iP];
			}
			else {
				SpeedCorOldb3 = potent[VZ][iB3];
			}
		}

		if (iE4 > -1) {
			if (!bE4) {
				SpeedCorOlde4 = feplus4 * potent[VX][iE4] + (1.0 - feplus4)*potent[VX][iP];
			}
			else {
				SpeedCorOlde4 = potent[VX][iE4];
			}
		}
		if (iN4 > -1) {
			if (!bN4) {
				SpeedCorOldn4 = fnplus4 * potent[VY][iN4] + (1.0 - fnplus4)*potent[VY][iP];
			}
			else {
				SpeedCorOldn4 = potent[VY][iN4];
			}
		}
		if (iT4 > -1) {
			if (!bT4) {
				SpeedCorOldt4 = ftplus4 * potent[VZ][iT4] + (1.0 - ftplus4)*potent[VZ][iP];
			}
			else {
				SpeedCorOldt4 = potent[VZ][iT4];
			}
		}
		if (iW4 > -1) {
			if (!bW4) {
				SpeedCorOldw4 = fwplus4 * potent[VX][iW4] + (1.0 - fwplus4)*potent[VX][iP];
			}
			else {
				SpeedCorOldw4 = potent[VX][iW4];
			}
		}
		if (iS4 > -1) {
			if (!bS4) {
				SpeedCorOlds4 = fsplus4 * potent[VY][iS4] + (1.0 - fsplus4)*potent[VY][iP];
			}
			else {
				SpeedCorOlds4 = potent[VY][iS4];
			}
		}
		if (iB4 > -1) {
			if (!bB4) {
				SpeedCorOldb4 = fbplus4 * potent[VZ][iB4] + (1.0 - fbplus4)*potent[VZ][iP];
			}
			else {
				SpeedCorOldb4 = potent[VZ][iB4];
			}
		}


		// mfold - значение массового потока с предыдущей итерации.
		Fe =  (rhoe*SpeedCorOlde+ rhoe2*SpeedCorOlde2 + rhoe3*SpeedCorOlde3 + rhoe4*SpeedCorOlde4)*dy*dz;
		Fn = (rhon*SpeedCorOldn + rhon2*SpeedCorOldn2 + rhon3*SpeedCorOldn3 + rhon4*SpeedCorOldn4)*dx*dz;
		Ft = (rhot*SpeedCorOldt + rhot2*SpeedCorOldt2 + rhot3*SpeedCorOldt3 + rhot4*SpeedCorOldt4)*dx*dy;
		Fw = (rhow*SpeedCorOldw + rhow2*SpeedCorOldw2 + rhow3*SpeedCorOldw3 + rhow4*SpeedCorOldw4)*dy*dz;
		Fs = (rhos*SpeedCorOlds + rhos2*SpeedCorOlds2 + rhos3*SpeedCorOlds3 + rhos4*SpeedCorOlds4)*dx*dz;
		Fb = (rhob*SpeedCorOldb + rhob2*SpeedCorOldb2 + rhob3*SpeedCorOldb3 + rhob4*SpeedCorOldb4)*dx*dy;
	
		//printf("dx=%e dy=%e rhot=%e rhob=%e SpeedCorOldt=%e SpeedCorOldb=%e Ft=%e Fb=%e \n",dx,dy,rhot,rhob, SpeedCorOldt, SpeedCorOldb, Ft, Fb);
		//printf("dx=%e dy=%e rhot2=%e rhob2=%e SpeedCorOldt2=%e SpeedCorOldb2=%e Ft=%e Fb=%e \n", dx, dy, rhot2, rhob2, SpeedCorOldt2, SpeedCorOldb2, Ft, Fb);
		//printf("dx=%e dy=%e rhot3=%e rhob3=%e SpeedCorOldt3=%e SpeedCorOldb3=%e Ft=%e Fb=%e \n", dx, dy, rhot3, rhob3, SpeedCorOldt3, SpeedCorOldb3, Ft, Fb);
		//printf("dx=%e dy=%e rhot4=%e rhob4=%e SpeedCorOldt4=%e SpeedCorOldb4=%e Ft=%e Fb=%e \n", dx, dy, rhot4, rhob4, SpeedCorOldt4, SpeedCorOldb4, Ft, Fb);
		//getchar();

		//if (fabs(Ft) + fabs(Fb) > 0.0) {
			//printf("non zero Ft=%e and Fb=%e. Ok.\n",Ft,Fb);
			//getchar();
		//}
	
		mfcurrentretune[ESIDE] = Fe;
		mfcurrentretune[NSIDE] = Fn;
		mfcurrentretune[TSIDE] = Ft;
		mfcurrentretune[WSIDE] = Fw;
		mfcurrentretune[SSIDE] = Fs;
		mfcurrentretune[BSIDE] = Fb;

		// Для прохождения проверки корректности необходимо
		// принудительно занулить скорость и поток скорости
		// на твердой сенке (нормальную к стенке компоненту).
		// Твердая стенка идентифицируентся маркером MCB == (ls + lw).

		if (iE >= maxelm) {
			// граничный узел
			integer inumber = iE - maxelm;
			if (sosedb[inumber].MCB == (ls + lw)) {
				mfcurrentretune[ESIDE] = 0.0;
			}
		}
				
		if (iW >= maxelm) {
			// граничный узел
			integer inumber = iW - maxelm;
			if (sosedb[inumber].MCB == (ls + lw)) {
				mfcurrentretune[WSIDE] = 0.0;
			}
		}

		if (iN >= maxelm) {
			// граничный узел
			integer inumber = iN - maxelm;
			if (sosedb[inumber].MCB == (ls + lw)) {
				mfcurrentretune[NSIDE] = 0.0;
			}
		}

		if (iS >= maxelm) {
			// граничный узел
			integer inumber = iS - maxelm;
			if (sosedb[inumber].MCB == (ls + lw)) {
				mfcurrentretune[SSIDE] = 0.0;
			}
		}

		if (iT >= maxelm) {
			// граничный узел
			integer inumber = iT - maxelm;
			if (sosedb[inumber].MCB == (ls + lw)) {
				mfcurrentretune[TSIDE] = 0.0;
			}
		}

		if (iB >= maxelm) {
			// граничный узел
			integer inumber = iB - maxelm;
			if (sosedb[inumber].MCB == (ls + lw)) {
				mfcurrentretune[BSIDE] = 0.0;
			}
		}


} // return_correct_mass_flux_only_interpolation


// возвращает значение скорости на грани с использованием монотонизирующей поправки Рхи-Чоу.
// Внимание : данная реализация подходит не только для стационарных задач, но и для нестационарных.
// Для нестационарных задач требуется хранить массовый поток на грани с предыдущего временного слоя.
// Основная сложность данной функции заключается в её правильном вызове (требуется правильно указать все параметры
// при вызове) из-за большого количества передаваемых параметров.
// реализовано 31 марта 2012 года.
doublereal calcFg(bool bG, doublereal fgplus, doublereal VG, doublereal VP, 
	        integer iP, integer G, doublereal alpha, doublereal rhog, 
			doublereal dS, doublereal dV, bool btimedepend, 
			doublereal** speedoldtimestep, doublereal* mfoldtimestep,
			doublereal dG, doublereal dP, doublereal dtimestep,
			doublereal RCh, integer** nvtx, ALICE_PARTITION** sosedi,
			integer maxelm, doublereal* pressure,  TOCHKA* pa, doublereal **diag_coef,
			bool bRhieChowi, bool bRhieChowb, bool bsimplelinearinterpol) {
	/*
	*  bG - является ли грань граничной ? 
	*  fgplus - интерполяционный фактор,
	*  VG, VP - значение соответствующей скорости в узлах iP и iG. 
	*  RCh (= 1.0 по умолчанию) - используется для уменьшения вклада поправки Рхи-Чоу.
	*  iP - номер рассматриваемого контрольного объёма,
	*  G - идентификатор одной из 6 граней,
	* alpha - коэффициент релаксации для скорости,
	* rhog - плотность на грани КО, 
	* dS - площадь грани контрольного объёма,
	* dV - объём контрольного объёма,
	* btimedepend- стационарный или нестационарный солвер ?
	* speedoldtimestep[VX,VY,VZ][iP,iG] - массив скоростей с предыдущего временного шага,
	* VGoldtimestep - скорость в узле iG с предыдущего временного слоя,
	* VPoldtimestep - скорость в узле iP с предыдущего временного слоя,
	* mfoldtimestep - массовый поток с предыдущего временного слоя на грани G ( массив из шести значений для любого значения G=E,W,N,S,T,B).
	* dG - диагональный коэффициент в матрице для закона сохранения импульса в узле iG,
	* dP - диагональный коэффициент в матрице для закона сохранения импульса в узле iP,
	* dtimestep - размер шага по времени,
	* pressure - поле давления,
	* diag_coef - диагональные коэффициенты в матрицах для компонент скорости которые были запомнены сразу после решения ур-ий сохранения импульса.
	* bRhieChowi (internal) - задействовать ли монотонизирующую поправку Рхи-Чоу для внутренних узлов расчётной области.
	* bRhieChowb (border) - задействовать ли монотонизирующую поправку Рхи-Чоу для граничных узлов расчётной области.
	* если bsimplelinearinterpol равно истине то нужно выполнить просто линейную интерполяцию скорости на грань без использования поправки Рхи-Чоу.
	* такая возможность требуется на первой итерации стационарного алгоритма.
	*/

	doublereal ug=0.0; // скорость на грани контрольного объёма.
	doublereal Fg=0.0; // массовый поток на грани контрольного объёма.

	if (bsimplelinearinterpol) {

		// самая обычная линейная интерполяция скорости на грань,
		// используется на первой итерации стационарного солвера.
        if (!bG) { 
	          ug=fgplus*VG+(1.0-fgplus)*VP;
	    }
	    else {
	          ug=VG;
	    }

		// конвективный поток через грань КО.
	    Fg=rhog*ug*dS;
	}
	else
	{

	// это делается лишь для того чтобы добавить дополнительную
	// нижнюю релаксацию так как это написано у I. Sezai.
	if (!bG) { 
	   ug=fgplus*VG+(1.0-fgplus)*VP;
	   if (bRhieChowi) {
	       ug+=RCh*ugRhieChow_internal(iP, G, alpha, nvtx, sosedi, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
	   }
	}
	else {
	   ug=VG;
	   if (bRhieChowb) {
	      ug+=RCh*ugRhieChow_internal_border(iP, G, alpha, nvtx, sosedi, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
	   }
	}

	// конвективный поток через грань КО.
	Fg=rhog*ug*dS;

	// Внимание здесь предполагается что плотность строго постоянна.
	if (btimedepend) {
		// в нестационарном случае к поправке Рхи-Чоу требуется сделать добавку.

		doublereal VGoldtimestep=0.0;
		doublereal VPoldtimestep=0.0;
		integer iG=sosedi[G][iP].iNODE1;
		switch (G) {
		 case ESIDE : case WSIDE : VPoldtimestep=speedoldtimestep[VX][iP];
			               VGoldtimestep=speedoldtimestep[VX][iG];
			      break;
		 case NSIDE : case SSIDE :VPoldtimestep=speedoldtimestep[VY][iP];
			               VGoldtimestep=speedoldtimestep[VY][iG];
			      break;
		 case TSIDE :  case BSIDE : VPoldtimestep=speedoldtimestep[VZ][iP];
			                VGoldtimestep=speedoldtimestep[VZ][iG];
			      break;
		 }

		// Скорость на грани КО с предыдущего временного слоя:
		doublereal ugold=fgplus*VGoldtimestep+(1.0-fgplus)*VPoldtimestep; 


		// диагональный коэффициент матрицы должен быть взят на границе контрольного объёма.
		if (!bG) {
			// внутренний контрольный объём
		    doublereal dc_e=dP*dG/(fgplus*dG+(1.0-fgplus)*dP); // диагональный коэффициент на грани
		    doublereal tau_e=(rhog*alpha*dV)/(dc_e); // псевдовремя (rhog - плотность на грани контрольного объёма.
			if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) tau_e/=(1.0-alpha);
		    Fg+=(tau_e/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
		}
		else {
			// граничный контрольный объём
			doublereal tau_e=(rhog*alpha*dV)/dG;
			if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) tau_e/=(1.0-alpha);
			Fg+=(tau_e/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
		}

	}
	}

	return Fg;

} // calcFg


// возвращает значение скорости на грани с использованием монотонизирующей поправки Рхи-Чоу.
// Внимание : данная реализация подходит не только для стационарных задач, но и для нестационарных.
// Для нестационарных задач требуется хранить массовый поток на грани с предыдущего временного слоя.
// Основная сложность данной функции заключается в её правильном вызове (требуется правильно указать все параметры
// при вызове) из-за большого количества передаваемых параметров.
// реализовано 21 июня 2012 года. (более простой вариант по сравнению с calcFg).
// Смотри Гаврилов Андрей опыт по ВГД (CFD).
doublereal calcFg2(bool bG, doublereal fgplus, 
	        integer iP, integer G, doublereal rhog, 
			doublereal dS, bool btimedepend, 
			doublereal** speedoldtimestep, doublereal* mfoldtimestep,
			doublereal dtimestep, doublereal RCh, ALICE_PARTITION** sosedi,
			bool bRhieChowi, bool bRhieChowb, bool bsimplelinearinterpol,
			doublereal* tau, doublereal drg, doublereal** potent) {
	/*
	*  bG - является ли грань граничной ? 
	*  fgplus - интерполяционный фактор,
	*  VG, VP - значение соответствующей скорости в узлах iP и iG. 
	*  RCh (= 1.0 по умолчанию) - используется для уменьшения вклада поправки Рхи-Чоу.
	*  iP - номер рассматриваемого контрольного объёма,
	*  G - идентификатор одной из 6 граней,
	* rhog - плотность на грани КО, 
	* dS - площадь грани контрольного объёма,
	* btimedepend- стационарный или нестационарный солвер ?
	* speedoldtimestep[VX,VY,VZ][iP,iG] - массив скоростей с предыдущего временного шага,
	* VGoldtimestep - скорость в узле iG с предыдущего временного слоя,
	* VPoldtimestep - скорость в узле iP с предыдущего временного слоя,
	* mfoldtimestep - массовый поток с предыдущего временного слоя на грани G ( массив из шести значений для любого значения G=E,W,N,S,T,B).
	* dtimestep - размер шага по времени,
	* potent - давление и градиент давления,
	* bRhieChowi (internal) - задействовать ли монотонизирующую поправку Рхи-Чоу для внутренних узлов расчётной области.
	* bRhieChowb (border) - задействовать ли монотонизирующую поправку Рхи-Чоу для граничных узлов расчётной области.
	* если bsimplelinearinterpol равно истине то нужно выполнить просто линейную интерполяцию скорости на грань без использования поправки Рхи-Чоу.
	* такая возможность требуется на первой итерации стационарного алгоритма.
	* tau - сглаженное псевдовремя.
	* drg - расстояние между контрольными объёмами iP и iG.
	*/

	integer iG=sosedi[G][iP].iNODE1;
	doublereal VG=0.0, VP=0.0;
	switch (G) {
		case ESIDE : case WSIDE : VG=potent[VX][iG]; VP=potent[VX][iP]; break;
		case NSIDE : case SSIDE :VG=potent[VY][iG]; VP=potent[VY][iP]; break;
		case TSIDE : case BSIDE : VG=potent[VZ][iG]; VP=potent[VZ][iP]; break;
	} // end switch
	

	doublereal ug=0.0; // скорость на грани контрольного объёма.
	doublereal Fg=0.0; // массовый поток на грани контрольного объёма.

	if (bsimplelinearinterpol) {

		// самая обычная линейная интерполяция скорости на грань,
		// используется на первой итерации стационарного солвера.
        if (!bG) { 
	          ug=fgplus*VG+(1.0-fgplus)*VP;
	    }
	    else {
	          ug=VG;
	    }

		// конвективный поток через грань КО.
	    Fg=rhog*ug*dS;
	}
	else
	{

		
        doublereal taug=tau[iP]*tau[iG]/(fgplus*tau[iG]+(1.0-fgplus)*tau[iP]); // псевдовремя на грани КО.
	    doublereal gradP=0.0, gradG=0.0;
	    switch(G) {
	       case ESIDE : gradP=potent[GRADXPRESS][iP]; gradG=potent[GRADXPRESS][sosedi[ESIDE][iP].iNODE1]; break;
	       case WSIDE : gradP=potent[GRADXPRESS][iP]; gradG=potent[GRADXPRESS][sosedi[WSIDE][iP].iNODE1]; break;
	       case NSIDE : gradP=potent[GRADYPRESS][iP]; gradG=potent[GRADYPRESS][sosedi[NSIDE][iP].iNODE1]; break;
	       case SSIDE :gradP=potent[GRADYPRESS][iP]; gradG=potent[GRADYPRESS][sosedi[SSIDE][iP].iNODE1]; break;
	       case TSIDE : gradP=potent[GRADZPRESS][iP]; gradG=potent[GRADZPRESS][sosedi[TSIDE][iP].iNODE1]; break;
	       case BSIDE : gradP=potent[GRADZPRESS][iP]; gradG=potent[GRADZPRESS][sosedi[BSIDE][iP].iNODE1]; break;
	    }


	    // это делается лишь для того чтобы добавить дополнительную
	    // нижнюю релаксацию так как это написано у I. Sezai.
	    if (!bG) { 
	       ug=fgplus*VG+(1.0-fgplus)*VP;
	       if (bRhieChowi) {
	           //ug+=RCh*ugRhieChow_internal(iP, G, alpha, nvtx, sosedi, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
		       ug+=RCh*(taug/rhog)*(-(potent[PRESS][iG]-potent[PRESS][iP])/drg+(1.0-fgplus)*gradP+fgplus*gradG); // Вклад поправки Рхи-Чоу
	       }
	    }
	    else {
	       ug=VG;
	       /*
	       // поправи для граничного узла нету. Т.е. считаем её нулевой.
	       if (bRhieChowb) {
	          set additional ug+=0.0;
	          //ug+=RCh*ugRhieChow_internal_border(iP, G, alpha, nvtx, sosedi, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
	       }
	       */
	    }

	    // конвективный поток через грань КО.
	    Fg=rhog*ug*dS;

 	    // Внимание здесь предполагается что плотность строго постоянна.
	    if (btimedepend) {
			// в нестационарном случае к поправке Рхи-Чоу требуется сделать добавку.

		    doublereal VGoldtimestep=0.0;
		    doublereal VPoldtimestep=0.0;

		    switch (G) {
		       case ESIDE : case WSIDE : VPoldtimestep=speedoldtimestep[VX][iP];
			                     VGoldtimestep=speedoldtimestep[VX][iG];
			                     break;
		       case NSIDE : case SSIDE :VPoldtimestep=speedoldtimestep[VY][iP];
			                     VGoldtimestep=speedoldtimestep[VY][iG];
			                     break;
		       case TSIDE :  case BSIDE : VPoldtimestep=speedoldtimestep[VZ][iP];
			                      VGoldtimestep=speedoldtimestep[VZ][iG];
			                     break;
		    }

		    // Скорость на грани КО с предыдущего временного слоя:
		    doublereal ugold=fgplus*VGoldtimestep+(1.0-fgplus)*VPoldtimestep; 


		    // диагональный коэффициент матрицы должен быть взят на границе контрольного объёма.
		    if (!bG) {
			    // внутренний контрольный объём
		        Fg+=(taug/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
		    }
		    else {
			    // граничный контрольный объём
			    Fg+=(taug/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
		    }

    	}
	}

	return Fg;

} // calcFg2

// возвращает значение скорости на грани с использованием монотонизирующей поправки Рхи-Чоу.
// Внимание : данная реализация подходит не только для стационарных задач, но и для нестационарных.
// Для нестационарных задач требуется хранить массовый поток на грани с предыдущего временного слоя.
// Основная сложность данной функции заключается в её правильном вызове (требуется правильно указать все параметры
// при вызове) из-за большого количества передаваемых параметров.
// реализовано 23 июня 2012 года. (более простой вариант по сравнению с calcFg).
// Смотри Гаврилов Андрей опыт по ВГД (CFD).
// Данная версия основана на трёх скалярных полях псевдовремени.
doublereal calcFg3(bool bG, doublereal fgplus, 
	        integer iP, integer G, doublereal rhog, 
			doublereal dS, bool btimedepend, 
			doublereal** speedoldtimestep, doublereal* mfoldtimestep,
			doublereal dtimestep, doublereal RCh, ALICE_PARTITION** sosedi,
			bool bRhieChowi, bool bRhieChowb, bool bsimplelinearinterpol,
			doublereal** tau, doublereal drg, doublereal** potent, integer iNODEid) {
	/*
	*  bG - является ли грань граничной ? 
	*  fgplus - интерполяционный фактор,
	*  VG, VP - значение соответствующей скорости в узлах iP и iG. 
	*  RCh (= 1.0 по умолчанию) - используется для уменьшения вклада поправки Рхи-Чоу.
	*  iP - номер рассматриваемого контрольного объёма,
	*  G - идентификатор одной из 6 граней,
	* rhog - плотность на грани КО, 
	* dS - площадь грани контрольного объёма,
	* btimedepend- стационарный или нестационарный солвер ?
	* speedoldtimestep[VX,VY,VZ][iP,iG] - массив скоростей с предыдущего временного шага,
	* VGoldtimestep - скорость в узле iG с предыдущего временного слоя,
	* VPoldtimestep - скорость в узле iP с предыдущего временного слоя,
	* mfoldtimestep - массовый поток с предыдущего временного слоя на грани G ( массив из шести значений для любого значения G=E,W,N,S,T,B).
	* dtimestep - размер шага по времени,
	* potent - давление и градиент давления,
	* bRhieChowi (internal) - задействовать ли монотонизирующую поправку Рхи-Чоу для внутренних узлов расчётной области.
	* bRhieChowb (border) - задействовать ли монотонизирующую поправку Рхи-Чоу для граничных узлов расчётной области.
	* если bsimplelinearinterpol равно истине то нужно выполнить просто линейную интерполяцию скорости на грань без использования поправки Рхи-Чоу.
	* такая возможность требуется на первой итерации стационарного алгоритма.
	* tau - сглаженное псевдовремя.
	* drg - расстояние между контрольными объёмами iP и iG.
	*/

	integer iG=sosedi[G][iP].iNODE1;
	switch (iNODEid) {
		case 1: iG = sosedi[G][iP].iNODE1; break;
		case 2: iG = sosedi[G][iP].iNODE2; break;
		case 3: iG = sosedi[G][iP].iNODE3; break;
		case 4: iG = sosedi[G][iP].iNODE4; break;
		default: iG = sosedi[G][iP].iNODE1;  break;
	}

	doublereal VG=0.0, VP=0.0;
	integer VGid=-1;
	switch (G) {
		case ESIDE : case WSIDE : VG=potent[VX][iG]; VP=potent[VX][iP]; VGid=VX; break;
		case NSIDE : case SSIDE : VG=potent[VY][iG]; VP=potent[VY][iP]; VGid=VY; break;
		case TSIDE : case BSIDE : VG=potent[VZ][iG]; VP=potent[VZ][iP]; VGid=VZ; break;
	} // end switch
	

	doublereal ug=0.0; // скорость на грани контрольного объёма.
	doublereal Fg=0.0; // массовый поток на грани контрольного объёма.

	if (bsimplelinearinterpol) {

		// самая обычная линейная интерполяция скорости на грань,
		// используется на первой итерации стационарного солвера.
        if (!bG) { 
	          ug=fgplus*VG+(1.0-fgplus)*VP;
	    }
	    else {
	          ug=VG;
	    }

		// конвективный поток через грань КО.
	    Fg=rhog*ug*dS;
	}
	else
	{

		// формулы справедливы и для граничных узлов тоже.
        doublereal taug=tau[VGid][iP]*tau[VGid][iG]/(fgplus*tau[VGid][iG]+(1.0-fgplus)*tau[VGid][iP]); // псевдовремя на грани КО.
		//28.07.2016
		if (bG) {
			taug = tau[VGid][iG];
		}
	    doublereal gradP=0.0, gradG=0.0;
		switch (iNODEid) {
		case 1:
			switch (G) {
			case ESIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[ESIDE][iP].iNODE1]; break;
			case WSIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[WSIDE][iP].iNODE1]; break;
			case NSIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[NSIDE][iP].iNODE1]; break;
			case SSIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[SSIDE][iP].iNODE1]; break;
			case TSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[TSIDE][iP].iNODE1]; break;
			case BSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[BSIDE][iP].iNODE1]; break;
			}
			break;
		case 2:
			switch (G) {
			case ESIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[ESIDE][iP].iNODE2]; break;
			case WSIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[WSIDE][iP].iNODE2]; break;
			case NSIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[NSIDE][iP].iNODE2]; break;
			case SSIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[SSIDE][iP].iNODE2]; break;
			case TSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[TSIDE][iP].iNODE2]; break;
			case BSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[BSIDE][iP].iNODE2]; break;
			}
			break;
		case 3:
			switch (G) {
			case ESIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[ESIDE][iP].iNODE3]; break;
			case WSIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[WSIDE][iP].iNODE3]; break;
			case NSIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[NSIDE][iP].iNODE3]; break;
			case SSIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[SSIDE][iP].iNODE3]; break;
			case TSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[TSIDE][iP].iNODE3]; break;
			case BSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[BSIDE][iP].iNODE3]; break;
			}
			break;
		case 4:
			switch (G) {
			case ESIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[ESIDE][iP].iNODE4]; break;
			case WSIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[WSIDE][iP].iNODE4]; break;
			case NSIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[NSIDE][iP].iNODE4]; break;
			case SSIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[SSIDE][iP].iNODE4]; break;
			case TSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[TSIDE][iP].iNODE4]; break;
			case BSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[BSIDE][iP].iNODE4]; break;
			}
			break;
		default:
			switch (G) {
			case ESIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[ESIDE][iP].iNODE1]; break;
			case WSIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][sosedi[WSIDE][iP].iNODE1]; break;
			case NSIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[NSIDE][iP].iNODE1]; break;
			case SSIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][sosedi[SSIDE][iP].iNODE1]; break;
			case TSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[TSIDE][iP].iNODE1]; break;
			case BSIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][sosedi[BSIDE][iP].iNODE1]; break;
			}
			break;
		}

	    // это делается лишь для того чтобы добавить дополнительную
	    // нижнюю релаксацию так как это написано у I. Sezai.
	    if (!bG) { 
	       ug=fgplus*VG+(1.0-fgplus)*VP;
	       if (bRhieChowi) {
	           //ug+=RCh*ugRhieChow_internal(iP, G, alpha, nvtx, sosedi, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу

			   
			   // Пример того как делать неправильно: 
			   //ug+=RCh*(taug/rhog)*(-(potent[PRESS][iG]-potent[PRESS][iP])/drg+(1.0-fgplus)*gradP+fgplus*gradG);
			   // Вышенаписанная строчка неправильная так как неправильно вычисляется градиент давления на грани potent[PRESS][iG]-potent[PRESS][iP].
			   // Запись (potent[PRESS][iG]-potent[PRESS][iP]) верна только для направлений E, N , T. НО НЕВЕРНА ДЛЯ WSIDE, S, B.

			   // Правильный вариант поправки Рхи-Чоу.
			   switch (G) {// Вклад поправки Рхи-Чоу
			      case ESIDE : case NSIDE : case TSIDE : ug+=RCh*(taug/rhog)*(-(potent[PRESS][iG]-potent[PRESS][iP])/drg+(1.0-fgplus)*gradP+fgplus*gradG); 
				                           break;
			      case WSIDE : case SSIDE :case BSIDE : ug+=RCh*(taug/rhog)*(-(potent[PRESS][iP]-potent[PRESS][iG])/drg+(1.0-fgplus)*gradP+fgplus*gradG);
				                           break;
			   }
			   
			   
	       }
	    }
	    else {
	       ug=VG;
	       
	       // поправки для граничного узла нету. Т.е. считаем её нулевой.
		   // но с другой стороны при тестировании обнаруживаем нефизичные осцилляции вблизи выходной границы.
	       if (bRhieChowb) {
	          //set additional ug+=0.0;
	          //ug+=RCh*ugRhieChow_internal_border(iP, G, alpha, nvtx, sosedi, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу

			   // в данном случае drg==0.5*dr; fgplus==0.5;

              // Правильный вариант поправки Рхи-Чоу.
			   switch (G) {// Вклад поправки Рхи-Чоу
			      case ESIDE : case NSIDE : case TSIDE : ug+=RCh*(taug/rhog)*(-(potent[PRESS][iG]-potent[PRESS][iP])/drg+(1.0-fgplus)*gradP+fgplus*gradG); 
				                           break;
			      case WSIDE : case SSIDE :case BSIDE : ug+=RCh*(taug/rhog)*(-(potent[PRESS][iP]-potent[PRESS][iG])/drg+(1.0-fgplus)*gradP+fgplus*gradG);
				                           break;
			   }

	       }
	       
	    }

	    // конвективный поток через грань КО.
	    Fg=rhog*ug*dS;

 	    // Внимание здесь предполагается что плотность строго постоянна.
	    if (btimedepend) {
			// в нестационарном случае к поправке Рхи-Чоу требуется сделать добавку.

		    doublereal VGoldtimestep=0.0;
		    doublereal VPoldtimestep=0.0;

		    switch (G) {
		       case ESIDE : case WSIDE : VPoldtimestep=speedoldtimestep[VX][iP];
			                     VGoldtimestep=speedoldtimestep[VX][iG];
			                     break;
		       case NSIDE : case SSIDE :VPoldtimestep=speedoldtimestep[VY][iP];
			                     VGoldtimestep=speedoldtimestep[VY][iG];
			                     break;
		       case TSIDE :  case BSIDE : VPoldtimestep=speedoldtimestep[VZ][iP];
			                      VGoldtimestep=speedoldtimestep[VZ][iG];
			                     break;
		    }

		    // Скорость на грани КО с предыдущего временного слоя:
		    doublereal ugold=fgplus*VGoldtimestep+(1.0-fgplus)*VPoldtimestep; 


		    // диагональный коэффициент матрицы должен быть взят на границе контрольного объёма.
		    if (!bG) {
				if (bRhieChowi) {
			       // внутренний контрольный объём
		           Fg+=(taug/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
				}
		    }
		    else {
				if (bRhieChowb) {
			       // граничный контрольный объём
			       Fg+=(taug/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
				}
		    }

    	}
	}

	return Fg;

} // calcFg3


// возвращает скорректированный массовый поток.
// скоректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
// скорректированной скорости и давления, а также сохранённых диагональных 
// коэффициентов матрицы СЛАУ для компонент скорости.
void return_calc_correct_mass_flux(integer iP, doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
	integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm, doublereal **diag_coef,
						doublereal* alpha, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, bool bsimplelinearinterpol,
						doublereal** SpeedCorOld, doublereal *mfold)
{

	// SpeedCorOld - скоректированная скорость на предыдущей итерации.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border.

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

	// плотность аппроксимируется средним гармоническим
	doublereal rhoe, rhow, rhon, rhos, rhot, rhob;
	doublereal rP, rE, rN, rT, rW, rS, rB;

    rP=prop[RHO][iP];
	if (!bE) rE=prop[RHO][iE]; else rE=prop_b[RHO][iE-maxelm];
    if (!bN) rN=prop[RHO][iN]; else rN=prop_b[RHO][iN-maxelm];
    if (!bT) rT=prop[RHO][iT]; else rT=prop_b[RHO][iT-maxelm];
	if (!bW) rW=prop[RHO][iW]; else rW=prop_b[RHO][iW-maxelm];
    if (!bS) rS=prop[RHO][iS]; else rS=prop_b[RHO][iS-maxelm];
    if (!bB) rB=prop[RHO][iB]; else rB=prop_b[RHO][iB-maxelm];

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (!bE) rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP); else rhoe=rE;
	if (!bW) rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP); else rhow=rW;
	if (!bN) rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP); else rhon=rN;
	if (!bS) rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP); else rhos=rS;
    if (!bT) rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP); else rhot=rT;
	if (!bB) rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP); else rhob=rB;


	doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 

	Fe=calcFg(bE, feplus, potent[VX][iE], potent[VX][iP],  
	          iP, ESIDE, alpha[VX], rhoe, 
		      dy*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VX][iE], diag_coef[VX][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb, bsimplelinearinterpol);

	Fw=calcFg(bW, fwplus, potent[VX][iW], potent[VX][iP],  
	          iP, WSIDE, alpha[VX], rhow, 
		      dy*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VX][iW], diag_coef[VX][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb, bsimplelinearinterpol);

	Fn=calcFg(bN, fnplus, potent[VY][iN], potent[VY][iP],  
	          iP, NSIDE, alpha[VY], rhon, 
	          dx*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VY][iN], diag_coef[VY][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb, bsimplelinearinterpol);

	Fs=calcFg(bS, fsplus, potent[VY][iS], potent[VY][iP],  
	          iP, SSIDE, alpha[VY], rhos, 
		      dx*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VY][iS], diag_coef[VY][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb, bsimplelinearinterpol);

    Ft=calcFg(bT, ftplus, potent[VZ][iT], potent[VZ][iP],  
	          iP, TSIDE, alpha[VZ], rhot, 
	          dx*dy, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VZ][iT], diag_coef[VZ][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb, bsimplelinearinterpol);

	Fb=calcFg(bB, fbplus, potent[VZ][iB], potent[VZ][iP],  
	          iP, BSIDE, alpha[VZ], rhob, 
		      dx*dy, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VZ][iB], diag_coef[VZ][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb, bsimplelinearinterpol);

	bool ISezai=true;
	if (ISezai) {

	    doublereal SpeedCorOlde, SpeedCorOldw,  SpeedCorOldn,  SpeedCorOlds,  SpeedCorOldt,  SpeedCorOldb; 
	    if (!bE) { 
		   SpeedCorOlde=feplus*SpeedCorOld[VX][iE]+(1.0-feplus)*SpeedCorOld[VX][iP];
	    }
	    else {
	       SpeedCorOlde=SpeedCorOld[VX][iE];
	    }
	    if (!bN) { 
		   SpeedCorOldn=fnplus*SpeedCorOld[VY][iN]+(1.0-fnplus)*SpeedCorOld[VY][iP];
	    }
	    else {
	       SpeedCorOldn=SpeedCorOld[VY][iN];
	    }
	    if (!bT) { 
		   SpeedCorOldt=ftplus*SpeedCorOld[VZ][iT]+(1.0-ftplus)*SpeedCorOld[VZ][iP];
	    }
	    else {
	       SpeedCorOldt=SpeedCorOld[VZ][iT];
	    }
	    if (!bW) { 
		   SpeedCorOldw=fwplus*SpeedCorOld[VX][iW]+(1.0-fwplus)*SpeedCorOld[VX][iP];
     	}
        else {
	       SpeedCorOldw=SpeedCorOld[VX][iW];
	    }
	    if (!bS) { 
		   SpeedCorOlds=fsplus*SpeedCorOld[VY][iS]+(1.0-fsplus)*SpeedCorOld[VY][iP];
	    }
	    else {
	       SpeedCorOlds=SpeedCorOld[VY][iS];
	    }
	    if (!bB) { 
		   SpeedCorOldb=fbplus*SpeedCorOld[VZ][iB]+(1.0-fbplus)*SpeedCorOld[VZ][iP];
	    }
	    else {
	       SpeedCorOldb=SpeedCorOld[VZ][iB];
	    }

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    //mfcurrentretune[ESIDE]=Fe+(1.0-alpha[VX])*(mfcurrentretune[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    //mfcurrentretune[NSIDE]=Fn+(1.0-alpha[VY])*(mfcurrentretune[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    //mfcurrentretune[TSIDE]=Ft+(1.0-alpha[VZ])*(mfcurrentretune[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    //mfcurrentretune[WSIDE]=Fw+(1.0-alpha[VX])*(mfcurrentretune[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    //mfcurrentretune[SSIDE]=Fs+(1.0-alpha[VY])*(mfcurrentretune[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    //mfcurrentretune[BSIDE]=Fb+(1.0-alpha[VZ])*(mfcurrentretune[BSIDE]-rhob*SpeedCorOldb*dx*dy);

	    // mfold - значение массового потока с предыдущей итерации.
        mfcurrentretune[ESIDE]=Fe+(1.0-alpha[VX])*(mfold[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    mfcurrentretune[NSIDE]=Fn+(1.0-alpha[VY])*(mfold[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    mfcurrentretune[TSIDE]=Ft+(1.0-alpha[VZ])*(mfold[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    mfcurrentretune[WSIDE]=Fw+(1.0-alpha[VX])*(mfold[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    mfcurrentretune[SSIDE]=Fs+(1.0-alpha[VY])*(mfold[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    mfcurrentretune[BSIDE]=Fb+(1.0-alpha[VZ])*(mfold[BSIDE]-rhob*SpeedCorOldb*dx*dy);
	}

	if (!ISezai) {
	   mfcurrentretune[ESIDE]=Fe;
	   mfcurrentretune[NSIDE]=Fn;
	   mfcurrentretune[TSIDE]=Ft;
	   mfcurrentretune[WSIDE]=Fw;
	   mfcurrentretune[SSIDE]=Fs;
	   mfcurrentretune[BSIDE]=Fb;
	}

} // return_correct_mass_flux

// возвращает скорректированный массовый поток.
// скоректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
// скорректированной скорости и давления, а также сохранённых диагональных 
// коэффициентов матрицы СЛАУ для компонент скорости.
// В данной реализации учитывается сглаженное псевдовремя tau.
// реализовано 22 июня 2012 года. Основывается на сглаженном псевдовремени tau.
void return_calc_correct_mass_flux2(integer iP, doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
	integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm, doublereal* alpha, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, bool bsimplelinearinterpol,
						doublereal** SpeedCorOld, doublereal *mfold, doublereal* tau)
{

	// SpeedCorOld - скоректированная скорость на предыдущей итерации.
	// tau - сглаженное псевдовремя.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border.

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

	// плотность аппроксимируется средним гармоническим
	doublereal rhoe, rhow, rhon, rhos, rhot, rhob;
	doublereal rP, rE, rN, rT, rW, rS, rB;

    rP=prop[RHO][iP];
	if (!bE) rE=prop[RHO][iE]; else rE=prop_b[RHO][iE-maxelm];
    if (!bN) rN=prop[RHO][iN]; else rN=prop_b[RHO][iN-maxelm];
    if (!bT) rT=prop[RHO][iT]; else rT=prop_b[RHO][iT-maxelm];
	if (!bW) rW=prop[RHO][iW]; else rW=prop_b[RHO][iW-maxelm];
    if (!bS) rS=prop[RHO][iS]; else rS=prop_b[RHO][iS-maxelm];
    if (!bB) rB=prop[RHO][iB]; else rB=prop_b[RHO][iB-maxelm];

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (!bE) rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP); else rhoe=rE;
	if (!bW) rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP); else rhow=rW;
	if (!bN) rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP); else rhon=rN;
	if (!bS) rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP); else rhos=rS;
    if (!bT) rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP); else rhot=rT;
	if (!bB) rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP); else rhob=rB;


	doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 

	Fe=calcFg2(bE, feplus, iP, ESIDE,  rhoe, 
			   dy*dz, btimedepend, 
			   speedoldtimestep,  mfoldtimestep,
		       dtimestep,	 RCh, sosedi, 
			   bRhieChowi, bRhieChowb, 
			   bsimplelinearinterpol, tau, dxe, potent);

	Fw=calcFg2(bW, fwplus, iP, WSIDE, rhow, 
			   dy*dz, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep,	RCh, sosedi, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dxw, potent);

	Fn=calcFg2(bN, fnplus, iP, NSIDE, rhon, 
			   dx*dz, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep, RCh, sosedi, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dyn, potent);

	Fs=calcFg2(bS, fsplus, iP, SSIDE, rhos, 
			   dx*dz, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep, RCh, sosedi, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dys, potent);

	Ft=calcFg2(bT, ftplus, iP, TSIDE, rhot, 
		       dx*dy, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep,	RCh, sosedi, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dzt, potent);

	Fb=calcFg2(bB, fbplus, iP, BSIDE, rhob, 
			   dx*dy, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep,	RCh, sosedi, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dzb, potent);

	bool ISezai=true;
	if (ISezai) {

	    doublereal SpeedCorOlde, SpeedCorOldw,  SpeedCorOldn,  SpeedCorOlds,  SpeedCorOldt,  SpeedCorOldb; 
	    if (!bE) { 
		   SpeedCorOlde=feplus*SpeedCorOld[VX][iE]+(1.0-feplus)*SpeedCorOld[VX][iP];
	    }
	    else {
	       SpeedCorOlde=SpeedCorOld[VX][iE];
	    }
	    if (!bN) { 
		   SpeedCorOldn=fnplus*SpeedCorOld[VY][iN]+(1.0-fnplus)*SpeedCorOld[VY][iP];
	    }
	    else {
	       SpeedCorOldn=SpeedCorOld[VY][iN];
	    }
	    if (!bT) { 
		   SpeedCorOldt=ftplus*SpeedCorOld[VZ][iT]+(1.0-ftplus)*SpeedCorOld[VZ][iP];
	    }
	    else {
	       SpeedCorOldt=SpeedCorOld[VZ][iT];
	    }
	    if (!bW) { 
		   SpeedCorOldw=fwplus*SpeedCorOld[VX][iW]+(1.0-fwplus)*SpeedCorOld[VX][iP];
     	}
        else {
	       SpeedCorOldw=SpeedCorOld[VX][iW];
	    }
	    if (!bS) { 
		   SpeedCorOlds=fsplus*SpeedCorOld[VY][iS]+(1.0-fsplus)*SpeedCorOld[VY][iP];
	    }
	    else {
	       SpeedCorOlds=SpeedCorOld[VY][iS];
	    }
	    if (!bB) { 
		   SpeedCorOldb=fbplus*SpeedCorOld[VZ][iB]+(1.0-fbplus)*SpeedCorOld[VZ][iP];
	    }
	    else {
	       SpeedCorOldb=SpeedCorOld[VZ][iB];
	    }

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	   // mfold - значение массового потока с предыдущей итерации.
        mfcurrentretune[ESIDE]=Fe+(1.0-alpha[VX])*(mfold[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    mfcurrentretune[NSIDE]=Fn+(1.0-alpha[VY])*(mfold[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    mfcurrentretune[TSIDE]=Ft+(1.0-alpha[VZ])*(mfold[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    mfcurrentretune[WSIDE]=Fw+(1.0-alpha[VX])*(mfold[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    mfcurrentretune[SSIDE]=Fs+(1.0-alpha[VY])*(mfold[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    mfcurrentretune[BSIDE]=Fb+(1.0-alpha[VZ])*(mfold[BSIDE]-rhob*SpeedCorOldb*dx*dy);
	}

	if (!ISezai) {
	   mfcurrentretune[ESIDE]=Fe;
	   mfcurrentretune[NSIDE]=Fn;
	   mfcurrentretune[TSIDE]=Ft;
	   mfcurrentretune[WSIDE]=Fw;
	   mfcurrentretune[SSIDE]=Fs;
	   mfcurrentretune[BSIDE]=Fb;
	}

} // return_correct_mass_flux2

// возвращает скорректированный массовый поток.
// скоректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
// скорректированной скорости и давления, а также сохранённых диагональных 
// коэффициентов матрицы СЛАУ для компонент скорости.
// В данной реализации учитывается сглаженное псевдовремя tau.
// реализовано 23 июня 2012 года. Основывается на сглаженном псевдовремени tau.
void return_calc_correct_mass_flux3(integer iP, doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
	integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm, doublereal* alpha, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, bool bsimplelinearinterpol,
						doublereal** SpeedCorOld, doublereal *mfold, doublereal** tau,
	                    BOUND* sosedb, integer *ilevel_alice, integer* ptr)
{

	// SpeedCorOld - скоректированная скорость на предыдущей итерации.
	// tau - сглаженное псевдовремя.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border.

	// iP - номер центрального контрольного объёма
	integer iE=-1, iN=-1, iT=-1, iW=-1, iS=-1, iB=-1; // номера соседних контрольных объёмов
	iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;


	// 26.09.2016 Добавок для АЛИС сетки.
	integer iE2=-1, iN2=-1, iT2=-1, iW2=-1, iS2=-1, iB2=-1; // номера соседних контрольных объёмов
	integer iE3=-1, iN3=-1, iT3=-1, iW3=-1, iS3=-1, iB3=-1; // номера соседних контрольных объёмов
	integer iE4=-1, iN4=-1, iT4=-1, iW4=-1, iS4=-1, iB4=-1; // номера соседних контрольных объёмов

										  // -1 если не используется и [0..maxelm+maxbound-1] если используется.

	iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
	iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;
	iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
	iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;
	iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
	iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


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

	//printf("%e %e %e %e %e %e\n",feplus, fwplus, fnplus, fsplus, ftplus, fbplus);
	//getchar();

	// плотность аппроксимируется средним гармоническим
	
	// плотность на грани КО аппроксимируется средним гармоническим	
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

	
	doublereal rhoe = 0.0, rhow = 0.0, rhon = 0.0, rhos = 0.0, rhot = 0.0, rhob = 0.0;
	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE > -1) {
		if (!bE) rhoe = rE * rP / (feplus*rE + (1.0 - feplus)*rP); else rhoe = rE; // проверено !
	}
	if (iW > -1) {
		if (!bW) rhow = rW * rP / (fwplus*rW + (1.0 - fwplus)*rP); else rhow = rW;
	}
	if (iN > -1) {
		if (!bN) rhon = rN * rP / (fnplus*rN + (1.0 - fnplus)*rP); else rhon = rN;
	}
	if (iS > -1) {
		if (!bS) rhos = rS * rP / (fsplus*rS + (1.0 - fsplus)*rP); else rhos = rS;
	}
	if (iT > -1) {
		if (!bT) rhot = rT * rP / (ftplus*rT + (1.0 - ftplus)*rP); else rhot = rT;
	}
	if (iB > -1) {
		if (!bB) rhob = rB * rP / (fbplus*rB + (1.0 - fbplus)*rP); else rhob = rB;
	}

	doublereal rhoe2 = 0.0, rhow2 = 0.0, rhon2 = 0.0, rhos2 = 0.0, rhot2 = 0.0, rhob2 = 0.0;
	doublereal rhoe3 = 0.0, rhow3 = 0.0, rhon3 = 0.0, rhos3 = 0.0, rhot3 = 0.0, rhob3 = 0.0;
	doublereal rhoe4 = 0.0, rhow4 = 0.0, rhon4 = 0.0, rhos4 = 0.0, rhot4 = 0.0, rhob4 = 0.0;

	if (iE2 > -1) {
		if (!bE2)  rhoe2 = rE2 * rP / (feplus2*rE2 + (1.0 - feplus2)*rP); else rhoe2 = rE2; // проверено !
	}
	if (iW2 > -1) {
		if (!bW2)  rhow2 = rW2 * rP / (fwplus2*rW2 + (1.0 - fwplus2)*rP); else rhow2 = rW2;
	}
	if (iN2 > -1) {
		if (!bN2) rhon2 = rN2 * rP / (fnplus2*rN2 + (1.0 - fnplus2)*rP); else rhon2 = rN2;
	}
	if (iS2 > -1) {
		if (!bS2)  rhos2 = rS2 * rP / (fsplus2*rS2 + (1.0 - fsplus2)*rP); else rhos2 = rS2;
	}
	if (iT2 > -1) {
		if (!bT2)  rhot2 = rT2 * rP / (ftplus2*rT2 + (1.0 - ftplus2)*rP); else rhot2 = rT2;
	}
	if (iB2 > -1) {
		if (!bB2) rhob2 = rB2 * rP / (fbplus2*rB2 + (1.0 - fbplus2)*rP); else rhob2 = rB2;
	}

	if (iE3 > -1) {
		if (!bE3) rhoe3 = rE3 * rP / (feplus3*rE3 + (1.0 - feplus3)*rP); else rhoe3 = rE3;
	}
	if (iW3 > -1) {
		if (!bW3) rhow3 = rW3 * rP / (fwplus3*rW3 + (1.0 - fwplus3)*rP); else rhow3 = rW3;
	}
	if (iN3 > -1) {
		if (!bN3) rhon3 = rN3 * rP / (fnplus3*rN3 + (1.0 - fnplus3)*rP); else rhon3 = rN3;
	}
	if (iS3 > -1) {
		if (!bS3) rhos3 = rS3 * rP / (fsplus3*rS3 + (1.0 - fsplus3)*rP); else rhos3 = rS3;
	}
	if (iT3 > -1) {
		if (!bT3) rhot3 = rT3 * rP / (ftplus3*rT3 + (1.0 - ftplus3)*rP); else rhot3 = rT3;
	}
	if (iB3 > -1) {
		if (!bB3) rhob3 = rB3 * rP / (fbplus3*rB3 + (1.0 - fbplus3)*rP); else rhob3 = rB3;
	}

	if (iE4 > -1) {
		if (!bE4) rhoe4 = rE4 * rP / (feplus4*rE4 + (1.0 - feplus4)*rP); else rhoe4 = rE4;
	}
	if (iW4 > -1) {
		if (!bW4) rhow4 = rW4 * rP / (fwplus4*rW4 + (1.0 - fwplus4)*rP); else rhow4 = rW4;
	}
	if (iN4 > -1) {
		if (!bN4) rhon4 = rN4 * rP / (fnplus4*rN4 + (1.0 - fnplus4)*rP); else rhon4 = rN4;
	}
	if (iS4 > -1) {
		if (!bS4) rhos4 = rS4 * rP / (fsplus4*rS4 + (1.0 - fsplus4)*rP); else rhos4 = rS4;
	}
	if (iT4 > -1) {
		if (!bT4) rhot4 = rT4 * rP / (ftplus4*rT4 + (1.0 - ftplus4)*rP); else rhot4 = rT4;
	}
	if (iB4 > -1) {
		if (!bB4) rhob4 = rB4 * rP / (fbplus4*rB4 + (1.0 - fbplus4)*rP); else rhob4 = rB4;
	}


	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.

	doublereal Fw = 0.0, Fe = 0.0, Fs = 0.0, Fn = 0.0, Ft = 0.0, Fb = 0.0;


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


		Fe = calcFg3(bE, feplus, iP, ESIDE, rhoe,
			dSqe, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe, potent, 1);
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

		Fw = calcFg3(bW, fwplus, iP, WSIDE, rhow,
			dSqw, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw, potent, 1);
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

		Fn = calcFg3(bN, fnplus, iP, NSIDE, rhon,
			dSqn, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn, potent, 1);
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

		Fs = calcFg3(bS, fsplus, iP, SSIDE, rhos,
			dSqs, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys, potent, 1);
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

		Ft = calcFg3(bT, ftplus, iP, TSIDE, rhot,
			dSqt, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt, potent, 1);
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

		Fb = calcFg3(bB, fbplus, iP, BSIDE, rhob,
			dSqb, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb, potent, 1);
	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal Fw2 = 0.0, Fe2 = 0.0, Fs2 = 0.0, Fn2 = 0.0, Ft2 = 0.0, Fb2 = 0.0;


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

		Fe2 = calcFg3(bE2, feplus2, iP, ESIDE, rhoe2,
			dSqe2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe2, potent, 2);
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

		Fw2 = calcFg3(bW2, fwplus2, iP, WSIDE, rhow2,
			dSqw2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw2, potent, 2);
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

		Fn2 = calcFg3(bN2, fnplus2, iP, NSIDE, rhon2,
			dSqn2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn2, potent, 2);
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

		Fs2 = calcFg3(bS2, fsplus2, iP, SSIDE, rhos2,
			dSqs2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys2, potent, 2);
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

		Ft2 = calcFg3(bT2, ftplus2, iP, TSIDE, rhot2,
			dSqt2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt2, potent, 2);
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

		Fb2 = calcFg3(bB2, fbplus2, iP, BSIDE, rhob2,
			dSqb2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb2, potent, 2);
	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal Fw3 = 0.0, Fe3 = 0.0, Fs3 = 0.0, Fn3 = 0.0, Ft3 = 0.0, Fb3 = 0.0;


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

		Fe3 = calcFg3(bE3, feplus3, iP, ESIDE, rhoe3,
			dSqe3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe3, potent, 3);
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

		Fw3 = calcFg3(bW3, fwplus3, iP, WSIDE, rhow3,
			dSqw3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw3, potent, 3);
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

		Fn3 = calcFg3(bN3, fnplus3, iP, NSIDE, rhon3,
			dSqn3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn3, potent, 3);
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

		Fs3 = calcFg3(bS3, fsplus3, iP, SSIDE, rhos3,
			dSqs3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys3, potent, 3);
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

		Ft3 = calcFg3(bT3, ftplus3, iP, TSIDE, rhot3,
			dSqt3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt3, potent, 3);
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

		Fb3 = calcFg3(bB3, fbplus3, iP, BSIDE, rhob3,
			dSqb3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb3, potent, 3);
	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.
	doublereal Fw4 = 0.0, Fe4 = 0.0, Fs4 = 0.0, Fn4 = 0.0, Ft4 = 0.0, Fb4 = 0.0;


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

		Fe4 = calcFg3(bE4, feplus4, iP, ESIDE, rhoe4,
			dSqe4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe4, potent, 4);
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

		Fw4 = calcFg3(bW4, fwplus4, iP, WSIDE, rhow4,
			dSqw4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw4, potent, 4);
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

		Fn4 = calcFg3(bN4, fnplus4, iP, NSIDE, rhon4,
			dSqn4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn4, potent, 4);
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

		Fs4 = calcFg3(bS4, fsplus4, iP, SSIDE, rhos4,
			dSqs4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys4, potent, 4);
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

		Ft4 = calcFg3(bT4, ftplus4, iP, TSIDE, rhot4,
			dSqt4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt4, potent, 4);
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

		Fb4 = calcFg3(bB4, fbplus4, iP, BSIDE, rhob4,
			dSqb4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb4, potent, 4);
	}

	doublereal Fe_sum = Fe + Fe2 + Fe3 + Fe4;
	doublereal Fw_sum = Fw + Fw2 + Fw3 + Fw4;
	doublereal Fn_sum = Fn + Fn2 + Fn3 + Fn4;
	doublereal Fs_sum = Fs + Fs2 + Fs3 + Fs4;
	doublereal Ft_sum = Ft + Ft2 + Ft3 + Ft4;
	doublereal Fb_sum = Fb + Fb2 + Fb3 + Fb4;



	bool ISezai=true; // true
	if (ISezai) {

		doublereal SpeedCorOlde = 0.0, SpeedCorOldw = 0.0, SpeedCorOldn = 0.0, SpeedCorOlds = 0.0, SpeedCorOldt = 0.0, SpeedCorOldb = 0.0;
		if (iE > -1) {
			if (!bE) {
				SpeedCorOlde = feplus * potent[VXCOR][iE] + (1.0 - feplus)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde = potent[VXCOR][iE];
			}
		}
		if (iN > -1) {
			if (!bN) {
				SpeedCorOldn = fnplus * potent[VYCOR][iN] + (1.0 - fnplus)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn = potent[VYCOR][iN];
			}
		}
		if (iT > -1) {
			if (!bT) {
				SpeedCorOldt = ftplus * potent[VZCOR][iT] + (1.0 - ftplus)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt = potent[VZCOR][iT];
			}
		}
		if (iW > -1) {
			if (!bW) {
				SpeedCorOldw = fwplus * potent[VXCOR][iW] + (1.0 - fwplus)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw = potent[VXCOR][iW];
			}
		}
		if (iS > -1) {
			if (!bS) {
				SpeedCorOlds = fsplus * potent[VYCOR][iS] + (1.0 - fsplus)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds = potent[VYCOR][iS];
			}
		}
		if (iB > -1) {
			if (!bB) {
				SpeedCorOldb = fbplus * potent[VZCOR][iB] + (1.0 - fbplus)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb = potent[VZCOR][iB];
			}
		}

		doublereal SpeedCorOlde2 = 0.0, SpeedCorOldw2 = 0.0, SpeedCorOldn2 = 0.0, SpeedCorOlds2 = 0.0, SpeedCorOldt2 = 0.0, SpeedCorOldb2 = 0.0;
		if (iE2 > -1) {
			if (!bE2) {
				SpeedCorOlde2 = feplus2 * potent[VXCOR][iE2] + (1.0 - feplus2)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde2 = potent[VXCOR][iE2];
			}
		}
		if (iN2 > -1) {
			if (!bN2) {
				SpeedCorOldn2 = fnplus2 * potent[VYCOR][iN2] + (1.0 - fnplus2)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn2 = potent[VYCOR][iN2];
			}
		}
		if (iT2 > -1) {
			if (!bT2) {
				SpeedCorOldt2 = ftplus2 * potent[VZCOR][iT2] + (1.0 - ftplus2)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt2 = potent[VZCOR][iT2];
			}
		}
		if (iW2 > -1) {
			if (!bW2) {
				SpeedCorOldw2 = fwplus2 * potent[VXCOR][iW2] + (1.0 - fwplus2)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw2 = potent[VXCOR][iW2];
			}
		}
		if (iS2 > -1) {
			if (!bS2) {
				SpeedCorOlds2 = fsplus2 * potent[VYCOR][iS2] + (1.0 - fsplus2)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds2 = potent[VYCOR][iS2];
			}
		}
		if (iB2 > -1) {
			if (!bB2) {
				SpeedCorOldb2 = fbplus2 * potent[VZCOR][iB2] + (1.0 - fbplus2)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb2 = potent[VZCOR][iB2];
			}
		}

		doublereal SpeedCorOlde3 = 0.0, SpeedCorOldw3 = 0.0, SpeedCorOldn3 = 0.0, SpeedCorOlds3 = 0.0, SpeedCorOldt3 = 0.0, SpeedCorOldb3 = 0.0;
		if (iE3 > -1) {
			if (!bE3) {
				SpeedCorOlde3 = feplus3 * potent[VXCOR][iE3] + (1.0 - feplus3)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde3 = potent[VXCOR][iE3];
			}
		}
		if (iN3 > -1) {
			if (!bN3) {
				SpeedCorOldn3 = fnplus3 * potent[VYCOR][iN3] + (1.0 - fnplus3)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn3 = potent[VYCOR][iN3];
			}
		}
		if (iT3 > -1) {
			if (!bT3) {
				SpeedCorOldt3 = ftplus3 * potent[VZCOR][iT3] + (1.0 - ftplus3)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt3 = potent[VZCOR][iT3];
			}
		}
		if (iW3 > -1) {
			if (!bW3) {
				SpeedCorOldw3 = fwplus3 * potent[VXCOR][iW3] + (1.0 - fwplus3)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw3 = potent[VXCOR][iW3];
			}
		}
		if (iS3 > -1) {
			if (!bS3) {
				SpeedCorOlds3 = fsplus3 * potent[VYCOR][iS3] + (1.0 - fsplus3)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds3 = potent[VYCOR][iS3];
			}
		}
		if (iB3 > -1) {
			if (!bB3) {
				SpeedCorOldb3 = fbplus3 * potent[VZCOR][iB3] + (1.0 - fbplus3)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb3 = potent[VZCOR][iB3];
			}
		}


		doublereal SpeedCorOlde4 = 0.0, SpeedCorOldw4 = 0.0, SpeedCorOldn4 = 0.0, SpeedCorOlds4 = 0.0, SpeedCorOldt4 = 0.0, SpeedCorOldb4 = 0.0;
		if (iE4 > -1) {
			if (!bE4) {
				SpeedCorOlde4 = feplus4 * potent[VXCOR][iE4] + (1.0 - feplus4)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde4 = potent[VXCOR][iE4];
			}
		}
		if (iN4 > -1) {
			if (!bN4) {
				SpeedCorOldn4 = fnplus4 * potent[VYCOR][iN4] + (1.0 - fnplus4)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn4 = potent[VYCOR][iN4];
			}
		}
		if (iT4 > -1) {
			if (!bT4) {
				SpeedCorOldt4 = ftplus4 * potent[VZCOR][iT4] + (1.0 - ftplus4)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt4 = potent[VZCOR][iT4];
			}
		}
		if (iW4 > -1) {
			if (!bW4) {
				SpeedCorOldw4 = fwplus4 * potent[VXCOR][iW4] + (1.0 - fwplus4)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw4 = potent[VXCOR][iW4];
			}
		}
		if (iS4 > -1) {
			if (!bS4) {
				SpeedCorOlds4 = fsplus4 * potent[VYCOR][iS4] + (1.0 - fsplus4)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds4 = potent[VYCOR][iS4];
			}
		}
		if (iB4 > -1) {
			if (!bB4) {
				SpeedCorOldb4 = fbplus4 * potent[VZCOR][iB4] + (1.0 - fbplus4)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb4 = potent[VZCOR][iB4];
			}
		}

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	   // mfold - значение массового потока с предыдущей итерации.
        mfcurrentretune[ESIDE]=Fe_sum+(1.0-alpha[VX])*(mfold[ESIDE]-rhoe*SpeedCorOlde*dSqe - rhoe2 * SpeedCorOlde2*dSqe2 - rhoe3 * SpeedCorOlde3*dSqe3 - rhoe4 * SpeedCorOlde4*dSqe4);
	    mfcurrentretune[NSIDE]=Fn_sum +(1.0-alpha[VY])*(mfold[NSIDE]-rhon*SpeedCorOldn*dSqn - rhon2 * SpeedCorOldn2*dSqn2 - rhon3 * SpeedCorOldn3*dSqn3 - rhon4 * SpeedCorOldn4*dSqn4);
	    mfcurrentretune[TSIDE]=Ft_sum +(1.0-alpha[VZ])*(mfold[TSIDE]-rhot*SpeedCorOldt*dSqt - rhot2 * SpeedCorOldt2*dSqt2 - rhot3 * SpeedCorOldt3*dSqt3 - rhot4 * SpeedCorOldt4*dSqt4);
	    mfcurrentretune[WSIDE]=Fw_sum +(1.0-alpha[VX])*(mfold[WSIDE]-rhow*SpeedCorOldw*dSqw - rhow2 * SpeedCorOldw2*dSqw2 - rhow3 * SpeedCorOldw3*dSqw3 - rhow4 * SpeedCorOldw4*dSqw4);
	    mfcurrentretune[SSIDE]=Fs_sum +(1.0-alpha[VY])*(mfold[SSIDE]-rhos*SpeedCorOlds*dSqs - rhos2 * SpeedCorOlds2*dSqs2 - rhos3 * SpeedCorOlds3*dSqs3 - rhos4 * SpeedCorOlds4*dSqs4);
	    mfcurrentretune[BSIDE]=Fb_sum +(1.0-alpha[VZ])*(mfold[BSIDE]-rhob*SpeedCorOldb*dSqb - rhob2 * SpeedCorOldb2*dSqb2 - rhob3 * SpeedCorOldb3*dSqb3 - rhob4 * SpeedCorOldb4*dSqb4);
	}

	if (!ISezai) {
	   mfcurrentretune[ESIDE]=Fe_sum;
	   mfcurrentretune[NSIDE]=Fn_sum;
	   mfcurrentretune[TSIDE]=Ft_sum;
	   mfcurrentretune[WSIDE]=Fw_sum;
	   mfcurrentretune[SSIDE]=Fs_sum;
	   mfcurrentretune[BSIDE]=Fb_sum;
	}

} // return_correct_mass_flux3



// Составляет матрицу для уравнения 
// поправки давления
void my_elmatr_quad_PAm(integer iP, equation3D** &sl, equation3D_bon** &slb,  
						doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
						integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm, doublereal **diag_coef,
						doublereal* alpha, doublereal dbeta, doublereal** &rhie_chow, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep) {

    // 31 марта 2012 из этой функции удалено много закоментированного кода. Его можно найти в backup`ах AliceFlowv0_07.

    /* btimedepend - стационарный (false), нестационарный (true).
	*  dtimestep - величина текущего шага по времени.
	*  
	*
	*  mfoldtimestep - массовый поток через каждую грань КО с предыдущего временного слоя. (неизменяемая в данной функции величина).
	*  mfcurrentretune - возвращаемый текущий массовый поток через каждую грань КО.
	*  КО имеет 6 граней : E,N,T,W,S,B.
	*
	* speedoldtimestep[VX or VY or VZ][iP : 0..maxelm+maxbound-1] - поле скорости с предыдущего временного слоя.
	*
	*/

	bool bhighorder=false;

	doublereal eps=admission; // для отделения вещественного нуля.
	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;
	sl[PAM][iP].iP=iP;
	sl[PAM][iP].iE=iE; sl[PAM][iP].iN=iN; 
	sl[PAM][iP].iS=iS; sl[PAM][iP].iW=iW;
    sl[PAM][iP].iT=iT; sl[PAM][iP].iB=iB;


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
		

	doublereal apue=1.0, apuw=1.0, apvn=1.0, apvs=1.0, apwt=1.0, apwb=1.0;
	if (!bE) apue=sl[VX][iE].ap*sl[VX][iP].ap/(feplus*sl[VX][iE].ap+(1-feplus)*sl[VX][iP].ap); else apue=slb[VX][iE-maxelm].aw;
	if (!bW) apuw=sl[VX][iW].ap*sl[VX][iP].ap/(fwplus*sl[VX][iW].ap+(1-fwplus)*sl[VX][iP].ap); else apuw=slb[VX][iW-maxelm].aw;
	if (!bN) apvn=sl[VY][iN].ap*sl[VY][iP].ap/(fnplus*sl[VY][iN].ap+(1-fnplus)*sl[VY][iP].ap); else apvn=slb[VY][iN-maxelm].aw;
	if (!bS) apvs=sl[VY][iS].ap*sl[VY][iP].ap/(fsplus*sl[VY][iS].ap+(1-fsplus)*sl[VY][iP].ap); else apvs=slb[VY][iS-maxelm].aw;
	if (!bT) apwt=sl[VZ][iT].ap*sl[VZ][iP].ap/(ftplus*sl[VZ][iT].ap+(1-ftplus)*sl[VZ][iP].ap); else apwt=slb[VZ][iT-maxelm].aw;
	if (!bB) apwb=sl[VZ][iB].ap*sl[VZ][iP].ap/(fbplus*sl[VZ][iB].ap+(1-fbplus)*sl[VZ][iP].ap); else apwb=slb[VZ][iB-maxelm].aw;

	doublereal de=dx, dw=dx, dn=dy, ds=dy, dt=dz, db=dz;
	// Так рекомендует делать С. Патанкар, надо только
	// помнить что это должно быть согласовано с этой-же
	// величиной при коррекции компонент скорости.
	//de=dy*dz/apue; dw=dy*dz/apuw; 
	//dn=dx*dz/apvn; ds=dx*dz/apvs;
	//dt=dx*dy/apwt; db=dx*dy/apwb;
	// Так рекомендуют делать в следующих статьях:
	//1. В.В.Винников, Д.Л.Ревизников
	// Применение декартовых сеток для решения уравнений Навье-Стокса
	// в областях с криволинейными границами. 
	// МАИ (ГТУ) Математическое Моделирование 2005г, том 17, номер 8, стр 15-30
	// 2. SIMPLE METHOD FOR THE SOLUTION OF INCOMPRESSIBLE FLOWS ON NON-STAGGERED GRIDS
	// I.Sezai - Eastern Mediterranean University. Revised January, 2011.
	// 3. Гаврилов Андрей sigma-flow.
	// Такие значения также должны быть согласованы с теми что используются при коррекции скоростей.
	if (iSIMPLE_alg==SIMPLE_Carretto) {
		// SIMPLE
		de=alpha[VX]*dy*dz/apue; dw=alpha[VX]*dy*dz/apuw; 
	    dn=alpha[VY]*dx*dz/apvn; ds=alpha[VY]*dx*dz/apvs;
	    dt=alpha[VZ]*dx*dy/apwt; db=alpha[VZ]*dx*dy/apwb;
	}
	if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC
		de=alpha[VX]*dy*dz/((1.0-alpha[VX])*apue); dw=alpha[VX]*dy*dz/((1.0-alpha[VX])*apuw); 
	    dn=alpha[VY]*dx*dz/((1.0-alpha[VY])*apvn); ds=alpha[VY]*dx*dz/((1.0-alpha[VY])*apvs);
	    dt=alpha[VZ]*dx*dy/((1.0-alpha[VZ])*apwt); db=alpha[VZ]*dx*dy/((1.0-alpha[VZ])*apwb);
	}

    // плотность аппроксимируется средним гармоническим
	doublereal rhoe, rhow, rhon, rhos, rhot, rhob;
	doublereal rP, rE, rN, rT, rW, rS, rB;

    rP=prop[RHO][iP];
	if (!bE) rE=prop[RHO][iE]; else rE=prop_b[RHO][iE-maxelm];
    if (!bN) rN=prop[RHO][iN]; else rN=prop_b[RHO][iN-maxelm];
    if (!bT) rT=prop[RHO][iT]; else rT=prop_b[RHO][iT-maxelm];
	if (!bW) rW=prop[RHO][iW]; else rW=prop_b[RHO][iW-maxelm];
    if (!bS) rS=prop[RHO][iS]; else rS=prop_b[RHO][iS-maxelm];
    if (!bB) rB=prop[RHO][iB]; else rB=prop_b[RHO][iB-maxelm];

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (!bE) rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP); else rhoe=rE;
	if (!bW) rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP); else rhow=rW;
	if (!bN) rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP); else rhon=rN;
	if (!bS) rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP); else rhos=rS;
    if (!bT) rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP); else rhot=rT;
	if (!bB) rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP); else rhob=rB;


    // 0.5 - параметр релаксации для скоростей (см Гаврилов Андрей)
	// уравнение поправки давления.

	doublereal De, Dw, Ds, Dn, Dt, Db; // диффузионный поток через грани КО.
	

	if (!bE) {
		if (bW) De=dbeta*rhoe*de*dy*dz;
		else De=rhoe*de*dy*dz;
	}
	else De=dbeta*rhoe*de*dy*dz;

	if (!bW) {
		if (bE) Dw=dbeta*rhow*dw*dy*dz;
		else Dw=rhow*dw*dy*dz; 
	}
	else Dw=dbeta*rhow*dw*dy*dz;

	if (!bN) {
		if (bS) Dn=dbeta*rhon*dn*dx*dz;
		else Dn=rhon*dn*dx*dz; 
	}
	else Dn=dbeta*rhon*dn*dx*dz;


	if (!bS) {
		if (bN) Ds=dbeta*rhos*ds*dx*dz;
		else Ds=rhos*ds*dx*dz;
	}
	else Ds=dbeta*rhos*ds*dx*dz;


	if (!bT) {
		if (bB) Dt=dbeta*rhot*dt*dx*dy;
		else Dt=rhot*dt*dx*dy;
	} 
	else Dt=dbeta*rhot*dt*dx*dy;


	if (!bB) {
		if (bT) Db=dbeta*rhob*db*dx*dy;
		else Db=rhob*db*dx*dy;
	}
	else Db=dbeta*rhob*db*dx*dy;

	// Число Пекле равно нулю.
	sl[PAM][iP].ae=De*1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[PAM][iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[PAM][iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[PAM][iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[PAM][iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[PAM][iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	sl[PAM][iP].ap=sl[PAM][iP].ae+sl[PAM][iP].aw+sl[PAM][iP].an+sl[PAM][iP].as+sl[PAM][iP].at+sl[PAM][iP].ab;

	doublereal baddDFLUX2=0.0;
	if (bhighorder) {
		// если bborder == false то узел строго внутренний.
		// если bborder == true то мы находимся вблизи граничного узла.
		bool bborder=false;
		doublereal myflux=0.0;
		myflux=De*(dxe*DFDXiP(potent[PAM], iP, ESIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iE]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
		myflux=Dw*(-dxw*DFDXiP(potent[PAM], iP, WSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iW]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dn*(dyn*DFDXiP(potent[PAM], iP, NSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iN]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
        myflux=Ds*(-dys*DFDXiP(potent[PAM], iP, SSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iS]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dt*(dzt*DFDXiP(potent[PAM], iP, TSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iT]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Db*(-dzb*DFDXiP(potent[PAM], iP, BSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iB]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
	}



    doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 
	     
	
	Fe=calcFg(bE, feplus, potent[VX][iE], potent[VX][iP],  
	          iP, ESIDE, alpha[VX], rhoe, 
		      dy*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VX][iE], diag_coef[VX][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fw=calcFg(bW, fwplus, potent[VX][iW], potent[VX][iP],  
	          iP, WSIDE, alpha[VX], rhow, 
		      dy*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VX][iW], diag_coef[VX][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fn=calcFg(bN, fnplus, potent[VY][iN], potent[VY][iP],  
	          iP, NSIDE, alpha[VY], rhon, 
		      dx*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VY][iN], diag_coef[VY][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fs=calcFg(bS, fsplus, potent[VY][iS], potent[VY][iP],  
	          iP, SSIDE, alpha[VY], rhos, 
		      dx*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VY][iS], diag_coef[VY][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

    Ft=calcFg(bT, ftplus, potent[VZ][iT], potent[VZ][iP],  
	          iP, TSIDE, alpha[VZ], rhot, 
		      dx*dy, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VZ][iT], diag_coef[VZ][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fb=calcFg(bB, fbplus, potent[VZ][iB], potent[VZ][iP],  
	          iP, BSIDE, alpha[VZ], rhob, 
		      dx*dy, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VZ][iB], diag_coef[VZ][iP], dtimestep,
		      RCh, nvtx, sosedi, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	
	bool ISezai=true;
	if (ISezai) {

		doublereal SpeedCorOlde, SpeedCorOldw,  SpeedCorOldn,  SpeedCorOlds,  SpeedCorOldt,  SpeedCorOldb; 
	    if (!bE) { 
		   SpeedCorOlde=feplus*potent[VXCOR][iE]+(1.0-feplus)*potent[VXCOR][iP];
	    }
	    else {
	       SpeedCorOlde=potent[VXCOR][iE];
	    }
	    if (!bN) { 
		   SpeedCorOldn=fnplus*potent[VYCOR][iN]+(1.0-fnplus)*potent[VYCOR][iP];
	    }
	    else {
	       SpeedCorOldn=potent[VYCOR][iN];
	    }
	    if (!bT) { 
		   SpeedCorOldt=ftplus*potent[VZCOR][iT]+(1.0-ftplus)*potent[VZCOR][iP];
	    }
	    else {
	       SpeedCorOldt=potent[VZCOR][iT];
	    }
	    if (!bW) { 
		   SpeedCorOldw=fwplus*potent[VXCOR][iW]+(1.0-fwplus)*potent[VXCOR][iP];
	    }
	    else {
	       SpeedCorOldw=potent[VXCOR][iW];
	    }
	    if (!bS) { 
		   SpeedCorOlds=fsplus*potent[VYCOR][iS]+(1.0-fsplus)*potent[VYCOR][iP];
	    }
	    else {
	       SpeedCorOlds=potent[VYCOR][iS];
    	}
	    if (!bB) { 
		   SpeedCorOldb=fbplus*potent[VZCOR][iB]+(1.0-fbplus)*potent[VZCOR][iP];
	    }
	    else {
	       SpeedCorOldb=potent[VZCOR][iB];
	    }

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    // дополнительная поправка осуществляется на основе скоректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    //mfcurrentretune[ESIDE]=Fe+(1.0-alpha[VX])*(mfcurrentretune[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    //mfcurrentretune[NSIDE]=Fn+(1.0-alpha[VY])*(mfcurrentretune[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    //mfcurrentretune[TSIDE]=Ft+(1.0-alpha[VZ])*(mfcurrentretune[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    //mfcurrentretune[WSIDE]=Fw+(1.0-alpha[VX])*(mfcurrentretune[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    //mfcurrentretune[SSIDE]=Fs+(1.0-alpha[VY])*(mfcurrentretune[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    //mfcurrentretune[BSIDE]=Fb+(1.0-alpha[VZ])*(mfcurrentretune[BSIDE]-rhob*SpeedCorOldb*dx*dy);

		// возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    // дополнительная поправка осуществляется на основе скоректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    Fe+=(1.0-alpha[VX])*(mfcurrentretune[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    Fn+=(1.0-alpha[VY])*(mfcurrentretune[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    Ft+=(1.0-alpha[VZ])*(mfcurrentretune[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    Fw+=(1.0-alpha[VX])*(mfcurrentretune[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    Fs+=(1.0-alpha[VY])*(mfcurrentretune[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    Fb+=(1.0-alpha[VZ])*(mfcurrentretune[BSIDE]-rhob*SpeedCorOldb*dx*dy);

	} 
	

	
	// возвращаем значение потока на грани КО.
	mfcurrentretune[ESIDE]=Fe;
	mfcurrentretune[NSIDE]=Fn;
	mfcurrentretune[TSIDE]=Ft;
	mfcurrentretune[WSIDE]=Fw;
	mfcurrentretune[SSIDE]=Fs;
	mfcurrentretune[BSIDE]=Fb;
	
	

	sl[PAM][iP].b=Fw-Fe+Fs-Fn+Fb-Ft+baddDFLUX2;


	// Симметризация СЛАУ :

	// для поправки давления получается эллиптическое уравнение с SPD матрицей.
	
	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

		// 1. Учёт условия Дирихле:
		if ((iE>=maxelm) && (fabs(slb[PAM][iE-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].ae*slb[PAM][iE-maxelm].b/slb[PAM][iE-maxelm].aw;
			sl[PAM][iP].ae=0.0;
			sl[PAM][iP].iE=-1; // не входит в матрицу СЛАУ.
		}
		if ((iW>=maxelm) && (fabs(slb[PAM][iW-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].aw*slb[PAM][iW-maxelm].b/slb[PAM][iW-maxelm].aw;
			sl[PAM][iP].aw=0.0;
			sl[PAM][iP].iW=-1; // не входит в матрицу СЛАУ.
		}
		if ((iN>=maxelm) && (fabs(slb[PAM][iN-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].an*slb[PAM][iN-maxelm].b/slb[PAM][iN-maxelm].aw;
			sl[PAM][iP].an=0.0;
			sl[PAM][iP].iN=-1; // не входит в матрицу СЛАУ.
		}
		if ((iS>=maxelm) && (fabs(slb[PAM][iS-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].as*slb[PAM][iS-maxelm].b/slb[PAM][iS-maxelm].aw;
			sl[PAM][iP].as=0.0;
			sl[PAM][iP].iS=-1; // не входит в матрицу СЛАУ.
		}
		if ((iT>=maxelm) && (fabs(slb[PAM][iT-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].at*slb[PAM][iT-maxelm].b/slb[PAM][iT-maxelm].aw;
			sl[PAM][iP].at=0.0;
			sl[PAM][iP].iT=-1; // не входит в матрицу СЛАУ.
		}
		if ((iB>=maxelm) && (fabs(slb[PAM][iB-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].ab*slb[PAM][iB-maxelm].b/slb[PAM][iB-maxelm].aw;
			sl[PAM][iP].ab=0.0;
			sl[PAM][iP].iB=-1; // не входит в матрицу СЛАУ.
		}


} // my_elmatr_quad_PAm

// Составляет матрицу для уравнения 
// поправки давления
// Данная сборка матрицы для поправки давления основана на сглаженном псевдовремени,
// так как рекомендует Гаврилов Андрей. Даже если сглаженное псевдовремя не поможет,
// всё равно этот способ проще и привлекательней т.к. содержит меньше деталей все детали спрятаны 
// в заранее вычисленный вектор tau (псевдовремя).
// реализовано 22 июня 2012 года.
void my_elmatr_quad_PAm2(integer iP, equation3D** &sl, equation3D_bon** &slb,  
						doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
						integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm,
						doublereal* alpha, doublereal dbeta, doublereal** &rhie_chow, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, doublereal* tau) {

    // 31 марта 2012 из этой функции удалено много закоментированного кода. Его можно найти в backup`ах AliceFlowv0_07.

	// tau - заранее вычисленное псевдовремя, оно служит коэффициентом диффузии в эллиптическом операторе, к тому-же оно
	// также входит в монотонизирующую поправку Рхи-Чоу (1983 года).

    /* btimedepend - стационарный (false), нестационарный (true).
	*  dtimestep - величина текущего шага по времени.
	*  
	*
	*  mfoldtimestep - массовый поток через каждую грань КО с предыдущего временного слоя. (неизменяемая в данной функции величина).
	*  mfcurrentretune - возвращаемый текущий массовый поток через каждую грань КО.
	*  КО имеет 6 граней : E,N,T,W,S,B.
	*
	* speedoldtimestep[VX or VY or VZ][iP : 0..maxelm+maxbound-1] - поле скорости с предыдущего временного слоя.
	*
	*/

	bool bhighorder=false;

	doublereal eps=admission; // для отделения вещественного нуля.
	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;
	sl[PAM][iP].iP=iP;
	sl[PAM][iP].iE=iE; sl[PAM][iP].iN=iN; 
	sl[PAM][iP].iS=iS; sl[PAM][iP].iW=iW;
    sl[PAM][iP].iT=iT; sl[PAM][iP].iB=iB;


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
	
	// Значение коэффициента диффузии (псевдовремени) на грани контрольного объёма.
	doublereal taue=1.0, tauw=1.0, taun=1.0, taus=1.0, taut=1.0, taub=1.0;
	if (!bE) taue=tau[iE]*tau[iP]/(feplus*tau[iE]+(1-feplus)*tau[iP]); else taue=tau[iE];
	if (!bW) tauw=tau[iW]*tau[iP]/(fwplus*tau[iW]+(1-fwplus)*tau[iP]); else tauw=tau[iW];
	if (!bN) taun=tau[iN]*tau[iP]/(fnplus*tau[iN]+(1-fnplus)*tau[iP]); else taun=tau[iN];
	if (!bS) taus=tau[iS]*tau[iP]/(fsplus*tau[iS]+(1-fsplus)*tau[iP]); else taus=tau[iS];
	if (!bT) taut=tau[iT]*tau[iP]/(ftplus*tau[iT]+(1-ftplus)*tau[iP]); else taut=tau[iT];
	if (!bB) taub=tau[iB]*tau[iP]/(fbplus*tau[iB]+(1-fbplus)*tau[iP]); else taub=tau[iB];

	
	/*
	doublereal de, dw, dn, ds, dt, db;
	// Так рекомендует делать С. Патанкар, надо только
	// помнить что это должно быть согласовано с этой-же
	// величиной при коррекции компонент скорости.
	//de=dy*dz/apue; dw=dy*dz/apuw; 
	//dn=dx*dz/apvn; ds=dx*dz/apvs;
	//dt=dx*dy/apwt; db=dx*dy/apwb;
	// Так рекомендуют делать в следующих статьях:
	//1. В.В.Винников, Д.Л.Ревизников
	// Применение декартовых сеток для решения уравнений Навье-Стокса
	// в областях с криволинейными границами. 
	// МАИ (ГТУ) Математическое Моделирование 2005г, том 17, номер 8, стр 15-30
	// 2. SIMPLE METHOD FOR THE SOLUTION OF INCOMPRESSIBLE FLOWS ON NON-STAGGERED GRIDS
	// I.Sezai - Eastern Mediterranean University. Revised January, 2011.
	// 3. Гаврилов Андрей sigma-flow.
	// Такие значения также должны быть согласованы с теми что используются при коррекции скоростей.
	if (iSIMPLE_alg==SIMPLE_Carretto) {
		// SIMPLE
		de=alpha[VX]*dy*dz/apue; dw=alpha[VX]*dy*dz/apuw; 
	    dn=alpha[VY]*dx*dz/apvn; ds=alpha[VY]*dx*dz/apvs;
	    dt=alpha[VZ]*dx*dy/apwt; db=alpha[VZ]*dx*dy/apwb;
	}
	if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC
		de=alpha[VX]*dy*dz/((1.0-alpha[VX])*apue); dw=alpha[VX]*dy*dz/((1.0-alpha[VX])*apuw); 
	    dn=alpha[VY]*dx*dz/((1.0-alpha[VY])*apvn); ds=alpha[VY]*dx*dz/((1.0-alpha[VY])*apvs);
	    dt=alpha[VZ]*dx*dy/((1.0-alpha[VZ])*apwt); db=alpha[VZ]*dx*dy/((1.0-alpha[VZ])*apwb);
	}
	*/
    // плотность аппроксимируется средним гармоническим
	doublereal rhoe, rhow, rhon, rhos, rhot, rhob;
	doublereal rP, rE, rN, rT, rW, rS, rB;

    rP=prop[RHO][iP];
	if (!bE) rE=prop[RHO][iE]; else rE=prop_b[RHO][iE-maxelm];
    if (!bN) rN=prop[RHO][iN]; else rN=prop_b[RHO][iN-maxelm];
    if (!bT) rT=prop[RHO][iT]; else rT=prop_b[RHO][iT-maxelm];
	if (!bW) rW=prop[RHO][iW]; else rW=prop_b[RHO][iW-maxelm];
    if (!bS) rS=prop[RHO][iS]; else rS=prop_b[RHO][iS-maxelm];
    if (!bB) rB=prop[RHO][iB]; else rB=prop_b[RHO][iB-maxelm];

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (!bE) rhoe=rE*rP/(feplus*rE+(1.0-feplus)*rP); else rhoe=rE;
	if (!bW) rhow=rW*rP/(fwplus*rW+(1.0-fwplus)*rP); else rhow=rW;
	if (!bN) rhon=rN*rP/(fnplus*rN+(1.0-fnplus)*rP); else rhon=rN;
	if (!bS) rhos=rS*rP/(fsplus*rS+(1.0-fsplus)*rP); else rhos=rS;
    if (!bT) rhot=rT*rP/(ftplus*rT+(1.0-ftplus)*rP); else rhot=rT;
	if (!bB) rhob=rB*rP/(fbplus*rB+(1.0-fbplus)*rP); else rhob=rB;


    // 0.8 - SIMPLEC параметр релаксации для скоростей (см Гаврилов Андрей)
	// уравнение поправки давления.

	doublereal De, Dw, Ds, Dn, Dt, Db; // диффузионный поток через грани КО.
	
	// Псевдовремя tau является коэффициентом диффузии 
	// в эллиптическом уравнении на поправку давления.

	if (!bE) {
		if (bW) De=dbeta*taue*dy*dz/dxe;
		else De=taue*dy*dz/dxe;
	}
	else De=dbeta*taue*dy*dz/dxe;

	if (!bW) {
		if (bE) Dw=dbeta*tauw*dy*dz/dxw;
		else Dw=tauw*dy*dz/dxw; 
	}
	else Dw=dbeta*tauw*dy*dz/dxw;

	if (!bN) {
		if (bS) Dn=dbeta*taun*dx*dz/dyn;
		else Dn=taun*dx*dz/dyn; 
	}
	else Dn=dbeta*taun*dx*dz/dyn;
	
	if (!bS) {
		if (bN) Ds=dbeta*taus*dx*dz/dys;
		else Ds=taus*dx*dz/dys;
	}
	else Ds=dbeta*taus*dx*dz/dys;


	if (!bT) {
		if (bB) Dt=dbeta*taut*dx*dy/dzt;
		else Dt=taut*dx*dy/dzt;
	} 
	else Dt=dbeta*taut*dx*dy/dzt;


	if (!bB) {
		if (bT) Db=dbeta*taub*dx*dy/dzb;
		else Db=taub*dx*dy/dzb;
	}
	else Db=dbeta*taub*dx*dy/dzb;

	
	// Число Пекле равно нулю.
	sl[PAM][iP].ae=De*1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[PAM][iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[PAM][iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[PAM][iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[PAM][iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[PAM][iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	sl[PAM][iP].ap=sl[PAM][iP].ae+sl[PAM][iP].aw+sl[PAM][iP].an+sl[PAM][iP].as+sl[PAM][iP].at+sl[PAM][iP].ab;

	doublereal baddDFLUX2=0.0;
	if (bhighorder) {
		// если bborder == false то узел строго внутренний.
		// если bborder == true то мы находимся вблизи граничного узла.
		bool bborder=false;
		doublereal myflux=0.0;
		myflux=De*(dxe*DFDXiP(potent[PAM], iP, ESIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iE]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
		myflux=Dw*(-dxw*DFDXiP(potent[PAM], iP, WSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iW]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dn*(dyn*DFDXiP(potent[PAM], iP, NSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iN]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Ds*(-dys*DFDXiP(potent[PAM], iP, SSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iS]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dt*(dzt*DFDXiP(potent[PAM], iP, TSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iT]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Db*(-dzb*DFDXiP(potent[PAM], iP, BSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iB]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
	}



    doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 
	     
	Fe=calcFg2(bE, feplus, iP, ESIDE, rhoe, 
			dy*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, sosedi, 
			bRhieChowi, bRhieChowb, false,
			tau, dxe, potent);

	Fw=calcFg2(bW, fwplus, iP, WSIDE, rhow, 
			dy*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, sosedi, 
			bRhieChowi, bRhieChowb, false,
			tau, dxw, potent);

	Fn=calcFg2(bN, fnplus, iP, NSIDE, rhon, 
			dx*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, sosedi, 
			bRhieChowi, bRhieChowb, false,
			tau, dyn, potent);

	Fs=calcFg2(bS, fsplus, iP, SSIDE, rhos, 
			dx*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, sosedi, 
			bRhieChowi, bRhieChowb, false,
			tau, dys, potent);

	Ft=calcFg2(bT, ftplus, iP, TSIDE, rhot, 
			dx*dy, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, sosedi, 
			bRhieChowi, bRhieChowb, false,
			tau, dzt, potent);

	Fb=calcFg2(bB, fbplus, iP, BSIDE, rhob, 
			dx*dy, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, sosedi, 
			bRhieChowi, bRhieChowb, false,
			tau, dzb, potent);
	
	bool ISezai=true;
	if (ISezai) {

		doublereal SpeedCorOlde, SpeedCorOldw,  SpeedCorOldn,  SpeedCorOlds,  SpeedCorOldt,  SpeedCorOldb; 
	    if (!bE) { 
		   SpeedCorOlde=feplus*potent[VXCOR][iE]+(1.0-feplus)*potent[VXCOR][iP];
	    }
	    else {
	       SpeedCorOlde=potent[VXCOR][iE];
	    }
	    if (!bN) { 
		   SpeedCorOldn=fnplus*potent[VYCOR][iN]+(1.0-fnplus)*potent[VYCOR][iP];
	    }
	    else {
	       SpeedCorOldn=potent[VYCOR][iN];
	    }
	    if (!bT) { 
		   SpeedCorOldt=ftplus*potent[VZCOR][iT]+(1.0-ftplus)*potent[VZCOR][iP];
	    }
	    else {
	       SpeedCorOldt=potent[VZCOR][iT];
	    }
	    if (!bW) { 
		   SpeedCorOldw=fwplus*potent[VXCOR][iW]+(1.0-fwplus)*potent[VXCOR][iP];
	    }
	    else {
	       SpeedCorOldw=potent[VXCOR][iW];
	    }
	    if (!bS) { 
		   SpeedCorOlds=fsplus*potent[VYCOR][iS]+(1.0-fsplus)*potent[VYCOR][iP];
	    }
	    else {
	       SpeedCorOlds=potent[VYCOR][iS];
    	}
	    if (!bB) { 
		   SpeedCorOldb=fbplus*potent[VZCOR][iB]+(1.0-fbplus)*potent[VZCOR][iP];
	    }
	    else {
	       SpeedCorOldb=potent[VZCOR][iB];
	    }

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    // дополнительная поправка осуществляется на основе скоректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    Fe+=(1.0-alpha[VX])*(mfcurrentretune[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    Fn+=(1.0-alpha[VY])*(mfcurrentretune[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    Ft+=(1.0-alpha[VZ])*(mfcurrentretune[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    Fw+=(1.0-alpha[VX])*(mfcurrentretune[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    Fs+=(1.0-alpha[VY])*(mfcurrentretune[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    Fb+=(1.0-alpha[VZ])*(mfcurrentretune[BSIDE]-rhob*SpeedCorOldb*dx*dy);

	} 
	

	
	// возвращаем значение потока на грани КО.
	mfcurrentretune[ESIDE]=Fe;
	mfcurrentretune[NSIDE]=Fn;
	mfcurrentretune[TSIDE]=Ft;
	mfcurrentretune[WSIDE]=Fw;
	mfcurrentretune[SSIDE]=Fs;
	mfcurrentretune[BSIDE]=Fb;
	
	

	sl[PAM][iP].b=Fw-Fe+Fs-Fn+Fb-Ft+baddDFLUX2;



	// Симметризация СЛАУ :

	// для поправки давления получается эллиптическое уравнение с SPD матрицей.
	
	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

		// 1. Учёт условия Дирихле:
		if ((iE>=maxelm) && (fabs(slb[PAM][iE-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].ae*slb[PAM][iE-maxelm].b/slb[PAM][iE-maxelm].aw;
			sl[PAM][iP].ae=0.0;
			sl[PAM][iP].iE=-1; // не входит в матрицу СЛАУ.
		}
		if ((iW>=maxelm) && (fabs(slb[PAM][iW-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].aw*slb[PAM][iW-maxelm].b/slb[PAM][iW-maxelm].aw;
			sl[PAM][iP].aw=0.0;
			sl[PAM][iP].iW=-1; // не входит в матрицу СЛАУ.
		}
		if ((iN>=maxelm) && (fabs(slb[PAM][iN-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].an*slb[PAM][iN-maxelm].b/slb[PAM][iN-maxelm].aw;
			sl[PAM][iP].an=0.0;
			sl[PAM][iP].iN=-1; // не входит в матрицу СЛАУ.
		}
		if ((iS>=maxelm) && (fabs(slb[PAM][iS-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].as*slb[PAM][iS-maxelm].b/slb[PAM][iS-maxelm].aw;
			sl[PAM][iP].as=0.0;
			sl[PAM][iP].iS=-1; // не входит в матрицу СЛАУ.
		}
		if ((iT>=maxelm) && (fabs(slb[PAM][iT-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].at*slb[PAM][iT-maxelm].b/slb[PAM][iT-maxelm].aw;
			sl[PAM][iP].at=0.0;
			sl[PAM][iP].iT=-1; // не входит в матрицу СЛАУ.
		}
		if ((iB>=maxelm) && (fabs(slb[PAM][iB-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].ab*slb[PAM][iB-maxelm].b/slb[PAM][iB-maxelm].aw;
			sl[PAM][iP].ab=0.0;
			sl[PAM][iP].iB=-1; // не входит в матрицу СЛАУ.
		}

		

} // my_elmatr_quad_PAm2

// Составляет матрицу для уравнения 
// поправки давления
// Данная сборка матрицы для поправки давления основана на сглаженном псевдовремени,
// так как рекомендует Гаврилов Андрей. Даже если сглаженное псевдовремя не поможет,
// всё равно этот способ проще и привлекательней т.к. содержит меньше деталей все детали спрятаны 
// в заранее вычисленный вектор tau (псевдовремя).
// реализовано 23 июня 2012 года.
void my_elmatr_quad_PAm3(integer iP, equation3D** &sl, equation3D_bon** &slb,  
						doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
						integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm,
						doublereal* alpha, doublereal dbeta, doublereal** &rhie_chow, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, doublereal** tau,
						bool bmyhighorder, bool bdeltapfinish, bool bRhieChowi, bool bRhieChowb,
	                    BOUND* sosedb, integer *ilevel_alice, integer* ptr) {

	// точность решения уравнения для поправки давления неудовлетворительная, поэтому будем решать несколько раз
	// каждый последующий раз увеличивая точность. 26 июня 2012 года.
	// bdeltapfinish==true. - финишное решение.

    // 31 марта 2012 из этой функции удалено много закоментированного кода. Его можно найти в backup`ах AliceFlowv0_07.

	// tau - заранее вычисленное псевдовремя, оно служит коэффициентом диффузии в эллиптическом операторе, к тому-же оно
	// также входит в монотонизирующую поправку Рхи-Чоу (1983 года).

    /* btimedepend - стационарный (false), нестационарный (true).
	*  dtimestep - величина текущего шага по времени.
	*  
	*
	*  mfoldtimestep - массовый поток через каждую грань КО с предыдущего временного слоя. (неизменяемая в данной функции величина).
	*  mfcurrentretune - возвращаемый текущий массовый поток через каждую грань КО.
	*  КО имеет 6 граней : E,N,T,W,S,B.
	*
	* speedoldtimestep[VX or VY or VZ][iP : 0..maxelm+maxbound-1] - поле скорости с предыдущего временного слоя.
	*
	*/

	// Возможно если bdeltapfinish==true то имеет смысл включить bhighorder.
	bool bhighorder=bmyhighorder;

	doublereal eps=admission; // для отделения вещественного нуля.
	// отключает или включает поправку Рхи-Чоу 1983г.
	//bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border теперь данные параметры передаются внутрь функции извне.

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=sosedi[ESIDE][iP].iNODE1; iN=sosedi[NSIDE][iP].iNODE1; iT=sosedi[TSIDE][iP].iNODE1;
	iW=sosedi[WSIDE][iP].iNODE1; iS=sosedi[SSIDE][iP].iNODE1; iB=sosedi[BSIDE][iP].iNODE1;
	sl[PAM][iP].iP=iP;
	sl[PAM][iP].iE=iE; sl[PAM][iP].iN=iN; 
	sl[PAM][iP].iS=iS; sl[PAM][iP].iW=iW;
    sl[PAM][iP].iT=iT; sl[PAM][iP].iB=iB;


	// 26.09.2016 Добавок для АЛИС сетки.
	integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
	integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
	integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов

										  // -1 если не используется и [0..maxelm+maxbound-1] если используется.

	iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
	iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;
	iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
	iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;
	iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
	iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;

	sl[PAM][iP].iE2 = iE2; sl[PAM][iP].iN2 = iN2; sl[PAM][iP].iT2 = iT2;
	sl[PAM][iP].iS2 = iS2; sl[PAM][iP].iW2 = iW2; sl[PAM][iP].iB2 = iB2;

	sl[PAM][iP].iE3 = iE3; sl[PAM][iP].iN3 = iN3; sl[PAM][iP].iT3 = iT3;
	sl[PAM][iP].iS3 = iS3; sl[PAM][iP].iW3 = iW3; sl[PAM][iP].iB3 = iB3;

	sl[PAM][iP].iE4 = iE4; sl[PAM][iP].iN4 = iN4; sl[PAM][iP].iT4 = iT4;
	sl[PAM][iP].iS4 = iS4; sl[PAM][iP].iW4 = iW4; sl[PAM][iP].iB4 = iB4;

	// Признак прсутствия связи.
	// От булевых флагов можно избавиться в целях экономии памяти ЭВМ.
	sl[PAM][iP].bE2 = false; sl[PAM][iP].bW2 = false; sl[PAM][iP].bS2 = false;
	sl[PAM][iP].bN2 = false; sl[PAM][iP].bB2 = false; sl[PAM][iP].bT2 = false;

	sl[PAM][iP].bE3 = false; sl[PAM][iP].bW3 = false; sl[PAM][iP].bS3 = false;
	sl[PAM][iP].bN3 = false; sl[PAM][iP].bB3 = false; sl[PAM][iP].bT3 = false;

	sl[PAM][iP].bE4 = false; sl[PAM][iP].bW4 = false; sl[PAM][iP].bS4 = false;
	sl[PAM][iP].bN4 = false; sl[PAM][iP].bB4 = false; sl[PAM][iP].bT4 = false;

	if (iE2 > -1) sl[PAM][iP].bE2 = true;
	if (iW2 > -1) sl[PAM][iP].bW2 = true;
	if (iN2 > -1) sl[PAM][iP].bN2 = true;
	if (iS2 > -1) sl[PAM][iP].bS2 = true;
	if (iT2 > -1) sl[PAM][iP].bT2 = true;
	if (iB2 > -1) sl[PAM][iP].bB2 = true;

	if (iE3 > -1) sl[PAM][iP].bE3 = true;
	if (iW3 > -1) sl[PAM][iP].bW3 = true;
	if (iN3 > -1) sl[PAM][iP].bN3 = true;
	if (iS3 > -1) sl[PAM][iP].bS3 = true;
	if (iT3 > -1) sl[PAM][iP].bT3 = true;
	if (iB3 > -1) sl[PAM][iP].bB3 = true;

	if (iE4 > -1) sl[PAM][iP].bE4 = true;
	if (iW4 > -1) sl[PAM][iP].bW4 = true;
	if (iN4 > -1) sl[PAM][iP].bN4 = true;
	if (iS4 > -1) sl[PAM][iP].bS4 = true;
	if (iT4 > -1) sl[PAM][iP].bT4 = true;
	if (iB4 > -1) sl[PAM][iP].bB4 = true;

	// если с одной из сторон граница расчётной области 
	// то переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

	if (iE > -1) if (iE>=maxelm) bE=true;
	if (iN > -1) if (iN>=maxelm) bN=true;
	if (iT > -1) if (iT>=maxelm) bT=true;
	if (iW > -1) if (iW>=maxelm) bW=true;
	if (iS > -1) if (iS>=maxelm) bS=true;
	if (iB > -1) if (iB>=maxelm) bB=true;


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

	// Значение коэффициента диффузии (псевдовремени) на грани контрольного объёма.
	// Учитывается что псевдовремя это три скалярных поля.
	doublereal taue=1.0, tauw=1.0, taun=1.0, taus=1.0, taut=1.0, taub=1.0;
	if (iE > -1) {
		if (!bE) taue = tau[VX][iE] * tau[VX][iP] / (feplus*tau[VX][iE] + (1.0 - feplus)*tau[VX][iP]); else taue = tau[VX][iE]; // проверено !
	}
	if (iW > -1) {
		if (!bW) tauw = tau[VX][iW] * tau[VX][iP] / (fwplus*tau[VX][iW] + (1.0 - fwplus)*tau[VX][iP]); else tauw = tau[VX][iW];
	}
	if (iN > -1) {
		if (!bN) taun = tau[VY][iN] * tau[VY][iP] / (fnplus*tau[VY][iN] + (1.0 - fnplus)*tau[VY][iP]); else taun = tau[VY][iN];
	}
	if (iS > -1) {
		if (!bS) taus = tau[VY][iS] * tau[VY][iP] / (fsplus*tau[VY][iS] + (1.0 - fsplus)*tau[VY][iP]); else taus = tau[VY][iS];
	}
	if (iT > -1) {
		if (!bT) taut = tau[VZ][iT] * tau[VZ][iP] / (ftplus*tau[VZ][iT] + (1.0 - ftplus)*tau[VZ][iP]); else taut = tau[VZ][iT];
	}
	if (iB > -1) {
		if (!bB) taub = tau[VZ][iB] * tau[VZ][iP] / (fbplus*tau[VZ][iB] + (1.0 - fbplus)*tau[VZ][iP]); else taub = tau[VZ][iB];
	}

	doublereal taue2 = 1.0, tauw2 = 1.0, taun2 = 1.0, taus2 = 1.0, taut2 = 1.0, taub2 = 1.0;
	if (iE2 > -1) {
		if (!bE2) taue2 = tau[VX][iE2] * tau[VX][iP] / (feplus2*tau[VX][iE2] + (1.0 - feplus2)*tau[VX][iP]); else taue2 = tau[VX][iE2]; // проверено !
	}
	if (iW2 > -1) {
		if (!bW2) tauw2 = tau[VX][iW2] * tau[VX][iP] / (fwplus2*tau[VX][iW2] + (1.0 - fwplus2)*tau[VX][iP]); else tauw2 = tau[VX][iW2];
	}
	if (iN2 > -1) {
		if (!bN2) taun2 = tau[VY][iN2] * tau[VY][iP] / (fnplus2*tau[VY][iN2] + (1.0 - fnplus2)*tau[VY][iP]); else taun2 = tau[VY][iN2];
	}
	if (iS2 > -1) {
		if (!bS2) taus2 = tau[VY][iS2] * tau[VY][iP] / (fsplus2*tau[VY][iS2] + (1.0 - fsplus2)*tau[VY][iP]); else taus2 = tau[VY][iS2];
	}
	if (iT2 > -1) {
		if (!bT2) taut2 = tau[VZ][iT2] * tau[VZ][iP] / (ftplus2*tau[VZ][iT2] + (1.0 - ftplus2)*tau[VZ][iP]); else taut2 = tau[VZ][iT2];
	}
	if (iB2 > -1) {
		if (!bB2) taub2 = tau[VZ][iB2] * tau[VZ][iP] / (fbplus2*tau[VZ][iB2] + (1.0 - fbplus2)*tau[VZ][iP]); else taub2 = tau[VZ][iB2];
	}

	doublereal taue3 = 1.0, tauw3 = 1.0, taun3 = 1.0, taus3 = 1.0, taut3 = 1.0, taub3 = 1.0;
	if (iE3 > -1) {
		if (!bE3) taue3 = tau[VX][iE3] * tau[VX][iP] / (feplus3*tau[VX][iE3] + (1.0 - feplus3)*tau[VX][iP]); else taue3 = tau[VX][iE3]; // проверено !
	}
	if (iW3 > -1) {
		if (!bW3) tauw3 = tau[VX][iW3] * tau[VX][iP] / (fwplus3*tau[VX][iW3] + (1.0 - fwplus3)*tau[VX][iP]); else tauw3 = tau[VX][iW3];
	}
	if (iN3 > -1) {
		if (!bN3) taun3 = tau[VY][iN3] * tau[VY][iP] / (fnplus3*tau[VY][iN3] + (1.0 - fnplus3)*tau[VY][iP]); else taun3 = tau[VY][iN3];
	}
	if (iS3 > -1) {
		if (!bS3) taus3 = tau[VY][iS3] * tau[VY][iP] / (fsplus3*tau[VY][iS3] + (1.0 - fsplus3)*tau[VY][iP]); else taus3 = tau[VY][iS3];
	}
	if (iT3 > -1) {
		if (!bT3) taut3 = tau[VZ][iT3] * tau[VZ][iP] / (ftplus3*tau[VZ][iT3] + (1.0 - ftplus3)*tau[VZ][iP]); else taut3 = tau[VZ][iT3];
	}
	if (iB3 > -1) {
		if (!bB3) taub3 = tau[VZ][iB3] * tau[VZ][iP] / (fbplus3*tau[VZ][iB3] + (1.0 - fbplus3)*tau[VZ][iP]); else taub3 = tau[VZ][iB3];
	}


	doublereal taue4 = 1.0, tauw4 = 1.0, taun4 = 1.0, taus4 = 1.0, taut4 = 1.0, taub4 = 1.0;
	if (iE4 > -1) {
		if (!bE4) taue4 = tau[VX][iE4] * tau[VX][iP] / (feplus4*tau[VX][iE4] + (1.0 - feplus4)*tau[VX][iP]); else taue4 = tau[VX][iE4]; // проверено !
	}
	if (iW4 > -1) {
		if (!bW4) tauw4 = tau[VX][iW4] * tau[VX][iP] / (fwplus4*tau[VX][iW4] + (1.0 - fwplus4)*tau[VX][iP]); else tauw4 = tau[VX][iW4];
	}
	if (iN4 > -1) {
		if (!bN4) taun4 = tau[VY][iN4] * tau[VY][iP] / (fnplus4*tau[VY][iN4] + (1.0 - fnplus4)*tau[VY][iP]); else taun4 = tau[VY][iN4];
	}
	if (iS4 > -1) {
		if (!bS4) taus4 = tau[VY][iS4] * tau[VY][iP] / (fsplus4*tau[VY][iS4] + (1.0 - fsplus4)*tau[VY][iP]); else taus4 = tau[VY][iS4];
	}
	if (iT4 > -1) {
		if (!bT4) taut4 = tau[VZ][iT4] * tau[VZ][iP] / (ftplus4*tau[VZ][iT4] + (1.0 - ftplus4)*tau[VZ][iP]); else taut4 = tau[VZ][iT4];
	}
	if (iB4 > -1) {
		if (!bB4) taub4 = tau[VZ][iB4] * tau[VZ][iP] / (fbplus4*tau[VZ][iB4] + (1.0 - fbplus4)*tau[VZ][iP]); else taub4 = tau[VZ][iB4];
	}


	/*
	doublereal de, dw, dn, ds, dt, db;
	// Так рекомендует делать С. Патанкар, надо только
	// помнить что это должно быть согласовано с этой-же
	// величиной при коррекции компонент скорости.
	//de=dy*dz/apue; dw=dy*dz/apuw; 
	//dn=dx*dz/apvn; ds=dx*dz/apvs;
	//dt=dx*dy/apwt; db=dx*dy/apwb;
	// Так рекомендуют делать в следующих статьях:
	//1. В.В.Винников, Д.Л.Ревизников
	// Применение декартовых сеток для решения уравнений Навье-Стокса
	// в областях с криволинейными границами. 
	// МАИ (ГТУ) Математическое Моделирование 2005г, том 17, номер 8, стр 15-30
	// 2. SIMPLE METHOD FOR THE SOLUTION OF INCOMPRESSIBLE FLOWS ON NON-STAGGERED GRIDS
	// I.Sezai - Eastern Mediterranean University. Revised January, 2011.
	// 3. Гаврилов Андрей sigma-flow.
	// Такие значения также должны быть согласованы с теми что используются при коррекции скоростей.
	if (iSIMPLE_alg==SIMPLE_Carretto) {
		// SIMPLE
		de=alpha[VX]*dy*dz/apue; dw=alpha[VX]*dy*dz/apuw; 
	    dn=alpha[VY]*dx*dz/apvn; ds=alpha[VY]*dx*dz/apvs;
	    dt=alpha[VZ]*dx*dy/apwt; db=alpha[VZ]*dx*dy/apwb;
	}
	if (iSIMPLE_alg==SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC
		de=alpha[VX]*dy*dz/((1.0-alpha[VX])*apue); dw=alpha[VX]*dy*dz/((1.0-alpha[VX])*apuw); 
	    dn=alpha[VY]*dx*dz/((1.0-alpha[VY])*apvn); ds=alpha[VY]*dx*dz/((1.0-alpha[VY])*apvs);
	    dt=alpha[VZ]*dx*dy/((1.0-alpha[VZ])*apwt); db=alpha[VZ]*dx*dy/((1.0-alpha[VZ])*apwb);
	}
	*/
    // плотность аппроксимируется средним гармоническим
	
	doublereal rP, rE=0.0, rN=0.0, rT=0.0, rW=0.0, rS=0.0, rB=0.0;

    rP=prop[RHO][iP];
	
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

	

	// В уравнении теплопроводности положительный эффект был обнаружен при выравнивании коэффициентов 
	// теплопроводности в граничном и ближайшем внутреннем узле. Зделаем здесь аналогичное действо.
	// Экспериментально было обнаружено, что данное изменение не влияет не на температуру ни на continity в данном случае.
	/* // пока мы от него отказались !
	if (!bE) rE=prop[RHO][iE]; else rE=rP;
    if (!bN) rN=prop[RHO][iN]; else rN=rP;
    if (!bT) rT=prop[RHO][iT]; else rT=rP;
	if (!bW) rW=prop[RHO][iW]; else rW=rP;
    if (!bS) rS=prop[RHO][iS]; else rS=rP;
    if (!bB) rB=prop[RHO][iB]; else rB=rP;
	*/

	doublereal rhoe = 0.0, rhow = 0.0, rhon = 0.0, rhos = 0.0, rhot = 0.0, rhob = 0.0;
	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE > -1) {
		if (!bE) rhoe = rE * rP / (feplus*rE + (1.0 - feplus)*rP); else rhoe = rE; // проверено !
	}
	if (iW > -1) {
		if (!bW) rhow = rW * rP / (fwplus*rW + (1.0 - fwplus)*rP); else rhow = rW;
	}
	if (iN > -1) {
		if (!bN) rhon = rN * rP / (fnplus*rN + (1.0 - fnplus)*rP); else rhon = rN;
	}
	if (iS > -1) {
		if (!bS) rhos = rS * rP / (fsplus*rS + (1.0 - fsplus)*rP); else rhos = rS;
	}
	if (iT > -1) {
		if (!bT) rhot = rT * rP / (ftplus*rT + (1.0 - ftplus)*rP); else rhot = rT;
	}
	if (iB > -1) {
		if (!bB) rhob = rB * rP / (fbplus*rB + (1.0 - fbplus)*rP); else rhob = rB;
	}

	doublereal rhoe2=0.0, rhow2=0.0, rhon2=0.0, rhos2=0.0, rhot2=0.0, rhob2=0.0;
	doublereal rhoe3=0.0, rhow3=0.0, rhon3=0.0, rhos3=0.0, rhot3=0.0, rhob3=0.0;
	doublereal rhoe4=0.0, rhow4=0.0, rhon4=0.0, rhos4=0.0, rhot4=0.0, rhob4=0.0;

	if (iE2 > -1) {
		if (!bE2)  rhoe2 = rE2 * rP / (feplus2*rE2 + (1.0 - feplus2)*rP); else rhoe2 = rE2; // проверено !
	}
	if (iW2 > -1) {
		if (!bW2)  rhow2 = rW2 * rP / (fwplus2*rW2 + (1.0 - fwplus2)*rP); else rhow2 = rW2;
	}
	if (iN2 > -1) {
		if (!bN2) rhon2 = rN2 * rP / (fnplus2*rN2 + (1.0 - fnplus2)*rP); else rhon2 = rN2;
	}
	if (iS2 > -1) {
		if (!bS2)  rhos2 = rS2 * rP / (fsplus2*rS2 + (1.0 - fsplus2)*rP); else rhos2 = rS2;
	}
	if (iT2 > -1) {
		if (!bT2)  rhot2 = rT2 * rP / (ftplus2*rT2 + (1.0 - ftplus2)*rP); else rhot2 = rT2;
	}
	if (iB2 > -1) {
		if (!bB2) rhob2 = rB2 * rP / (fbplus2*rB2 + (1.0 - fbplus2)*rP); else rhob2 = rB2;
	}

	if (iE3 > -1) {
		if (!bE3) rhoe3 = rE3 * rP / (feplus3*rE3 + (1.0 - feplus3)*rP); else rhoe3 = rE3;
	}
	if (iW3 > -1) {
		if (!bW3) rhow3 = rW3 * rP / (fwplus3*rW3 + (1.0 - fwplus3)*rP); else rhow3 = rW3;
	}
	if (iN3 > -1) {
		if (!bN3) rhon3 = rN3 * rP / (fnplus3*rN3 + (1.0 - fnplus3)*rP); else rhon3 = rN3;
	}
	if (iS3 > -1) {
		if (!bS3) rhos3 = rS3 * rP / (fsplus3*rS3 + (1.0 - fsplus3)*rP); else rhos3 = rS3;
	}
	if (iT3 > -1) {
		if (!bT3) rhot3 = rT3 * rP / (ftplus3*rT3 + (1.0 - ftplus3)*rP); else rhot3 = rT3;
	}
	if (iB3 > -1) {
		if (!bB3) rhob3 = rB3 * rP / (fbplus3*rB3 + (1.0 - fbplus3)*rP); else rhob3 = rB3;
	}

	if (iE4 > -1) {
		if (!bE4) rhoe4 = rE4 * rP / (feplus4*rE4 + (1.0 - feplus4)*rP); else rhoe4 = rE4;
	}
	if (iW4 > -1) {
		if (!bW4) rhow4 = rW4 * rP / (fwplus4*rW4 + (1.0 - fwplus4)*rP); else rhow4 = rW4;
	}
	if (iN4 > -1) {
		if (!bN4) rhon4 = rN4 * rP / (fnplus4*rN4 + (1.0 - fnplus4)*rP); else rhon4 = rN4;
	}
	if (iS4 > -1) {
		if (!bS4) rhos4 = rS4 * rP / (fsplus4*rS4 + (1.0 - fsplus4)*rP); else rhos4 = rS4;
	}
	if (iT4 > -1) {
		if (!bT4) rhot4 = rT4 * rP / (ftplus4*rT4 + (1.0 - ftplus4)*rP); else rhot4 = rT4;
	}
	if (iB4 > -1) {
		if (!bB4) rhob4 = rB4 * rP / (fbplus4*rB4 + (1.0 - fbplus4)*rP); else rhob4 = rB4;
	}


    // 0.8 - SIMPLEC параметр релаксации для скоростей (см Гаврилов Андрей)
	// уравнение поправки давления.

	const doublereal ZeroDiffusion = 1.0e-30;
	// Диффузионная составляющая потока:
	doublereal De = ZeroDiffusion, Dw = ZeroDiffusion, Dn = ZeroDiffusion, Ds = ZeroDiffusion, Dt = ZeroDiffusion, Db = ZeroDiffusion; // инициализация
	doublereal De2 = ZeroDiffusion, Dw2 = ZeroDiffusion, Dn2 = ZeroDiffusion, Ds2 = ZeroDiffusion, Dt2 = ZeroDiffusion, Db2 = ZeroDiffusion; // инициализация
	doublereal De3 = ZeroDiffusion, Dw3 = ZeroDiffusion, Dn3 = ZeroDiffusion, Ds3 = ZeroDiffusion, Dt3 = ZeroDiffusion, Db3 = ZeroDiffusion; // инициализация
	doublereal De4 = ZeroDiffusion, Dw4 = ZeroDiffusion, Dn4 = ZeroDiffusion, Ds4 = ZeroDiffusion, Dt4 = ZeroDiffusion, Db4 = ZeroDiffusion; // инициализация
		
	// Псевдовремя tau является коэффициентом диффузии 
	// в эллиптическом уравнении на поправку давления.
	if (iE > -1) {
		if (!bE) {
			if (bW) De = dbeta * taue*dy*dz / dxe;
			else De = taue * dy*dz / dxe;
		}
		else De = dbeta * taue*dy*dz / dxe;
	}

	if (iW > -1) {
		if (!bW) {
			if (bE) Dw = dbeta * tauw*dy*dz / dxw;
			else Dw = tauw * dy*dz / dxw;
		}
		else Dw = dbeta * tauw*dy*dz / dxw;
	}

	if (iN > -1) {
		if (!bN) {
			if (bS) Dn = dbeta * taun*dx*dz / dyn;
			else Dn = taun * dx*dz / dyn;
		}
		else Dn = dbeta * taun*dx*dz / dyn;
	}
	
	if (iS > -1) {
		if (!bS) {
			if (bN) Ds = dbeta * taus*dx*dz / dys;
			else Ds = taus * dx*dz / dys;
		}
		else Ds = dbeta * taus*dx*dz / dys;
	}

	if (iT > -1) {
		if (!bT) {
			if (bB) Dt = dbeta * taut*dx*dy / dzt;
			else Dt = taut * dx*dy / dzt;
		}
		else Dt = dbeta * taut*dx*dy / dzt;
	}

	if (iB > -1) {
		if (!bB) {
			if (bT) Db = dbeta * taub*dx*dy / dzb;
			else Db = taub * dx*dy / dzb;
		}
		else Db = dbeta * taub*dx*dy / dzb;
	}
	
	if (iE2 > -1) {
		if (!bE2) {
			if (bW2) De2 = dbeta * taue2*dy*dz / dxe2;
			else De2 = taue2 * dy*dz / dxe2;
		}
		else De2 = dbeta * taue2*dy*dz / dxe2;
	}

	if (iW2 > -1) {
		if (!bW2) {
			if (bE2) Dw2 = dbeta * tauw2*dy*dz / dxw2;
			else Dw2 = tauw2 * dy*dz / dxw2;
		}
		else Dw2 = dbeta * tauw2*dy*dz / dxw2;
	}

	if (iN2 > -1) {
		if (!bN2) {
			if (bS2) Dn2 = dbeta * taun2*dx*dz / dyn2;
			else Dn2 = taun2 * dx*dz / dyn2;
		}
		else Dn2 = dbeta * taun2*dx*dz / dyn2;
	}

	if (iS2 > -1) {
		if (!bS2) {
			if (bN2) Ds2 = dbeta * taus2*dx*dz / dys2;
			else Ds2 = taus2 * dx*dz / dys2;
		}
		else Ds2 = dbeta * taus2*dx*dz / dys2;
	}

	if (iT2 > -1) {
		if (!bT2) {
			if (bB2) Dt2 = dbeta * taut2*dx*dy / dzt2;
			else Dt2 = taut2 * dx*dy / dzt2;
		}
		else Dt2 = dbeta * taut2*dx*dy / dzt2;
	}

	if (iB2 > -1) {
		if (!bB2) {
			if (bT2) Db2 = dbeta * taub2*dx*dy / dzb2;
			else Db2 = taub2 * dx*dy / dzb2;
		}
		else Db2 = dbeta * taub2*dx*dy / dzb2;
	}

	if (iE3 > -1) {
		if (!bE3) {
			if (bW3) De3 = dbeta * taue3*dy*dz / dxe3;
			else De3 = taue3 * dy*dz / dxe3;
		}
		else De3 = dbeta * taue3*dy*dz / dxe3;
	}

	if (iW3 > -1) {
		if (!bW3) {
			if (bE3) Dw3 = dbeta * tauw3*dy*dz / dxw3;
			else Dw3 = tauw3 * dy*dz / dxw3;
		}
		else Dw3 = dbeta * tauw3*dy*dz / dxw3;
	}

	if (iN3 > -1) {
		if (!bN3) {
			if (bS3) Dn3 = dbeta * taun3*dx*dz / dyn3;
			else Dn3 = taun3 * dx*dz / dyn3;
		}
		else Dn3 = dbeta * taun3*dx*dz / dyn3;
	}

	if (iS3 > -1) {
		if (!bS3) {
			if (bN3) Ds3 = dbeta * taus3*dx*dz / dys3;
			else Ds3 = taus3 * dx*dz / dys3;
		}
		else Ds3 = dbeta * taus3*dx*dz / dys3;
	}

	if (iT3 > -1) {
		if (!bT3) {
			if (bB3) Dt = dbeta * taut3*dx*dy / dzt3;
			else Dt3 = taut3 * dx*dy / dzt3;
		}
		else Dt3 = dbeta * taut3*dx*dy / dzt3;
	}

	if (iB3 > -1) {
		if (!bB3) {
			if (bT3) Db3 = dbeta * taub3*dx*dy / dzb3;
			else Db3 = taub3 * dx*dy / dzb3;
		}
		else Db3 = dbeta * taub3*dx*dy / dzb3;
	}

	if (iE4 > -1) {
		if (!bE4) {
			if (bW4) De4 = dbeta * taue4*dy*dz / dxe4;
			else De4 = taue4 * dy*dz / dxe4;
		}
		else De4 = dbeta * taue4*dy*dz / dxe4;
	}

	if (iW4 > -1) {
		if (!bW4) {
			if (bE4) Dw4 = dbeta * tauw4*dy*dz / dxw4;
			else Dw4 = tauw4 * dy*dz / dxw4;
		}
		else Dw4 = dbeta * tauw4*dy*dz / dxw4;
	}

	if (iN4 > -1) {
		if (!bN4) {
			if (bS4) Dn4 = dbeta * taun4*dx*dz / dyn4;
			else Dn4 = taun4 * dx*dz / dyn4;
		}
		else Dn4 = dbeta * taun4*dx*dz / dyn4;
	}

	if (iS4 > -1) {
		if (!bS4) {
			if (bN4) Ds4 = dbeta * taus4*dx*dz / dys4;
			else Ds4 = taus4 * dx*dz / dys4;
		}
		else Ds4 = dbeta * taus4*dx*dz / dys4;
	}

	if (iT4 > -1) {
		if (!bT4) {
			if (bB4) Dt4 = dbeta * taut4*dx*dy / dzt4;
			else Dt4 = taut4 * dx*dy / dzt4;
		}
		else Dt4 = dbeta * taut4*dx*dy / dzt4;
	}

	if (iB4 > -1) {
		if (!bB4) {
			if (bT4) Db4 = dbeta * taub4*dx*dy / dzb4;
			else Db4 = taub4 * dx*dy / dzb4;
		}
		else Db4 = dbeta * taub4*dx*dy / dzb4;
	}


	// Число Пекле равно нулю.
	sl[PAM][iP].ae=De*1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[PAM][iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[PAM][iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[PAM][iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[PAM][iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[PAM][iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	sl[PAM][iP].ae2 = De2 * 1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[PAM][iP].aw2 = Dw2; // *fD(0.0, EXP2, true, fwplus2); равна строго 1.0;
	sl[PAM][iP].an2 = Dn2; // *fD(0.0, EXP2, true, fnplus2); равна строго 1.0;
	sl[PAM][iP].as2 = Ds2; // *fD(0.0, EXP2, true, fsplus2); равна строго 1.0;
	sl[PAM][iP].at2 = Dt2; // *fD(0.0, EXP2, true, ftplus2); равна строго 1.0;
	sl[PAM][iP].ab2 = Db2; // *fD(0.0, EXP2, true, fbplus2); равна строго 1.0;

	sl[PAM][iP].ae3 = De3 * 1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus3); равна строго 1.0;
	sl[PAM][iP].aw3 = Dw3; // *fD(0.0, EXP2, true, fwplus3); равна строго 1.0;
	sl[PAM][iP].an3 = Dn3; // *fD(0.0, EXP2, true, fnplus3); равна строго 1.0;
	sl[PAM][iP].as3 = Ds3; // *fD(0.0, EXP2, true, fsplus3); равна строго 1.0;
	sl[PAM][iP].at3 = Dt3; // *fD(0.0, EXP2, true, ftplus3); равна строго 1.0;
	sl[PAM][iP].ab3 = Db3; // *fD(0.0, EXP2, true, fbplus3); равна строго 1.0;

	sl[PAM][iP].ae4 = De4 * 1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus4); равна строго 1.0;
	sl[PAM][iP].aw4 = Dw4; // *fD(0.0, EXP2, true, fwplus4); равна строго 1.0;
	sl[PAM][iP].an4 = Dn4; // *fD(0.0, EXP2, true, fnplus4); равна строго 1.0;
	sl[PAM][iP].as4 = Ds4; // *fD(0.0, EXP2, true, fsplus4); равна строго 1.0;
	sl[PAM][iP].at4 = Dt4; // *fD(0.0, EXP2, true, ftplus4); равна строго 1.0;
	sl[PAM][iP].ab4 = Db4; // *fD(0.0, EXP2, true, fbplus4); равна строго 1.0;

	sl[PAM][iP].ap=sl[PAM][iP].ae+sl[PAM][iP].aw+sl[PAM][iP].an+sl[PAM][iP].as+sl[PAM][iP].at+sl[PAM][iP].ab+
		sl[PAM][iP].ae2 + sl[PAM][iP].aw2 + sl[PAM][iP].an2 + sl[PAM][iP].as2 + sl[PAM][iP].at2 + sl[PAM][iP].ab2+
		sl[PAM][iP].ae3 + sl[PAM][iP].aw3 + sl[PAM][iP].an3 + sl[PAM][iP].as3 + sl[PAM][iP].at3 + sl[PAM][iP].ab3+
		sl[PAM][iP].ae4 + sl[PAM][iP].aw4 + sl[PAM][iP].an4 + sl[PAM][iP].as4 + sl[PAM][iP].at4 + sl[PAM][iP].ab4;

	doublereal baddDFLUX2=0.0;
	if (0&&bhighorder) {
		// На АЛИС сетке это использовать не представляется возможным.

		// если bborder == false то узел строго внутренний.
		// если bborder == true то мы находимся вблизи граничного узла.
		bool bborder=false;
		doublereal myflux=0.0;
		myflux=De*(dxe*DFDXiP(potent[PAM], iP, ESIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iE]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
		myflux=Dw*(-dxw*DFDXiP(potent[PAM], iP, WSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iW]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dn*(dyn*DFDXiP(potent[PAM], iP, NSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iN]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Ds*(-dys*DFDXiP(potent[PAM], iP, SSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iS]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dt*(dzt*DFDXiP(potent[PAM], iP, TSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iT]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Db*(-dzb*DFDXiP(potent[PAM], iP, BSIDE, sosedi, maxelm, nvtx, pa, bborder)-(potent[PAM][iB]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
	}


	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.

    doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 

	
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
		

		Fe = calcFg3(bE, feplus, iP, ESIDE, rhoe,
			dSqe, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe, potent, 1);
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

		Fw = calcFg3(bW, fwplus, iP, WSIDE, rhow,
			dSqw, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw, potent, 1);
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

		Fn = calcFg3(bN, fnplus, iP, NSIDE, rhon,
			dSqn, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn, potent, 1);
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

		Fs = calcFg3(bS, fsplus, iP, SSIDE, rhos,
			dSqs, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys, potent, 1);
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

		Ft = calcFg3(bT, ftplus, iP, TSIDE, rhot,
			dSqt, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt, potent, 1);
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

		Fb = calcFg3(bB, fbplus, iP, BSIDE, rhob,
			dSqb, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb, potent, 1);
	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal Fw2 = 0.0, Fe2 = 0.0, Fs2 = 0.0, Fn2 = 0.0, Ft2 = 0.0, Fb2 = 0.0;

	
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

		Fe2 = calcFg3(bE2, feplus2, iP, ESIDE, rhoe2,
			dSqe2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe2, potent, 2);
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

		Fw2 = calcFg3(bW2, fwplus2, iP, WSIDE, rhow2,
			dSqw2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw2, potent, 2);
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

		Fn2 = calcFg3(bN2, fnplus2, iP, NSIDE, rhon2,
			dSqn2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn2, potent, 2);
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

		Fs2 = calcFg3(bS2, fsplus2, iP, SSIDE, rhos2,
			dSqs2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys2, potent, 2);
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

		Ft2 = calcFg3(bT2, ftplus2, iP, TSIDE, rhot2,
			dSqt2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt2, potent, 2);
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

		Fb2 = calcFg3(bB2, fbplus2, iP, BSIDE, rhob2,
			dSqb2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb2, potent, 2);
	}
	

	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal Fw3 = 0.0, Fe3 = 0.0, Fs3 = 0.0, Fn3 = 0.0, Ft3 = 0.0, Fb3 = 0.0;
	
	
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

		Fe3 = calcFg3(bE3, feplus3, iP, ESIDE, rhoe3,
			dSqe3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe3, potent, 3);
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

		Fw3 = calcFg3(bW3, fwplus3, iP, WSIDE, rhow3,
			dSqw3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw3, potent, 3);
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

		Fn3 = calcFg3(bN3, fnplus3, iP, NSIDE, rhon3,
			dSqn3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn3, potent, 3);
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

		Fs3 = calcFg3(bS3, fsplus3, iP, SSIDE, rhos3,
			dSqs3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys3, potent, 3);
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

		Ft3 = calcFg3(bT3, ftplus3, iP, TSIDE, rhot3,
			dSqt3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt3, potent, 3);
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

		Fb3 = calcFg3(bB3, fbplus3, iP, BSIDE, rhob3,
			dSqb3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb3, potent, 3);
	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.
	doublereal Fw4 = 0.0, Fe4 = 0.0, Fs4 = 0.0, Fn4 = 0.0, Ft4 = 0.0, Fb4 = 0.0;

	
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

		Fe4 = calcFg3(bE4, feplus4, iP, ESIDE, rhoe4,
			dSqe4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxe4, potent, 4);
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

		Fw4 = calcFg3(bW4, fwplus4, iP, WSIDE, rhow4,
			dSqw4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dxw4, potent, 4);
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

		Fn4 = calcFg3(bN4, fnplus4, iP, NSIDE, rhon4,
			dSqn4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dyn4, potent, 4);
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

		Fs4 = calcFg3(bS4, fsplus4, iP, SSIDE, rhos4,
			dSqs4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dys4, potent, 4);
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

		Ft4 = calcFg3(bT4, ftplus4, iP, TSIDE, rhot4,
			dSqt4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzt4, potent, 4);
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

		Fb4 = calcFg3(bB4, fbplus4, iP, BSIDE, rhob4,
			dSqb4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, sosedi,
			bRhieChowi, bRhieChowb, false,
			tau, dzb4, potent, 4);
	}

	doublereal Fe_sum = Fe + Fe2 + Fe3 + Fe4;
	doublereal Fw_sum = Fw + Fw2 + Fw3 + Fw4;
	doublereal Fn_sum = Fn + Fn2 + Fn3 + Fn4;
	doublereal Fs_sum = Fs + Fs2 + Fs3 + Fs4;
	doublereal Ft_sum = Ft + Ft2 + Ft3 + Ft4;
	doublereal Fb_sum = Fb + Fb2 + Fb3 + Fb4;


	bool ISezai=true;
	if (ISezai) {

		doublereal SpeedCorOlde=0.0, SpeedCorOldw=0.0,  SpeedCorOldn=0.0,  SpeedCorOlds=0.0,  SpeedCorOldt=0.0,  SpeedCorOldb=0.0; 
		if (iE > -1) {
			if (!bE) {
				SpeedCorOlde = feplus * potent[VXCOR][iE] + (1.0 - feplus)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde = potent[VXCOR][iE];
			}
		}
		if (iN > -1) {
			if (!bN) {
				SpeedCorOldn = fnplus * potent[VYCOR][iN] + (1.0 - fnplus)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn = potent[VYCOR][iN];
			}
		}
		if (iT > -1) {
			if (!bT) {
				SpeedCorOldt = ftplus * potent[VZCOR][iT] + (1.0 - ftplus)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt = potent[VZCOR][iT];
			}
		}
		if (iW > -1) {
			if (!bW) {
				SpeedCorOldw = fwplus * potent[VXCOR][iW] + (1.0 - fwplus)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw = potent[VXCOR][iW];
			}
		}
		if (iS > -1) {
			if (!bS) {
				SpeedCorOlds = fsplus * potent[VYCOR][iS] + (1.0 - fsplus)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds = potent[VYCOR][iS];
			}
		}
		if (iB > -1) {
			if (!bB) {
				SpeedCorOldb = fbplus * potent[VZCOR][iB] + (1.0 - fbplus)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb = potent[VZCOR][iB];
			}
		}

		doublereal SpeedCorOlde2 = 0.0, SpeedCorOldw2 = 0.0, SpeedCorOldn2 = 0.0, SpeedCorOlds2 = 0.0, SpeedCorOldt2 = 0.0, SpeedCorOldb2 = 0.0;
		if (iE2 > -1) {
			if (!bE2) {
				SpeedCorOlde2 = feplus2 * potent[VXCOR][iE2] + (1.0 - feplus2)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde2 = potent[VXCOR][iE2];
			}
		}
		if (iN2 > -1) {
			if (!bN2) {
				SpeedCorOldn2 = fnplus2 * potent[VYCOR][iN2] + (1.0 - fnplus2)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn2 = potent[VYCOR][iN2];
			}
		}
		if (iT2 > -1) {
			if (!bT2) {
				SpeedCorOldt2 = ftplus2 * potent[VZCOR][iT2] + (1.0 - ftplus2)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt2 = potent[VZCOR][iT2];
			}
		}
		if (iW2 > -1) {
			if (!bW2) {
				SpeedCorOldw2 = fwplus2 * potent[VXCOR][iW2] + (1.0 - fwplus2)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw2 = potent[VXCOR][iW2];
			}
		}
		if (iS2 > -1) {
			if (!bS2) {
				SpeedCorOlds2 = fsplus2 * potent[VYCOR][iS2] + (1.0 - fsplus2)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds2 = potent[VYCOR][iS2];
			}
		}
		if (iB2 > -1) {
			if (!bB2) {
				SpeedCorOldb2 = fbplus2 * potent[VZCOR][iB2] + (1.0 - fbplus2)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb2 = potent[VZCOR][iB2];
			}
		}

		doublereal SpeedCorOlde3 = 0.0, SpeedCorOldw3 = 0.0, SpeedCorOldn3 = 0.0, SpeedCorOlds3 = 0.0, SpeedCorOldt3 = 0.0, SpeedCorOldb3 = 0.0;
		if (iE3 > -1) {
			if (!bE3) {
				SpeedCorOlde3 = feplus3 * potent[VXCOR][iE3] + (1.0 - feplus3)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde3 = potent[VXCOR][iE3];
			}
		}
		if (iN3 > -1) {
			if (!bN3) {
				SpeedCorOldn3 = fnplus3 * potent[VYCOR][iN3] + (1.0 - fnplus3)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn3 = potent[VYCOR][iN3];
			}
		}
		if (iT3 > -1) {
			if (!bT3) {
				SpeedCorOldt3 = ftplus3 * potent[VZCOR][iT3] + (1.0 - ftplus3)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt3 = potent[VZCOR][iT3];
			}
		}
		if (iW3 > -1) {
			if (!bW3) {
				SpeedCorOldw3 = fwplus3 * potent[VXCOR][iW3] + (1.0 - fwplus3)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw3 = potent[VXCOR][iW3];
			}
		}
		if (iS3 > -1) {
			if (!bS3) {
				SpeedCorOlds3 = fsplus3 * potent[VYCOR][iS3] + (1.0 - fsplus3)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds3 = potent[VYCOR][iS3];
			}
		}
		if (iB3 > -1) {
			if (!bB3) {
				SpeedCorOldb3 = fbplus3 * potent[VZCOR][iB3] + (1.0 - fbplus3)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb3 = potent[VZCOR][iB3];
			}
		}


		doublereal SpeedCorOlde4 = 0.0, SpeedCorOldw4 = 0.0, SpeedCorOldn4 = 0.0, SpeedCorOlds4 = 0.0, SpeedCorOldt4 = 0.0, SpeedCorOldb4 = 0.0;
		if (iE4 > -1) {
			if (!bE4) {
				SpeedCorOlde4 = feplus4 * potent[VXCOR][iE4] + (1.0 - feplus4)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOlde4 = potent[VXCOR][iE4];
			}
		}
		if (iN4 > -1) {
			if (!bN4) {
				SpeedCorOldn4 = fnplus4 * potent[VYCOR][iN4] + (1.0 - fnplus4)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOldn4 = potent[VYCOR][iN4];
			}
		}
		if (iT4 > -1) {
			if (!bT4) {
				SpeedCorOldt4 = ftplus4 * potent[VZCOR][iT4] + (1.0 - ftplus4)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldt4 = potent[VZCOR][iT4];
			}
		}
		if (iW4 > -1) {
			if (!bW4) {
				SpeedCorOldw4 = fwplus4 * potent[VXCOR][iW4] + (1.0 - fwplus4)*potent[VXCOR][iP];
			}
			else {
				SpeedCorOldw4 = potent[VXCOR][iW4];
			}
		}
		if (iS4 > -1) {
			if (!bS4) {
				SpeedCorOlds4 = fsplus4 * potent[VYCOR][iS4] + (1.0 - fsplus4)*potent[VYCOR][iP];
			}
			else {
				SpeedCorOlds4 = potent[VYCOR][iS4];
			}
		}
		if (iB4 > -1) {
			if (!bB4) {
				SpeedCorOldb4 = fbplus4 * potent[VZCOR][iB4] + (1.0 - fbplus4)*potent[VZCOR][iP];
			}
			else {
				SpeedCorOldb4 = potent[VZCOR][iB4];
			}
		}

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    // дополнительная поправка осуществляется на основе скоректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    Fe_sum +=(1.0-alpha[VX])*(mfcurrentretune[ESIDE] - rhoe * SpeedCorOlde*dSqe - rhoe2 * SpeedCorOlde2*dSqe2 - rhoe3 * SpeedCorOlde3*dSqe3 - rhoe4 * SpeedCorOlde4*dSqe4);
	    Fn_sum +=(1.0-alpha[VY])*(mfcurrentretune[NSIDE] - rhon * SpeedCorOldn*dSqn - rhon2 * SpeedCorOldn2*dSqn2 - rhon3 * SpeedCorOldn3*dSqn3 - rhon4 * SpeedCorOldn4*dSqn4);
	    Ft_sum +=(1.0-alpha[VZ])*(mfcurrentretune[TSIDE] - rhot * SpeedCorOldt*dSqt - rhot2 * SpeedCorOldt2*dSqt2 - rhot3 * SpeedCorOldt3*dSqt3 - rhot4 * SpeedCorOldt4*dSqt4);
	    Fw_sum +=(1.0-alpha[VX])*(mfcurrentretune[WSIDE] - rhow * SpeedCorOldw*dSqw - rhow2 * SpeedCorOldw2*dSqw2 - rhow3 * SpeedCorOldw3*dSqw3 - rhow4 * SpeedCorOldw4*dSqw4);
	    Fs_sum +=(1.0-alpha[VY])*(mfcurrentretune[SSIDE] - rhos * SpeedCorOlds*dSqs - rhos2 * SpeedCorOlds2*dSqs2 - rhos3 * SpeedCorOlds3*dSqs3 - rhos4 * SpeedCorOlds4*dSqs4);
	    Fb_sum +=(1.0-alpha[VZ])*(mfcurrentretune[BSIDE] - rhob * SpeedCorOldb*dSqb - rhob2 * SpeedCorOldb2*dSqb2 - rhob3 * SpeedCorOldb3*dSqb3 - rhob4 * SpeedCorOldb4*dSqb4);

	} 
	

	
	// возвращаем значение потока на грани КО.
	if (bdeltapfinish==true) {
	   // запоминаем финишный поток через грань.
	   mfcurrentretune[ESIDE]=Fe_sum;
	   mfcurrentretune[NSIDE]=Fn_sum;
	   mfcurrentretune[TSIDE]=Ft_sum;
	   mfcurrentretune[WSIDE]=Fw_sum;
	   mfcurrentretune[SSIDE]=Fs_sum;
	   mfcurrentretune[BSIDE]=Fb_sum;
	}
	
	

	sl[PAM][iP].b=(Fw_sum -Fe_sum +Fs_sum -Fn_sum +Fb_sum -Ft_sum + baddDFLUX2);
	//sl[PAM][iP].b = Fw_sum - Fe_sum + Fs_sum - Fn_sum + Fb_sum - Ft_sum; // 20.07.2016

	
	// Симметризация СЛАУ :
	/*
	// для поправки давления получается эллиптическое уравнение с SPD матрицей.
	
	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

		// 1. Учёт условия Дирихле:
		if ((iE>=maxelm) && (fabs(slb[PAM][iE-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].ae*slb[PAM][iE-maxelm].b/slb[PAM][iE-maxelm].aw;
			sl[PAM][iP].ae=0.0;
			sl[PAM][iP].iE=-1; // не входит в матрицу СЛАУ.
		}
		if ((iW>=maxelm) && (fabs(slb[PAM][iW-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].aw*slb[PAM][iW-maxelm].b/slb[PAM][iW-maxelm].aw;
			sl[PAM][iP].aw=0.0;
			sl[PAM][iP].iW=-1; // не входит в матрицу СЛАУ.
		}
		if ((iN>=maxelm) && (fabs(slb[PAM][iN-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].an*slb[PAM][iN-maxelm].b/slb[PAM][iN-maxelm].aw;
			sl[PAM][iP].an=0.0;
			sl[PAM][iP].iN=-1; // не входит в матрицу СЛАУ.
		}
		if ((iS>=maxelm) && (fabs(slb[PAM][iS-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].as*slb[PAM][iS-maxelm].b/slb[PAM][iS-maxelm].aw;
			sl[PAM][iP].as=0.0;
			sl[PAM][iP].iS=-1; // не входит в матрицу СЛАУ.
		}
		if ((iT>=maxelm) && (fabs(slb[PAM][iT-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].at*slb[PAM][iT-maxelm].b/slb[PAM][iT-maxelm].aw;
			sl[PAM][iP].at=0.0;
			sl[PAM][iP].iT=-1; // не входит в матрицу СЛАУ.
		}
		if ((iB>=maxelm) && (fabs(slb[PAM][iB-maxelm].ai)<eps)) {
			sl[PAM][iP].b+=sl[PAM][iP].ab*slb[PAM][iB-maxelm].b/slb[PAM][iB-maxelm].aw;
			sl[PAM][iP].ab=0.0;
			sl[PAM][iP].iB=-1; // не входит в матрицу СЛАУ.
		}
	
	*/

} // my_elmatr_quad_PAm3

// Цель данной функции скоректировать поправку давления на границе
// так чтобы она отвечала изменинию поправки внутри расчётной области, а 
// именно интерполировать поправку из центральных узлов на границу с помощью 
// квадратичной интерполяции.
void free_pressure(FLOW &f) {
	const integer LINEAR_INTERPOL=0;
	const integer QUAD_INTERPOL=1;
	const integer NO_INTERPOL = 2;

	integer imy_interpol = QUAD_INTERPOL;// QUAD_INTERPOL; // LINEAR_INTERPOL QUAD_INTERPOL

		for (integer iP=0; iP<f.maxelm; iP++) {
			// iP - номер центрального контрольного объёма
	        integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	        iE=f.sosedi[ESIDE][iP].iNODE1; iN=f.sosedi[NSIDE][iP].iNODE1; iT=f.sosedi[TSIDE][iP].iNODE1;
	        iW=f.sosedi[WSIDE][iP].iNODE1; iS=f.sosedi[SSIDE][iP].iNODE1; iB=f.sosedi[BSIDE][iP].iNODE1;
	        

            // если с одной из сторон граница расчётной области 
	        // то переменная равна true
	        bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

            if (iE>=f.maxelm) bE=true;
	        if (iN>=f.maxelm) bN=true;
	        if (iT>=f.maxelm) bT=true;
            if (iW>=f.maxelm) bW=true;
	        if (iS>=f.maxelm) bS=true;
	        if (iB>=f.maxelm) bB=true;

	        // вычисление размеров текущего контрольного объёма:
	        doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
            volume3D(iP, f.nvtx, f.pa, dx, dy, dz);
	    

    	
			// Линейная  и квадратичная интерполяции.
            if (bW) {
				// узел W граничный
				TOCHKA pp,pb,pbb;
				integer iEE=f.sosedi[ESIDE][iE].iNODE1;
				switch (imy_interpol) {
				case LINEAR_INTERPOL : 
		                   center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                   center_cord3D(iE, f.nvtx, f.pa, pb,ESIDE);
		                   f.potent[PAM][iW]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iE], pp.x, pb.x, pp.x-0.5*dx);
						   break;
				case QUAD_INTERPOL :
					        center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                    center_cord3D(iE, f.nvtx, f.pa, pb,ESIDE);
						    center_cord3D(iEE, f.nvtx, f.pa, pbb,EE);
							f.potent[PAM][iW]=my_quadratic_interpolation('-', f.potent[PAM][iEE], f.potent[PAM][iE], f.potent[PAM][iP],
	                            pbb.x , pb.x, pp.x, pp.x-0.5*dx);
					       break;
				case NO_INTERPOL: f.potent[PAM][iW] = f.potent[PAM][iP];
					break;
				default :
					// узел W граничный
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iE, f.nvtx, f.pa, pb,ESIDE);
		            f.potent[PAM][iW]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iE], pp.x, pb.x, pp.x-0.5*dx);
					break;
				}
	        }
	        if (bE) {
		        // узел E граничный
				TOCHKA pp,pb,pbb;
				integer iWW=f.sosedi[WSIDE][iW].iNODE1;
				switch (imy_interpol) {
				case LINEAR_INTERPOL :
		                  center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                  center_cord3D(iW, f.nvtx, f.pa, pb,WSIDE);
		                  f.potent[PAM][iE]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iW], pp.x, pb.x, pp.x+0.5*dx);
						  break;
				case QUAD_INTERPOL :
					      center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                  center_cord3D(iW, f.nvtx, f.pa, pb,WSIDE);
						  center_cord3D(iWW, f.nvtx, f.pa, pbb,WW);
						  
						  f.potent[PAM][iE]=my_quadratic_interpolation('+', f.potent[PAM][iWW], f.potent[PAM][iW], f.potent[PAM][iP], pbb.x, pb.x, pp.x, pp.x+0.5*dx);	  
					break;
				case NO_INTERPOL: f.potent[PAM][iE] = f.potent[PAM][iP];
					break;
				default : 
					center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iW, f.nvtx, f.pa, pb,WSIDE);
		            f.potent[PAM][iE]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iW], pp.x, pb.x, pp.x+0.5*dx);  
					break;
				}
	        }
	        if (bS) {
		        // узел S граничный
				TOCHKA pp,pb,pbb;
				integer iNN=f.sosedi[NSIDE][iN].iNODE1;
				switch (imy_interpol) {
                case LINEAR_INTERPOL :
		                 center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                 center_cord3D(iN, f.nvtx, f.pa, pb,NSIDE);
		                 f.potent[PAM][iS]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iN], pp.y, pb.y, pp.y-0.5*dy);
				         break;
				case QUAD_INTERPOL :
		                 center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                 center_cord3D(iN, f.nvtx, f.pa, pb,NSIDE);
						 center_cord3D(iNN, f.nvtx, f.pa, pbb,NN);
						 f.potent[PAM][iS]=my_quadratic_interpolation('-', f.potent[PAM][iNN], f.potent[PAM][iN], f.potent[PAM][iP],
	                            pbb.y, pb.y, pp.y, pp.y-0.5*dy);
					break;
				case NO_INTERPOL: f.potent[PAM][iS] = f.potent[PAM][iP];
					break;
				default :
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iN, f.nvtx, f.pa, pb,NSIDE);
		            f.potent[PAM][iS]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iN], pp.y, pb.y, pp.y-0.5*dy);
					break;
				}
	        }
	        if (bN) {
		        // узел N граничный
				TOCHKA pp,pb,pbb;
				integer iSS=f.sosedi[SSIDE][iS].iNODE1;
				switch (imy_interpol) {
					case LINEAR_INTERPOL :
		                 center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                 center_cord3D(iS, f.nvtx, f.pa, pb,SSIDE);
		                 f.potent[PAM][iN]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iS], pp.y, pb.y, pp.y+0.5*dy);
				    break;
					case QUAD_INTERPOL :
		                center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                center_cord3D(iS, f.nvtx, f.pa, pb,SSIDE);
						center_cord3D(iSS, f.nvtx, f.pa, pbb,SS);
						f.potent[PAM][iN]=my_quadratic_interpolation('+', f.potent[PAM][iSS], f.potent[PAM][iS], f.potent[PAM][iP],
	                            pbb.y , pb.y, pp.y, pp.y+0.5*dy);
					break;
					case NO_INTERPOL: f.potent[PAM][iN] = f.potent[PAM][iP];
						break;
					default :
						center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                center_cord3D(iS, f.nvtx, f.pa, pb,SSIDE);
		                f.potent[PAM][iN]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iS], pp.y, pb.y, pp.y+0.5*dy);
					break;
				}
	        }
	       if (bB) {
		        // узел B граничный
			   TOCHKA pp,pb,pbb;
			   integer iTT=f.sosedi[TSIDE][iT].iNODE1;
				switch (imy_interpol) {
		        case LINEAR_INTERPOL :
		               center_cord3D(iP, f.nvtx, f.pa, pp,100);
		               center_cord3D(iT, f.nvtx, f.pa, pb,TSIDE);
		               f.potent[PAM][iB]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iT], pp.z, pb.z, pp.z-0.5*dz);
				   break;
				case QUAD_INTERPOL :
					center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iT, f.nvtx, f.pa, pb,TSIDE);
					center_cord3D(iTT, f.nvtx, f.pa, pbb,TTSIDE);
					f.potent[PAM][iB]=my_quadratic_interpolation('-', f.potent[PAM][iTT], f.potent[PAM][iT], f.potent[PAM][iP],
	                            pbb.z, pb.z, pp.z, pp.z-0.5*dz);

					break;
				case NO_INTERPOL: f.potent[PAM][iB] = f.potent[PAM][iP];
					break;
				default :
					center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iT, f.nvtx, f.pa, pb,TSIDE);
		            f.potent[PAM][iB]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iT], pp.z, pb.z, pp.z-0.5*dz);
					break;
				}
	       }
	       if (bT) { 
		        // узел T граничный
		        TOCHKA pp,pb, pbb;
				integer iBB=f.sosedi[BSIDE][iB].iNODE1;
				switch (imy_interpol) {
				case LINEAR_INTERPOL :
		              center_cord3D(iP, f.nvtx, f.pa, pp,100);
		              center_cord3D(iB, f.nvtx, f.pa, pb,BSIDE);
		              f.potent[PAM][iT]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iB], pp.z, pb.z, pp.z+0.5*dz);
					  break;
                case QUAD_INTERPOL :
					  center_cord3D(iP, f.nvtx, f.pa, pp,100);
		              center_cord3D(iB, f.nvtx, f.pa, pb,BSIDE);
					  center_cord3D(iBB, f.nvtx, f.pa, pbb,BB);
					  f.potent[PAM][iT]=my_quadratic_interpolation('+', f.potent[PAM][iBB], f.potent[PAM][iB], f.potent[PAM][iP],
	                            pbb.z , pb.z, pp.z, pp.z+0.5*dz);
					break;
				case NO_INTERPOL: f.potent[PAM][iT] = f.potent[PAM][iP];
					break;
				default :
					  center_cord3D(iP, f.nvtx, f.pa, pp,100);
		              center_cord3D(iB, f.nvtx, f.pa, pb,BSIDE);
		              f.potent[PAM][iT]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iB], pp.z, pb.z, pp.z+0.5*dz);
					break;
				}
	       }
	
		
	}
} // free_pressure

// полезное действие этой функции сомнительно.
void correctpressureoutlet(FLOW &f, integer lw, integer ls, WALL* w) {
	// экстраполируем давление изнутри расчётной области
	// на границу.
	free_pressure(f);

	// Минимальное значение давления должно достигаться на выходной границе расчётной области.

	doublereal pamout=0.0;
	for (integer iP=0; iP<f.maxelm; iP++) {
			// iP - номер центрального контрольного объёма
	        integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	        iE=f.sosedi[ESIDE][iP].iNODE1; iN=f.sosedi[NSIDE][iP].iNODE1; iT=f.sosedi[TSIDE][iP].iNODE1;
	        iW=f.sosedi[WSIDE][iP].iNODE1; iS=f.sosedi[SSIDE][iP].iNODE1; iB=f.sosedi[BSIDE][iP].iNODE1;
	        

            // если с одной из сторон граница расчётной области 
	        // то переменная равна true
	        bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

            if (iE>=f.maxelm) bE=true;
	        if (iN>=f.maxelm) bN=true;
	        if (iT>=f.maxelm) bT=true;
            if (iW>=f.maxelm) bW=true;
	        if (iS>=f.maxelm) bS=true;
	        if (iB>=f.maxelm) bB=true;



	        integer inumber;

			if (bW) {
		       // узел W граничный
				inumber=iW-f.maxelm;
		       if ((f.sosedb[inumber].MCB<(ls+lw)) && (f.sosedb[inumber].MCB>=ls) && (w[f.sosedb[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iW],pamout);
			   }
	        }
	        if (bE) {
		        // узел E граничный
		        inumber=iE-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls+lw)) && (f.sosedb[inumber].MCB>=ls) && (w[f.sosedb[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iE],pamout);
			   }
	        }
	        if (bS) {
		        // узел S граничный
		        inumber=iS-f.maxelm;
                if ((f.sosedb[inumber].MCB<(ls+lw)) && (f.sosedb[inumber].MCB>=ls) && (w[f.sosedb[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iS],pamout);
			   }
	        }
	        if (bN) {
		        // узел N граничный
		        inumber=iN-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls+lw)) && (f.sosedb[inumber].MCB>=ls) && (w[f.sosedb[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iN],pamout);
			    }
	        }
	       if (bB) {
		        // узел B граничный
		        inumber=iB-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls+lw)) && (f.sosedb[inumber].MCB>=ls) && (w[f.sosedb[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iB],pamout);
			   }
	       }
	       if (bT) { 
		        // узел T граничный
		        inumber=iT-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls+lw)) && (f.sosedb[inumber].MCB>=ls) && (w[f.sosedb[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iT],pamout);
			   }
	       }
	
		
	}

	for (integer i=0; i<(f.maxelm+f.maxbound); i++) f.potent[PAM][i]-=pamout;

} // correctpressureoutlet

// Интерполяция скорости изнутри расчётной области на границу.
// Это избавляет от однородных условий Неймана.
void interpolatevel(FLOW &f, integer lw, integer ls, WALL* w, integer iVar) {
	// экстраполируем давление изнутри расчётной области
	// на границу.
	free_pressure(f);

	
	for (integer iP=0; iP<f.maxelm; iP++) {
			// iP - номер центрального контрольного объёма
	        integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	        iE=f.sosedi[ESIDE][iP].iNODE1; iN=f.sosedi[NSIDE][iP].iNODE1; iT=f.sosedi[TSIDE][iP].iNODE1;
	        iW=f.sosedi[WSIDE][iP].iNODE1; iS=f.sosedi[SSIDE][iP].iNODE1; iB=f.sosedi[BSIDE][iP].iNODE1;
	        

            // если с одной из сторон граница расчётной области 
	        // то переменная равна true
	        bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

            if (iE>=f.maxelm) bE=true;
	        if (iN>=f.maxelm) bN=true;
	        if (iT>=f.maxelm) bT=true;
            if (iW>=f.maxelm) bW=true;
	        if (iS>=f.maxelm) bS=true;
	        if (iB>=f.maxelm) bB=true;

			 // вычисление размеров текущего контрольного объёма:
	        doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
            volume3D(iP, f.nvtx, f.pa, dx, dy, dz);

	        integer inumber;

			if (bW) {
		       // узел W граничный
				inumber=iW-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && ((w[f.sosedb[inumber].MCB - ls].bpressure || w[f.sosedb[inumber].MCB - ls].bopening))) {
				   // Линейная интерполяция
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iE, f.nvtx, f.pa, pb,ESIDE);
		           f.potent[iVar][iW]=my_linear_interpolation('-', f.potent[iVar][iP], f.potent[iVar][iE], pp.x, pb.x, pp.x-0.5*dx);
			   }
	        }
	        if (bE) {
		        // узел E граничный
		        inumber=iE-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && ((w[f.sosedb[inumber].MCB - ls].bpressure || w[f.sosedb[inumber].MCB - ls].bopening))) {
					// Линейная интерполяция
				    TOCHKA pp,pb;
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iW, f.nvtx, f.pa, pb,WSIDE);
		            f.potent[iVar][iE]=my_linear_interpolation('+', f.potent[iVar][iP], f.potent[iVar][iW], pp.x, pb.x, pp.x+0.5*dx);
			   }
	        }
	        if (bS) {
		        // узел S граничный
		        inumber=iS-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && ((w[f.sosedb[inumber].MCB - ls].bpressure || w[f.sosedb[inumber].MCB - ls].bopening))) {
				    TOCHKA pp,pb;
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iN, f.nvtx, f.pa, pb,NSIDE);
		            f.potent[iVar][iS]=my_linear_interpolation('-', f.potent[iVar][iP], f.potent[iVar][iN], pp.y, pb.y, pp.y-0.5*dy);
			   }
	        }
	        if (bN) {
		        // узел N граничный
		        inumber=iN-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && ((w[f.sosedb[inumber].MCB - ls].bpressure || w[f.sosedb[inumber].MCB - ls].bopening))) {
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iS, f.nvtx, f.pa, pb,SSIDE);
		           f.potent[iVar][iN]=my_linear_interpolation('+', f.potent[iVar][iP], f.potent[iVar][iS], pp.y, pb.y, pp.y+0.5*dy);
			    }
	        }
	       if (bB) {
		        // узел B граничный
		        inumber=iB-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && ((w[f.sosedb[inumber].MCB - ls].bpressure || w[f.sosedb[inumber].MCB - ls].bopening))) {
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iT, f.nvtx, f.pa, pb,TSIDE);
		           f.potent[iVar][iB]=my_linear_interpolation('-', f.potent[iVar][iP], f.potent[iVar][iT], pp.z, pb.z, pp.z-0.5*dz);
			   }
	       }
	       if (bT) { 
		        // узел T граничный
		        inumber=iT-f.maxelm;
				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && ((w[f.sosedb[inumber].MCB - ls].bpressure ||  w[f.sosedb[inumber].MCB - ls].bopening))) {
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iB, f.nvtx, f.pa, pb,BSIDE);
		           f.potent[iVar][iT]=my_linear_interpolation('+', f.potent[iVar][iP], f.potent[iVar][iB], pp.z, pb.z, pp.z+0.5*dz);
			   }
	       }
	
		
	}

} // interpolatevel

#endif
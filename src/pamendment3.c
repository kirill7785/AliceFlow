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
#include "my_elmatr_quad_f3d.cpp"
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
							BOUND* &border_neighbor, int** nvtx, bool bPfix,
							doublereal dbeta, TOCHKA* pa, doublereal** potent,
							float** prop, float** prop_b, doublereal* alpha,
							integer ls, integer lw, WALL* w, bool bDirichlet, 
							int*** neighbors_for_the_internal_node, doublereal **diag_coef, doublereal RCh,
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
		   slb[PAM][inumber].iW=border_neighbor[inumber].iB; 

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

		if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure)) {
            
			if (!bFReeStyle) {
				doublereal rCe=0.0; // значение давления на выходной границе
			    integer ioutflowcondition=BERNULLI; // PRESSUREOUTLET // BERNULLI И.Ю.Чумаков
			    doublereal kineticenergy=0.0; // кинетическая энергия потока: 0.5*rho*Vel!2.
			    doublereal rsign=1.0; // >0 скорость направлена из расчётной области наружу, если <0 то внутрь и используется соотношение Бернулли.
			
			    if (bDirichlet) {

			    	switch (ioutflowcondition) {
				     case PRESSUREOUTLET: // поправка давления фиксированна и равна нулю:
				                      // простейший способ - традиционно используется во многих расчётах.
				                      // Но только не при наличии рециркуляционных зон на границе.
		                              slb[PAM][inumber].aw=1.0;
		                              slb[PAM][inumber].ai=0.0;
		                              slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                              slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                              slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					                  break;
			      	 case BERNULLI: // Рассматривается выходная граница потока и условие для поправки давления для неё.
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
	                      switch (border_neighbor[inumber].Norm) {
						  case E_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case W_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case N_SIDE: if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case S_SIDE:if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case T_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case B_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]<0.0) {
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
		                     slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                     slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                     slb[PAM][inumber].iW=border_neighbor[inumber].iB; 
						  }
						  else {
							  kineticenergy=0.5*prop_b[RHO][border_neighbor[inumber].iB-maxelm];
							  kineticenergy*=(potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]*potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]+
								              potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]*potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]+
											  potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]*potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]);
							  slb[PAM][inumber].aw=1.0;
		                      slb[PAM][inumber].ai=0.0;
		                      slb[PAM][inumber].b=rCe-kineticenergy-potent[PRESS][border_neighbor[inumber].iB]; // Соотношение Бернулли
		                      slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                      slb[PAM][inumber].iW=border_neighbor[inumber].iB; 
						  }
					      break;
				     default: // поправка давления фиксированна и равна нулю:
				          // простейший способ - традиционно используется во многих расчётах.
				          // Но только не при наличии рециркуляционных зон на границе.
		                  slb[PAM][inumber].aw=1.0;
		                  slb[PAM][inumber].ai=0.0;
		                  slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                  slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                  slb[PAM][inumber].iW=border_neighbor[inumber].iB;
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
	        switch (border_neighbor[inumber].Norm) {
		        case E_SIDE:
			    
				    dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_X_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_X_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI] + fiplus*prop[RHO][border_neighbor[inumber].iII]);// проверено !
					aUPi=(sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_X_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_X_COMPONENT]);

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

					//slb[PAM][inumber].b+=rhoi*dS*(potent[VX][border_neighbor[inumber].iB]-potent[VX][border_neighbor[inumber].iI]);
					//printf("iII=%d, iI=%d, iB=%d, maxelm=%d",border_neighbor[inumber].iII,border_neighbor[inumber].iI,border_neighbor[inumber].iB,maxelm);
					//getchar(); // debug внутренняя нормаль
					  //FgRhieChow=rFgRhieChow_internal_border(border_neighbor[inumber].iI, WSIDE, prop_b[RHO][border_neighbor[inumber].iB-maxelm], alpha[VX], nvtx, neighbors_for_the_internal_node, maxelm, potent[PRESS], pa, diag_coef);
					//slb[PAM][inumber].b+=prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*0.5*(potent[VX][border_neighbor[inumber].iB]+potent[VX][border_neighbor[inumber].iI]);
					  //slb[PAM][inumber].b+=prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*potent[VX][border_neighbor[inumber].iB];
					  //slb[PAM][inumber].b+=RCh*FgRhieChow;
					  // данная добавка по-видимому приводит к расходимости.
					   //printf("RCh=%e, Fe=%e\n",RCh*FgRhieChow,slb[PAM][inumber].b+=prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*potent[VX][border_neighbor[inumber].iB]);
					   // getchar();
				    break;			
			
		        case N_SIDE:

                    dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Y_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Y_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI] + fiplus*prop[RHO][border_neighbor[inumber].iII]);// проверено !
					aUPi=(sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Y_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Y_COMPONENT]);

				    break;

			    case T_SIDE: 

                    dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Z_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI] + fiplus*prop[RHO][border_neighbor[inumber].iII]);// проверено !
					aUPi=(sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Z_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					
				    break;

			    case W_SIDE:

                    dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_X_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_X_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]);// проверено !
					aUPi=(sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_X_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_X_COMPONENT]);
					
					break;

		        case S_SIDE:

                    dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Y_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Y_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI] + fiplus*prop[RHO][border_neighbor[inumber].iII]);// проверено !
					aUPi=(sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap);// проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Y_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Y_COMPONENT]);
					
				   break;

		        case B_SIDE: 

                    dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Z_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Z_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					
				    break;
	        } // switch

			//*/
            integer j,l,xitem,k;
			// сортировка по возрастанию
			for (j=0; j<5; j++) {
				k=j; xitem=border_neighbor[inumber].iW[j];
				for (l=j+1; l<6; l++) {
					if (border_neighbor[inumber].iW[l] < xitem) {
						k=l; xitem=border_neighbor[inumber].iW[k];
					}
				}
                border_neighbor[inumber].iW[k]=border_neighbor[inumber].iW[j];
				border_neighbor[inumber].iW[j]=xitem;
			}

            j=0; l=0;
			while (border_neighbor[inumber].iW[j]==(-1)) j++;

			if (j<6) { slb[PAM][inumber].iW1=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW2=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW3=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW4=border_neighbor[inumber].iW[j++]; l++; } 

			switch (l) {
				case 0: slb[PAM][inumber].iW1=-1;
		                 slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
		                 break;
				case 1: slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 2: slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 3: slb[PAM][inumber].iW4=-1;
						 break;
			}

			
				}
			}

	    }
		 else if ((!bDirichlet) && !((inumber == (maxbound - 1)) && bPfix) && !((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure )))
	    {  
			// Не условие Дирихле, Не фиксация давления в точке, Не выходная граница
			//printf("neiman!\n"); getchar(); // debug

			/*
            // поправка давления фиксированно равна нулю:
		    slb[PAM][inumber].aw=1.0;
		    slb[PAM][inumber].ai=0.0;
		    slb[PAM][inumber].b=0.0;
			slb[PAM][inumber].iI=-1;//border_neighbor[inumber].iI; //  присутствует в матрице
		    slb[PAM][inumber].iW=border_neighbor[inumber].iB;
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
	        switch (border_neighbor[inumber].Norm) {
		        case E_SIDE:
			    
				    dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_X_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_X_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_X_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_X_COMPONENT]);

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

					//slb[PAM][inumber].b+=rhoi*dS*(potent[VX][border_neighbor[inumber].iB]-potent[VX][border_neighbor[inumber].iI]);
					//printf("iII=%d, iI=%d, iB=%d, maxelm=%d",border_neighbor[inumber].iII,border_neighbor[inumber].iI,border_neighbor[inumber].iB,maxelm);
					//getchar(); // debug внутренняя нормаль
					  //FgRhieChow=rFgRhieChow_internal_border(border_neighbor[inumber].iI, WSIDE, prop_b[RHO][border_neighbor[inumber].iB-maxelm], alpha[VX], nvtx, neighbors_for_the_internal_node, maxelm, potent[PRESS], pa, diag_coef);
					//slb[PAM][inumber].b+=prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*0.5*(potent[VX][border_neighbor[inumber].iB]+potent[VX][border_neighbor[inumber].iI]);
					  //slb[PAM][inumber].b+=prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*potent[VX][border_neighbor[inumber].iB];
					  //slb[PAM][inumber].b+=RCh*FgRhieChow;
					  // данная добавка по-видимому приводит к расходимости.
					   //printf("RCh=%e, Fe=%e\n",RCh*FgRhieChow,slb[PAM][inumber].b+=prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*potent[VX][border_neighbor[inumber].iB]);
					   // getchar();
				    break;			
			
		        case N_SIDE:

                    dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Y_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Y_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Y_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Y_COMPONENT]);

				    break;

			    case T_SIDE: 

                    dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Z_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Z_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					
				    break;

			    case W_SIDE:

                    dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_X_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_X_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI] + fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_X_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_X_COMPONENT]);
					
					break;

		        case S_SIDE:

                    dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Y_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Y_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Y_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Y_COMPONENT]);
					
				   break;

		        case B_SIDE: 

                    dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
					slb[PAM][inumber].ai=dbeta*alpha[VELOCITY_Z_COMPONENT]*prop_b[RHO][border_neighbor[inumber].iB-maxelm]*dS*dS/slb[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB-maxelm].aw;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].ai/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					rhoi=(prop[RHO][border_neighbor[inumber].iI]*prop[RHO][border_neighbor[inumber].iII]);
					rhoi=rhoi/((1.0-fiplus)*prop[RHO][border_neighbor[inumber].iI]+ fiplus*prop[RHO][border_neighbor[inumber].iII]); // проверено !
					aUPi=(sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap);
					aUPi=aUPi/((1.0-fiplus)*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI].ap + fiplus*sl[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII].ap); // проверено !
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*alpha[VELOCITY_Z_COMPONENT]*rhoi*dS*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/aUPi;
					if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) slb[PAM][inumber].b/=(1.0-alpha[VELOCITY_Z_COMPONENT]);
					
				    break;
	        } // switch

			//*/
            integer j,l,xitem,k;
			// сортировка по возрастанию
			for (j=0; j<5; j++) {
				k=j; xitem=border_neighbor[inumber].iW[j];
				for (l=j+1; l<6; l++) {
					if (border_neighbor[inumber].iW[l] < xitem) {
						k=l; xitem=border_neighbor[inumber].iW[k];
					}
				}
                border_neighbor[inumber].iW[k]=border_neighbor[inumber].iW[j];
				border_neighbor[inumber].iW[j]=xitem;
			}

            j=0; l=0;
			while (border_neighbor[inumber].iW[j]==(-1)) j++;

			if (j<6) { slb[PAM][inumber].iW1=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW2=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW3=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW4=border_neighbor[inumber].iW[j++]; l++; } 

			switch (l) {
				case 0: slb[PAM][inumber].iW1=-1;
		                 slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
		                 break;
				case 1: slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 2: slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 3: slb[PAM][inumber].iW4=-1;
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
							BOUND* &border_neighbor, int** nvtx, bool bPfix,
							doublereal dbeta, TOCHKA* pa, doublereal** potent,
							float** prop, float** prop_b, doublereal* alpha,
							integer ls, integer lw, WALL* w, bool bDirichlet, 
							int*** neighbors_for_the_internal_node, doublereal RCh,
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
		   slb[PAM][inumber].iW=border_neighbor[inumber].iB; 

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

		if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure )) {
            
			if (!bFReeStyle) {
				doublereal rCe=0.0; // значение давления на выходной границе
			    integer ioutflowcondition=PRESSUREOUTLET; // PRESSUREOUTLET // BERNULLI И.Ю.Чумаков
			    doublereal kineticenergy=0.0; // кинетическая энергия потока: 0.5*rho*Vel!2.
			    doublereal rsign=1.0; // >0 скорость направлена из расчётной области наружу, если <0 то внутрь и используется соотношение Бернулли.
			
			    if (bDirichlet) {

			    	switch (ioutflowcondition) {
				     case PRESSUREOUTLET: // поправка давления фиксированна и равна нулю:
				                      // простейший способ - традиционно используется во многих расчётах.
				                      // Но только не при наличии рециркуляционных зон на границе.
		                              slb[PAM][inumber].aw=1.0;
		                              slb[PAM][inumber].ai=0.0;
		                              slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                              slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                              slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					                  break;
			      	 case BERNULLI: // Рассматривается выходная граница потока и условие для поправки давления для неё.
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
	                      switch (border_neighbor[inumber].Norm) {
						  case E_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case W_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case N_SIDE: if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case S_SIDE:if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case T_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case B_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]<0.0) {
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
		                     slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                     slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                     slb[PAM][inumber].iW=border_neighbor[inumber].iB; 
						  }
						  else {
							  kineticenergy=0.5*prop_b[RHO][border_neighbor[inumber].iB-maxelm];
							  kineticenergy*=(potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]*potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]+
								              potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]*potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]+
											  potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]*potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]);
							  slb[PAM][inumber].aw=1.0;
		                      slb[PAM][inumber].ai=0.0;
		                      slb[PAM][inumber].b=rCe-kineticenergy-potent[PRESS][border_neighbor[inumber].iB]; // Соотношение Бернулли
		                      slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                      slb[PAM][inumber].iW=border_neighbor[inumber].iB; 
						  }
					      break;
				     default: // поправка давления фиксированна и равна нулю:
				          // простейший способ - традиционно используется во многих расчётах.
				          // Но только не при наличии рециркуляционных зон на границе.
		                  slb[PAM][inumber].aw=1.0;
		                  slb[PAM][inumber].ai=0.0;
		                  slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                  slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                  slb[PAM][inumber].iW=border_neighbor[inumber].iB;
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
	                switch (border_neighbor[inumber].Norm) {
		              case E_SIDE:
			    
				          dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					      dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
					      dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					      deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					      deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;

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
			
		              case N_SIDE:

                          dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					      dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					      deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					      deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+ fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					

				          break;

			          case T_SIDE: 

                          dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					      dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					      deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					      deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
				          break;

			          case W_SIDE:

                          dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					      dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y; 
					      dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					      deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					      deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
					      break;

		              case S_SIDE:

                          dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					      dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					      deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					      deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
					
				          break;

		              case B_SIDE: 

                          dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					      dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					      dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани

					      slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					      slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					      slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					      slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					      deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					      deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                          fiplus=0.5*dl/deltal;

					      taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					      taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					      
			        
                          // правая часть:  
					      slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					      
					
				          break;
	                } // switch

			        //*/
                    integer j,l,xitem,k;
			        // сортировка по возрастанию
			        for (j=0; j<5; j++) {
				        k=j; xitem=border_neighbor[inumber].iW[j];
				        for (l=j+1; l<6; l++) {
					        if (border_neighbor[inumber].iW[l] < xitem) {
						       k=l; xitem=border_neighbor[inumber].iW[k];
					        }
				        }
                        border_neighbor[inumber].iW[k]=border_neighbor[inumber].iW[j];
				        border_neighbor[inumber].iW[j]=xitem;
			        }

                    j=0; l=0;
			        while (border_neighbor[inumber].iW[j]==(-1)) j++;

			        if (j<6) { slb[PAM][inumber].iW1=border_neighbor[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW2=border_neighbor[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW3=border_neighbor[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW4=border_neighbor[inumber].iW[j++]; l++; } 

			        switch (l) {
				       case 0: slb[PAM][inumber].iW1=-1;
		                        slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
		                        break;
				       case 1: slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 2: slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 3: slb[PAM][inumber].iW4=-1;
						        break;
			        }

			
				}
			}

	    }
	     else if ((!bDirichlet) && !((inumber==(maxbound-1)) && bPfix) && !((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure) )
	    {  
			// Не условие Дирихле, Не фиксация давления в точке, Не выходная граница
			//printf("neiman!\n"); getchar(); // debug


			/*
            // поправка давления фиксированно равна нулю:
		    slb[PAM][inumber].aw=1.0;
		    slb[PAM][inumber].ai=0.0;
		    slb[PAM][inumber].b=0.0;
			slb[PAM][inumber].iI=-1;//border_neighbor[inumber].iI; //  присутствует в матрице
		    slb[PAM][inumber].iW=border_neighbor[inumber].iB;
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
	        switch (border_neighbor[inumber].Norm) {
		        case E_SIDE:
			    
				    dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;

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
			
		        case N_SIDE:

                    dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;

				    break;

			    case T_SIDE: 

                    dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
								        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
				    break;

			    case W_SIDE:

                    dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
					dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
					deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                    fiplus=0.5*dl/deltal;

					taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
					break;

		        case S_SIDE:

                    dl=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
					deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                    fiplus=0.5*dl/deltal;

					taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
				   break;

		        case B_SIDE: 

                    dl=pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
					dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
					dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани

					slb[PAM][inumber].ai=2.0*dbeta*tau[border_neighbor[inumber].iB]*dS/dl;
					slb[PAM][inumber].iI=border_neighbor[inumber].iI;
					slb[PAM][inumber].aw=slb[PAM][inumber].ai;
					slb[PAM][inumber].iW=border_neighbor[inumber].iB;
					
					deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
					deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                    fiplus=0.5*dl/deltal;

					taui=(tau[border_neighbor[inumber].iI]*tau[border_neighbor[inumber].iII]);
					taui=taui/((1.0-fiplus)*tau[border_neighbor[inumber].iI]+fiplus*tau[border_neighbor[inumber].iII]);
					
			        
                    // правая часть:  
					slb[PAM][inumber].b=(dbeta-1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI]-potent[PAM][border_neighbor[inumber].iII])/deltal;
					
				    break;
	        } // switch

			//*/
            integer j,l,xitem,k;
			// сортировка по возрастанию
			for (j=0; j<5; j++) {
				k=j; xitem=border_neighbor[inumber].iW[j];
				for (l=j+1; l<6; l++) {
					if (border_neighbor[inumber].iW[l] < xitem) {
						k=l; xitem=border_neighbor[inumber].iW[k];
					}
				}
                border_neighbor[inumber].iW[k]=border_neighbor[inumber].iW[j];
				border_neighbor[inumber].iW[j]=xitem;
			}

            j=0; l=0;
			while (border_neighbor[inumber].iW[j]==(-1)) j++;

			if (j<6) { slb[PAM][inumber].iW1=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW2=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW3=border_neighbor[inumber].iW[j++]; l++; }
			if (j<6) { slb[PAM][inumber].iW4=border_neighbor[inumber].iW[j++]; l++; } 

			switch (l) {
				case 0: slb[PAM][inumber].iW1=-1;
		                 slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
		                 break;
				case 1: slb[PAM][inumber].iW2=-1;
		                 slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 2: slb[PAM][inumber].iW3=-1;
		                 slb[PAM][inumber].iW4=-1;
						 break;
				case 3: slb[PAM][inumber].iW4=-1;
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
							BOUND* &border_neighbor, int** nvtx, bool bPfix,
							doublereal dbeta, TOCHKA* pa, doublereal** potent,
							float** prop, float** prop_b, doublereal* alpha,
							integer ls, integer lw, WALL* w, bool bDirichlet, 
							int*** neighbors_for_the_internal_node, doublereal RCh,
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
		   slb[PAM][inumber].iW=border_neighbor[inumber].iB; 

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

	if ((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB<(ls + lw)) && (w[border_neighbor[inumber].MCB - ls].bpressure || w[border_neighbor[inumber].MCB - ls].bopening)) {
            
			if (!bFReeStyle) {
				


				doublereal rCe=w[border_neighbor[inumber].MCB-ls].P; // значение давления на выходной границе
				//printf("rCe=%e %d \n",rCe,border_neighbor[inumber].MCB-ls);
				//getchar();

				integer ioutflowcondition =  BERNULLI; // PRESSUREOUTLET // BERNULLI И.Ю.Чумаков
			    doublereal kineticenergy=0.0; // кинетическая энергия потока: 0.5*rho*Vel!2.
			    doublereal rsign=1.0; // >0 скорость направлена из расчётной области наружу, если <0 то внутрь и используется соотношение Бернулли.
			
			    if (bDirichlet) {

					bcheck = true;

			    	switch (ioutflowcondition) {
				     case PRESSUREOUTLET: // поправка давления фиксированна и равна нулю:
				                      // простейший способ - традиционно используется во многих расчётах.
				                      // Но только не при наличии рециркуляционных зон на границе.
		                              slb[PAM][inumber].aw=1.0;
		                              slb[PAM][inumber].ai=0.0;
		                              slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                              slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                              slb[PAM][inumber].iW=border_neighbor[inumber].iB;
									 // printf("%d\n",slb[PAM][inumber].iW); getchar();// больше либо равно maxelm
					                  break;
			      	 case BERNULLI: // Рассматривается выходная граница потока и условие для поправки давления для неё.
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
	                      switch (border_neighbor[inumber].Norm) {
						  case E_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case W_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case N_SIDE: if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case S_SIDE:if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB]<0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case T_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]>0.0) {
							          // жидкость втекает внутрь расчётной области через выходную границу.
							          rsign=-1.0;
									  breversedflow=true;
								   } else rsign=1.0;
							       break;
						  case B_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB]<0.0) {
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
		                     slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                     slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                     slb[PAM][inumber].iW=border_neighbor[inumber].iB; 
						  }
						  else {
							  // Жидкость проникает внутрь расчётной области.

							  kineticenergy=0.5*prop_b[RHO][border_neighbor[inumber].iB-maxelm];
							  // VComponentCOR вместо VComponent не повлияло.
							  kineticenergy*=(potent[VXCOR][border_neighbor[inumber].iB]*potent[VXCOR][border_neighbor[inumber].iB]+
								              potent[VYCOR][border_neighbor[inumber].iB]*potent[VYCOR][border_neighbor[inumber].iB]+
											  potent[VZCOR][border_neighbor[inumber].iB]*potent[VZCOR][border_neighbor[inumber].iB]);
							  slb[PAM][inumber].aw=1.0;
		                      slb[PAM][inumber].ai=0.0;
		                      slb[PAM][inumber].b=rCe-kineticenergy-potent[PRESS][border_neighbor[inumber].iB]; // Соотношение Бернулли
		                      slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                      slb[PAM][inumber].iW=border_neighbor[inumber].iB; 
						  }
					      break;
				     default: // поправка давления фиксированна и равна нулю:
				          // простейший способ - традиционно используется во многих расчётах.
				          // Но только не при наличии рециркуляционных зон на границе.
		                  slb[PAM][inumber].aw=1.0;
		                  slb[PAM][inumber].ai=0.0;
		                  slb[PAM][inumber].b=rCe-potent[PRESS][border_neighbor[inumber].iB];
		                  slb[PAM][inumber].iI=-1; // не присутствует в матрице
		                  slb[PAM][inumber].iW=border_neighbor[inumber].iB;
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
					 if (w[border_neighbor[inumber].MCB - ls].bopening) {
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
							 switch (border_neighbor[inumber].Norm) {
							 case E_SIDE:  if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB] > 0.0) {
								        // жидкость втекает внутрь расчётной области через выходную границу.
								        rsign = -1.0;
										// Повторно не увеличиваем счетчик КО с рециркуляцией.
								        // breversedflow = true;
							         }
									 else rsign = 1.0;
									 break;
							 case W_SIDE: if (potent[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB] < 0.0) {
								          // жидкость втекает внутрь расчётной области через выходную границу.
								          rsign = -1.0;
										  // Повторно не увеличиваем счетчик КО с рециркуляцией.
								          //breversedflow = true;
							         }
									 else rsign = 1.0;
									 break;
							 case N_SIDE: if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB] > 0.0) {
								             // жидкость втекает внутрь расчётной области через выходную границу.
								             rsign = -1.0;
											 // Повторно не увеличиваем счетчик КО с рециркуляцией.
								             //breversedflow = true;
							             }
										 else rsign = 1.0;
										 break;
							 case S_SIDE:if (potent[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB] < 0.0) {
								            // жидкость втекает внутрь расчётной области через выходную границу.
								            rsign = -1.0;
											// Повторно не увеличиваем счетчик КО с рециркуляцией.
								            //breversedflow = true;
							            }
										else rsign = 1.0;
										break;
							 case T_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB] > 0.0) {
								             // жидкость втекает внутрь расчётной области через выходную границу.
								             rsign = -1.0;
											 // Повторно не увеличиваем счетчик КО с рециркуляцией.				             
											 //breversedflow = true;
							             }
										 else rsign = 1.0;
										 break;
							 case B_SIDE: if (potent[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB] < 0.0) {
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
								 slb[PAM][inumber].iI = border_neighbor[inumber].iI;
								 slb[PAM][inumber].iW = border_neighbor[inumber].iB;
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

					if (b_on_adaptive_local_refinement_mesh) {
						// Т.к. отсутствует информация о дальней связи iII.
						slb[PAM][inumber].ai = 1.0;
						slb[PAM][inumber].iI = border_neighbor[inumber].iI;
						slb[PAM][inumber].aw = 1.0;
						slb[PAM][inumber].iW = border_neighbor[inumber].iB;
						slb[PAM][inumber].b = 0.0;
					}
					else {

						// внутренняя нормаль
						switch (border_neighbor[inumber].Norm) {
						case E_SIDE:

							dl = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
							dS = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[1][border_neighbor[inumber].iI] - 1].y;
							dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

							slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
							slb[PAM][inumber].iI = border_neighbor[inumber].iI;
							slb[PAM][inumber].aw = slb[PAM][inumber].ai;
							slb[PAM][inumber].iW = border_neighbor[inumber].iB;
							if (!b_prosto) {


								deltal = 0.5*(pa[nvtx[1][border_neighbor[inumber].iII] - 1].x + pa[nvtx[0][border_neighbor[inumber].iII] - 1].x);
								deltal -= 0.5*(pa[nvtx[1][border_neighbor[inumber].iI] - 1].x + pa[nvtx[0][border_neighbor[inumber].iI] - 1].x);
								fiplus = 0.5*dl / deltal;

								taui = (tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);
								taui = taui / ((1.0 - fiplus)*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]); // проверено !


								// правая часть:  
								slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
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

						case N_SIDE:

							dl = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
							dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
							dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

							slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
							slb[PAM][inumber].iI = border_neighbor[inumber].iI;
							slb[PAM][inumber].aw = slb[PAM][inumber].ai;
							slb[PAM][inumber].iW = border_neighbor[inumber].iB;
							if (!b_prosto) {

								deltal = 0.5*(pa[nvtx[2][border_neighbor[inumber].iII] - 1].y + pa[nvtx[0][border_neighbor[inumber].iII] - 1].y);
								deltal -= 0.5*(pa[nvtx[2][border_neighbor[inumber].iI] - 1].y + pa[nvtx[0][border_neighbor[inumber].iI] - 1].y);
								fiplus = 0.5*dl / deltal;

								taui = (tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);
								taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]); // проверено !


								// правая часть:  
								slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
							}
							else slb[PAM][inumber].b = 0.0;

							break;

						case T_SIDE:

							dl = pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z;
							dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
							dS *= (pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y); // площадь грани

							slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
							slb[PAM][inumber].iI = border_neighbor[inumber].iI;
							slb[PAM][inumber].aw = slb[PAM][inumber].ai;
							slb[PAM][inumber].iW = border_neighbor[inumber].iB;
							if (!b_prosto) {

								deltal = 0.5*(pa[nvtx[4][border_neighbor[inumber].iII] - 1].z + pa[nvtx[0][border_neighbor[inumber].iII] - 1].z);
								deltal -= 0.5*(pa[nvtx[4][border_neighbor[inumber].iI] - 1].z + pa[nvtx[0][border_neighbor[inumber].iI] - 1].z);
								fiplus = 0.5*dl / deltal;

								taui = (tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);
								taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);// проверено !


								// правая часть:  
								slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
							}
							else slb[PAM][inumber].b = 0.0;

							break;

						case W_SIDE:

							dl = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
							dS = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
							dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

							slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
							slb[PAM][inumber].iI = border_neighbor[inumber].iI;
							slb[PAM][inumber].aw = slb[PAM][inumber].ai;
							slb[PAM][inumber].iW = border_neighbor[inumber].iB;

							if (!b_prosto) {

								deltal = -0.5*(pa[nvtx[1][border_neighbor[inumber].iII] - 1].x + pa[nvtx[0][border_neighbor[inumber].iII] - 1].x);
								deltal += 0.5*(pa[nvtx[1][border_neighbor[inumber].iI] - 1].x + pa[nvtx[0][border_neighbor[inumber].iI] - 1].x);
								fiplus = 0.5*dl / deltal;

								taui = (tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);
								taui = taui / ((1.0 - fiplus)*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);// проверено !


								// правая часть:  
								slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
							}
							else slb[PAM][inumber].b = 0.0;

							break;

						case S_SIDE:

							dl = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
							dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
							dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

							slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
							slb[PAM][inumber].iI = border_neighbor[inumber].iI;
							slb[PAM][inumber].aw = slb[PAM][inumber].ai;
							slb[PAM][inumber].iW = border_neighbor[inumber].iB;
							if (!b_prosto) {

								deltal = -0.5*(pa[nvtx[2][border_neighbor[inumber].iII] - 1].y + pa[nvtx[0][border_neighbor[inumber].iII] - 1].y);
								deltal += 0.5*(pa[nvtx[2][border_neighbor[inumber].iI] - 1].y + pa[nvtx[0][border_neighbor[inumber].iI] - 1].y);
								fiplus = 0.5*dl / deltal;

								taui = (tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);
								taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);// проверено !


								// правая часть:  
								slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
							}
							else slb[PAM][inumber].b = 0.0;

							break;

						case B_SIDE:

							dl = pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z;
							dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
							dS *= (pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y); // площадь грани

							slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
							slb[PAM][inumber].iI = border_neighbor[inumber].iI;
							slb[PAM][inumber].aw = slb[PAM][inumber].ai;
							slb[PAM][inumber].iW = border_neighbor[inumber].iB;

							if (!b_prosto) {

								deltal = -0.5*(pa[nvtx[4][border_neighbor[inumber].iII] - 1].z + pa[nvtx[0][border_neighbor[inumber].iII] - 1].z);
								deltal += 0.5*(pa[nvtx[4][border_neighbor[inumber].iI] - 1].z + pa[nvtx[0][border_neighbor[inumber].iI] - 1].z);
								fiplus = 0.5*dl / deltal;

								taui = (tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);
								taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);// проверено !


								// правая часть:  
								slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;

							}
							else slb[PAM][inumber].b = 0.0;

							break;
						} // switch
					}

			        //*/
                    integer j,l,xitem,k;
			        // сортировка по возрастанию
			        for (j=0; j<5; j++) {
				        k=j; xitem=border_neighbor[inumber].iW[j];
				        for (l=j+1; l<6; l++) {
					        if (border_neighbor[inumber].iW[l] < xitem) {
						       k=l; xitem=border_neighbor[inumber].iW[k];
					        }
				        }
                        border_neighbor[inumber].iW[k]=border_neighbor[inumber].iW[j];
				        border_neighbor[inumber].iW[j]=xitem;
			        }

                    j=0; l=0;
			        while (border_neighbor[inumber].iW[j]==(-1)) j++;

			        if (j<6) { slb[PAM][inumber].iW1=border_neighbor[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW2=border_neighbor[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW3=border_neighbor[inumber].iW[j++]; l++; }
			        if (j<6) { slb[PAM][inumber].iW4=border_neighbor[inumber].iW[j++]; l++; } 

			        switch (l) {
				       case 0: slb[PAM][inumber].iW1=-1;
		                        slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
		                        break;
				       case 1: slb[PAM][inumber].iW2=-1;
		                        slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 2: slb[PAM][inumber].iW3=-1;
		                        slb[PAM][inumber].iW4=-1;
						        break;
				       case 3: slb[PAM][inumber].iW4=-1;
						        break;
			        }

			
				}
			}

	    }
	     else if ((!bDirichlet) /*&& !((inumber==(maxbound-1)) && bPfix && (!bFReeStyle))*/ && !((border_neighbor[inumber].MCB>=ls) && (border_neighbor[inumber].MCB<(ls+lw)) && w[border_neighbor[inumber].MCB-ls].bpressure) )
	    {  
			// Не условие Дирихле, Не фиксация давления в точке, Не выходная граница
			//printf("neiman!\n"); getchar(); // debug

		

			/*
            // поправка давления фиксированно равна нулю:
		    slb[PAM][inumber].aw=1.0;
		    slb[PAM][inumber].ai=0.0;
		    slb[PAM][inumber].b=0.0;
			slb[PAM][inumber].iI=-1;//border_neighbor[inumber].iI; //  присутствует в матрице
		    slb[PAM][inumber].iW=border_neighbor[inumber].iB;
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

			if (b_on_adaptive_local_refinement_mesh) {
				slb[PAM][inumber].ai = 1.0;
				slb[PAM][inumber].iI = border_neighbor[inumber].iI;
				slb[PAM][inumber].aw = 1.0;
				slb[PAM][inumber].iW = border_neighbor[inumber].iB;
				slb[PAM][inumber].b = 0.0;
			}
			else {

				// внутренняя нормаль
				switch (border_neighbor[inumber].Norm) {
				case E_SIDE:

					dl = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
					dS = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[1][border_neighbor[inumber].iI] - 1].y;
					dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;

					if (!b_prosto) {

						deltal = 0.5*(pa[nvtx[1][border_neighbor[inumber].iII] - 1].x + pa[nvtx[0][border_neighbor[inumber].iII] - 1].x);
						deltal -= 0.5*(pa[nvtx[1][border_neighbor[inumber].iI] - 1].x + pa[nvtx[0][border_neighbor[inumber].iI] - 1].x);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;

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

				case N_SIDE:

					dl = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
					dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
					dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;

					if (!b_prosto) {


						deltal = 0.5*(pa[nvtx[2][border_neighbor[inumber].iII] - 1].y + pa[nvtx[0][border_neighbor[inumber].iII] - 1].y);
						deltal -= 0.5*(pa[nvtx[2][border_neighbor[inumber].iI] - 1].y + pa[nvtx[0][border_neighbor[inumber].iI] - 1].y);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;

				case T_SIDE:

					dl = pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z;
					dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
					dS *= (pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;

					if (!b_prosto) {


						deltal = 0.5*(pa[nvtx[4][border_neighbor[inumber].iII] - 1].z + pa[nvtx[0][border_neighbor[inumber].iII] - 1].z);
						deltal -= 0.5*(pa[nvtx[4][border_neighbor[inumber].iI] - 1].z + pa[nvtx[0][border_neighbor[inumber].iI] - 1].z);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);// проверено !

						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;

				case W_SIDE:

					dl = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
					dS = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
					dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;

					if (!b_prosto) {

						deltal = -0.5*(pa[nvtx[1][border_neighbor[inumber].iII] - 1].x + pa[nvtx[0][border_neighbor[inumber].iII] - 1].x);
						deltal += 0.5*(pa[nvtx[1][border_neighbor[inumber].iI] - 1].x + pa[nvtx[0][border_neighbor[inumber].iI] - 1].x);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_X_COMPONENT][border_neighbor[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;

				case S_SIDE:

					dl = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
					dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
					dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;

					if (!b_prosto) {


						deltal = -0.5*(pa[nvtx[2][border_neighbor[inumber].iII] - 1].y + pa[nvtx[0][border_neighbor[inumber].iII] - 1].y);
						deltal += 0.5*(pa[nvtx[2][border_neighbor[inumber].iI] - 1].y + pa[nvtx[0][border_neighbor[inumber].iI] - 1].y);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Y_COMPONENT][border_neighbor[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;

				case B_SIDE:

					dl = pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z;
					dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
					dS *= (pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y); // площадь грани

					slb[PAM][inumber].ai = 2.0*dbeta*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iB] * dS / dl;
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = slb[PAM][inumber].ai;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;
					if (!b_prosto) {


						deltal = -0.5*(pa[nvtx[4][border_neighbor[inumber].iII] - 1].z + pa[nvtx[0][border_neighbor[inumber].iII] - 1].z);
						deltal += 0.5*(pa[nvtx[4][border_neighbor[inumber].iI] - 1].z + pa[nvtx[0][border_neighbor[inumber].iI] - 1].z);
						fiplus = 0.5*dl / deltal;

						taui = (tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);
						taui = taui / ((1.0 - fiplus)*tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iI] + fiplus * tau[VELOCITY_Z_COMPONENT][border_neighbor[inumber].iII]);// проверено !


						// правая часть:  
						slb[PAM][inumber].b = (dbeta - 1.0)*taui*dS*(potent[PAM][border_neighbor[inumber].iI] - potent[PAM][border_neighbor[inumber].iII]) / deltal;
					}
					else slb[PAM][inumber].b = 0.0;

					break;
				} // switch
			}

			//*/
			integer j, l, xitem, k;
			// сортировка по возрастанию
			for (j = 0; j < 5; j++) {
				k = j; xitem = border_neighbor[inumber].iW[j];
				for (l = j + 1; l < 6; l++) {
					if (border_neighbor[inumber].iW[l] < xitem) {
						k = l; xitem = border_neighbor[inumber].iW[k];
					}
				}
				border_neighbor[inumber].iW[k] = border_neighbor[inumber].iW[j];
				border_neighbor[inumber].iW[j] = xitem;
			}

			j = 0; l = 0;
			while (border_neighbor[inumber].iW[j] == (-1)) j++;

			if (j < 6) { slb[PAM][inumber].iW1 = border_neighbor[inumber].iW[j++]; l++; }
			if (j < 6) { slb[PAM][inumber].iW2 = border_neighbor[inumber].iW[j++]; l++; }
			if (j < 6) { slb[PAM][inumber].iW3 = border_neighbor[inumber].iW[j++]; l++; }
			if (j < 6) { slb[PAM][inumber].iW4 = border_neighbor[inumber].iW[j++]; l++; }

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
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
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
				if (((border_neighbor[inumber].MCB >= ls) && (border_neighbor[inumber].MCB < (ls + lw)) && w[border_neighbor[inumber].MCB - ls].bopening)) {
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
					slb[PAM][inumber].iI = border_neighbor[inumber].iI;
					slb[PAM][inumber].aw = 1.0;
					slb[PAM][inumber].iW = border_neighbor[inumber].iB;
					slb[PAM][inumber].b = 0.0;
					*/
				}
			}
	}


		if (slb[PAM][inumber].ai != slb[PAM][inumber].ai) {
			printf("ai!=ai assemble bug. inumber=%lld ai=%e\n", inumber, slb[PAM][inumber].ai);
			printf("PAM \n");
			system("pause");
		}
		if (slb[PAM][inumber].aw != slb[PAM][inumber].aw) {
			printf("aw!=aw assemble bug. inumber=%lld aw=%e\n", inumber, slb[PAM][inumber].aw);
			printf("PAM \n");
			system("pause");
		}
		if (slb[PAM][inumber].b != slb[PAM][inumber].b) {
			printf("b!=b assemble bug. inumber=%lld b=%e\n", inumber, slb[PAM][inumber].b);
			printf("PAM \n");
			system("pause");
		}

		if (slb[PAM][inumber].ai<0.) {
			printf("Negative diffusion coefficients in PAM equation boundary condition ai. \n");
			system("PAUSE");
		}
		if (slb[PAM][inumber].aw<0.) {
			printf("Negative diffusion coefficients in PAM equation boundary condition aw. \n");
			system("PAUSE");
		}
		if (slb[PAM][inumber].b<0.) {
		//	printf("Negative rthdsd in PAM equation boundary condition b. \n");
			//system("PAUSE");
		}

} // my_elmatr_quad_PAm_bon3


  // возвращает скорректированный массовый поток.
  // скорректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
  // скорректированной скорости или просто скорости на основе простейшей интерполяции.
// Это нужно для отдельного решения уравнения конвекции-диффузии где задана пользовательская скорость,
// никакой поправки Рхи-Чоу просто интерполяция.
// 26.03.2017
void return_calc_correct_mass_flux_only_interpolation(integer iP, doublereal** potent,
	TOCHKA* pa, float** prop, float** prop_b,
	int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, 
	doublereal* &mfcurrentretune, BOUND* &border_neighbor, integer &ls, integer &lw,
	integer* ilevel_alice, int* ptr)
{

	// SpeedCorOld - скорректированная скорость на предыдущей итерации.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	//bool bRhieChowi = true, bRhieChowb = false; // i - internal, b - border.

	// iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = neighbors_for_the_internal_node[B_SIDE][0][iP];


	// 26.09.2016 Добавок для АЛИС сетки.
	integer iE2=-1, iN2=-1, iT2=-1, iW2=-1, iS2=-1, iB2=-1; // номера соседних контрольных объёмов
	integer iE3=-1, iN3=-1, iT3=-1, iW3=-1, iS3=-1, iB3=-1; // номера соседних контрольных объёмов
	integer iE4=-1, iN4=-1, iT4=-1, iW4=-1, iS4=-1, iB4=-1; // номера соседних контрольных объёмов

	// -1 если не используется и [0..maxelm+maxbound-1] если используется.
	if (b_on_adaptive_local_refinement_mesh) {

		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}

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

	if (b_on_adaptive_local_refinement_mesh) {

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

	}

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


	doublereal feplus2=0.5, fwplus2=0.5, fnplus2=0.5, fsplus2=0.5, ftplus2=0.5, fbplus2=0.5;
	if (b_on_adaptive_local_refinement_mesh) {
		// x-direction
		feplus2 = 0.5 * dx / dxe2;
		fwplus2 = 0.5 * dx / dxw2;
		// y-direction
		fnplus2 = 0.5 * dy / dyn2;
		fsplus2 = 0.5 * dy / dys2;
		// z-direction
		ftplus2 = 0.5 * dz / dzt2;
		fbplus2 = 0.5 * dz / dzb2;
	}
	doublereal feplus3=0.5, fwplus3=0.5, fnplus3=0.5, fsplus3=0.5, ftplus3=0.5, fbplus3=0.5;

	if (b_on_adaptive_local_refinement_mesh) {
		// x-direction
		feplus3 = 0.5 * dx / dxe3;
		fwplus3 = 0.5 * dx / dxw3;
		// y-direction
		fnplus3 = 0.5 * dy / dyn3;
		fsplus3 = 0.5 * dy / dys3;
		// z-direction
		ftplus3 = 0.5 * dz / dzt3;
		fbplus3 = 0.5 * dz / dzb3;
	}

	doublereal feplus4=0.5, fwplus4=0.5, fnplus4=0.5, fsplus4=0.5, ftplus4=0.5, fbplus4=0.5;

	if (b_on_adaptive_local_refinement_mesh) {
		// x-direction
		feplus4 = 0.5 * dx / dxe4;
		fwplus4 = 0.5 * dx / dxw4;
		// y-direction
		fnplus4 = 0.5 * dy / dyn4;
		fsplus4 = 0.5 * dy / dys4;
		// z-direction
		ftplus4 = 0.5 * dz / dzt4;
		fbplus4 = 0.5 * dz / dzb4;
	}

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

	if (b_on_adaptive_local_refinement_mesh) {

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
	}

	doublereal  rE3 = 0.0, rN3 = 0.0, rT3 = 0.0, rW3 = 0.0, rS3 = 0.0, rB3 = 0.0;

	if (b_on_adaptive_local_refinement_mesh) {

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
	}

	doublereal  rE4 = 0.0, rN4 = 0.0, rT4 = 0.0, rW4 = 0.0, rS4 = 0.0, rB4 = 0.0;

	if (b_on_adaptive_local_refinement_mesh) {

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

	if (b_on_adaptive_local_refinement_mesh) {

		// интерполяция плотности сделана так, чтобы выполнялись 
		// предельные соотношения.
		if (iE2 > -1) {
			if (!bE2) rhoe2 = rE2 * rP / (feplus2 * rE2 + (1.0 - feplus2) * rP); else rhoe2 = rE2;
		}
		if (iW2 > -1) {
			if (!bW2) rhow2 = rW2 * rP / (fwplus2 * rW2 + (1.0 - fwplus2) * rP); else rhow2 = rW2;
		}
		if (iN2 > -1) {
			if (!bN2) rhon2 = rN2 * rP / (fnplus2 * rN2 + (1.0 - fnplus2) * rP); else rhon2 = rN2;
		}
		if (iS2 > -1) {
			if (!bS2) rhos2 = rS2 * rP / (fsplus2 * rS2 + (1.0 - fsplus2) * rP); else rhos2 = rS2;
		}
		if (iT2 > -1) {
			if (!bT2) rhot2 = rT2 * rP / (ftplus2 * rT2 + (1.0 - ftplus2) * rP); else rhot2 = rT2;
		}
		if (iB2 > -1) {
			if (!bB2) rhob2 = rB2 * rP / (fbplus2 * rB2 + (1.0 - fbplus2) * rP); else rhob2 = rB2;
		}

		// интерполяция плотности сделана так, чтобы выполнялись 
		// предельные соотношения.
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
				SpeedCorOlde = feplus * potent[VELOCITY_X_COMPONENT][iE] + (1.0 - feplus)*potent[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOlde = potent[VELOCITY_X_COMPONENT][iE];
			}
		}
		if (iN > -1) {
			if (!bN) {
				SpeedCorOldn = fnplus * potent[VELOCITY_Y_COMPONENT][iN] + (1.0 - fnplus)*potent[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOldn = potent[VELOCITY_Y_COMPONENT][iN];
			}
		}
		if (iT > -1) {
			if (!bT) {
				SpeedCorOldt = ftplus * potent[VELOCITY_Z_COMPONENT][iT] + (1.0 - ftplus)*potent[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldt = potent[VELOCITY_Z_COMPONENT][iT];
			}
		}
		if (iW > -1) {
			if (!bW) {
				SpeedCorOldw = fwplus * potent[VELOCITY_X_COMPONENT][iW] + (1.0 - fwplus)*potent[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOldw = potent[VELOCITY_X_COMPONENT][iW];
			}
		}
		if (iS > -1) {
			if (!bS) {
				SpeedCorOlds = fsplus * potent[VELOCITY_Y_COMPONENT][iS] + (1.0 - fsplus)*potent[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOlds = potent[VELOCITY_Y_COMPONENT][iS];
			}
		}
		if (iB > -1) {
			if (!bB) {
				SpeedCorOldb = fbplus * potent[VELOCITY_Z_COMPONENT][iB] + (1.0 - fbplus)*potent[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldb = potent[VELOCITY_Z_COMPONENT][iB];
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

		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				if (!bE2) {
					SpeedCorOlde2 = feplus2 * potent[VELOCITY_X_COMPONENT][iE2] + (1.0 - feplus2) * potent[VELOCITY_X_COMPONENT][iP];
				}
				else {
					SpeedCorOlde2 = potent[VELOCITY_X_COMPONENT][iE2];
				}
			}
			if (iN2 > -1) {
				if (!bN2) {
					SpeedCorOldn2 = fnplus2 * potent[VELOCITY_Y_COMPONENT][iN2] + (1.0 - fnplus2) * potent[VELOCITY_Y_COMPONENT][iP];
				}
				else {
					SpeedCorOldn2 = potent[VELOCITY_Y_COMPONENT][iN2];
				}
			}
			if (iT2 > -1) {
				if (!bT2) {
					SpeedCorOldt2 = ftplus2 * potent[VELOCITY_Z_COMPONENT][iT2] + (1.0 - ftplus2) * potent[VELOCITY_Z_COMPONENT][iP];
				}
				else {
					SpeedCorOldt2 = potent[VELOCITY_Z_COMPONENT][iT2];
				}
			}
			if (iW2 > -1) {
				if (!bW2) {
					SpeedCorOldw2 = fwplus2 * potent[VELOCITY_X_COMPONENT][iW2] + (1.0 - fwplus2) * potent[VELOCITY_X_COMPONENT][iP];
				}
				else {
					SpeedCorOldw2 = potent[VELOCITY_X_COMPONENT][iW2];
				}
			}
			if (iS2 > -1) {
				if (!bS2) {
					SpeedCorOlds2 = fsplus2 * potent[VELOCITY_Y_COMPONENT][iS2] + (1.0 - fsplus2) * potent[VELOCITY_Y_COMPONENT][iP];
				}
				else {
					SpeedCorOlds2 = potent[VELOCITY_Y_COMPONENT][iS2];
				}
			}
			if (iB2 > -1) {
				if (!bB2) {
					SpeedCorOldb2 = fbplus2 * potent[VELOCITY_Z_COMPONENT][iB2] + (1.0 - fbplus2) * potent[VELOCITY_Z_COMPONENT][iP];
				}
				else {
					SpeedCorOldb2 = potent[VELOCITY_Z_COMPONENT][iB2];
				}
			}

			if (iE3 > -1) {
				if (!bE3) {
					SpeedCorOlde3 = feplus3 * potent[VELOCITY_X_COMPONENT][iE3] + (1.0 - feplus3) * potent[VELOCITY_X_COMPONENT][iP];
				}
				else {
					SpeedCorOlde3 = potent[VELOCITY_X_COMPONENT][iE3];
				}
			}
			if (iN3 > -1) {
				if (!bN3) {
					SpeedCorOldn3 = fnplus3 * potent[VELOCITY_Y_COMPONENT][iN3] + (1.0 - fnplus3) * potent[VELOCITY_Y_COMPONENT][iP];
				}
				else {
					SpeedCorOldn3 = potent[VELOCITY_Y_COMPONENT][iN3];
				}
			}
			if (iT3 > -1) {
				if (!bT3) {
					SpeedCorOldt3 = ftplus3 * potent[VELOCITY_Z_COMPONENT][iT3] + (1.0 - ftplus3) * potent[VELOCITY_Z_COMPONENT][iP];
				}
				else {
					SpeedCorOldt3 = potent[VELOCITY_Z_COMPONENT][iT3];
				}
			}
			if (iW3 > -1) {
				if (!bW3) {
					SpeedCorOldw3 = fwplus3 * potent[VELOCITY_X_COMPONENT][iW3] + (1.0 - fwplus3) * potent[VELOCITY_X_COMPONENT][iP];
				}
				else {
					SpeedCorOldw3 = potent[VELOCITY_X_COMPONENT][iW3];
				}
			}
			if (iS3 > -1) {
				if (!bS3) {
					SpeedCorOlds3 = fsplus3 * potent[VELOCITY_Y_COMPONENT][iS3] + (1.0 - fsplus3) * potent[VELOCITY_Y_COMPONENT][iP];
				}
				else {
					SpeedCorOlds3 = potent[VELOCITY_Y_COMPONENT][iS3];
				}
			}
			if (iB3 > -1) {
				if (!bB3) {
					SpeedCorOldb3 = fbplus3 * potent[VELOCITY_Z_COMPONENT][iB3] + (1.0 - fbplus3) * potent[VELOCITY_Z_COMPONENT][iP];
				}
				else {
					SpeedCorOldb3 = potent[VELOCITY_Z_COMPONENT][iB3];
				}
			}

			if (iE4 > -1) {
				if (!bE4) {
					SpeedCorOlde4 = feplus4 * potent[VELOCITY_X_COMPONENT][iE4] + (1.0 - feplus4) * potent[VELOCITY_X_COMPONENT][iP];
				}
				else {
					SpeedCorOlde4 = potent[VELOCITY_X_COMPONENT][iE4];
				}
			}
			if (iN4 > -1) {
				if (!bN4) {
					SpeedCorOldn4 = fnplus4 * potent[VELOCITY_Y_COMPONENT][iN4] + (1.0 - fnplus4) * potent[VELOCITY_Y_COMPONENT][iP];
				}
				else {
					SpeedCorOldn4 = potent[VELOCITY_Y_COMPONENT][iN4];
				}
			}
			if (iT4 > -1) {
				if (!bT4) {
					SpeedCorOldt4 = ftplus4 * potent[VELOCITY_Z_COMPONENT][iT4] + (1.0 - ftplus4) * potent[VELOCITY_Z_COMPONENT][iP];
				}
				else {
					SpeedCorOldt4 = potent[VELOCITY_Z_COMPONENT][iT4];
				}
			}
			if (iW4 > -1) {
				if (!bW4) {
					SpeedCorOldw4 = fwplus4 * potent[VELOCITY_X_COMPONENT][iW4] + (1.0 - fwplus4) * potent[VELOCITY_X_COMPONENT][iP];
				}
				else {
					SpeedCorOldw4 = potent[VELOCITY_X_COMPONENT][iW4];
				}
			}
			if (iS4 > -1) {
				if (!bS4) {
					SpeedCorOlds4 = fsplus4 * potent[VELOCITY_Y_COMPONENT][iS4] + (1.0 - fsplus4) * potent[VELOCITY_Y_COMPONENT][iP];
				}
				else {
					SpeedCorOlds4 = potent[VELOCITY_Y_COMPONENT][iS4];
				}
			}
			if (iB4 > -1) {
				if (!bB4) {
					SpeedCorOldb4 = fbplus4 * potent[VELOCITY_Z_COMPONENT][iB4] + (1.0 - fbplus4) * potent[VELOCITY_Z_COMPONENT][iP];
				}
				else {
					SpeedCorOldb4 = potent[VELOCITY_Z_COMPONENT][iB4];
				}
			}
		}
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

		}

		doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
		
		if (b_on_adaptive_local_refinement_mesh) {

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

		}

		doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.
		
		if (b_on_adaptive_local_refinement_mesh) {

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

		}

		// mfold - значение массового потока с предыдущей итерации.
		//Fe =  (rhoe*SpeedCorOlde+ rhoe2*SpeedCorOlde2 + rhoe3*SpeedCorOlde3 + rhoe4*SpeedCorOlde4)*dy*dz;
		//Fn = (rhon*SpeedCorOldn + rhon2*SpeedCorOldn2 + rhon3*SpeedCorOldn3 + rhon4*SpeedCorOldn4)*dx*dz;
		//Ft = (rhot*SpeedCorOldt + rhot2*SpeedCorOldt2 + rhot3*SpeedCorOldt3 + rhot4*SpeedCorOldt4)*dx*dy;
		//Fw = (rhow*SpeedCorOldw + rhow2*SpeedCorOldw2 + rhow3*SpeedCorOldw3 + rhow4*SpeedCorOldw4)*dy*dz;
		//Fs = (rhos*SpeedCorOlds + rhos2*SpeedCorOlds2 + rhos3*SpeedCorOlds3 + rhos4*SpeedCorOlds4)*dx*dz;
		//Fb = (rhob*SpeedCorOldb + rhob2*SpeedCorOldb2 + rhob3*SpeedCorOldb3 + rhob4*SpeedCorOldb4)*dx*dy;

		Fe = (rhoe * SpeedCorOlde * dSqe + rhoe2 * SpeedCorOlde2 * dSqe2 + rhoe3 * SpeedCorOlde3 * dSqe3 + rhoe4 * SpeedCorOlde4 * dSqe4);
		Fn = (rhon * SpeedCorOldn * dSqn + rhon2 * SpeedCorOldn2 * dSqn2 + rhon3 * SpeedCorOldn3 * dSqn3 + rhon4 * SpeedCorOldn4 * dSqn4);
		Ft = (rhot * SpeedCorOldt * dSqt + rhot2 * SpeedCorOldt2 * dSqt2 + rhot3 * SpeedCorOldt3 * dSqt3 + rhot4 * SpeedCorOldt4 * dSqt4);
		Fw = (rhow * SpeedCorOldw * dSqw + rhow2 * SpeedCorOldw2 * dSqw2 + rhow3 * SpeedCorOldw3 * dSqw3 + rhow4 * SpeedCorOldw4 * dSqw4);
		Fs = (rhos * SpeedCorOlds * dSqs + rhos2 * SpeedCorOlds2 * dSqs2 + rhos3 * SpeedCorOlds3 * dSqs3 + rhos4 * SpeedCorOlds4 * dSqs4);
		Fb = (rhob * SpeedCorOldb * dSqb + rhob2 * SpeedCorOldb2 * dSqb2 + rhob3 * SpeedCorOldb3 * dSqb3 + rhob4 * SpeedCorOldb4 * dSqb4);
	
		//printf("dx=%e dy=%e rhot=%e rhob=%e SpeedCorOldt=%e SpeedCorOldb=%e Ft=%e Fb=%e \n",dx,dy,rhot,rhob, SpeedCorOldt, SpeedCorOldb, Ft, Fb);
		//printf("dx=%e dy=%e rhot2=%e rhob2=%e SpeedCorOldt2=%e SpeedCorOldb2=%e Ft=%e Fb=%e \n", dx, dy, rhot2, rhob2, SpeedCorOldt2, SpeedCorOldb2, Ft, Fb);
		//printf("dx=%e dy=%e rhot3=%e rhob3=%e SpeedCorOldt3=%e SpeedCorOldb3=%e Ft=%e Fb=%e \n", dx, dy, rhot3, rhob3, SpeedCorOldt3, SpeedCorOldb3, Ft, Fb);
		//printf("dx=%e dy=%e rhot4=%e rhob4=%e SpeedCorOldt4=%e SpeedCorOldb4=%e Ft=%e Fb=%e \n", dx, dy, rhot4, rhob4, SpeedCorOldt4, SpeedCorOldb4, Ft, Fb);
		//getchar();

		//if (fabs(Ft) + fabs(Fb) > 0.0) {
			//printf("non zero Ft=%e and Fb=%e. Ok.\n",Ft,Fb);
			//getchar();
		//}
	
		mfcurrentretune[E_SIDE] = Fe;
		mfcurrentretune[N_SIDE] = Fn;
		mfcurrentretune[T_SIDE] = Ft;
		mfcurrentretune[W_SIDE] = Fw;
		mfcurrentretune[S_SIDE] = Fs;
		mfcurrentretune[B_SIDE] = Fb;

		// Для прохождения проверки корректности необходимо
		// принудительно занулить скорость и поток скорости
		// на твердой сенке (нормальную к стенке компоненту).
		// Твердая стенка идентифицируентся маркером MCB == (ls + lw).

		if (iE >= maxelm) {
			// граничный узел
			integer inumber = iE - maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				mfcurrentretune[E_SIDE] = 0.0;
			}
		}
				
		if (iW >= maxelm) {
			// граничный узел
			integer inumber = iW - maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				mfcurrentretune[W_SIDE] = 0.0;
			}
		}

		if (iN >= maxelm) {
			// граничный узел
			integer inumber = iN - maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				mfcurrentretune[N_SIDE] = 0.0;
			}
		}

		if (iS >= maxelm) {
			// граничный узел
			integer inumber = iS - maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				mfcurrentretune[S_SIDE] = 0.0;
			}
		}

		if (iT >= maxelm) {
			// граничный узел
			integer inumber = iT - maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				mfcurrentretune[T_SIDE] = 0.0;
			}
		}

		if (iB >= maxelm) {
			// граничный узел
			integer inumber = iB - maxelm;
			if (border_neighbor[inumber].MCB == (ls + lw)) {
				mfcurrentretune[B_SIDE] = 0.0;
			}
		}


} // return_correct_mass_flux_only_interpolation


// возвращает значение скорости на грани с использованием монотонизирующей поправки Рхи-Чоу.
// Внимание: данная реализация подходит не только для стационарных задач, но и для нестационарных.
// Для нестационарных задач требуется хранить массовый поток на грани с предыдущего временного слоя.
// Основная сложность данной функции заключается в её правильном вызове (требуется правильно указать все параметры
// при вызове) из-за большого количества передаваемых параметров.
// реализовано 31 марта 2012 года.
doublereal calcFg(bool bG, doublereal fgplus, doublereal VG, doublereal VP, 
	        integer iP, integer G, doublereal alpha, doublereal rhog, 
			doublereal dS, doublereal dV, bool btimedepend, 
			doublereal** speedoldtimestep, doublereal* mfoldtimestep,
			doublereal dG, doublereal dP, doublereal dtimestep,
			doublereal RCh, int** nvtx, int*** neighbors_for_the_internal_node,
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
	       ug+=RCh*ugRhieChow_internal(iP, G, alpha, nvtx, neighbors_for_the_internal_node, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
	   }
	}
	else {
	   ug=VG;
	   if (bRhieChowb) {
	      ug+=RCh*ugRhieChow_internal_border(iP, G, alpha, nvtx, neighbors_for_the_internal_node, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
	   }
	}

	// конвективный поток через грань КО.
	Fg=rhog*ug*dS;

	// Внимание здесь предполагается что плотность строго постоянна.
	if (btimedepend) {
		// в нестационарном случае к поправке Рхи-Чоу требуется сделать добавку.

		doublereal VGoldtimestep=0.0;
		doublereal VPoldtimestep=0.0;
		integer iG=neighbors_for_the_internal_node[G][0][iP];
		switch (G) {
		 case E_SIDE: case W_SIDE: VPoldtimestep=speedoldtimestep[VELOCITY_X_COMPONENT][iP];
			               VGoldtimestep=speedoldtimestep[VELOCITY_X_COMPONENT][iG];
			      break;
		 case N_SIDE: case S_SIDE:VPoldtimestep=speedoldtimestep[VELOCITY_Y_COMPONENT][iP];
			               VGoldtimestep=speedoldtimestep[VELOCITY_Y_COMPONENT][iG];
			      break;
		 case T_SIDE:  case B_SIDE: VPoldtimestep=speedoldtimestep[VELOCITY_Z_COMPONENT][iP];
			                VGoldtimestep=speedoldtimestep[VELOCITY_Z_COMPONENT][iG];
			      break;
		 }

		// Скорость на грани КО с предыдущего временного слоя:
		doublereal ugold=fgplus*VGoldtimestep+(1.0-fgplus)*VPoldtimestep; 


		// диагональный коэффициент матрицы должен быть взят на границе контрольного объёма.
		if (!bG) {
			// внутренний контрольный объём
		    doublereal dc_e=dP*dG/(fgplus*dG+(1.0-fgplus)*dP); // диагональный коэффициент на грани
		    doublereal tau_e=(rhog*alpha*dV)/(dc_e); // псевдовремя (rhog - плотность на грани контрольного объёма.
			if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) tau_e/=(1.0-alpha);
		    Fg+=(tau_e/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
		}
		else {
			// граничный контрольный объём
			doublereal tau_e=(rhog*alpha*dV)/dG;
			if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) tau_e/=(1.0-alpha);
			Fg+=(tau_e/dtimestep)*(mfoldtimestep[G]-rhog*ugold*dS);
		}

	}
	}

	return Fg;

} // calcFg


// возвращает значение скорости на грани с использованием монотонизирующей поправки Рхи-Чоу.
// Внимание: данная реализация подходит не только для стационарных задач, но и для нестационарных.
// Для нестационарных задач требуется хранить массовый поток на грани с предыдущего временного слоя.
// Основная сложность данной функции заключается в её правильном вызове (требуется правильно указать все параметры
// при вызове) из-за большого количества передаваемых параметров.
// реализовано 21 июня 2012 года. (более простой вариант по сравнению с calcFg).
// Смотри Гаврилов Андрей опыт по ВГД (CFD).
doublereal calcFg2(bool bG, doublereal fgplus, 
	        integer iP, integer G, doublereal rhog, 
			doublereal dS, bool btimedepend, 
			doublereal** speedoldtimestep, doublereal* mfoldtimestep,
			doublereal dtimestep, doublereal RCh, int*** neighbors_for_the_internal_node,
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

	integer iG=neighbors_for_the_internal_node[G][0][iP];
	doublereal VG=0.0, VP=0.0;
	switch (G) {
		case E_SIDE: case W_SIDE: VG=potent[VELOCITY_X_COMPONENT][iG]; VP=potent[VELOCITY_X_COMPONENT][iP]; break;
		case N_SIDE: case S_SIDE:VG=potent[VELOCITY_Y_COMPONENT][iG]; VP=potent[VELOCITY_Y_COMPONENT][iP]; break;
		case T_SIDE: case B_SIDE: VG=potent[VELOCITY_Z_COMPONENT][iG]; VP=potent[VELOCITY_Z_COMPONENT][iP]; break;
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
	       case E_SIDE: gradP=potent[GRADXPRESS][iP]; gradG=potent[GRADXPRESS][neighbors_for_the_internal_node[E_SIDE][0][iP]]; break;
	       case W_SIDE: gradP=potent[GRADXPRESS][iP]; gradG=potent[GRADXPRESS][neighbors_for_the_internal_node[W_SIDE][0][iP]]; break;
	       case N_SIDE: gradP=potent[GRADYPRESS][iP]; gradG=potent[GRADYPRESS][neighbors_for_the_internal_node[N_SIDE][0][iP]]; break;
	       case S_SIDE:gradP=potent[GRADYPRESS][iP]; gradG=potent[GRADYPRESS][neighbors_for_the_internal_node[S_SIDE][0][iP]]; break;
	       case T_SIDE: gradP=potent[GRADZPRESS][iP]; gradG=potent[GRADZPRESS][neighbors_for_the_internal_node[T_SIDE][0][iP]]; break;
	       case B_SIDE: gradP=potent[GRADZPRESS][iP]; gradG=potent[GRADZPRESS][neighbors_for_the_internal_node[B_SIDE][0][iP]]; break;
	    }


	    // это делается лишь для того чтобы добавить дополнительную
	    // нижнюю релаксацию так как это написано у I. Sezai.
	    if (!bG) { 
	       ug=fgplus*VG+(1.0-fgplus)*VP;
	       if (bRhieChowi) {
	           //ug+=RCh*ugRhieChow_internal(iP, G, alpha, nvtx, neighbors_for_the_internal_node, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
		       ug+=RCh*(taug/rhog)*(-(potent[PRESS][iG]-potent[PRESS][iP])/drg+(1.0-fgplus)*gradP+fgplus*gradG); // Вклад поправки Рхи-Чоу
	       }
	    }
	    else {
	       ug=VG;
	       /*
	       // поправи для граничного узла нету. Т.е. считаем её нулевой.
	       if (bRhieChowb) {
	          set additional ug+=0.0;
	          //ug+=RCh*ugRhieChow_internal_border(iP, G, alpha, nvtx, neighbors_for_the_internal_node, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу
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
		       case E_SIDE: case W_SIDE: VPoldtimestep=speedoldtimestep[VELOCITY_X_COMPONENT][iP];
			                     VGoldtimestep=speedoldtimestep[VELOCITY_X_COMPONENT][iG];
			                     break;
		       case N_SIDE: case S_SIDE:VPoldtimestep=speedoldtimestep[VELOCITY_Y_COMPONENT][iP];
			                     VGoldtimestep=speedoldtimestep[VELOCITY_Y_COMPONENT][iG];
			                     break;
		       case T_SIDE:  case B_SIDE: VPoldtimestep=speedoldtimestep[VELOCITY_Z_COMPONENT][iP];
			                      VGoldtimestep=speedoldtimestep[VELOCITY_Z_COMPONENT][iG];
			                     break;
		    }

		    // Скорость на грани КО с предыдущего временного слоя:
		    doublereal ugold=fgplus*VGoldtimestep+(1.0-fgplus)*VPoldtimestep; 


		    // диагональный коэффициент матрицы должен быть взят на границе контрольного объёма.
			// одинаково работает и для внутреннего и для граничного КО.
			Fg += (taug / dtimestep)*(mfoldtimestep[G] - rhog*ugold*dS);
		    

    	}
	}

	return Fg;

} // calcFg2

// возвращает значение скорости на грани с использованием монотонизирующей поправки Рхи-Чоу.
// Внимание: данная реализация подходит не только для стационарных задач, но и для нестационарных.
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
			doublereal dtimestep, doublereal RCh, int*** neighbors_for_the_internal_node,
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

	integer iG=neighbors_for_the_internal_node[G][0][iP];
	switch (iNODEid) {
		case 1: iG = neighbors_for_the_internal_node[G][0][iP]; break;
		case 2: iG = neighbors_for_the_internal_node[G][1][iP]; break;
		case 3: iG = neighbors_for_the_internal_node[G][2][iP]; break;
		case 4: iG = neighbors_for_the_internal_node[G][3][iP]; break;
		default: iG = neighbors_for_the_internal_node[G][0][iP];  break;
	}

	doublereal Fg = 0.0; // массовый поток на грани контрольного объёма.

	if (iG > -1) {

		doublereal VG = 0.0, VP = 0.0;
		integer VGid = -1;
		switch (G) {
		case E_SIDE: case W_SIDE: VG = potent[VELOCITY_X_COMPONENT][iG]; VP = potent[VELOCITY_X_COMPONENT][iP]; VGid = VELOCITY_X_COMPONENT; break;
		case N_SIDE: case S_SIDE: VG = potent[VELOCITY_Y_COMPONENT][iG]; VP = potent[VELOCITY_Y_COMPONENT][iP]; VGid = VELOCITY_Y_COMPONENT; break;
		case T_SIDE: case B_SIDE: VG = potent[VELOCITY_Z_COMPONENT][iG]; VP = potent[VELOCITY_Z_COMPONENT][iP]; VGid = VELOCITY_Z_COMPONENT; break;
		} // end switch


		if (VG != VG) {
			printf("ERROR !!! function calcFg3: VG=%e VP=%e iG=%lld iP=%lld G=%lld Vid=%lld\n",VG, VP, iG, iP, G, VGid);
			// Возможно вы забыли задать выходную границу потока в своей модели.
			printf("ATTANTION!!! Perhaps you forgot to set the output boundary of the flow in your model...\n");
			system("pause");
		}

		doublereal ug = 0.0; // скорость на грани контрольного объёма.
		

		if (bsimplelinearinterpol) {

			// самая обычная линейная интерполяция скорости на грань,
			// используется на первой итерации стационарного солвера.
			if (!bG) {
				ug = fgplus * VG + (1.0 - fgplus)*VP;
			}
			else {
				ug = VG;
			}

			// конвективный поток через грань КО.
			Fg = rhog * ug*dS;
		}
		else
		{

			// формулы справедливы и для граничных узлов тоже.
			doublereal taug = tau[VGid][iP] * tau[VGid][iG] / (fgplus*tau[VGid][iG] + (1.0 - fgplus)*tau[VGid][iP]); // псевдовремя на грани КО.
			//28.07.2016
			if (bG) {
				taug = tau[VGid][iG];
			}
			doublereal gradP = 0.0, gradG = 0.0;
			switch (iNODEid) {
			case 1:
				switch (G) {
				case E_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[E_SIDE][0][iP]]; break;
				case W_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[W_SIDE][0][iP]]; break;
				case N_SIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[N_SIDE][0][iP]]; break;
				case S_SIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[S_SIDE][0][iP]]; break;
				case T_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[T_SIDE][0][iP]]; break;
				case B_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[B_SIDE][0][iP]]; break;
				}
				break;
			case 2:
				switch (G) {
				case E_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[E_SIDE][1][iP]]; break;
				case W_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[W_SIDE][1][iP]]; break;
				case N_SIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[N_SIDE][1][iP]]; break;
				case S_SIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[S_SIDE][1][iP]]; break;
				case T_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[T_SIDE][1][iP]]; break;
				case B_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[B_SIDE][1][iP]]; break;
				}
				break;
			case 3:
				switch (G) {
				case E_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[E_SIDE][2][iP]]; break;
				case W_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[W_SIDE][2][iP]]; break;
				case N_SIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[N_SIDE][2][iP]]; break;
				case S_SIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[S_SIDE][2][iP]]; break;
				case T_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[T_SIDE][2][iP]]; break;
				case B_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[B_SIDE][2][iP]]; break;
				}
				break;
			case 4:
				switch (G) {
				case E_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[E_SIDE][3][iP]]; break;
				case W_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[W_SIDE][3][iP]]; break;
				case N_SIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[N_SIDE][3][iP]]; break;
				case S_SIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[S_SIDE][3][iP]]; break;
				case T_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[T_SIDE][3][iP]]; break;
				case B_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[B_SIDE][3][iP]]; break;
				}
				break;
			default:
				switch (G) {
				case E_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[E_SIDE][0][iP]]; break;
				case W_SIDE: gradP = potent[GRADXPRESS][iP]; gradG = potent[GRADXPRESS][neighbors_for_the_internal_node[W_SIDE][0][iP]]; break;
				case N_SIDE: gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[N_SIDE][0][iP]]; break;
				case S_SIDE:gradP = potent[GRADYPRESS][iP]; gradG = potent[GRADYPRESS][neighbors_for_the_internal_node[S_SIDE][0][iP]]; break;
				case T_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[T_SIDE][0][iP]]; break;
				case B_SIDE: gradP = potent[GRADZPRESS][iP]; gradG = potent[GRADZPRESS][neighbors_for_the_internal_node[B_SIDE][0][iP]]; break;
				}
				break;
			}

			

			// это делается лишь для того чтобы добавить дополнительную
			// нижнюю релаксацию так как это написано у I. Sezai.
			if (!bG) {
				ug = fgplus * VG + (1.0 - fgplus)*VP;
				if (ug != ug) {
					printf("ug!=ug 3\n");
					printf("fgplus=%e VG=%e VP=%e\n",fgplus, VG, VP);
					system("pause");
				}

				if (bRhieChowi) {
					//ug+=RCh*ugRhieChow_internal(iP, G, alpha, nvtx, neighbors_for_the_internal_node, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу


					// Пример того как делать неправильно: 
					//ug+=RCh*(taug/rhog)*(-(potent[PRESS][iG]-potent[PRESS][iP])/drg+(1.0-fgplus)*gradP+fgplus*gradG);
					// Вышенаписанная строчка неправильная так как неправильно вычисляется градиент давления на грани potent[PRESS][iG]-potent[PRESS][iP].
					// Запись (potent[PRESS][iG]-potent[PRESS][iP]) верна только для направлений E, N , T. НО НЕВЕРНА ДЛЯ WSIDE, S, B.

					// Правильный вариант поправки Рхи-Чоу.
					switch (G) {// Вклад поправки Рхи-Чоу
					case E_SIDE: case N_SIDE: case T_SIDE:
						ug += RCh * (taug / rhog)*(-(potent[PRESS][iG] - potent[PRESS][iP]) / drg + (1.0 - fgplus)*gradP + fgplus * gradG);						
						break;
					case W_SIDE: case S_SIDE:case B_SIDE:
						ug += RCh * (taug / rhog)*(-(potent[PRESS][iP] - potent[PRESS][iG]) / drg + (1.0 - fgplus)*gradP + fgplus * gradG);
						break;
					}


					if (ug != ug) {
						printf("ug!=ug 4\n");
						system("pause");
					}
				}
			}
			else {
				ug = VG;
				if (ug != ug) {
					printf("ug!=ug 1\n");
					system("pause");
				}

				// поправки для граничного узла нету. Т.е. считаем её нулевой.
				// но с другой стороны при тестировании обнаруживаем нефизичные осцилляции вблизи выходной границы.
				if (bRhieChowb) {
					//set additional ug+=0.0;
					//ug+=RCh*ugRhieChow_internal_border(iP, G, alpha, nvtx, neighbors_for_the_internal_node, maxelm, pressure, pa, diag_coef); // Вклад поправки Рхи-Чоу

					 // в данном случае drg==0.5*dr; fgplus==0.5;

					// Правильный вариант поправки Рхи-Чоу.
					switch (G) {// Вклад поправки Рхи-Чоу
					case E_SIDE: case N_SIDE: case T_SIDE: 
						ug += RCh * (taug / rhog)*(-(potent[PRESS][iG] - potent[PRESS][iP]) / drg + (1.0 - fgplus)*gradP + fgplus * gradG);
						break;
					case W_SIDE: case S_SIDE:case B_SIDE:
						ug += RCh * (taug / rhog)*(-(potent[PRESS][iP] - potent[PRESS][iG]) / drg + (1.0 - fgplus)*gradP + fgplus * gradG);
						break;
					}

				}
				if (ug != ug) {
					printf("ug!=ug 2\n");
					system("pause");
				}
			}

			// конвективный поток через грань КО.
			Fg = rhog * ug * dS;
			if (Fg != Fg) {
				printf("Fg3!=Fg3 calcFg3\n");
				printf("rhog=%e ug=%e dS=%e\n",rhog,ug,dS);
				system("pause");
			}

			// Внимание здесь предполагается что плотность строго постоянна.
			if (btimedepend) {
				// в нестационарном случае к поправке Рхи-Чоу требуется сделать добавку.

				doublereal VGoldtimestep = 0.0;
				doublereal VPoldtimestep = 0.0;

				switch (G) {
				case E_SIDE: case W_SIDE: VPoldtimestep = speedoldtimestep[VELOCITY_X_COMPONENT][iP];
					VGoldtimestep = speedoldtimestep[VELOCITY_X_COMPONENT][iG];
					break;
				case N_SIDE: case S_SIDE:VPoldtimestep = speedoldtimestep[VELOCITY_Y_COMPONENT][iP];
					VGoldtimestep = speedoldtimestep[VELOCITY_Y_COMPONENT][iG];
					break;
				case T_SIDE:  case B_SIDE: VPoldtimestep = speedoldtimestep[VELOCITY_Z_COMPONENT][iP];
					VGoldtimestep = speedoldtimestep[VELOCITY_Z_COMPONENT][iG];
					break;
				}

				// Скорость на грани КО с предыдущего временного слоя:
				doublereal ugold = fgplus * VGoldtimestep + (1.0 - fgplus)*VPoldtimestep;
				if (Fg != Fg) {
					printf("Fg2!=Fg2 calcFg3\n");
					system("pause");
				}

				// диагональный коэффициент матрицы должен быть взят на границе контрольного объёма.
				if (!bG) {
					if (bRhieChowi) {
						// внутренний контрольный объём
						Fg += (taug / dtimestep)*(mfoldtimestep[G] - rhog * ugold*dS);
					}
				}
				else {
					if (bRhieChowb) {
						// граничный контрольный объём
						Fg += (taug / dtimestep)*(mfoldtimestep[G] - rhog * ugold*dS);
					}
				}

			}
		}
		if (Fg != Fg) {
			printf("Fg!=Fg calcFg3\n");
			system("pause");
		}

	}
	else {
		Fg = 0.0;
	}

	if (Fg != Fg) {
		// не число.
		return 0.0;
	}
	else {
		return Fg;
	}

} // calcFg3


// возвращает скорректированный массовый поток.
// скорректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
// скорректированной скорости и давления, а также сохранённых диагональных 
// коэффициентов матрицы СЛАУ для компонент скорости.
void return_calc_correct_mass_flux(integer iP, doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
	int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, doublereal **diag_coef,
						doublereal* alpha, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, bool bsimplelinearinterpol,
						doublereal** SpeedCorOld, doublereal *mfold,
	                    BOUND* border_neighbor, integer *ilevel_alice, int* ptr)
{

	// SpeedCorOld - скорректированная скорость на предыдущей итерации.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border.

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

	doublereal dxe3 = 0.5*dx, dxw3 = 0.5*dx, dyn3 = 0.5*dy, dys3 = 0.5*dy, dzt3 = 0.5*dz, dzb3 = 0.5*dz;
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

	doublereal dxe4 = 0.5*dx, dxw4 = 0.5*dx, dyn4 = 0.5*dy, dys4 = 0.5*dy, dzt4 = 0.5*dz, dzb4 = 0.5*dz;
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

	// плотность аппроксимируется средним гармоническим
	doublereal rhoe=0.0, rhow = 0.0, rhon = 0.0, rhos = 0.0, rhot = 0.0, rhob = 0.0;
	doublereal rP, rE, rN, rT, rW, rS, rB;

    rP=prop[RHO][iP];
	rE = rP; rN = rP; rT = rP;
	rW = rP; rS = rP; rB = rP;

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

	doublereal rhoe2 = 0.0, rhow2 = 0.0, rhon2 = 0.0, rhos2 = 0.0, rhot2 = 0.0, rhob2 = 0.0;
	doublereal rE2, rN2, rT2, rW2, rS2, rB2;

	rE2 = rP; rN2 = rP; rT2 = rP;
	rW2 = rP; rS2 = rP; rB2 = rP;

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

	doublereal rhoe3 = 0.0, rhow3 = 0.0, rhon3 = 0.0, rhos3 = 0.0, rhot3 = 0.0, rhob3 = 0.0;
	doublereal rE3, rN3, rT3, rW3, rS3, rB3;

	rE3 = rP; rN3 = rP; rT3 = rP;
	rW3 = rP; rS3 = rP; rB3 = rP;

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


	doublereal rhoe4 = 0.0, rhow4 = 0.0, rhon4 = 0.0, rhos4 = 0.0, rhot4 = 0.0, rhob4 = 0.0;
	doublereal rE4, rN4, rT4, rW4, rS4, rB4;

	rE4 = rP; rN4 = rP; rT4 = rP;
	rW4 = rP; rS4 = rP; rB4 = rP;

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
		if (!bE) rhoe = rE*rP / (feplus*rE + (1.0 - feplus)*rP); else rhoe = rE;
	}
	if (iW > -1) {
		if (!bW) rhow = rW*rP / (fwplus*rW + (1.0 - fwplus)*rP); else rhow = rW;
	}
	if (iN > -1) {
		if (!bN) rhon = rN*rP / (fnplus*rN + (1.0 - fnplus)*rP); else rhon = rN;
	}
	if (iS > -1) {
		if (!bS) rhos = rS*rP / (fsplus*rS + (1.0 - fsplus)*rP); else rhos = rS;
	}
	if (iT > -1) {
		if (!bT) rhot = rT*rP / (ftplus*rT + (1.0 - ftplus)*rP); else rhot = rT;
	}
	if (iB > -1) {
		if (!bB) rhob = rB*rP / (fbplus*rB + (1.0 - fbplus)*rP); else rhob = rB;
	}

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE2 > -1) {
		if (!bE2) rhoe2 = rE2*rP / (feplus2*rE2 + (1.0 - feplus2)*rP); else rhoe2 = rE2;
	}
	if (iW2 > -1) {
		if (!bW2) rhow2 = rW2*rP / (fwplus2*rW2 + (1.0 - fwplus2)*rP); else rhow2 = rW2;
	}
	if (iN2 > -1) {
		if (!bN2) rhon2 = rN2*rP / (fnplus2*rN2 + (1.0 - fnplus2)*rP); else rhon2 = rN2;
	}
	if (iS2 > -1) {
		if (!bS2) rhos2 = rS2*rP / (fsplus2*rS2 + (1.0 - fsplus2)*rP); else rhos2 = rS2;
	}
	if (iT2 > -1) {
		if (!bT2) rhot2 = rT2*rP / (ftplus2*rT2 + (1.0 - ftplus2)*rP); else rhot2 = rT2;
	}
	if (iB2 > -1) {
		if (!bB2) rhob2 = rB2*rP / (fbplus2*rB2 + (1.0 - fbplus2)*rP); else rhob2 = rB2;
	}


	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE3 > -1) {
		if (!bE3) rhoe3 = rE3*rP / (feplus3*rE3 + (1.0 - feplus3)*rP); else rhoe3 = rE3;
	}
	if (iW3 > -1) {
		if (!bW3) rhow3 = rW3*rP / (fwplus3*rW3 + (1.0 - fwplus3)*rP); else rhow3 = rW3;
	}
	if (iN3 > -1) {
		if (!bN3) rhon3 = rN3*rP / (fnplus3*rN3 + (1.0 - fnplus3)*rP); else rhon3 = rN3;
	}
	if (iS3 > -1) {
		if (!bS3) rhos3 = rS3*rP / (fsplus3*rS3 + (1.0 - fsplus3)*rP); else rhos3 = rS3;
	}
	if (iT3 > -1) {
		if (!bT3) rhot3 = rT3*rP / (ftplus3*rT3 + (1.0 - ftplus3)*rP); else rhot3 = rT3;
	}
	if (iB3 > -1) {
		if (!bB3) rhob3 = rB3*rP / (fbplus3*rB3 + (1.0 - fbplus3)*rP); else rhob3 = rB3;
	}

	// интерполяция плотности сделана так, чтобы выполнялись 
	// предельные соотношения.
	if (iE4 > -1) {
		if (!bE4) rhoe4 = rE4*rP / (feplus4*rE4 + (1.0 - feplus4)*rP); else rhoe4 = rE4;
	}
	if (iW4 > -1) {
		if (!bW4) rhow4 = rW4*rP / (fwplus4*rW4 + (1.0 - fwplus4)*rP); else rhow4 = rW4;
	}
	if (iN4 > -1) {
		if (!bN4) rhon4 = rN4*rP / (fnplus4*rN4 + (1.0 - fnplus4)*rP); else rhon4 = rN4;
	}
	if (iS4 > -1) {
		if (!bS4) rhos4 = rS4*rP / (fsplus4*rS4 + (1.0 - fsplus4)*rP); else rhos4 = rS4;
	}
	if (iT4 > -1) {
		if (!bT4) rhot4 = rT4*rP / (ftplus4*rT4 + (1.0 - ftplus4)*rP); else rhot4 = rT4;
	}
	if (iB4 > -1) {
		if (!bB4) rhob4 = rB4*rP / (fbplus4*rB4 + (1.0 - fbplus4)*rP); else rhob4 = rB4;
	}

	doublereal Fw = 0.0, Fe = 0.0, Fs = 0.0, Fn = 0.0, Ft = 0.0, Fb = 0.0;
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

		Fe = calcFg(bE, feplus, potent[VELOCITY_X_COMPONENT][iE], potent[VELOCITY_X_COMPONENT][iP],
			iP, E_SIDE, alpha[VELOCITY_X_COMPONENT], rhoe,
			dSqe, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iE], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);

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

		Fw = calcFg(bW, fwplus, potent[VELOCITY_X_COMPONENT][iW], potent[VELOCITY_X_COMPONENT][iP],
			iP, W_SIDE, alpha[VELOCITY_X_COMPONENT], rhow,
			dSqw, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iW], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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


		Fn = calcFg(bN, fnplus, potent[VELOCITY_Y_COMPONENT][iN], potent[VELOCITY_Y_COMPONENT][iP],
			iP, N_SIDE, alpha[VELOCITY_Y_COMPONENT], rhon,
			dSqn, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iN], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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


		Fs = calcFg(bS, fsplus, potent[VELOCITY_Y_COMPONENT][iS], potent[VELOCITY_Y_COMPONENT][iP],
			iP, S_SIDE, alpha[VELOCITY_Y_COMPONENT], rhos,
			dSqs, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iS], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);

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

		Ft = calcFg(bT, ftplus, potent[VELOCITY_Z_COMPONENT][iT], potent[VELOCITY_Z_COMPONENT][iP],
			iP, T_SIDE, alpha[VELOCITY_Z_COMPONENT], rhot,
			dSqt, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iT], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fb = calcFg(bB, fbplus, potent[VELOCITY_Z_COMPONENT][iB], potent[VELOCITY_Z_COMPONENT][iP],
			iP, B_SIDE, alpha[VELOCITY_Z_COMPONENT], rhob,
			dSqb, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iB], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal Fw2 = 0.0, Fe2 = 0.0, Fs2 = 0.0, Fn2 = 0.0, Ft2 = 0.0, Fb2 = 0.0;


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

		Fe2 = calcFg(bE2, feplus2, potent[VELOCITY_X_COMPONENT][iE2], potent[VELOCITY_X_COMPONENT][iP],
			iP, E_SIDE, alpha[VELOCITY_X_COMPONENT], rhoe2,
			dSqe2, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iE2], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fw2 = calcFg(bW2, fwplus2, potent[VELOCITY_X_COMPONENT][iW2], potent[VELOCITY_X_COMPONENT][iP],
			iP, W_SIDE, alpha[VELOCITY_X_COMPONENT], rhow2,
			dSqw2, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iW2], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fn2 = calcFg(bN2, fnplus2, potent[VELOCITY_Y_COMPONENT][iN2], potent[VELOCITY_Y_COMPONENT][iP],
			iP, N_SIDE, alpha[VELOCITY_Y_COMPONENT], rhon2,
			dSqn2, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iN2], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fs2 = calcFg(bS2, fsplus2, potent[VELOCITY_Y_COMPONENT][iS2], potent[VELOCITY_Y_COMPONENT][iP],
			iP, S_SIDE, alpha[VELOCITY_Y_COMPONENT], rhos2,
			dSqs2, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iS2], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Ft2 = calcFg(bT2, ftplus2, potent[VELOCITY_Z_COMPONENT][iT2], potent[VELOCITY_Z_COMPONENT][iP],
			iP, T_SIDE, alpha[VELOCITY_Z_COMPONENT], rhot2,
			dSqt2, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iT2], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fb2 = calcFg(bB2, fbplus2, potent[VELOCITY_Z_COMPONENT][iB2], potent[VELOCITY_Z_COMPONENT][iP],
			iP, B_SIDE, alpha[VELOCITY_Z_COMPONENT], rhob2,
			dSqb2, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iB2], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal Fw3 = 0.0, Fe3 = 0.0, Fs3 = 0.0, Fn3 = 0.0, Ft3 = 0.0, Fb3 = 0.0;


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

		Fe3 = calcFg(bE3, feplus3, potent[VELOCITY_X_COMPONENT][iE3], potent[VELOCITY_X_COMPONENT][iP],
			iP, E_SIDE, alpha[VELOCITY_X_COMPONENT], rhoe3,
			dSqe3, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iE3], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fw3 = calcFg(bW3, fwplus3, potent[VELOCITY_X_COMPONENT][iW3], potent[VELOCITY_X_COMPONENT][iP],
			iP, W_SIDE, alpha[VELOCITY_X_COMPONENT], rhow3,
			dSqw3, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iW3], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fn3 = calcFg(bN3, fnplus3, potent[VELOCITY_Y_COMPONENT][iN3], potent[VELOCITY_Y_COMPONENT][iP],
			iP, N_SIDE, alpha[VELOCITY_Y_COMPONENT], rhon3,
			dSqn3, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iN3], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fs3 = calcFg(bS3, fsplus3, potent[VELOCITY_Y_COMPONENT][iS3], potent[VELOCITY_Y_COMPONENT][iP],
			iP, S_SIDE, alpha[VELOCITY_Y_COMPONENT], rhos3,
			dSqs3, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iS3], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Ft3 = calcFg(bT3, ftplus3, potent[VELOCITY_Z_COMPONENT][iT3], potent[VELOCITY_Z_COMPONENT][iP],
			iP, T_SIDE, alpha[VELOCITY_Z_COMPONENT], rhot3,
			dSqt3, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iT3], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fb3 = calcFg(bB3, fbplus3, potent[VELOCITY_Z_COMPONENT][iB3], potent[VELOCITY_Z_COMPONENT][iP],
			iP, B_SIDE, alpha[VELOCITY_Z_COMPONENT], rhob3,
			dSqb3, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iB3], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.
	doublereal Fw4 = 0.0, Fe4 = 0.0, Fs4 = 0.0, Fn4 = 0.0, Ft4 = 0.0, Fb4 = 0.0;


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

		Fe4 = calcFg(bE4, feplus4, potent[VELOCITY_X_COMPONENT][iE4], potent[VELOCITY_X_COMPONENT][iP],
			iP, E_SIDE, alpha[VELOCITY_X_COMPONENT], rhoe4,
			dSqe4, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iE4], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fw4 = calcFg(bW4, fwplus4, potent[VELOCITY_X_COMPONENT][iW4], potent[VELOCITY_X_COMPONENT][iP],
			iP, W_SIDE, alpha[VELOCITY_X_COMPONENT], rhow4,
			dSqw4, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_X_COMPONENT][iW4], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fn4 = calcFg(bN4, fnplus4, potent[VELOCITY_Y_COMPONENT][iN4], potent[VELOCITY_Y_COMPONENT][iP],
			iP, N_SIDE, alpha[VELOCITY_Y_COMPONENT], rhon4,
			dSqn4, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iN4], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fs4 = calcFg(bS4, fsplus4, potent[VELOCITY_Y_COMPONENT][iS4], potent[VELOCITY_Y_COMPONENT][iP],
			iP, S_SIDE, alpha[VELOCITY_Y_COMPONENT], rhos4,
			dSqs4, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Y_COMPONENT][iS4], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Ft4 = calcFg(bT4, ftplus4, potent[VELOCITY_Z_COMPONENT][iT4], potent[VELOCITY_Z_COMPONENT][iP],
			iP, T_SIDE, alpha[VELOCITY_Z_COMPONENT], rhot4,
			dSqt4, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iT4], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
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

		Fb4 = calcFg(bB4, fbplus4, potent[VELOCITY_Z_COMPONENT][iB4], potent[VELOCITY_Z_COMPONENT][iP],
			iP, B_SIDE, alpha[VELOCITY_Z_COMPONENT], rhob4,
			dSqb4, dx*dy*dz, btimedepend,
			speedoldtimestep, mfoldtimestep,
			diag_coef[VELOCITY_Z_COMPONENT][iB4], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
			RCh, nvtx, neighbors_for_the_internal_node, maxelm,
			potent[PRESS], pa, diag_coef,
			bRhieChowi, bRhieChowb, bsimplelinearinterpol);
	}

	doublereal Fe_sum = Fe + Fe2 + Fe3 + Fe4;
	doublereal Fw_sum = Fw + Fw2 + Fw3 + Fw4;
	doublereal Fn_sum = Fn + Fn2 + Fn3 + Fn4;
	doublereal Fs_sum = Fs + Fs2 + Fs3 + Fs4;
	doublereal Ft_sum = Ft + Ft2 + Ft3 + Ft4;
	doublereal Fb_sum = Fb + Fb2 + Fb3 + Fb4;


	

	bool ISezai=true;
	if (ISezai) {

	    doublereal SpeedCorOlde = 0.0, SpeedCorOldw = 0.0,  SpeedCorOldn = 0.0,  SpeedCorOlds = 0.0,  SpeedCorOldt = 0.0,  SpeedCorOldb = 0.0;
		if (iE > -1) {
			if (!bE) {
				SpeedCorOlde = feplus*SpeedCorOld[VELOCITY_X_COMPONENT][iE] + (1.0 - feplus)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOlde = SpeedCorOld[VELOCITY_X_COMPONENT][iE];
			}
		}
		if (iN > -1) {
			if (!bN) {
				SpeedCorOldn = fnplus*SpeedCorOld[VELOCITY_Y_COMPONENT][iN] + (1.0 - fnplus)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOldn = SpeedCorOld[VELOCITY_Y_COMPONENT][iN];
			}
		}
		if (iT > -1) {
			if (!bT) {
				SpeedCorOldt = ftplus*SpeedCorOld[VELOCITY_Z_COMPONENT][iT] + (1.0 - ftplus)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldt = SpeedCorOld[VELOCITY_Z_COMPONENT][iT];
			}
		}
		if (iW > -1) {
			if (!bW) {
				SpeedCorOldw = fwplus*SpeedCorOld[VELOCITY_X_COMPONENT][iW] + (1.0 - fwplus)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOldw = SpeedCorOld[VELOCITY_X_COMPONENT][iW];
			}
		}
		if (iS > -1) {
			if (!bS) {
				SpeedCorOlds = fsplus*SpeedCorOld[VELOCITY_Y_COMPONENT][iS] + (1.0 - fsplus)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOlds = SpeedCorOld[VELOCITY_Y_COMPONENT][iS];
			}
		}
		if (iB > -1) {
			if (!bB) {
				SpeedCorOldb = fbplus*SpeedCorOld[VELOCITY_Z_COMPONENT][iB] + (1.0 - fbplus)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldb = SpeedCorOld[VELOCITY_Z_COMPONENT][iB];
			}
		}

		doublereal SpeedCorOlde2 = 0.0, SpeedCorOldw2 = 0.0, SpeedCorOldn2 = 0.0, SpeedCorOlds2 = 0.0, SpeedCorOldt2 = 0.0, SpeedCorOldb2 = 0.0;
		if (iE2 > -1) {
			if (!bE2) {
				SpeedCorOlde2 = feplus2*SpeedCorOld[VELOCITY_X_COMPONENT][iE2] + (1.0 - feplus2)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOlde2 = SpeedCorOld[VELOCITY_X_COMPONENT][iE2];
			}
		}
		if (iN2 > -1) {
			if (!bN2) {
				SpeedCorOldn2 = fnplus2*SpeedCorOld[VELOCITY_Y_COMPONENT][iN2] + (1.0 - fnplus2)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOldn2 = SpeedCorOld[VELOCITY_Y_COMPONENT][iN2];
			}
		}
		if (iT2 > -1) {
			if (!bT2) {
				SpeedCorOldt2 = ftplus2*SpeedCorOld[VELOCITY_Z_COMPONENT][iT2] + (1.0 - ftplus2)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldt2 = SpeedCorOld[VELOCITY_Z_COMPONENT][iT2];
			}
		}
		if (iW2 > -1) {
			if (!bW2) {
				SpeedCorOldw2 = fwplus2*SpeedCorOld[VELOCITY_X_COMPONENT][iW2] + (1.0 - fwplus2)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOldw2 = SpeedCorOld[VELOCITY_X_COMPONENT][iW2];
			}
		}
		if (iS2 > -1) {
			if (!bS2) {
				SpeedCorOlds2 = fsplus2*SpeedCorOld[VELOCITY_Y_COMPONENT][iS2] + (1.0 - fsplus2)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOlds2 = SpeedCorOld[VELOCITY_Y_COMPONENT][iS2];
			}
		}
		if (iB2 > -1) {
			if (!bB2) {
				SpeedCorOldb2 = fbplus2*SpeedCorOld[VELOCITY_Z_COMPONENT][iB2] + (1.0 - fbplus2)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldb2 = SpeedCorOld[VELOCITY_Z_COMPONENT][iB2];
			}
		}

		doublereal SpeedCorOlde3 = 0.0, SpeedCorOldw3 = 0.0, SpeedCorOldn3 = 0.0, SpeedCorOlds3 = 0.0, SpeedCorOldt3 = 0.0, SpeedCorOldb3 = 0.0;
		if (iE3 > -1) {
			if (!bE3) {
				SpeedCorOlde3 = feplus3*SpeedCorOld[VELOCITY_X_COMPONENT][iE3] + (1.0 - feplus3)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOlde3 = SpeedCorOld[VELOCITY_X_COMPONENT][iE3];
			}
		}
		if (iN3 > -1) {
			if (!bN3) {
				SpeedCorOldn3 = fnplus3*SpeedCorOld[VELOCITY_Y_COMPONENT][iN3] + (1.0 - fnplus3)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOldn3 = SpeedCorOld[VELOCITY_Y_COMPONENT][iN3];
			}
		}
		if (iT3 > -1) {
			if (!bT3) {
				SpeedCorOldt3 = ftplus3*SpeedCorOld[VELOCITY_Z_COMPONENT][iT3] + (1.0 - ftplus3)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldt3 = SpeedCorOld[VELOCITY_Z_COMPONENT][iT3];
			}
		}
		if (iW3 > -1) {
			if (!bW3) {
				SpeedCorOldw3 = fwplus3*SpeedCorOld[VELOCITY_X_COMPONENT][iW3] + (1.0 - fwplus3)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOldw3 = SpeedCorOld[VELOCITY_X_COMPONENT][iW3];
			}
		}
		if (iS3 > -1) {
			if (!bS3) {
				SpeedCorOlds3 = fsplus3*SpeedCorOld[VELOCITY_Y_COMPONENT][iS3] + (1.0 - fsplus3)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOlds3 = SpeedCorOld[VELOCITY_Y_COMPONENT][iS3];
			}
		}
		if (iB3 > -1) {
			if (!bB3) {
				SpeedCorOldb3 = fbplus3*SpeedCorOld[VELOCITY_Z_COMPONENT][iB3] + (1.0 - fbplus3)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldb3 = SpeedCorOld[VELOCITY_Z_COMPONENT][iB3];
			}
		}


		doublereal SpeedCorOlde4 = 0.0, SpeedCorOldw4 = 0.0, SpeedCorOldn4 = 0.0, SpeedCorOlds4 = 0.0, SpeedCorOldt4 = 0.0, SpeedCorOldb4 = 0.0;
		if (iE4 > -1) {
			if (!bE4) {
				SpeedCorOlde4 = feplus4*SpeedCorOld[VELOCITY_X_COMPONENT][iE4] + (1.0 - feplus4)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOlde4 = SpeedCorOld[VELOCITY_X_COMPONENT][iE4];
			}
		}
		if (iN4 > -1) {
			if (!bN4) {
				SpeedCorOldn4 = fnplus4*SpeedCorOld[VELOCITY_Y_COMPONENT][iN4] + (1.0 - fnplus4)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOldn4 = SpeedCorOld[VELOCITY_Y_COMPONENT][iN4];
			}
		}
		if (iT4 > -1) {
			if (!bT4) {
				SpeedCorOldt4 = ftplus4*SpeedCorOld[VELOCITY_Z_COMPONENT][iT4] + (1.0 - ftplus4)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldt4 = SpeedCorOld[VELOCITY_Z_COMPONENT][iT4];
			}
		}
		if (iW4 > -1) {
			if (!bW4) {
				SpeedCorOldw4 = fwplus4*SpeedCorOld[VELOCITY_X_COMPONENT][iW4] + (1.0 - fwplus4)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
			}
			else {
				SpeedCorOldw4 = SpeedCorOld[VELOCITY_X_COMPONENT][iW4];
			}
		}
		if (iS4 > -1) {
			if (!bS4) {
				SpeedCorOlds4 = fsplus4*SpeedCorOld[VELOCITY_Y_COMPONENT][iS4] + (1.0 - fsplus4)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
			}
			else {
				SpeedCorOlds4 = SpeedCorOld[VELOCITY_Y_COMPONENT][iS4];
			}
		}
		if (iB4 > -1) {
			if (!bB4) {
				SpeedCorOldb4 = fbplus4*SpeedCorOld[VELOCITY_Z_COMPONENT][iB4] + (1.0 - fbplus4)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
			}
			else {
				SpeedCorOldb4 = SpeedCorOld[VELOCITY_Z_COMPONENT][iB4];
			}
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
        //mfcurrentretune[ESIDE]=Fe+(1.0-alpha[VX])*(mfold[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    //mfcurrentretune[NSIDE]=Fn +(1.0-alpha[VY])*(mfold[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    //mfcurrentretune[TSIDE]=Ft +(1.0-alpha[VZ])*(mfold[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    //mfcurrentretune[WSIDE]=Fw +(1.0-alpha[VX])*(mfold[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    //mfcurrentretune[SSIDE]=Fs +(1.0-alpha[VY])*(mfold[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    //mfcurrentretune[BSIDE]=Fb +(1.0-alpha[VZ])*(mfold[BSIDE]-rhob*SpeedCorOldb*dx*dy);

		// возвращаем значение потока на грани КО.
		// С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
		// mfold - значение массового потока с предыдущей итерации.
		mfcurrentretune[E_SIDE] = Fe_sum + (1.0 - alpha[VELOCITY_X_COMPONENT])*(mfold[E_SIDE] - rhoe*SpeedCorOlde*dSqe - rhoe2 * SpeedCorOlde2*dSqe2 - rhoe3 * SpeedCorOlde3*dSqe3 - rhoe4 * SpeedCorOlde4*dSqe4);
		mfcurrentretune[N_SIDE] = Fn_sum + (1.0 - alpha[VELOCITY_Y_COMPONENT])*(mfold[N_SIDE] - rhon*SpeedCorOldn*dSqn - rhon2 * SpeedCorOldn2*dSqn2 - rhon3 * SpeedCorOldn3*dSqn3 - rhon4 * SpeedCorOldn4*dSqn4);
		mfcurrentretune[T_SIDE] = Ft_sum + (1.0 - alpha[VELOCITY_Z_COMPONENT])*(mfold[T_SIDE] - rhot*SpeedCorOldt*dSqt - rhot2 * SpeedCorOldt2*dSqt2 - rhot3 * SpeedCorOldt3*dSqt3 - rhot4 * SpeedCorOldt4*dSqt4);
		mfcurrentretune[W_SIDE] = Fw_sum + (1.0 - alpha[VELOCITY_X_COMPONENT])*(mfold[W_SIDE] - rhow*SpeedCorOldw*dSqw - rhow2 * SpeedCorOldw2*dSqw2 - rhow3 * SpeedCorOldw3*dSqw3 - rhow4 * SpeedCorOldw4*dSqw4);
		mfcurrentretune[S_SIDE] = Fs_sum + (1.0 - alpha[VELOCITY_Y_COMPONENT])*(mfold[S_SIDE] - rhos*SpeedCorOlds*dSqs - rhos2 * SpeedCorOlds2*dSqs2 - rhos3 * SpeedCorOlds3*dSqs3 - rhos4 * SpeedCorOlds4*dSqs4);
		mfcurrentretune[B_SIDE] = Fb_sum + (1.0 - alpha[VELOCITY_Z_COMPONENT])*(mfold[B_SIDE] - rhob*SpeedCorOldb*dSqb - rhob2 * SpeedCorOldb2*dSqb2 - rhob3 * SpeedCorOldb3*dSqb3 - rhob4 * SpeedCorOldb4*dSqb4);
	}

	if (!ISezai) {
	   //mfcurrentretune[ESIDE]=Fe;
	   //mfcurrentretune[NSIDE]=Fn;
	   //mfcurrentretune[TSIDE]=Ft;
	   //mfcurrentretune[WSIDE]=Fw;
	   //mfcurrentretune[SSIDE]=Fs;
	   //mfcurrentretune[BSIDE]=Fb;
		mfcurrentretune[E_SIDE] = Fe_sum;
		mfcurrentretune[N_SIDE] = Fn_sum;
		mfcurrentretune[T_SIDE] = Ft_sum;
		mfcurrentretune[W_SIDE] = Fw_sum;
		mfcurrentretune[S_SIDE] = Fs_sum;
		mfcurrentretune[B_SIDE] = Fb_sum;
	}

} // return_correct_mass_flux

// возвращает скорректированный массовый поток.
// скорректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
// скорректированной скорости и давления, а также сохранённых диагональных 
// коэффициентов матрицы СЛАУ для компонент скорости.
// В данной реализации учитывается сглаженное псевдовремя tau.
// реализовано 22 июня 2012 года. Основывается на сглаженном псевдовремени tau.
void return_calc_correct_mass_flux2(integer iP, doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
	int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, doublereal* alpha, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, bool bsimplelinearinterpol,
						doublereal** SpeedCorOld, doublereal *mfold, doublereal* tau)
{

	// SpeedCorOld - скорректированная скорость на предыдущей итерации.
	// tau - сглаженное псевдовремя.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border.

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

	Fe=calcFg2(bE, feplus, iP, E_SIDE,  rhoe, 
			   dy*dz, btimedepend, 
			   speedoldtimestep,  mfoldtimestep,
		       dtimestep,	 RCh, neighbors_for_the_internal_node, 
			   bRhieChowi, bRhieChowb, 
			   bsimplelinearinterpol, tau, dxe, potent);

	Fw=calcFg2(bW, fwplus, iP, W_SIDE, rhow, 
			   dy*dz, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep,	RCh, neighbors_for_the_internal_node, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dxw, potent);

	Fn=calcFg2(bN, fnplus, iP, N_SIDE, rhon, 
			   dx*dz, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep, RCh, neighbors_for_the_internal_node, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dyn, potent);

	Fs=calcFg2(bS, fsplus, iP, S_SIDE, rhos, 
			   dx*dz, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep, RCh, neighbors_for_the_internal_node, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dys, potent);

	Ft=calcFg2(bT, ftplus, iP, T_SIDE, rhot, 
		       dx*dy, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep,	RCh, neighbors_for_the_internal_node, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dzt, potent);

	Fb=calcFg2(bB, fbplus, iP, B_SIDE, rhob, 
			   dx*dy, btimedepend, 
			   speedoldtimestep, mfoldtimestep,
			   dtimestep,	RCh, neighbors_for_the_internal_node, 
			   bRhieChowi, bRhieChowb, bsimplelinearinterpol,
			   tau, dzb, potent);

	bool ISezai=true;
	if (ISezai) {

	    doublereal SpeedCorOlde, SpeedCorOldw,  SpeedCorOldn,  SpeedCorOlds,  SpeedCorOldt,  SpeedCorOldb; 
	    if (!bE) { 
		   SpeedCorOlde=feplus*SpeedCorOld[VELOCITY_X_COMPONENT][iE]+(1.0-feplus)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
	    }
	    else {
	       SpeedCorOlde=SpeedCorOld[VELOCITY_X_COMPONENT][iE];
	    }
	    if (!bN) { 
		   SpeedCorOldn=fnplus*SpeedCorOld[VELOCITY_Y_COMPONENT][iN]+(1.0-fnplus)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
	    }
	    else {
	       SpeedCorOldn=SpeedCorOld[VELOCITY_Y_COMPONENT][iN];
	    }
	    if (!bT) { 
		   SpeedCorOldt=ftplus*SpeedCorOld[VELOCITY_Z_COMPONENT][iT]+(1.0-ftplus)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
	    }
	    else {
	       SpeedCorOldt=SpeedCorOld[VELOCITY_Z_COMPONENT][iT];
	    }
	    if (!bW) { 
		   SpeedCorOldw=fwplus*SpeedCorOld[VELOCITY_X_COMPONENT][iW]+(1.0-fwplus)*SpeedCorOld[VELOCITY_X_COMPONENT][iP];
     	}
        else {
	       SpeedCorOldw=SpeedCorOld[VELOCITY_X_COMPONENT][iW];
	    }
	    if (!bS) { 
		   SpeedCorOlds=fsplus*SpeedCorOld[VELOCITY_Y_COMPONENT][iS]+(1.0-fsplus)*SpeedCorOld[VELOCITY_Y_COMPONENT][iP];
	    }
	    else {
	       SpeedCorOlds=SpeedCorOld[VELOCITY_Y_COMPONENT][iS];
	    }
	    if (!bB) { 
		   SpeedCorOldb=fbplus*SpeedCorOld[VELOCITY_Z_COMPONENT][iB]+(1.0-fbplus)*SpeedCorOld[VELOCITY_Z_COMPONENT][iP];
	    }
	    else {
	       SpeedCorOldb=SpeedCorOld[VELOCITY_Z_COMPONENT][iB];
	    }

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	   // mfold - значение массового потока с предыдущей итерации.
        mfcurrentretune[E_SIDE]=Fe+(1.0-alpha[VELOCITY_X_COMPONENT])*(mfold[E_SIDE]-rhoe*SpeedCorOlde*dy*dz);
	    mfcurrentretune[N_SIDE]=Fn+(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfold[N_SIDE]-rhon*SpeedCorOldn*dx*dz);
	    mfcurrentretune[T_SIDE]=Ft+(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfold[T_SIDE]-rhot*SpeedCorOldt*dx*dy);
	    mfcurrentretune[W_SIDE]=Fw+(1.0-alpha[VELOCITY_X_COMPONENT])*(mfold[W_SIDE]-rhow*SpeedCorOldw*dy*dz);
	    mfcurrentretune[S_SIDE]=Fs+(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfold[S_SIDE]-rhos*SpeedCorOlds*dx*dz);
	    mfcurrentretune[B_SIDE]=Fb+(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfold[B_SIDE]-rhob*SpeedCorOldb*dx*dy);
	}

	if (!ISezai) {
	   mfcurrentretune[E_SIDE]=Fe;
	   mfcurrentretune[N_SIDE]=Fn;
	   mfcurrentretune[T_SIDE]=Ft;
	   mfcurrentretune[W_SIDE]=Fw;
	   mfcurrentretune[S_SIDE]=Fs;
	   mfcurrentretune[B_SIDE]=Fb;
	}

} // return_correct_mass_flux2

// возвращает скорректированный массовый поток.
// скорректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
// скорректированной скорости и давления, а также сохранённых диагональных 
// коэффициентов матрицы СЛАУ для компонент скорости.
// В данной реализации учитывается сглаженное псевдовремя tau.
// реализовано 23 июня 2012 года. Основывается на сглаженном псевдовремени tau.
void return_calc_correct_mass_flux3(integer iP, doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
	int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, doublereal* alpha, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, bool bsimplelinearinterpol,
						doublereal** SpeedCorOld, doublereal *mfold, doublereal** tau,
	                    BOUND* border_neighbor, integer *ilevel_alice, int* ptr)
{

	// SpeedCorOld - скорректированная скорость на предыдущей итерации.
	// tau - сглаженное псевдовремя.

	// Если bsimplelinearinterpol равен true то выполняется простая линейная интерполяция скорости на грань контрольного объёма.

	// По-видимому имеет смысл включать поправку Рхи-Чоу только во внутренней грани,
	// для граничной грани скорость задана (из граничных условий) и по-видимому не 
	// требуется применять к ней монотонизирующую поправку.

	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border.

	// iP - номер центрального контрольного объёма
	integer iE=-1, iN=-1, iT=-1, iW=-1, iS=-1, iB=-1; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];


	// 26.09.2016 Добавок для АЛИС сетки.
	integer iE2=-1, iN2=-1, iT2=-1, iW2=-1, iS2=-1, iB2=-1; // номера соседних контрольных объёмов
	integer iE3=-1, iN3=-1, iT3=-1, iW3=-1, iS3=-1, iB3=-1; // номера соседних контрольных объёмов
	integer iE4=-1, iN4=-1, iT4=-1, iW4=-1, iS4=-1, iB4=-1; // номера соседних контрольных объёмов

										  // -1 если не используется и [0..maxelm+maxbound-1] если используется.

	iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
	iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
	iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
	iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
	iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
	iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


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


	

	doublereal Fw = 0.0, Fe = 0.0, Fs = 0.0, Fn = 0.0, Ft = 0.0, Fb = 0.0;

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


		Fe = calcFg3(bE, feplus, iP, E_SIDE, rhoe,
			dSqe, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxe, potent, 1);
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

		Fw = calcFg3(bW, fwplus, iP, W_SIDE, rhow,
			dSqw, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxw, potent, 1);
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

		Fn = calcFg3(bN, fnplus, iP, N_SIDE, rhon,
			dSqn, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dyn, potent, 1);
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

		Fs = calcFg3(bS, fsplus, iP, S_SIDE, rhos,
			dSqs, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dys, potent, 1);
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

		Ft = calcFg3(bT, ftplus, iP, T_SIDE, rhot,
			dSqt, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzt, potent, 1);
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

		Fb = calcFg3(bB, fbplus, iP, B_SIDE, rhob,
			dSqb, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzb, potent, 1);
	}

	doublereal dSqe2 = 0.0, dSqw2 = 0.0, dSqn2 = 0.0, dSqs2 = 0.0, dSqt2 = 0.0, dSqb2 = 0.0; // площадь грани.
	doublereal Fw2 = 0.0, Fe2 = 0.0, Fs2 = 0.0, Fn2 = 0.0, Ft2 = 0.0, Fb2 = 0.0;


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

		Fe2 = calcFg3(bE2, feplus2, iP, E_SIDE, rhoe2,
			dSqe2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxe2, potent, 2);
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

		Fw2 = calcFg3(bW2, fwplus2, iP, W_SIDE, rhow2,
			dSqw2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxw2, potent, 2);
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

		Fn2 = calcFg3(bN2, fnplus2, iP, N_SIDE, rhon2,
			dSqn2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dyn2, potent, 2);
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

		Fs2 = calcFg3(bS2, fsplus2, iP, S_SIDE, rhos2,
			dSqs2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dys2, potent, 2);
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

		Ft2 = calcFg3(bT2, ftplus2, iP, T_SIDE, rhot2,
			dSqt2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzt2, potent, 2);
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

		Fb2 = calcFg3(bB2, fbplus2, iP, B_SIDE, rhob2,
			dSqb2, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzb2, potent, 2);
	}


	doublereal dSqe3 = 0.0, dSqw3 = 0.0, dSqn3 = 0.0, dSqs3 = 0.0, dSqt3 = 0.0, dSqb3 = 0.0; // площадь грани.
	doublereal Fw3 = 0.0, Fe3 = 0.0, Fs3 = 0.0, Fn3 = 0.0, Ft3 = 0.0, Fb3 = 0.0;


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

		Fe3 = calcFg3(bE3, feplus3, iP, E_SIDE, rhoe3,
			dSqe3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxe3, potent, 3);
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

		Fw3 = calcFg3(bW3, fwplus3, iP, W_SIDE, rhow3,
			dSqw3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxw3, potent, 3);
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

		Fn3 = calcFg3(bN3, fnplus3, iP, N_SIDE, rhon3,
			dSqn3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dyn3, potent, 3);
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

		Fs3 = calcFg3(bS3, fsplus3, iP, S_SIDE, rhos3,
			dSqs3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dys3, potent, 3);
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

		Ft3 = calcFg3(bT3, ftplus3, iP, T_SIDE, rhot3,
			dSqt3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzt3, potent, 3);
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

		Fb3 = calcFg3(bB3, fbplus3, iP, B_SIDE, rhob3,
			dSqb3, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzb3, potent, 3);
	}

	doublereal dSqe4 = 0.0, dSqw4 = 0.0, dSqn4 = 0.0, dSqs4 = 0.0, dSqt4 = 0.0, dSqb4 = 0.0; // площадь грани.
	doublereal Fw4 = 0.0, Fe4 = 0.0, Fs4 = 0.0, Fn4 = 0.0, Ft4 = 0.0, Fb4 = 0.0;


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

		Fe4 = calcFg3(bE4, feplus4, iP, E_SIDE, rhoe4,
			dSqe4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxe4, potent, 4);
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

		Fw4 = calcFg3(bW4, fwplus4, iP, W_SIDE, rhow4,
			dSqw4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxw4, potent, 4);
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

		Fn4 = calcFg3(bN4, fnplus4, iP, N_SIDE, rhon4,
			dSqn4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dyn4, potent, 4);
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

		Fs4 = calcFg3(bS4, fsplus4, iP, S_SIDE, rhos4,
			dSqs4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dys4, potent, 4);
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

		Ft4 = calcFg3(bT4, ftplus4, iP, T_SIDE, rhot4,
			dSqt4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzt4, potent, 4);
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

		Fb4 = calcFg3(bB4, fbplus4, iP, B_SIDE, rhob4,
			dSqb4, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
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
        mfcurrentretune[E_SIDE]=Fe_sum+(1.0-alpha[VELOCITY_X_COMPONENT])*(mfold[E_SIDE]-rhoe*SpeedCorOlde*dSqe - rhoe2 * SpeedCorOlde2*dSqe2 - rhoe3 * SpeedCorOlde3*dSqe3 - rhoe4 * SpeedCorOlde4*dSqe4);
	    mfcurrentretune[N_SIDE]=Fn_sum +(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfold[N_SIDE]-rhon*SpeedCorOldn*dSqn - rhon2 * SpeedCorOldn2*dSqn2 - rhon3 * SpeedCorOldn3*dSqn3 - rhon4 * SpeedCorOldn4*dSqn4);
	    mfcurrentretune[T_SIDE]=Ft_sum +(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfold[T_SIDE]-rhot*SpeedCorOldt*dSqt - rhot2 * SpeedCorOldt2*dSqt2 - rhot3 * SpeedCorOldt3*dSqt3 - rhot4 * SpeedCorOldt4*dSqt4);
	    mfcurrentretune[W_SIDE]=Fw_sum +(1.0-alpha[VELOCITY_X_COMPONENT])*(mfold[W_SIDE]-rhow*SpeedCorOldw*dSqw - rhow2 * SpeedCorOldw2*dSqw2 - rhow3 * SpeedCorOldw3*dSqw3 - rhow4 * SpeedCorOldw4*dSqw4);
	    mfcurrentretune[S_SIDE]=Fs_sum +(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfold[S_SIDE]-rhos*SpeedCorOlds*dSqs - rhos2 * SpeedCorOlds2*dSqs2 - rhos3 * SpeedCorOlds3*dSqs3 - rhos4 * SpeedCorOlds4*dSqs4);
	    mfcurrentretune[B_SIDE]=Fb_sum +(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfold[B_SIDE]-rhob*SpeedCorOldb*dSqb - rhob2 * SpeedCorOldb2*dSqb2 - rhob3 * SpeedCorOldb3*dSqb3 - rhob4 * SpeedCorOldb4*dSqb4);
	}

	if (!ISezai) {
	   mfcurrentretune[E_SIDE]=Fe_sum;
	   mfcurrentretune[N_SIDE]=Fn_sum;
	   mfcurrentretune[T_SIDE]=Ft_sum;
	   mfcurrentretune[W_SIDE]=Fw_sum;
	   mfcurrentretune[S_SIDE]=Fs_sum;
	   mfcurrentretune[B_SIDE]=Fb_sum;
	}

} // return_correct_mass_flux3



// Составляет матрицу для уравнения 
// поправки давления
void my_elmatr_quad_PAm(integer iP, equation3D** &sl, equation3D_bon** &slb,  
						doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
						int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, doublereal **diag_coef,
						doublereal* alpha, doublereal dbeta, doublereal** &rhie_chow, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep) {

    // 31 марта 2012 из этой функции удалено много закомментированного кода. Его можно найти в backup`ах AliceFlowv0_07.

    /* btimedepend - стационарный (false), нестационарный (true).
	*  dtimestep - величина текущего шага по времени.
	*  
	*
	*  mfoldtimestep - массовый поток через каждую грань КО с предыдущего временного слоя. (неизменяемая в данной функции величина).
	*  mfcurrentretune - возвращаемый текущий массовый поток через каждую грань КО.
	*  КО имеет 6 граней: E,N,T,W,S,B.
	*
	* speedoldtimestep[VX or VY or VZ][iP: 0..maxelm+maxbound-1] - поле скорости с предыдущего временного слоя.
	*
	*/

	bool bhighorder=false;

	doublereal eps=admission; // для отделения вещественного нуля.
	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
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
	if (!bE) apue=sl[VELOCITY_X_COMPONENT][iE].ap*sl[VELOCITY_X_COMPONENT][iP].ap/(feplus*sl[VELOCITY_X_COMPONENT][iE].ap+(1-feplus)*sl[VELOCITY_X_COMPONENT][iP].ap); else apue=slb[VELOCITY_X_COMPONENT][iE-maxelm].aw;
	if (!bW) apuw=sl[VELOCITY_X_COMPONENT][iW].ap*sl[VELOCITY_X_COMPONENT][iP].ap/(fwplus*sl[VELOCITY_X_COMPONENT][iW].ap+(1-fwplus)*sl[VELOCITY_X_COMPONENT][iP].ap); else apuw=slb[VELOCITY_X_COMPONENT][iW-maxelm].aw;
	if (!bN) apvn=sl[VELOCITY_Y_COMPONENT][iN].ap*sl[VELOCITY_Y_COMPONENT][iP].ap/(fnplus*sl[VELOCITY_Y_COMPONENT][iN].ap+(1-fnplus)*sl[VELOCITY_Y_COMPONENT][iP].ap); else apvn=slb[VELOCITY_Y_COMPONENT][iN-maxelm].aw;
	if (!bS) apvs=sl[VELOCITY_Y_COMPONENT][iS].ap*sl[VELOCITY_Y_COMPONENT][iP].ap/(fsplus*sl[VELOCITY_Y_COMPONENT][iS].ap+(1-fsplus)*sl[VELOCITY_Y_COMPONENT][iP].ap); else apvs=slb[VELOCITY_Y_COMPONENT][iS-maxelm].aw;
	if (!bT) apwt=sl[VELOCITY_Z_COMPONENT][iT].ap*sl[VELOCITY_Z_COMPONENT][iP].ap/(ftplus*sl[VELOCITY_Z_COMPONENT][iT].ap+(1-ftplus)*sl[VELOCITY_Z_COMPONENT][iP].ap); else apwt=slb[VELOCITY_Z_COMPONENT][iT-maxelm].aw;
	if (!bB) apwb=sl[VELOCITY_Z_COMPONENT][iB].ap*sl[VELOCITY_Z_COMPONENT][iP].ap/(fbplus*sl[VELOCITY_Z_COMPONENT][iB].ap+(1-fbplus)*sl[VELOCITY_Z_COMPONENT][iP].ap); else apwb=slb[VELOCITY_Z_COMPONENT][iB-maxelm].aw;

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
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
		// SIMPLE
		de=alpha[VELOCITY_X_COMPONENT]*dy*dz/apue; dw=alpha[VELOCITY_X_COMPONENT]*dy*dz/apuw; 
	    dn=alpha[VELOCITY_Y_COMPONENT]*dx*dz/apvn; ds=alpha[VELOCITY_Y_COMPONENT]*dx*dz/apvs;
	    dt=alpha[VELOCITY_Z_COMPONENT]*dx*dy/apwt; db=alpha[VELOCITY_Z_COMPONENT]*dx*dy/apwb;
	}
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC
		de=alpha[VELOCITY_X_COMPONENT]*dy*dz/((1.0-alpha[VELOCITY_X_COMPONENT])*apue); dw=alpha[VELOCITY_X_COMPONENT]*dy*dz/((1.0-alpha[VELOCITY_X_COMPONENT])*apuw); 
	    dn=alpha[VELOCITY_Y_COMPONENT]*dx*dz/((1.0-alpha[VELOCITY_Y_COMPONENT])*apvn); ds=alpha[VELOCITY_Y_COMPONENT]*dx*dz/((1.0-alpha[VELOCITY_Y_COMPONENT])*apvs);
	    dt=alpha[VELOCITY_Z_COMPONENT]*dx*dy/((1.0-alpha[VELOCITY_Z_COMPONENT])*apwt); db=alpha[VELOCITY_Z_COMPONENT]*dx*dy/((1.0-alpha[VELOCITY_Z_COMPONENT])*apwb);
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
		// если bborder   то мы находимся вблизи граничного узла.
		bool bborder=false;
		doublereal myflux=0.0;
		myflux=De*(dxe*DFDXiP(potent[PAM], iP, E_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iE]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
		myflux=Dw*(-dxw*DFDXiP(potent[PAM], iP, W_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iW]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dn*(dyn*DFDXiP(potent[PAM], iP, N_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iN]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
        myflux=Ds*(-dys*DFDXiP(potent[PAM], iP, S_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iS]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dt*(dzt*DFDXiP(potent[PAM], iP, T_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iT]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Db*(-dzb*DFDXiP(potent[PAM], iP, B_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iB]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
	}



    doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 
	     
	
	Fe=calcFg(bE, feplus, potent[VELOCITY_X_COMPONENT][iE], potent[VELOCITY_X_COMPONENT][iP],  
	          iP, E_SIDE, alpha[VELOCITY_X_COMPONENT], rhoe, 
		      dy*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VELOCITY_X_COMPONENT][iE], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
		      RCh, nvtx, neighbors_for_the_internal_node, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fw=calcFg(bW, fwplus, potent[VELOCITY_X_COMPONENT][iW], potent[VELOCITY_X_COMPONENT][iP],  
	          iP, W_SIDE, alpha[VELOCITY_X_COMPONENT], rhow, 
		      dy*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VELOCITY_X_COMPONENT][iW], diag_coef[VELOCITY_X_COMPONENT][iP], dtimestep,
		      RCh, nvtx, neighbors_for_the_internal_node, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fn=calcFg(bN, fnplus, potent[VELOCITY_Y_COMPONENT][iN], potent[VELOCITY_Y_COMPONENT][iP],  
	          iP, N_SIDE, alpha[VELOCITY_Y_COMPONENT], rhon, 
		      dx*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VELOCITY_Y_COMPONENT][iN], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
		      RCh, nvtx, neighbors_for_the_internal_node, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fs=calcFg(bS, fsplus, potent[VELOCITY_Y_COMPONENT][iS], potent[VELOCITY_Y_COMPONENT][iP],  
	          iP, S_SIDE, alpha[VELOCITY_Y_COMPONENT], rhos, 
		      dx*dz, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VELOCITY_Y_COMPONENT][iS], diag_coef[VELOCITY_Y_COMPONENT][iP], dtimestep,
		      RCh, nvtx, neighbors_for_the_internal_node, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

    Ft=calcFg(bT, ftplus, potent[VELOCITY_Z_COMPONENT][iT], potent[VELOCITY_Z_COMPONENT][iP],  
	          iP, T_SIDE, alpha[VELOCITY_Z_COMPONENT], rhot, 
		      dx*dy, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VELOCITY_Z_COMPONENT][iT], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
		      RCh, nvtx, neighbors_for_the_internal_node, maxelm, 
			  potent[PRESS],  pa, diag_coef,
		      bRhieChowi, bRhieChowb,false);

	Fb=calcFg(bB, fbplus, potent[VELOCITY_Z_COMPONENT][iB], potent[VELOCITY_Z_COMPONENT][iP],  
	          iP, B_SIDE, alpha[VELOCITY_Z_COMPONENT], rhob, 
		      dx*dy, dx*dy*dz, btimedepend, 
		      speedoldtimestep, mfoldtimestep,
		      diag_coef[VELOCITY_Z_COMPONENT][iB], diag_coef[VELOCITY_Z_COMPONENT][iP], dtimestep,
		      RCh, nvtx, neighbors_for_the_internal_node, maxelm, 
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
	    // дополнительная поправка осуществляется на основе скорректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    //mfcurrentretune[ESIDE]=Fe+(1.0-alpha[VX])*(mfcurrentretune[ESIDE]-rhoe*SpeedCorOlde*dy*dz);
	    //mfcurrentretune[NSIDE]=Fn+(1.0-alpha[VY])*(mfcurrentretune[NSIDE]-rhon*SpeedCorOldn*dx*dz);
	    //mfcurrentretune[TSIDE]=Ft+(1.0-alpha[VZ])*(mfcurrentretune[TSIDE]-rhot*SpeedCorOldt*dx*dy);
	    //mfcurrentretune[WSIDE]=Fw+(1.0-alpha[VX])*(mfcurrentretune[WSIDE]-rhow*SpeedCorOldw*dy*dz);
	    //mfcurrentretune[SSIDE]=Fs+(1.0-alpha[VY])*(mfcurrentretune[SSIDE]-rhos*SpeedCorOlds*dx*dz);
	    //mfcurrentretune[BSIDE]=Fb+(1.0-alpha[VZ])*(mfcurrentretune[BSIDE]-rhob*SpeedCorOldb*dx*dy);

		// возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    // дополнительная поправка осуществляется на основе скорректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    Fe+=(1.0-alpha[VELOCITY_X_COMPONENT])*(mfcurrentretune[E_SIDE]-rhoe*SpeedCorOlde*dy*dz);
	    Fn+=(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfcurrentretune[N_SIDE]-rhon*SpeedCorOldn*dx*dz);
	    Ft+=(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfcurrentretune[T_SIDE]-rhot*SpeedCorOldt*dx*dy);
	    Fw+=(1.0-alpha[VELOCITY_X_COMPONENT])*(mfcurrentretune[W_SIDE]-rhow*SpeedCorOldw*dy*dz);
	    Fs+=(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfcurrentretune[S_SIDE]-rhos*SpeedCorOlds*dx*dz);
	    Fb+=(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfcurrentretune[B_SIDE]-rhob*SpeedCorOldb*dx*dy);

	} 
	

	
	// возвращаем значение потока на грани КО.
	mfcurrentretune[E_SIDE]=Fe;
	mfcurrentretune[N_SIDE]=Fn;
	mfcurrentretune[T_SIDE]=Ft;
	mfcurrentretune[W_SIDE]=Fw;
	mfcurrentretune[S_SIDE]=Fs;
	mfcurrentretune[B_SIDE]=Fb;
	
	

	sl[PAM][iP].b=Fw-Fe+Fs-Fn+Fb-Ft+baddDFLUX2;


	// Симметризация СЛАУ:

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
						doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
						int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm,
						doublereal* alpha, doublereal dbeta, doublereal** &rhie_chow, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, doublereal* tau) {

    // 31 марта 2012 из этой функции удалено много закомментированного кода. Его можно найти в backup`ах AliceFlowv0_07.

	// tau - заранее вычисленное псевдовремя, оно служит коэффициентом диффузии в эллиптическом операторе, к тому-же оно
	// также входит в монотонизирующую поправку Рхи-Чоу (1983 года).

    /* btimedepend - стационарный (false), нестационарный (true).
	*  dtimestep - величина текущего шага по времени.
	*  
	*
	*  mfoldtimestep - массовый поток через каждую грань КО с предыдущего временного слоя. (неизменяемая в данной функции величина).
	*  mfcurrentretune - возвращаемый текущий массовый поток через каждую грань КО.
	*  КО имеет 6 граней: E,N,T,W,S,B.
	*
	* speedoldtimestep[VX or VY or VZ][iP: 0..maxelm+maxbound-1] - поле скорости с предыдущего временного слоя.
	*
	*/

	bool bhighorder=false;

	doublereal eps=admission; // для отделения вещественного нуля.
	// отключает или включает поправку Рхи-Чоу 1983г.
	bool bRhieChowi=true, bRhieChowb=false; // i - internal, b - border

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
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
		// если bborder   то мы находимся вблизи граничного узла.
		bool bborder=false;
		doublereal myflux=0.0;
		myflux=De*(dxe*DFDXiP(potent[PAM], iP, E_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iE]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
		myflux=Dw*(-dxw*DFDXiP(potent[PAM], iP, W_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iW]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dn*(dyn*DFDXiP(potent[PAM], iP, N_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iN]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Ds*(-dys*DFDXiP(potent[PAM], iP, S_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iS]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dt*(dzt*DFDXiP(potent[PAM], iP, T_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iT]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Db*(-dzb*DFDXiP(potent[PAM], iP, B_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iB]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
	}



    doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 
	     
	Fe=calcFg2(bE, feplus, iP, E_SIDE, rhoe, 
			dy*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, neighbors_for_the_internal_node, 
			bRhieChowi, bRhieChowb, false,
			tau, dxe, potent);

	Fw=calcFg2(bW, fwplus, iP, W_SIDE, rhow, 
			dy*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, neighbors_for_the_internal_node, 
			bRhieChowi, bRhieChowb, false,
			tau, dxw, potent);

	Fn=calcFg2(bN, fnplus, iP, N_SIDE, rhon, 
			dx*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, neighbors_for_the_internal_node, 
			bRhieChowi, bRhieChowb, false,
			tau, dyn, potent);

	Fs=calcFg2(bS, fsplus, iP, S_SIDE, rhos, 
			dx*dz, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, neighbors_for_the_internal_node, 
			bRhieChowi, bRhieChowb, false,
			tau, dys, potent);

	Ft=calcFg2(bT, ftplus, iP, T_SIDE, rhot, 
			dx*dy, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, neighbors_for_the_internal_node, 
			bRhieChowi, bRhieChowb, false,
			tau, dzt, potent);

	Fb=calcFg2(bB, fbplus, iP, B_SIDE, rhob, 
			dx*dy, btimedepend, 
			speedoldtimestep, mfoldtimestep,
			dtimestep,	RCh, neighbors_for_the_internal_node, 
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
			SpeedCorOldt = ftplus*potent[VZCOR][iT] + (1.0 - ftplus)*potent[VZCOR][iP];
		}
		else {
			SpeedCorOldt = potent[VZCOR][iT];
		}
		if (!bW) {
			SpeedCorOldw = fwplus*potent[VXCOR][iW] + (1.0 - fwplus)*potent[VXCOR][iP];
		}
		else {
			SpeedCorOldw = potent[VXCOR][iW];
		}
		if (!bS) {
			SpeedCorOlds = fsplus*potent[VYCOR][iS] + (1.0 - fsplus)*potent[VYCOR][iP];
		}
		else {
			SpeedCorOlds = potent[VYCOR][iS];
		}
		if (!bB) {
			SpeedCorOldb = fbplus*potent[VZCOR][iB] + (1.0 - fbplus)*potent[VZCOR][iP];
		}
		else {
			SpeedCorOldb = potent[VZCOR][iB];
		}

		// возвращаем значение потока на грани КО.
		// С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
		// дополнительная поправка осуществляется на основе скорректированной скорости.
		// Вообще говоря поле плотности также должно быть с предыдущей итерации.
		Fe += (1.0 - alpha[VELOCITY_X_COMPONENT])*(mfcurrentretune[E_SIDE] - rhoe*SpeedCorOlde*dy*dz);
		Fn += (1.0 - alpha[VELOCITY_Y_COMPONENT])*(mfcurrentretune[N_SIDE] - rhon*SpeedCorOldn*dx*dz);
		Ft += (1.0 - alpha[VELOCITY_Z_COMPONENT])*(mfcurrentretune[T_SIDE] - rhot*SpeedCorOldt*dx*dy);
		Fw += (1.0 - alpha[VELOCITY_X_COMPONENT])*(mfcurrentretune[W_SIDE] - rhow*SpeedCorOldw*dy*dz);
		Fs += (1.0 - alpha[VELOCITY_Y_COMPONENT])*(mfcurrentretune[S_SIDE] - rhos*SpeedCorOlds*dx*dz);
		Fb += (1.0 - alpha[VELOCITY_Z_COMPONENT])*(mfcurrentretune[B_SIDE] - rhob*SpeedCorOldb*dx*dy);

	}



	// возвращаем значение потока на грани КО.
	mfcurrentretune[E_SIDE] = Fe;
	mfcurrentretune[N_SIDE] = Fn;
	mfcurrentretune[T_SIDE] = Ft;
	mfcurrentretune[W_SIDE] = Fw;
	mfcurrentretune[S_SIDE] = Fs;
	mfcurrentretune[B_SIDE] = Fb;



	sl[PAM][iP].b = Fw - Fe + Fs - Fn + Fb - Ft + baddDFLUX2;



	// Симметризация СЛАУ:

	// для поправки давления получается эллиптическое уравнение с SPD матрицей.

	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

		// 1. Учёт условия Дирихле:
	if ((iE >= maxelm) && (fabs(slb[PAM][iE - maxelm].ai) < eps)) {
		sl[PAM][iP].b += sl[PAM][iP].ae*slb[PAM][iE - maxelm].b / slb[PAM][iE - maxelm].aw;
		sl[PAM][iP].ae = 0.0;
		sl[PAM][iP].iE = -1; // не входит в матрицу СЛАУ.
	}
	if ((iW >= maxelm) && (fabs(slb[PAM][iW - maxelm].ai) < eps)) {
		sl[PAM][iP].b += sl[PAM][iP].aw*slb[PAM][iW - maxelm].b / slb[PAM][iW - maxelm].aw;
		sl[PAM][iP].aw = 0.0;
		sl[PAM][iP].iW = -1; // не входит в матрицу СЛАУ.
	}
	if ((iN >= maxelm) && (fabs(slb[PAM][iN - maxelm].ai) < eps)) {
		sl[PAM][iP].b += sl[PAM][iP].an*slb[PAM][iN - maxelm].b / slb[PAM][iN - maxelm].aw;
		sl[PAM][iP].an = 0.0;
		sl[PAM][iP].iN = -1; // не входит в матрицу СЛАУ.
	}
	if ((iS >= maxelm) && (fabs(slb[PAM][iS - maxelm].ai) < eps)) {
		sl[PAM][iP].b += sl[PAM][iP].as*slb[PAM][iS - maxelm].b / slb[PAM][iS - maxelm].aw;
		sl[PAM][iP].as = 0.0;
		sl[PAM][iP].iS = -1; // не входит в матрицу СЛАУ.
	}
	if ((iT >= maxelm) && (fabs(slb[PAM][iT - maxelm].ai) < eps)) {
		sl[PAM][iP].b += sl[PAM][iP].at*slb[PAM][iT - maxelm].b / slb[PAM][iT - maxelm].aw;
		sl[PAM][iP].at = 0.0;
		sl[PAM][iP].iT = -1; // не входит в матрицу СЛАУ.
	}
	if ((iB >= maxelm) && (fabs(slb[PAM][iB - maxelm].ai) < eps)) {
		sl[PAM][iP].b += sl[PAM][iP].ab*slb[PAM][iB - maxelm].b / slb[PAM][iB - maxelm].aw;
		sl[PAM][iP].ab = 0.0;
		sl[PAM][iP].iB = -1; // не входит в матрицу СЛАУ.
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
						doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
						int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm,
						doublereal* alpha, doublereal dbeta, doublereal** &rhie_chow, doublereal RCh,
						bool btimedepend, doublereal dtimestep, doublereal* mfoldtimestep,
						doublereal* &mfcurrentretune, doublereal** speedoldtimestep, doublereal** tau,
						bool bmyhighorder, bool bdeltapfinish, bool bRhieChowi, bool bRhieChowb,
	                    BOUND* border_neighbor, integer *ilevel_alice, int* ptr, integer maxbound,
	                    doublereal** sumanb) {

	if (iP > maxelm + maxbound) {
		printf("my_elmatr_quad_PAm3 iP=%lld maxelm=%lld maxbound=%lld\n", iP, maxelm, maxbound);
		system("pause");
	}

	// точность решения уравнения для поправки давления неудовлетворительная, поэтому будем решать несколько раз
	// каждый последующий раз увеличивая точность. 26 июня 2012 года.
	// bdeltapfinish==true. - финишное решение.

    // 31 марта 2012 из этой функции удалено много закомментированного кода. Его можно найти в backup`ах AliceFlowv0_07.

	// tau - заранее вычисленное псевдовремя, оно служит коэффициентом диффузии в эллиптическом операторе, к тому-же оно
	// также входит в монотонизирующую поправку Рхи-Чоу (1983 года).

    /* btimedepend - стационарный (false), нестационарный (true).
	*  dtimestep - величина текущего шага по времени.
	*  
	*
	*  mfoldtimestep - массовый поток через каждую грань КО с предыдущего временного слоя. (неизменяемая в данной функции величина).
	*  mfcurrentretune - возвращаемый текущий массовый поток через каждую грань КО.
	*  КО имеет 6 граней: E,N,T,W,S,B.
	*
	* speedoldtimestep[VX or VY or VZ][iP: 0..maxelm+maxbound-1] - поле скорости с предыдущего временного слоя.
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
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
	sl[PAM][iP].iP=iP;
	sl[PAM][iP].iE=iE; sl[PAM][iP].iN=iN; 
	sl[PAM][iP].iS=iS; sl[PAM][iP].iW=iW;
    sl[PAM][iP].iT=iT; sl[PAM][iP].iB=iB;


	// 26.09.2016 Добавок для АЛИС сетки.
	integer iE2 = -1, iN2 = -1, iT2 = -1, iW2 = -1, iS2 = -1, iB2 = -1; // номера соседних контрольных объёмов
	integer iE3 = -1, iN3 = -1, iT3 = -1, iW3 = -1, iS3 = -1, iB3 = -1; // номера соседних контрольных объёмов
	integer iE4 = -1, iN4 = -1, iT4 = -1, iW4 = -1, iS4 = -1, iB4 = -1; // номера соседних контрольных объёмов

										  // -1 если не используется и [0..maxelm+maxbound-1] если используется.
	if (b_on_adaptive_local_refinement_mesh) {
		iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
		iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];
		iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
		iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];
		iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
		iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
	}
	/*
	if (iE2 > -1 || iW2 > -1 || iN2 > -1 || iS2 > -1 || iT2 > -1 || iB2 > -1 || iE3 > -1 || iW3 > -1 || iN3 > -1 || iS3 > -1 || iT3 > -1 || iB3 > -1 || iE4 > -1 || iW4 > -1 || iN4 > -1 || iS4 > -1 || iT4 > -1 || iB4 > -1) {
		printf("first popal\n");
		getchar();
	}
	else {
		printf("incorrect\n");
	}
	*/

	sl[PAM][iP].iE2 = iE2; sl[PAM][iP].iN2 = iN2; sl[PAM][iP].iT2 = iT2;
	sl[PAM][iP].iS2 = iS2; sl[PAM][iP].iW2 = iW2; sl[PAM][iP].iB2 = iB2;

	sl[PAM][iP].iE3 = iE3; sl[PAM][iP].iN3 = iN3; sl[PAM][iP].iT3 = iT3;
	sl[PAM][iP].iS3 = iS3; sl[PAM][iP].iW3 = iW3; sl[PAM][iP].iB3 = iB3;

	sl[PAM][iP].iE4 = iE4; sl[PAM][iP].iN4 = iN4; sl[PAM][iP].iT4 = iT4;
	sl[PAM][iP].iS4 = iS4; sl[PAM][iP].iW4 = iW4; sl[PAM][iP].iB4 = iB4;

	// Инициализирующее обнуление.
	sl[PAM][iP].ae = 0.0;
	sl[PAM][iP].aw = 0.0;
	sl[PAM][iP].an = 0.0;
	sl[PAM][iP].as = 0.0;
	sl[PAM][iP].at = 0.0;
	sl[PAM][iP].ab = 0.0;

	sl[PAM][iP].ae2 = 0.0;
	sl[PAM][iP].aw2 = 0.0;
	sl[PAM][iP].an2 = 0.0;
	sl[PAM][iP].as2 = 0.0;
	sl[PAM][iP].at2 = 0.0;
	sl[PAM][iP].ab2 = 0.0;

	sl[PAM][iP].ae3 = 0.0;
	sl[PAM][iP].aw3 = 0.0;
	sl[PAM][iP].an3 = 0.0;
	sl[PAM][iP].as3 = 0.0;
	sl[PAM][iP].at3 = 0.0;
	sl[PAM][iP].ab3 = 0.0;

	sl[PAM][iP].ae4 = 0.0;
	sl[PAM][iP].aw4 = 0.0;
	sl[PAM][iP].an4 = 0.0;
	sl[PAM][iP].as4 = 0.0;
	sl[PAM][iP].at4 = 0.0;
	sl[PAM][iP].ab4 = 0.0;
	
	// Признак присутствия связи.
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
	doublereal feplus3, fwplus3, fnplus3, fsplus3, ftplus3, fbplus3;
	doublereal feplus4, fwplus4, fnplus4, fsplus4, ftplus4, fbplus4;

	if (b_on_adaptive_local_refinement_mesh) {

		// x-direction
		feplus2 = 0.5 * dx / dxe2;
		fwplus2 = 0.5 * dx / dxw2;
		// y-direction
		fnplus2 = 0.5 * dy / dyn2;
		fsplus2 = 0.5 * dy / dys2;
		// z-direction
		ftplus2 = 0.5 * dz / dzt2;
		fbplus2 = 0.5 * dz / dzb2;
		
		// x-direction
		feplus3 = 0.5 * dx / dxe3;
		fwplus3 = 0.5 * dx / dxw3;
		// y-direction
		fnplus3 = 0.5 * dy / dyn3;
		fsplus3 = 0.5 * dy / dys3;
		// z-direction
		ftplus3 = 0.5 * dz / dzt3;
		fbplus3 = 0.5 * dz / dzb3;
		
		// x-direction
		feplus4 = 0.5 * dx / dxe4;
		fwplus4 = 0.5 * dx / dxw4;
		// y-direction
		fnplus4 = 0.5 * dy / dyn4;
		fsplus4 = 0.5 * dy / dys4;
		// z-direction
		ftplus4 = 0.5 * dz / dzt4;
		fbplus4 = 0.5 * dz / dzb4;
	}
	// Значение коэффициента диффузии (псевдовремени) на грани контрольного объёма.
	// Учитывается что псевдовремя это три скалярных поля.
	doublereal taue=1.0, tauw=1.0, taun=1.0, taus=1.0, taut=1.0, taub=1.0;
	if (iE > -1) {
		if (!bE) taue = tau[VELOCITY_X_COMPONENT][iE] * tau[VELOCITY_X_COMPONENT][iP] / (feplus*tau[VELOCITY_X_COMPONENT][iE] + (1.0 - feplus)*tau[VELOCITY_X_COMPONENT][iP]); else taue = tau[VELOCITY_X_COMPONENT][iE]; // проверено !
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

	doublereal taue2 = 1.0, tauw2 = 1.0, taun2 = 1.0, taus2 = 1.0, taut2 = 1.0, taub2 = 1.0;
	doublereal taue3 = 1.0, tauw3 = 1.0, taun3 = 1.0, taus3 = 1.0, taut3 = 1.0, taub3 = 1.0;
	doublereal taue4 = 1.0, tauw4 = 1.0, taun4 = 1.0, taus4 = 1.0, taut4 = 1.0, taub4 = 1.0;

	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 > -1) {
			if (!bE2) taue2 = tau[VELOCITY_X_COMPONENT][iE2] * tau[VELOCITY_X_COMPONENT][iP] / (feplus2 * tau[VELOCITY_X_COMPONENT][iE2] + (1.0 - feplus2) * tau[VELOCITY_X_COMPONENT][iP]); else taue2 = tau[VELOCITY_X_COMPONENT][iE2]; // проверено !
		}
		if (iW2 > -1) {
			if (!bW2) tauw2 = tau[VELOCITY_X_COMPONENT][iW2] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus2 * tau[VELOCITY_X_COMPONENT][iW2] + (1.0 - fwplus2) * tau[VELOCITY_X_COMPONENT][iP]); else tauw2 = tau[VELOCITY_X_COMPONENT][iW2];
		}
		if (iN2 > -1) {
			if (!bN2) taun2 = tau[VELOCITY_Y_COMPONENT][iN2] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus2 * tau[VELOCITY_Y_COMPONENT][iN2] + (1.0 - fnplus2) * tau[VELOCITY_Y_COMPONENT][iP]); else taun2 = tau[VELOCITY_Y_COMPONENT][iN2];
		}
		if (iS2 > -1) {
			if (!bS2) taus2 = tau[VELOCITY_Y_COMPONENT][iS2] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus2 * tau[VELOCITY_Y_COMPONENT][iS2] + (1.0 - fsplus2) * tau[VELOCITY_Y_COMPONENT][iP]); else taus2 = tau[VELOCITY_Y_COMPONENT][iS2];
		}
		if (iT2 > -1) {
			if (!bT2) taut2 = tau[VELOCITY_Z_COMPONENT][iT2] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus2 * tau[VELOCITY_Z_COMPONENT][iT2] + (1.0 - ftplus2) * tau[VELOCITY_Z_COMPONENT][iP]); else taut2 = tau[VELOCITY_Z_COMPONENT][iT2];
		}
		if (iB2 > -1) {
			if (!bB2) taub2 = tau[VELOCITY_Z_COMPONENT][iB2] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus2 * tau[VELOCITY_Z_COMPONENT][iB2] + (1.0 - fbplus2) * tau[VELOCITY_Z_COMPONENT][iP]); else taub2 = tau[VELOCITY_Z_COMPONENT][iB2];
		}

		
		if (iE3 > -1) {
			if (!bE3) taue3 = tau[VELOCITY_X_COMPONENT][iE3] * tau[VELOCITY_X_COMPONENT][iP] / (feplus3 * tau[VELOCITY_X_COMPONENT][iE3] + (1.0 - feplus3) * tau[VELOCITY_X_COMPONENT][iP]); else taue3 = tau[VELOCITY_X_COMPONENT][iE3]; // проверено !
		}
		if (iW3 > -1) {
			if (!bW3) tauw3 = tau[VELOCITY_X_COMPONENT][iW3] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus3 * tau[VELOCITY_X_COMPONENT][iW3] + (1.0 - fwplus3) * tau[VELOCITY_X_COMPONENT][iP]); else tauw3 = tau[VELOCITY_X_COMPONENT][iW3];
		}
		if (iN3 > -1) {
			if (!bN3) taun3 = tau[VELOCITY_Y_COMPONENT][iN3] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus3 * tau[VELOCITY_Y_COMPONENT][iN3] + (1.0 - fnplus3) * tau[VELOCITY_Y_COMPONENT][iP]); else taun3 = tau[VELOCITY_Y_COMPONENT][iN3];
		}
		if (iS3 > -1) {
			if (!bS3) taus3 = tau[VELOCITY_Y_COMPONENT][iS3] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus3 * tau[VELOCITY_Y_COMPONENT][iS3] + (1.0 - fsplus3) * tau[VELOCITY_Y_COMPONENT][iP]); else taus3 = tau[VELOCITY_Y_COMPONENT][iS3];
		}
		if (iT3 > -1) {
			if (!bT3) taut3 = tau[VELOCITY_Z_COMPONENT][iT3] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus3 * tau[VELOCITY_Z_COMPONENT][iT3] + (1.0 - ftplus3) * tau[VELOCITY_Z_COMPONENT][iP]); else taut3 = tau[VELOCITY_Z_COMPONENT][iT3];
		}
		if (iB3 > -1) {
			if (!bB3) taub3 = tau[VELOCITY_Z_COMPONENT][iB3] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus3 * tau[VELOCITY_Z_COMPONENT][iB3] + (1.0 - fbplus3) * tau[VELOCITY_Z_COMPONENT][iP]); else taub3 = tau[VELOCITY_Z_COMPONENT][iB3];
		}

		
		if (iE4 > -1) {
			if (!bE4) taue4 = tau[VELOCITY_X_COMPONENT][iE4] * tau[VELOCITY_X_COMPONENT][iP] / (feplus4 * tau[VELOCITY_X_COMPONENT][iE4] + (1.0 - feplus4) * tau[VELOCITY_X_COMPONENT][iP]); else taue4 = tau[VELOCITY_X_COMPONENT][iE4]; // проверено !
		}
		if (iW4 > -1) {
			if (!bW4) tauw4 = tau[VELOCITY_X_COMPONENT][iW4] * tau[VELOCITY_X_COMPONENT][iP] / (fwplus4 * tau[VELOCITY_X_COMPONENT][iW4] + (1.0 - fwplus4) * tau[VELOCITY_X_COMPONENT][iP]); else tauw4 = tau[VELOCITY_X_COMPONENT][iW4];
		}
		if (iN4 > -1) {
			if (!bN4) taun4 = tau[VELOCITY_Y_COMPONENT][iN4] * tau[VELOCITY_Y_COMPONENT][iP] / (fnplus4 * tau[VELOCITY_Y_COMPONENT][iN4] + (1.0 - fnplus4) * tau[VELOCITY_Y_COMPONENT][iP]); else taun4 = tau[VELOCITY_Y_COMPONENT][iN4];
		}
		if (iS4 > -1) {
			if (!bS4) taus4 = tau[VELOCITY_Y_COMPONENT][iS4] * tau[VELOCITY_Y_COMPONENT][iP] / (fsplus4 * tau[VELOCITY_Y_COMPONENT][iS4] + (1.0 - fsplus4) * tau[VELOCITY_Y_COMPONENT][iP]); else taus4 = tau[VELOCITY_Y_COMPONENT][iS4];
		}
		if (iT4 > -1) {
			if (!bT4) taut4 = tau[VELOCITY_Z_COMPONENT][iT4] * tau[VELOCITY_Z_COMPONENT][iP] / (ftplus4 * tau[VELOCITY_Z_COMPONENT][iT4] + (1.0 - ftplus4) * tau[VELOCITY_Z_COMPONENT][iP]); else taut4 = tau[VELOCITY_Z_COMPONENT][iT4];
		}
		if (iB4 > -1) {
			if (!bB4) taub4 = tau[VELOCITY_Z_COMPONENT][iB4] * tau[VELOCITY_Z_COMPONENT][iP] / (fbplus4 * tau[VELOCITY_Z_COMPONENT][iB4] + (1.0 - fbplus4) * tau[VELOCITY_Z_COMPONENT][iP]); else taub4 = tau[VELOCITY_Z_COMPONENT][iB4];
		}
	}

	/*
	doublereal ae, aw, as, an, at, ab;
	if (iE > -1) {
		if (!bE) ae = sumanb[VX][iE] * sumanb[VX][iP] / (feplus * sumanb[VX][iE] + (1.0 - feplus) * sumanb[VX][iP]); else ae = sumanb[VX][iE]; // проверено !
	}
	if (iW > -1) {
		if (!bW) aw = sumanb[VX][iW] * sumanb[VX][iP] / (fwplus * sumanb[VX][iW] + (1.0 - fwplus) * sumanb[VX][iP]); else aw = sumanb[VX][iW]; // проверено !
	}
	if (iN > -1) {
		if (!bN) an = sumanb[VY][iN] * sumanb[VY][iP] / (fnplus * sumanb[VY][iN] + (1.0 - fnplus) * sumanb[VY][iP]); else an = sumanb[VY][iN]; // проверено !
	}
	if (iS > -1) {
		if (!bS) as = sumanb[VY][iS] * sumanb[VY][iP] / (fsplus * sumanb[VY][iS] + (1.0 - fsplus) * sumanb[VY][iP]); else as = sumanb[VY][iS]; // проверено !
	}
	if (iT > -1) {
		if (!bT) at = sumanb[VZ][iT] * sumanb[VZ][iP] / (ftplus * sumanb[VZ][iT] + (1.0 - ftplus) * sumanb[VZ][iP]); else at = sumanb[VZ][iT]; // проверено !
	}
	if (iB > -1) {
		if (!bB) ab = sumanb[VZ][iB] * sumanb[VZ][iP] / (fbplus * sumanb[VZ][iB] + (1.0 - fbplus) * sumanb[VZ][iP]); else ab = sumanb[VZ][iB]; // проверено !
	}

	doublereal ae2, aw2, as2, an2, at2, ab2;
	if (iE2 > -1) {
		if (!bE2) ae2 = sumanb[VX][iE2] * sumanb[VX][iP] / (feplus2 * sumanb[VX][iE2] + (1.0 - feplus2) * sumanb[VX][iP]); else ae2 = sumanb[VX][iE2]; // проверено !
	}
	if (iW2 > -1) {
		if (!bW2) aw2 = sumanb[VX][iW2] * sumanb[VX][iP] / (fwplus2 * sumanb[VX][iW2] + (1.0 - fwplus2) * sumanb[VX][iP]); else aw2 = sumanb[VX][iW2]; // проверено !
	}
	if (iN2 > -1) {
		if (!bN2) an2 = sumanb[VY][iN2] * sumanb[VY][iP] / (fnplus2 * sumanb[VY][iN2] + (1.0 - fnplus2) * sumanb[VY][iP]); else an2 = sumanb[VY][iN2]; // проверено !
	}
	if (iS2 > -1) {
		if (!bS2) as2 = sumanb[VY][iS2] * sumanb[VY][iP] / (fsplus2 * sumanb[VY][iS2] + (1.0 - fsplus2) * sumanb[VY][iP]); else as2 = sumanb[VY][iS2]; // проверено !
	}
	if (iT2 > -1) {
		if (!bT2) at2 = sumanb[VZ][iT2] * sumanb[VZ][iP] / (ftplus2 * sumanb[VZ][iT2] + (1.0 - ftplus2) * sumanb[VZ][iP]); else at2 = sumanb[VZ][iT2]; // проверено !
	}
	if (iB2 > -1) {
		if (!bB2) ab2 = sumanb[VZ][iB2] * sumanb[VZ][iP] / (fbplus2 * sumanb[VZ][iB2] + (1.0 - fbplus2) * sumanb[VZ][iP]); else ab2 = sumanb[VZ][iB2]; // проверено !
	}

	doublereal ae3, aw3, as3, an3, at3, ab3;
	if (iE3 > -1) {
		if (!bE3) ae3 = sumanb[VX][iE3] * sumanb[VX][iP] / (feplus3 * sumanb[VX][iE3] + (1.0 - feplus3) * sumanb[VX][iP]); else ae3 = sumanb[VX][iE3]; // проверено !
	}
	if (iW3 > -1) {
		if (!bW3) aw3 = sumanb[VX][iW3] * sumanb[VX][iP] / (fwplus3 * sumanb[VX][iW3] + (1.0 - fwplus3) * sumanb[VX][iP]); else aw3 = sumanb[VX][iW3]; // проверено !
	}
	if (iN3 > -1) {
		if (!bN3) an3 = sumanb[VY][iN3] * sumanb[VY][iP] / (fnplus3 * sumanb[VY][iN3] + (1.0 - fnplus3) * sumanb[VY][iP]); else an3 = sumanb[VY][iN3]; // проверено !
	}
	if (iS3 > -1) {
		if (!bS3) as3 = sumanb[VY][iS3] * sumanb[VY][iP] / (fsplus3 * sumanb[VY][iS3] + (1.0 - fsplus3) * sumanb[VY][iP]); else as3 = sumanb[VY][iS3]; // проверено !
	}
	if (iT3 > -1) {
		if (!bT3) at3 = sumanb[VZ][iT3] * sumanb[VZ][iP] / (ftplus3 * sumanb[VZ][iT3] + (1.0 - ftplus3) * sumanb[VZ][iP]); else at3 = sumanb[VZ][iT3]; // проверено !
	}
	if (iB3 > -1) {
		if (!bB3) ab3 = sumanb[VZ][iB3] * sumanb[VZ][iP] / (fbplus3 * sumanb[VZ][iB3] + (1.0 - fbplus3) * sumanb[VZ][iP]); else ab3 = sumanb[VZ][iB3]; // проверено !
	}

	doublereal ae4, aw4, as4, an4, at4, ab4;
	if (iE4 > -1) {
		if (!bE4) ae4 = sumanb[VX][iE4] * sumanb[VX][iP] / (feplus4 * sumanb[VX][iE4] + (1.0 - feplus4) * sumanb[VX][iP]); else ae4 = sumanb[VX][iE4]; // проверено !
	}
	if (iW4 > -1) {
		if (!bW4) aw4 = sumanb[VX][iW4] * sumanb[VX][iP] / (fwplus4 * sumanb[VX][iW4] + (1.0 - fwplus4) * sumanb[VX][iP]); else aw4 = sumanb[VX][iW4]; // проверено !
	}
	if (iN4 > -1) {
		if (!bN4) an4 = sumanb[VY][iN4] * sumanb[VY][iP] / (fnplus4 * sumanb[VY][iN4] + (1.0 - fnplus4) * sumanb[VY][iP]); else an4 = sumanb[VY][iN4]; // проверено !
	}
	if (iS4 > -1) {
		if (!bS4) as4 = sumanb[VY][iS4] * sumanb[VY][iP] / (fsplus4 * sumanb[VY][iS4] + (1.0 - fsplus4) * sumanb[VY][iP]); else as4 = sumanb[VY][iS4]; // проверено !
	}
	if (iT4 > -1) {
		if (!bT4) at4 = sumanb[VZ][iT4] * sumanb[VZ][iP] / (ftplus4 * sumanb[VZ][iT4] + (1.0 - ftplus4) * sumanb[VZ][iP]); else at4 = sumanb[VZ][iT4]; // проверено !
	}
	if (iB4 > -1) {
		if (!bB4) ab4 = sumanb[VZ][iB4] * sumanb[VZ][iP] / (fbplus4 * sumanb[VZ][iB4] + (1.0 - fbplus4) * sumanb[VZ][iP]); else ab4 = sumanb[VZ][iB4]; // проверено !
	}
	*/

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
	doublereal  rE3 = 0.0, rN3 = 0.0, rT3 = 0.0, rW3 = 0.0, rS3 = 0.0, rB3 = 0.0;
	doublereal  rE4 = 0.0, rN4 = 0.0, rT4 = 0.0, rW4 = 0.0, rS4 = 0.0, rB4 = 0.0;

	if (b_on_adaptive_local_refinement_mesh) {

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

	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 > -1) {
			if (!bE2)  rhoe2 = rE2 * rP / (feplus2 * rE2 + (1.0 - feplus2) * rP); else rhoe2 = rE2; // проверено !
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

	}

	doublereal dSqe = 0.0, dSqw = 0.0, dSqn = 0.0, dSqs = 0.0, dSqt = 0.0, dSqb = 0.0; // площадь грани.

	if (iE > -1) {

		dSqe = dy * dz;

		if (bE) {
			// граничный узел.
			dSqe = border_neighbor[iE - maxelm].dS;
			//printf("bon square %e\n", dSqe);// все верно
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

	}

	//dSqe = dy * dz; dSqw = dy * dz; dSqn = dx * dz; dSqs = dx * dz; dSqt = dx * dy; dSqb = dx * dy; // площадь грани.
	//dSqe2 = dy * dz; dSqw2 = dy * dz; dSqn2 = dx * dz; dSqs2 = dx * dz; dSqt2 = dx * dy; dSqb2 = dx * dy; // площадь грани.
	//dSqe3 = dy * dz; dSqw3 = dy * dz; dSqn3 = dx * dz; dSqs3 = dx * dz; dSqt3 = dx * dy; dSqb3 = dx * dy; // площадь грани.
	//dSqe4 = dy * dz; dSqw4 = dy * dz; dSqn4 = dx * dz; dSqs4 = dx * dz; dSqt4 = dx * dy; dSqb4 = dx * dy; // площадь грани.

    // 0.8 - SIMPLEC параметр релаксации для скоростей (см Гаврилов Андрей)
	// уравнение поправки давления.

	const doublereal ZeroDiffusion = 0.0;// 1.0e-30;
	// Диффузионная составляющая потока:
	doublereal De = ZeroDiffusion, Dw = ZeroDiffusion, Dn = ZeroDiffusion, Ds = ZeroDiffusion, Dt = ZeroDiffusion, Db = ZeroDiffusion; // инициализация
	doublereal De2 = ZeroDiffusion, Dw2 = ZeroDiffusion, Dn2 = ZeroDiffusion, Ds2 = ZeroDiffusion, Dt2 = ZeroDiffusion, Db2 = ZeroDiffusion; // инициализация
	doublereal De3 = ZeroDiffusion, Dw3 = ZeroDiffusion, Dn3 = ZeroDiffusion, Ds3 = ZeroDiffusion, Dt3 = ZeroDiffusion, Db3 = ZeroDiffusion; // инициализация
	doublereal De4 = ZeroDiffusion, Dw4 = ZeroDiffusion, Dn4 = ZeroDiffusion, Ds4 = ZeroDiffusion, Dt4 = ZeroDiffusion, Db4 = ZeroDiffusion; // инициализация
		
	// Псевдовремя tau является коэффициентом диффузии 
	// в эллиптическом уравнении на поправку давления.

	if (iE > -1) {
		if (!bE) {
			if (bW) De = dbeta * taue*dSqe / dxe;
			else De = taue * dSqe / dxe;
		}
		else De = dbeta * taue*dSqe / dxe;
	}

	if (iW > -1) {
		if (!bW) {
			if (bE) Dw = dbeta * tauw*dSqw / dxw;
			else Dw = tauw * dSqw / dxw;
		}
		else Dw = dbeta * tauw*dSqw / dxw;
	}

	if (iN > -1) {
		if (!bN) {
			if (bS) Dn = dbeta * taun*dSqn / dyn;
			else Dn = taun * dSqn / dyn;
		}
		else Dn = dbeta * taun*dSqn / dyn;
	}
	
	if (iS > -1) {
		if (!bS) {
			if (bN) Ds = dbeta * taus*dSqs / dys;
			else Ds = taus * dSqs / dys;
		}
		else Ds = dbeta * taus*dSqs / dys;
	}

	if (iT > -1) {
		if (!bT) {
			if (bB) Dt = dbeta * taut*dSqt / dzt;
			else Dt = taut * dSqt / dzt;
		}
		else Dt = dbeta * taut*dSqt / dzt;
	}

	if (iB > -1) {
		if (!bB) {
			if (bT) Db = dbeta * taub*dSqb / dzb;
			else Db = taub * dSqb / dzb;
		}
		else Db = dbeta * taub*dSqb / dzb;
	}
	
	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 > -1) {
			if (!bE2) {
				if (bW2) De2 = dbeta * taue2 * dSqe2 / dxe2;
				else De2 = taue2 * dSqe2 / dxe2;
			}
			else De2 = dbeta * taue2 * dSqe2 / dxe2;
		}

		if (iW2 > -1) {
			if (!bW2) {
				if (bE2) Dw2 = dbeta * tauw2 * dSqw2 / dxw2;
				else Dw2 = tauw2 * dSqw2 / dxw2;
			}
			else Dw2 = dbeta * tauw2 * dSqw2 / dxw2;
		}

		if (iN2 > -1) {
			if (!bN2) {
				if (bS2) Dn2 = dbeta * taun2 * dSqn2 / dyn2;
				else Dn2 = taun2 * dSqn2 / dyn2;
			}
			else Dn2 = dbeta * taun2 * dSqn2 / dyn2;
		}

		if (iS2 > -1) {
			if (!bS2) {
				if (bN2) Ds2 = dbeta * taus2 * dSqs2 / dys2;
				else Ds2 = taus2 * dSqs2 / dys2;
			}
			else Ds2 = dbeta * taus2 * dSqs2 / dys2;
		}

		if (iT2 > -1) {
			if (!bT2) {
				if (bB2) Dt2 = dbeta * taut2 * dSqt2 / dzt2;
				else Dt2 = taut2 * dSqt2 / dzt2;
			}
			else Dt2 = dbeta * taut2 * dSqt2 / dzt2;
		}

		if (iB2 > -1) {
			if (!bB2) {
				if (bT2) Db2 = dbeta * taub2 * dSqb2 / dzb2;
				else Db2 = taub2 * dSqb2 / dzb2;
			}
			else Db2 = dbeta * taub2 * dSqb2 / dzb2;
		}

		if (iE3 > -1) {
			if (!bE3) {
				if (bW3) De3 = dbeta * taue3 * dSqe3 / dxe3;
				else De3 = taue3 * dSqe3 / dxe3;
			}
			else De3 = dbeta * taue3 * dSqe3 / dxe3;
		}

		if (iW3 > -1) {
			if (!bW3) {
				if (bE3) Dw3 = dbeta * tauw3 * dSqw3 / dxw3;
				else Dw3 = tauw3 * dSqw3 / dxw3;
			}
			else Dw3 = dbeta * tauw3 * dSqw3 / dxw3;
		}

		if (iN3 > -1) {
			if (!bN3) {
				if (bS3) Dn3 = dbeta * taun3 * dSqn3 / dyn3;
				else Dn3 = taun3 * dSqn3 / dyn3;
			}
			else Dn3 = dbeta * taun3 * dSqn3 / dyn3;
		}

		if (iS3 > -1) {
			if (!bS3) {
				if (bN3) Ds3 = dbeta * taus3 * dSqs3 / dys3;
				else Ds3 = taus3 * dSqs3 / dys3;
			}
			else Ds3 = dbeta * taus3 * dSqs3 / dys3;
		}

		if (iT3 > -1) {
			if (!bT3) {
				if (bB3) Dt = dbeta * taut3 * dSqt3 / dzt3;
				else Dt3 = taut3 * dSqt3 / dzt3;
			}
			else Dt3 = dbeta * taut3 * dSqt3 / dzt3;
		}

		if (iB3 > -1) {
			if (!bB3) {
				if (bT3) Db3 = dbeta * taub3 * dSqb3 / dzb3;
				else Db3 = taub3 * dSqb3 / dzb3;
			}
			else Db3 = dbeta * taub3 * dSqb3 / dzb3;
		}

		if (iE4 > -1) {
			if (!bE4) {
				if (bW4) De4 = dbeta * taue4 * dSqe4 / dxe4;
				else De4 = taue4 * dSqe4 / dxe4;
			}
			else De4 = dbeta * taue4 * dSqe4 / dxe4;
		}

		if (iW4 > -1) {
			if (!bW4) {
				if (bE4) Dw4 = dbeta * tauw4 * dSqw4 / dxw4;
				else Dw4 = tauw4 * dSqw4 / dxw4;
			}
			else Dw4 = dbeta * tauw4 * dSqw4 / dxw4;
		}

		if (iN4 > -1) {
			if (!bN4) {
				if (bS4) Dn4 = dbeta * taun4 * dSqn4 / dyn4;
				else Dn4 = taun4 * dSqn4 / dyn4;
			}
			else Dn4 = dbeta * taun4 * dSqn4 / dyn4;
		}

		if (iS4 > -1) {
			if (!bS4) {
				if (bN4) Ds4 = dbeta * taus4 * dSqs4 / dys4;
				else Ds4 = taus4 * dSqs4 / dys4;
			}
			else Ds4 = dbeta * taus4 * dSqs4 / dys4;
		}

		if (iT4 > -1) {
			if (!bT4) {
				if (bB4) Dt4 = dbeta * taut4 * dSqt4 / dzt4;
				else Dt4 = taut4 * dSqt4 / dzt4;
			}
			else Dt4 = dbeta * taut4 * dSqt4 / dzt4;
		}

		if (iB4 > -1) {
			if (!bB4) {
				if (bT4) Db4 = dbeta * taub4 * dSqb4 / dzb4;
				else Db4 = taub4 * dSqb4 / dzb4;
			}
			else Db4 = dbeta * taub4 * dSqb4 / dzb4;
		}
	}

/*
	// Псевдовремя tau является коэффициентом диффузии 
	// в эллиптическом уравнении на поправку давления.
if (iE > -1) {
	if (!bE) {
		if (bW) De = dbeta *rhoe* alpha[VX] * dSqe * dSqe / (ae);
		else De = rhoe * alpha[VX] * dSqe * dSqe / (ae);
	}
	else De = dbeta * rhoe * alpha[VX] * dSqe * dSqe / (ae);
}

if (iW > -1) {
	if (!bW) {
		if (bE) Dw = dbeta * rhow * alpha[VX] * dSqw * dSqw / (aw);
		else Dw = rhow * alpha[VX] * dSqw * dSqw / (aw);
	}
	else Dw = dbeta * rhow * alpha[VX] * dSqw * dSqw / (aw);
}

if (iN > -1) {
	if (!bN) {
		if (bS) Dn = dbeta * rhon * alpha[VY] * dSqn * dSqn / (an);
		else Dn = rhon * alpha[VY] * dSqn * dSqn / (an);
	}
	else Dn = dbeta * rhon * alpha[VY] * dSqn * dSqn / (an);
}

if (iS > -1) {
	if (!bS) {
		if (bN) Ds = dbeta * rhos * alpha[VY] * dSqs * dSqs / (as);
		else Ds = rhos * alpha[VY] * dSqs * dSqs / (as);
	}
	else Ds = dbeta * rhos * alpha[VY] * dSqs * dSqs / (as);
}

if (iT > -1) {
	if (!bT) {
		if (bB) Dt = dbeta * rhot * alpha[VZ] * dSqt * dSqt / (at);
		else Dt = rhot * alpha[VZ] * dSqt * dSqt / (at);
	}
	else Dt = dbeta * rhot * alpha[VZ] * dSqt * dSqt / (at);
}

if (iB > -1) {
	if (!bB) {
		if (bT) Db = dbeta * rhob * alpha[VZ] * dSqb * dSqb / (ab);
		else Db = taub * rhob * alpha[VZ] * dSqb * dSqb / (ab);
	}
	else Db = dbeta * rhob * alpha[VZ] * dSqb * dSqb / (ab);
}

if (iE2 > -1) {
	if (!bE2) {
		if (bW2) De2 = dbeta * rhoe2 * alpha[VX] * dSqe2 * dSqe2 / (ae2);
		else De2 = rhoe2 * alpha[VX] * dSqe2 * dSqe2 / (ae2);
	}
	else De2 = dbeta * rhoe2 * alpha[VX] * dSqe2 * dSqe2 / (ae2);
}

if (iW2 > -1) {
	if (!bW2) {
		if (bE2) Dw2 = dbeta * rhow2 * alpha[VX] * dSqw2 * dSqw2 / (aw2);
		else Dw2 = rhow2 * alpha[VX] * dSqw2 * dSqw2 / (aw2);
	}
	else Dw2 = dbeta * rhow2 * alpha[VX] * dSqw2 * dSqw2 / (aw2);
}

if (iN2 > -1) {
	if (!bN2) {
		if (bS2) Dn2 = dbeta * rhon2 * alpha[VY] * dSqn2 * dSqn2 / (an2);
		else Dn2 = rhon2 * alpha[VY] * dSqn2 * dSqn2 / (an2);
	}
	else Dn2 = dbeta * rhon2 * alpha[VY] * dSqn2 * dSqn2 / (an2);
}

if (iS2 > -1) {
	if (!bS2) {
		if (bN2) Ds2 = dbeta * rhos2 * alpha[VY] * dSqs2 * dSqs2 / (as2);
		else Ds2 = rhos2 * alpha[VY] * dSqs2 * dSqs2 / (as2);
	}
	else Ds2 = dbeta * rhos2 * alpha[VY] * dSqs2 * dSqs2 / (as2);
}

if (iT2 > -1) {
	if (!bT2) {
		if (bB2) Dt2 = dbeta * rhot2 * alpha[VZ] * dSqt2 * dSqt2 / (at2);
		else Dt2 = rhot2 * alpha[VZ] * dSqt2 * dSqt2 / (at2);
	}
	else Dt2 = dbeta * rhot2 * alpha[VZ] * dSqt2 * dSqt2 / (at2);
}

if (iB2 > -1) {
	if (!bB2) {
		if (bT2) Db2 = dbeta * rhob2 * alpha[VZ] * dSqb2 * dSqb2 / (ab2);
		else Db2 = rhob2 * alpha[VZ] * dSqb2 * dSqb2 / (ab2);
	}
	else Db2 = dbeta * rhob2 * alpha[VZ] * dSqb2 * dSqb2 / (ab2);
}

if (iE3 > -1) {
	if (!bE3) {
		if (bW3) De3 = dbeta * rhoe3 * alpha[VX] * dSqe3 * dSqe3 / (ae3);
		else De3 = rhoe3 * alpha[VX] * dSqe3 * dSqe3 / (ae3);
	}
	else De3 = dbeta * rhoe3 * alpha[VX] * dSqe3 * dSqe3 / (ae3);
}

if (iW3 > -1) {
	if (!bW3) {
		if (bE3) Dw3 = dbeta * rhow3 * alpha[VX] * dSqw3 * dSqw3 / (aw3);
		else Dw3 = rhow3 * alpha[VX] * dSqw3 * dSqw3 / (aw3);
	}
	else Dw3 = dbeta * rhow3 * alpha[VX] * dSqw3 * dSqw3 / (aw3);
}

if (iN3 > -1) {
	if (!bN3) {
		if (bS3) Dn3 = dbeta * rhon3 * alpha[VY] * dSqn3 * dSqn3 / (an3);
		else Dn3 = rhon3 * alpha[VY] * dSqn3 * dSqn3 / (an3);
	}
	else Dn3 = dbeta * rhon3 * alpha[VY] * dSqn3 * dSqn3 / (an3);
}

if (iS3 > -1) {
	if (!bS3) {
		if (bN3) Ds3 = dbeta * rhos3 * alpha[VY] * dSqs3 * dSqs3 / (as3);
		else Ds3 = rhos3 * alpha[VY] * dSqs3 * dSqs3 / (as3);
	}
	else Ds3 = dbeta * rhos3 * alpha[VY] * dSqs3 * dSqs3 / (as3);
}

if (iT3 > -1) {
	if (!bT3) {
		if (bB3) Dt = dbeta * rhot3 * alpha[VZ] * dSqt3 * dSqt3 / (at3);
		else Dt3 = rhot3 * alpha[VZ] * dSqt3 * dSqt3  / (at3);
	}
	else Dt3 = dbeta * rhot3 * alpha[VZ] * dSqt3 * dSqt3 / (at3);
}

if (iB3 > -1) {
	if (!bB3) {
		if (bT3) Db3 = dbeta * rhob3 * alpha[VZ] * dSqb3 * dSqb3 / (ab3);
		else Db3 = rhob3 * alpha[VZ] * dSqb3 * dSqb3 / (ab3);
	}
	else Db3 = dbeta * rhob3 * alpha[VZ] * dSqb3 * dSqb3 / (ab3);
}

if (iE4 > -1) {
	if (!bE4) {
		if (bW4) De4 = dbeta * rhoe4 * alpha[VX] * dSqe4 * dSqe4 / (ae4);
		else De4 = rhoe4 * alpha[VX] * dSqe4 * dSqe4 / (ae4);
	}
	else De4 = dbeta * rhoe4 * alpha[VX] * dSqe4 * dSqe4 / (ae4);
}

if (iW4 > -1) {
	if (!bW4) {
		if (bE4) Dw4 = dbeta * rhow4 * alpha[VX] * dSqw4 * dSqw4 / (aw4);
		else Dw4 = rhow4 * alpha[VX] * dSqw4 * dSqw4 / (aw4);
	}
	else Dw4 = dbeta * rhow4 * alpha[VX] * dSqw4 * dSqw4 / (aw4);
}

if (iN4 > -1) {
	if (!bN4) {
		if (bS4) Dn4 = dbeta * rhon4 * alpha[VY] * dSqn4 * dSqn4 / (an4);
		else Dn4 = rhon4 * alpha[VY] * dSqn4 * dSqn4 / (an4);
	}
	else Dn4 = dbeta * rhon4 * alpha[VY] * dSqn4 * dSqn4 / (an4);
}

if (iS4 > -1) {
	if (!bS4) {
		if (bN4) Ds4 = dbeta * rhos4 * alpha[VY] * dSqs4 * dSqs4 / (as4);
		else Ds4 = rhos4 * alpha[VY] * dSqs4 * dSqs4 / (as4);
	}
	else Ds4 = dbeta * rhos4 * alpha[VY] * dSqs4 * dSqs4 / (as4);
}

if (iT4 > -1) {
	if (!bT4) {
		if (bB4) Dt4 = dbeta * rhot4 * alpha[VZ] * dSqt4 * dSqt4 / (at4);
		else Dt4 = rhot4 * alpha[VZ] * dSqt4 * dSqt4 / (at4);
	}
	else Dt4 = dbeta * rhot4 * alpha[VZ] * dSqt4 * dSqt4 / (at4);
}

if (iB4 > -1) {
	if (!bB4) {
		if (bT4) Db4 = dbeta * rhob4 * alpha[VZ] * dSqb4 * dSqb4 / (ab4);
		else Db4 = rhob4 * alpha[VZ] * dSqb4 * dSqb4 / (ab4);
	}
	else Db4 = dbeta * rhob4 * alpha[VZ] * dSqb4 * dSqb4 / (ab4);
}
*/

	if ((De < 0.) || (Dw < 0.) || (Dn < 0.) || (Ds < 0.) || (Dt < 0.) || (Db < 0.)) {
		printf("Negative diffusion coefficients in PAM equation. \n");
		system("PAUSE");
	}

	if (b_on_adaptive_local_refinement_mesh) {

		if ((De2 < 0.) || (Dw2 < 0.) || (Dn2 < 0.) || (Ds2 < 0.) || (Dt2 < 0.) || (Db2 < 0.)) {
			printf("Negative 2 diffusion coefficients in PAM equation. \n");
			system("PAUSE");
		}
		if ((De3 < 0.) || (Dw3 < 0.) || (Dn3 < 0.) || (Ds3 < 0.) || (Dt3 < 0.) || (Db3 < 0.)) {
			printf("Negative 3 diffusion coefficients in PAM equation. \n");
			system("PAUSE");
		}
		if ((De4 < 0.) || (Dw4 < 0.) || (Dn4 < 0.) || (Ds4 < 0.) || (Dt4 < 0.) || (Db4 < 0.)) {
			printf("Negative 4 diffusion coefficients in PAM equation. \n");
			system("PAUSE");
		}
	}

	// На каждой грани все коэффициенты диффузии участвующие в дискретизации
	// равны одному и тому-же среднему значению.
	if (0) {
		AVERAGE_DIFFUSION(iE, iE2, iE3, iE4, De, De2, De3, De4);
		AVERAGE_DIFFUSION(iW, iW2, iW3, iW4, Dw, Dw2, Dw3, Dw4);
		AVERAGE_DIFFUSION(iN, iN2, iN3, iN4, Dn, Dn2, Dn3, Dn4);
		AVERAGE_DIFFUSION(iS, iS2, iS3, iS4, Ds, Ds2, Ds3, Ds4);
		AVERAGE_DIFFUSION(iT, iT2, iT3, iT4, Dt, Dt2, Dt3, Dt4);
		AVERAGE_DIFFUSION(iB, iB2, iB3, iB4, Db, Db2, Db3, Db4);
	}

	// Число Пекле равно нулю.
	sl[PAM][iP].ae=De*1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[PAM][iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[PAM][iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[PAM][iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[PAM][iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[PAM][iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	if (b_on_adaptive_local_refinement_mesh) {

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

		sl[PAM][iP].ap = sl[PAM][iP].ae + sl[PAM][iP].aw + sl[PAM][iP].an + sl[PAM][iP].as + sl[PAM][iP].at + sl[PAM][iP].ab +
			sl[PAM][iP].ae2 + sl[PAM][iP].aw2 + sl[PAM][iP].an2 + sl[PAM][iP].as2 + sl[PAM][iP].at2 + sl[PAM][iP].ab2 +
			sl[PAM][iP].ae3 + sl[PAM][iP].aw3 + sl[PAM][iP].an3 + sl[PAM][iP].as3 + sl[PAM][iP].at3 + sl[PAM][iP].ab3 +
			sl[PAM][iP].ae4 + sl[PAM][iP].aw4 + sl[PAM][iP].an4 + sl[PAM][iP].as4 + sl[PAM][iP].at4 + sl[PAM][iP].ab4;

	}
	else {

		sl[PAM][iP].ap = sl[PAM][iP].ae + sl[PAM][iP].aw + sl[PAM][iP].an + sl[PAM][iP].as + sl[PAM][iP].at + sl[PAM][iP].ab;
	}

	doublereal baddDFLUX2=0.0;
	if (0&&bhighorder) {
		// На АЛИС сетке это использовать не представляется возможным.

		// если bborder == false то узел строго внутренний.
		// если bborder   то мы находимся вблизи граничного узла.
		bool bborder=false;
		doublereal myflux=0.0;
		myflux=De*(dxe*DFDXiP(potent[PAM], iP, E_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iE]-potent[PAM][iP]));
		baddDFLUX2+=myflux;
		myflux=Dw*(-dxw*DFDXiP(potent[PAM], iP, W_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iW]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dn*(dyn*DFDXiP(potent[PAM], iP, N_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iN]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Ds*(-dys*DFDXiP(potent[PAM], iP, S_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iS]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Dt*(dzt*DFDXiP(potent[PAM], iP, T_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iT]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
		myflux=Db*(-dzb*DFDXiP(potent[PAM], iP, B_SIDE, neighbors_for_the_internal_node, maxelm, nvtx, pa, bborder)-(potent[PAM][iB]-potent[PAM][iP]));
	    baddDFLUX2+=myflux;
	}


	

    doublereal Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0, Ft=0.0, Fb=0.0; 


	
	if (iE > -1) {

		
		

		Fe = calcFg3(bE, feplus, iP, E_SIDE, rhoe,
			dSqe, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxe, potent, 1);

		if (Fe != Fe) {
			printf("Fe!=Fe\n");
			system("pause");
		}
	}

	
	if (iW > -1) {		

		Fw = calcFg3(bW, fwplus, iP, W_SIDE, rhow,
			dSqw, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dxw, potent, 1);
	}

	
	if (iN > -1) {		

		Fn = calcFg3(bN, fnplus, iP, N_SIDE, rhon,
			dSqn, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dyn, potent, 1);
	}

	
	if (iS > -1) {


		Fs = calcFg3(bS, fsplus, iP, S_SIDE, rhos,
			dSqs, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dys, potent, 1);
	}

	
	if (iT > -1) {

		

		Ft = calcFg3(bT, ftplus, iP, T_SIDE, rhot,
			dSqt, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzt, potent, 1);
	}

	
	if (iB > -1) {

		

		Fb = calcFg3(bB, fbplus, iP, B_SIDE, rhob,
			dSqb, btimedepend,
			speedoldtimestep, mfoldtimestep,
			dtimestep, RCh, neighbors_for_the_internal_node,
			bRhieChowi, bRhieChowb, false,
			tau, dzb, potent, 1);
	}

	
	doublereal Fw2 = 0.0, Fe2 = 0.0, Fs2 = 0.0, Fn2 = 0.0, Ft2 = 0.0, Fb2 = 0.0;
	doublereal Fw3 = 0.0, Fe3 = 0.0, Fs3 = 0.0, Fn3 = 0.0, Ft3 = 0.0, Fb3 = 0.0;
	doublereal Fw4 = 0.0, Fe4 = 0.0, Fs4 = 0.0, Fn4 = 0.0, Ft4 = 0.0, Fb4 = 0.0;

	if (b_on_adaptive_local_refinement_mesh) {

		if (iE2 > -1) {



			Fe2 = calcFg3(bE2, feplus2, iP, E_SIDE, rhoe2,
				dSqe2, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dxe2, potent, 2);
		}


		if (iW2 > -1) {


			Fw2 = calcFg3(bW2, fwplus2, iP, W_SIDE, rhow2,
				dSqw2, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dxw2, potent, 2);
		}


		if (iN2 > -1) {


			Fn2 = calcFg3(bN2, fnplus2, iP, N_SIDE, rhon2,
				dSqn2, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dyn2, potent, 2);
		}


		if (iS2 > -1) {



			Fs2 = calcFg3(bS2, fsplus2, iP, S_SIDE, rhos2,
				dSqs2, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dys2, potent, 2);
		}


		if (iT2 > -1) {



			Ft2 = calcFg3(bT2, ftplus2, iP, T_SIDE, rhot2,
				dSqt2, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dzt2, potent, 2);
		}


		if (iB2 > -1) {



			Fb2 = calcFg3(bB2, fbplus2, iP, B_SIDE, rhob2,
				dSqb2, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dzb2, potent, 2);
		}



		

		if (iE3 > -1) {



			Fe3 = calcFg3(bE3, feplus3, iP, E_SIDE, rhoe3,
				dSqe3, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dxe3, potent, 3);
		}


		if (iW3 > -1) {



			Fw3 = calcFg3(bW3, fwplus3, iP, W_SIDE, rhow3,
				dSqw3, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dxw3, potent, 3);
		}


		if (iN3 > -1) {


			Fn3 = calcFg3(bN3, fnplus3, iP, N_SIDE, rhon3,
				dSqn3, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dyn3, potent, 3);
		}


		if (iS3 > -1) {



			Fs3 = calcFg3(bS3, fsplus3, iP, S_SIDE, rhos3,
				dSqs3, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dys3, potent, 3);
		}


		if (iT3 > -1) {



			Ft3 = calcFg3(bT3, ftplus3, iP, T_SIDE, rhot3,
				dSqt3, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dzt3, potent, 3);
		}


		if (iB3 > -1) {



			Fb3 = calcFg3(bB3, fbplus3, iP, B_SIDE, rhob3,
				dSqb3, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dzb3, potent, 3);
		}


		


		if (iE4 > -1) {



			Fe4 = calcFg3(bE4, feplus4, iP, E_SIDE, rhoe4,
				dSqe4, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dxe4, potent, 4);
		}


		if (iW4 > -1) {


			Fw4 = calcFg3(bW4, fwplus4, iP, W_SIDE, rhow4,
				dSqw4, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dxw4, potent, 4);
		}


		if (iN4 > -1) {



			Fn4 = calcFg3(bN4, fnplus4, iP, N_SIDE, rhon4,
				dSqn4, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dyn4, potent, 4);
		}


		if (iS4 > -1) {


			Fs4 = calcFg3(bS4, fsplus4, iP, S_SIDE, rhos4,
				dSqs4, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dys4, potent, 4);
		}


		if (iT4 > -1) {



			Ft4 = calcFg3(bT4, ftplus4, iP, T_SIDE, rhot4,
				dSqt4, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dzt4, potent, 4);
		}


		if (iB4 > -1) {


			Fb4 = calcFg3(bB4, fbplus4, iP, B_SIDE, rhob4,
				dSqb4, btimedepend,
				speedoldtimestep, mfoldtimestep,
				dtimestep, RCh, neighbors_for_the_internal_node,
				bRhieChowi, bRhieChowb, false,
				tau, dzb4, potent, 4);
		}

	}

	doublereal Fe_sum = Fe + Fe2 + Fe3 + Fe4;
	doublereal Fw_sum = Fw + Fw2 + Fw3 + Fw4;
	doublereal Fn_sum = Fn + Fn2 + Fn3 + Fn4;
	doublereal Fs_sum = Fs + Fs2 + Fs3 + Fs4;
	doublereal Ft_sum = Ft + Ft2 + Ft3 + Ft4;
	doublereal Fb_sum = Fb + Fb2 + Fb3 + Fb4;

	if (Fe_sum != Fe_sum) {
		printf("Fe=%e Fe2=%e Fe3=%e Fe4=%e\n",Fe,Fe2,Fe3,Fe4);
		system("pause");
	}
	if (Fw_sum != Fw_sum) {
		printf("Fw=%e Fw2=%e Fw3=%e Fw4=%e\n", Fw, Fw2, Fw3, Fw4);
		system("pause");
	}
	if (Fn_sum != Fn_sum) {
		printf("Fn=%e Fn2=%e Fn3=%e Fn4=%e\n", Fn, Fn2, Fn3, Fn4);
		system("pause");
	}

	if (Fs_sum != Fs_sum) {
		printf("Fs=%e Fs2=%e Fs3=%e Fs4=%e\n", Fs, Fs2, Fs3, Fs4);
		system("pause");
	}
	if (Ft_sum != Ft_sum) {
		printf("Ft=%e Ft2=%e Ft3=%e Ft4=%e\n", Ft, Ft2, Ft3, Ft4);
		system("pause");
	}
	if (Fb_sum != Fb_sum) {
		printf("Fb=%e Fb2=%e Fb3=%e Fb4=%e\n", Fb, Fb2, Fb3, Fb4);
		system("pause");
	}

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
		doublereal SpeedCorOlde3 = 0.0, SpeedCorOldw3 = 0.0, SpeedCorOldn3 = 0.0, SpeedCorOlds3 = 0.0, SpeedCorOldt3 = 0.0, SpeedCorOldb3 = 0.0;
		doublereal SpeedCorOlde4 = 0.0, SpeedCorOldw4 = 0.0, SpeedCorOldn4 = 0.0, SpeedCorOlds4 = 0.0, SpeedCorOldt4 = 0.0, SpeedCorOldb4 = 0.0;


		if (b_on_adaptive_local_refinement_mesh) {

			if (iE2 > -1) {
				if (!bE2) {
					SpeedCorOlde2 = feplus2 * potent[VXCOR][iE2] + (1.0 - feplus2) * potent[VXCOR][iP];
				}
				else {
					SpeedCorOlde2 = potent[VXCOR][iE2];
				}
			}
			if (iN2 > -1) {
				if (!bN2) {
					SpeedCorOldn2 = fnplus2 * potent[VYCOR][iN2] + (1.0 - fnplus2) * potent[VYCOR][iP];
				}
				else {
					SpeedCorOldn2 = potent[VYCOR][iN2];
				}
			}
			if (iT2 > -1) {
				if (!bT2) {
					SpeedCorOldt2 = ftplus2 * potent[VZCOR][iT2] + (1.0 - ftplus2) * potent[VZCOR][iP];
				}
				else {
					SpeedCorOldt2 = potent[VZCOR][iT2];
				}
			}
			if (iW2 > -1) {
				if (!bW2) {
					SpeedCorOldw2 = fwplus2 * potent[VXCOR][iW2] + (1.0 - fwplus2) * potent[VXCOR][iP];
				}
				else {
					SpeedCorOldw2 = potent[VXCOR][iW2];
				}
			}
			if (iS2 > -1) {
				if (!bS2) {
					SpeedCorOlds2 = fsplus2 * potent[VYCOR][iS2] + (1.0 - fsplus2) * potent[VYCOR][iP];
				}
				else {
					SpeedCorOlds2 = potent[VYCOR][iS2];
				}
			}
			if (iB2 > -1) {
				if (!bB2) {
					SpeedCorOldb2 = fbplus2 * potent[VZCOR][iB2] + (1.0 - fbplus2) * potent[VZCOR][iP];
				}
				else {
					SpeedCorOldb2 = potent[VZCOR][iB2];
				}
			}

			
			if (iE3 > -1) {
				if (!bE3) {
					SpeedCorOlde3 = feplus3 * potent[VXCOR][iE3] + (1.0 - feplus3) * potent[VXCOR][iP];
				}
				else {
					SpeedCorOlde3 = potent[VXCOR][iE3];
				}
			}
			if (iN3 > -1) {
				if (!bN3) {
					SpeedCorOldn3 = fnplus3 * potent[VYCOR][iN3] + (1.0 - fnplus3) * potent[VYCOR][iP];
				}
				else {
					SpeedCorOldn3 = potent[VYCOR][iN3];
				}
			}
			if (iT3 > -1) {
				if (!bT3) {
					SpeedCorOldt3 = ftplus3 * potent[VZCOR][iT3] + (1.0 - ftplus3) * potent[VZCOR][iP];
				}
				else {
					SpeedCorOldt3 = potent[VZCOR][iT3];
				}
			}
			if (iW3 > -1) {
				if (!bW3) {
					SpeedCorOldw3 = fwplus3 * potent[VXCOR][iW3] + (1.0 - fwplus3) * potent[VXCOR][iP];
				}
				else {
					SpeedCorOldw3 = potent[VXCOR][iW3];
				}
			}
			if (iS3 > -1) {
				if (!bS3) {
					SpeedCorOlds3 = fsplus3 * potent[VYCOR][iS3] + (1.0 - fsplus3) * potent[VYCOR][iP];
				}
				else {
					SpeedCorOlds3 = potent[VYCOR][iS3];
				}
			}
			if (iB3 > -1) {
				if (!bB3) {
					SpeedCorOldb3 = fbplus3 * potent[VZCOR][iB3] + (1.0 - fbplus3) * potent[VZCOR][iP];
				}
				else {
					SpeedCorOldb3 = potent[VZCOR][iB3];
				}
			}


			
			if (iE4 > -1) {
				if (!bE4) {
					SpeedCorOlde4 = feplus4 * potent[VXCOR][iE4] + (1.0 - feplus4) * potent[VXCOR][iP];
				}
				else {
					SpeedCorOlde4 = potent[VXCOR][iE4];
				}
			}
			if (iN4 > -1) {
				if (!bN4) {
					SpeedCorOldn4 = fnplus4 * potent[VYCOR][iN4] + (1.0 - fnplus4) * potent[VYCOR][iP];
				}
				else {
					SpeedCorOldn4 = potent[VYCOR][iN4];
				}
			}
			if (iT4 > -1) {
				if (!bT4) {
					SpeedCorOldt4 = ftplus4 * potent[VZCOR][iT4] + (1.0 - ftplus4) * potent[VZCOR][iP];
				}
				else {
					SpeedCorOldt4 = potent[VZCOR][iT4];
				}
			}
			if (iW4 > -1) {
				if (!bW4) {
					SpeedCorOldw4 = fwplus4 * potent[VXCOR][iW4] + (1.0 - fwplus4) * potent[VXCOR][iP];
				}
				else {
					SpeedCorOldw4 = potent[VXCOR][iW4];
				}
			}
			if (iS4 > -1) {
				if (!bS4) {
					SpeedCorOlds4 = fsplus4 * potent[VYCOR][iS4] + (1.0 - fsplus4) * potent[VYCOR][iP];
				}
				else {
					SpeedCorOlds4 = potent[VYCOR][iS4];
				}
			}
			if (iB4 > -1) {
				if (!bB4) {
					SpeedCorOldb4 = fbplus4 * potent[VZCOR][iB4] + (1.0 - fbplus4) * potent[VZCOR][iP];
				}
				else {
					SpeedCorOldb4 = potent[VZCOR][iB4];
				}
			}

		}

	    // возвращаем значение потока на грани КО.
	    // С включённой поправкой (дополнительная нижняя релаксация) из статьи I. Sezai. !!!
	    // дополнительная поправка осуществляется на основе скорректированной скорости.
	    // Вообще говоря поле плотности также должно быть с предыдущей итерации.
	    Fe_sum +=(1.0-alpha[VELOCITY_X_COMPONENT])*(mfcurrentretune[E_SIDE] - rhoe * SpeedCorOlde*dSqe - rhoe2 * SpeedCorOlde2*dSqe2 - rhoe3 * SpeedCorOlde3*dSqe3 - rhoe4 * SpeedCorOlde4*dSqe4);
	    Fn_sum +=(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfcurrentretune[N_SIDE] - rhon * SpeedCorOldn*dSqn - rhon2 * SpeedCorOldn2*dSqn2 - rhon3 * SpeedCorOldn3*dSqn3 - rhon4 * SpeedCorOldn4*dSqn4);
	    Ft_sum +=(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfcurrentretune[T_SIDE] - rhot * SpeedCorOldt*dSqt - rhot2 * SpeedCorOldt2*dSqt2 - rhot3 * SpeedCorOldt3*dSqt3 - rhot4 * SpeedCorOldt4*dSqt4);
	    Fw_sum +=(1.0-alpha[VELOCITY_X_COMPONENT])*(mfcurrentretune[W_SIDE] - rhow * SpeedCorOldw*dSqw - rhow2 * SpeedCorOldw2*dSqw2 - rhow3 * SpeedCorOldw3*dSqw3 - rhow4 * SpeedCorOldw4*dSqw4);
	    Fs_sum +=(1.0-alpha[VELOCITY_Y_COMPONENT])*(mfcurrentretune[S_SIDE] - rhos * SpeedCorOlds*dSqs - rhos2 * SpeedCorOlds2*dSqs2 - rhos3 * SpeedCorOlds3*dSqs3 - rhos4 * SpeedCorOlds4*dSqs4);
	    Fb_sum +=(1.0-alpha[VELOCITY_Z_COMPONENT])*(mfcurrentretune[B_SIDE] - rhob * SpeedCorOldb*dSqb - rhob2 * SpeedCorOldb2*dSqb2 - rhob3 * SpeedCorOldb3*dSqb3 - rhob4 * SpeedCorOldb4*dSqb4);

	} 
	

	
	// возвращаем значение потока на грани КО.
	if (bdeltapfinish==true) {
	   // запоминаем финишный поток через грань.
	   mfcurrentretune[E_SIDE]=Fe_sum;
	   mfcurrentretune[N_SIDE]=Fn_sum;
	   mfcurrentretune[T_SIDE]=Ft_sum;
	   mfcurrentretune[W_SIDE]=Fw_sum;
	   mfcurrentretune[S_SIDE]=Fs_sum;
	   mfcurrentretune[B_SIDE]=Fb_sum;
	}
	
	

	sl[PAM][iP].b=(Fw_sum -Fe_sum +Fs_sum -Fn_sum +Fb_sum -Ft_sum + baddDFLUX2);
	if (sl[PAM][iP].b != sl[PAM][iP].b) {
		printf("Fw_sum=%e Fe_sum=%e Fs_sum=%e Fn_sum=%e Fb_sum=%e Ft_sum=%e baddDFLUX2=%e\n", Fw_sum, Fe_sum, Fs_sum, Fn_sum, Fb_sum, Ft_sum, baddDFLUX2);
		system("PAUSE");
		exit(1);
	}
	//sl[PAM][iP].b = Fw_sum - Fe_sum + Fs_sum - Fn_sum + Fb_sum - Ft_sum; // 20.07.2016
	/*
	//02.01.2018 Вроде Ok.
	printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e\n", sl[PAM][iP].ae, sl[PAM][iP].aw, sl[PAM][iP].an, sl[PAM][iP].as, sl[PAM][iP].at, sl[PAM][iP].ab);
	printf("ae2=%e aw2=%e an2=%e as2=%e at2=%e ab2=%e\n", sl[PAM][iP].ae2, sl[PAM][iP].aw2, sl[PAM][iP].an2, sl[PAM][iP].as2, sl[PAM][iP].at2, sl[PAM][iP].ab2);
	printf("ae3=%e aw3=%e an3=%e as3=%e at3=%e ab3=%e\n", sl[PAM][iP].ae3, sl[PAM][iP].aw3, sl[PAM][iP].an3, sl[PAM][iP].as3, sl[PAM][iP].at3, sl[PAM][iP].ab3);
	printf("ae4=%e aw4=%e an4=%e as4=%e at4=%e ab4=%e\n", sl[PAM][iP].ae4, sl[PAM][iP].aw4, sl[PAM][iP].an4, sl[PAM][iP].as4, sl[PAM][iP].at4, sl[PAM][iP].ab4);
	printf("ap=%e b=%e\n", sl[PAM][iP].ap, sl[PAM][iP].b);
	*/
	//getchar();
	
	// Симметризация СЛАУ:
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

	if (fabs(sl[PAM][iP].ap) < 1.0e-300) {		
		printf("Error sl[PAM][%lld].ap=%e < 1.0e-300 PAM in function my_elmatr_quad_PAm3\n", iP, sl[PAM][iP].ap);
		printf("sl[PAM][%lld].ae=%e \n", iP, sl[PAM][iP].ae);
		printf("sl[PAM][%lld].aw=%e \n", iP, sl[PAM][iP].aw);
		printf("sl[PAM][%lld].an=%e \n", iP, sl[PAM][iP].an);
		printf("sl[PAM][%lld].as=%e \n", iP, sl[PAM][iP].as);
		printf("sl[PAM][%lld].at=%e \n", iP, sl[PAM][iP].at);
		printf("sl[PAM][%lld].ab=%e \n", iP, sl[PAM][iP].ab);
		if (b_on_adaptive_local_refinement_mesh) {
			printf("sl[PAM][%lld].ae2=%e \n", iP, sl[PAM][iP].ae2);
			printf("sl[PAM][%lld].aw2=%e \n", iP, sl[PAM][iP].aw2);
			printf("sl[PAM][%lld].an2=%e \n", iP, sl[PAM][iP].an2);
			printf("sl[PAM][%lld].as2=%e \n", iP, sl[PAM][iP].as2);
			printf("sl[PAM][%lld].at2=%e \n", iP, sl[PAM][iP].at2);
			printf("sl[PAM][%lld].ab2=%e \n", iP, sl[PAM][iP].ab2);
			printf("sl[PAM][%lld].ae3=%e \n", iP, sl[PAM][iP].ae3);
			printf("sl[PAM][%lld].aw3=%e \n", iP, sl[PAM][iP].aw3);
			printf("sl[PAM][%lld].an3=%e \n", iP, sl[PAM][iP].an3);
			printf("sl[PAM][%lld].as3=%e \n", iP, sl[PAM][iP].as3);
			printf("sl[PAM][%lld].at3=%e \n", iP, sl[PAM][iP].at3);
			printf("sl[PAM][%lld].ab3=%e \n", iP, sl[PAM][iP].ab3);
			printf("sl[PAM][%lld].ae4=%e \n", iP, sl[PAM][iP].ae4);
			printf("sl[PAM][%lld].aw4=%e \n", iP, sl[PAM][iP].aw4);
			printf("sl[PAM][%lld].an4=%e \n", iP, sl[PAM][iP].an4);
			printf("sl[PAM][%lld].as4=%e \n", iP, sl[PAM][iP].as4);
			printf("sl[PAM][%lld].at4=%e \n", iP, sl[PAM][iP].at4);
			printf("sl[PAM][%lld].ab4=%e \n", iP, sl[PAM][iP].ab4);
		}
		if (sl[PAM][iP].b != sl[PAM][iP].b) {
			printf("Zero ap in velocity component.\n");
		}
		if (fabs(sl[PAM][iP].ap) > 1.0e-301) {
			// Проведем ренормировку.
			sl[PAM][iP].ae /= sl[PAM][iP].ap;
			sl[PAM][iP].aw /= sl[PAM][iP].ap;
			sl[PAM][iP].an /= sl[PAM][iP].ap;
			sl[PAM][iP].as /= sl[PAM][iP].ap;
			sl[PAM][iP].at /= sl[PAM][iP].ap;
			sl[PAM][iP].ab /= sl[PAM][iP].ap;

			if (b_on_adaptive_local_refinement_mesh) {
				sl[PAM][iP].ae2 /= sl[PAM][iP].ap;
				sl[PAM][iP].aw2 /= sl[PAM][iP].ap;
				sl[PAM][iP].an2 /= sl[PAM][iP].ap;
				sl[PAM][iP].as2 /= sl[PAM][iP].ap;
				sl[PAM][iP].at2 /= sl[PAM][iP].ap;
				sl[PAM][iP].ab2 /= sl[PAM][iP].ap;
				sl[PAM][iP].ae3 /= sl[PAM][iP].ap;
				sl[PAM][iP].aw3 /= sl[PAM][iP].ap;
				sl[PAM][iP].an3 /= sl[PAM][iP].ap;
				sl[PAM][iP].as3 /= sl[PAM][iP].ap;
				sl[PAM][iP].at3 /= sl[PAM][iP].ap;
				sl[PAM][iP].ab3 /= sl[PAM][iP].ap;
				sl[PAM][iP].ae4 /= sl[PAM][iP].ap;
				sl[PAM][iP].aw4 /= sl[PAM][iP].ap;
				sl[PAM][iP].an4 /= sl[PAM][iP].ap;
				sl[PAM][iP].as4 /= sl[PAM][iP].ap;
				sl[PAM][iP].at4 /= sl[PAM][iP].ap;
				sl[PAM][iP].ab4 /= sl[PAM][iP].ap;
			}

			sl[PAM][iP].b/= sl[PAM][iP].ap;
			sl[PAM][iP].ap = 1.0;
		}
		else {
			sl[PAM][iP].ap = 1.0;
			system("pause");
		}
		
	}



	if (sl[PAM][iP].ap != sl[PAM][iP].ap) {
		printf("ap!=ap assemble bug. iP=%lld ap=%e\n", iP, sl[PAM][iP].ap);
		printf("PAM \n");
		system("pause");
	}
	if (sl[PAM][iP].ae != sl[PAM][iP].ae) {
		printf("ae!=ae assemble bug\n");
		printf("PAM \n");
		system("pause");
	}
	if (sl[PAM][iP].aw != sl[PAM][iP].aw) {
		printf("aw!=aw assemble bug\n");
		printf("PAM \n");
		system("pause");
	}
	if (sl[PAM][iP].an != sl[PAM][iP].an) {
		printf("an!=an assemble bug\n");
		printf("PAM \n");
		system("pause");
	}
	if (sl[PAM][iP].as != sl[PAM][iP].as) {
		printf("as!=as assemble bug\n");
		printf("PAM \n");
		system("pause");
	}
	if (sl[PAM][iP].at != sl[PAM][iP].at) {
		printf("at!=at assemble bug\n");
		printf("PAM \n");
		system("pause");
	}
	if (sl[PAM][iP].ab != sl[PAM][iP].ab) {
		printf("ab!=ab assemble bug\n");
		printf("PAM \n");
		system("pause");
	}
	if (b_on_adaptive_local_refinement_mesh) {

		if (sl[PAM][iP].ae2 != sl[PAM][iP].ae2) {
			printf("ae2!=ae2 assemble bug %e %e\n", sl[PAM][iP].ae2, sl[PAM][iP].ae2);
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].aw2 != sl[PAM][iP].aw2) {
			printf("aw2!=aw2 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].an2 != sl[PAM][iP].an2) {
			printf("an2!=an2 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].as2 != sl[PAM][iP].as2) {
			printf("as2!=as2 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].at2 != sl[PAM][iP].at2) {
			printf("at2!=at2 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].ab2 != sl[PAM][iP].ab2) {
			printf("ab2!=ab2 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].ae3 != sl[PAM][iP].ae3) {
			printf("ae3!=ae3 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].aw3 != sl[PAM][iP].aw3) {
			printf("aw3!=aw3 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].an3 != sl[PAM][iP].an3) {
			printf("an3!=an3 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].as3 != sl[PAM][iP].as3) {
			printf("as3!=as3 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].at3 != sl[PAM][iP].at3) {
			printf("at3!=at3 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].ab3 != sl[PAM][iP].ab3) {
			printf("ab3!=ab3 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].ae4 != sl[PAM][iP].ae4) {
			printf("ae4!=ae4 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].aw4 != sl[PAM][iP].aw4) {
			printf("aw4!=aw4 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].an4 != sl[PAM][iP].an4) {
			printf("an4!=an4 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].as4 != sl[PAM][iP].as4) {
			printf("as4!=as4 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].at4 != sl[PAM][iP].at4) {
			printf("at4!=at4 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
		if (sl[PAM][iP].ab4 != sl[PAM][iP].ab4) {
			printf("ab4!=ab4 assemble bug\n");
			printf("PAM \n");
			system("pause");
		}
	}

} // my_elmatr_quad_PAm3

// Цель данной функции скорректировать поправку давления на границе
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
	        iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP];
	        iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];
	        

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
				integer iEE=f.neighbors_for_the_internal_node[E_SIDE][0][iE];
				switch (imy_interpol) {
				case LINEAR_INTERPOL: 
		                   center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                   center_cord3D(iE, f.nvtx, f.pa, pb,E_SIDE);
		                   f.potent[PAM][iW]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iE], pp.x, pb.x, pp.x-0.5*dx);
						   break;
				case QUAD_INTERPOL:
					        center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                    center_cord3D(iE, f.nvtx, f.pa, pb,E_SIDE);
						    center_cord3D(iEE, f.nvtx, f.pa, pbb,EE_SIDE);
							f.potent[PAM][iW]=my_quadratic_interpolation('-', f.potent[PAM][iEE], f.potent[PAM][iE], f.potent[PAM][iP],
	                            pbb.x , pb.x, pp.x, pp.x-0.5*dx);
					       break;
				case NO_INTERPOL: f.potent[PAM][iW] = f.potent[PAM][iP];
					break;
				default:
					// узел W граничный
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iE, f.nvtx, f.pa, pb,E_SIDE);
		            f.potent[PAM][iW]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iE], pp.x, pb.x, pp.x-0.5*dx);
					break;
				}
	        }
	        if (bE) {
		        // узел E граничный
				TOCHKA pp,pb,pbb;
				integer iWW=f.neighbors_for_the_internal_node[W_SIDE][0][iW];
				switch (imy_interpol) {
				case LINEAR_INTERPOL:
		                  center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                  center_cord3D(iW, f.nvtx, f.pa, pb,W_SIDE);
		                  f.potent[PAM][iE]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iW], pp.x, pb.x, pp.x+0.5*dx);
						  break;
				case QUAD_INTERPOL:
					      center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                  center_cord3D(iW, f.nvtx, f.pa, pb,W_SIDE);
						  center_cord3D(iWW, f.nvtx, f.pa, pbb,WW_SIDE);
						  
						  f.potent[PAM][iE]=my_quadratic_interpolation('+', f.potent[PAM][iWW], f.potent[PAM][iW], f.potent[PAM][iP], pbb.x, pb.x, pp.x, pp.x+0.5*dx);	  
					break;
				case NO_INTERPOL: f.potent[PAM][iE] = f.potent[PAM][iP];
					break;
				default: 
					center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iW, f.nvtx, f.pa, pb,W_SIDE);
		            f.potent[PAM][iE]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iW], pp.x, pb.x, pp.x+0.5*dx);  
					break;
				}
	        }
	        if (bS) {
		        // узел S граничный
				TOCHKA pp,pb,pbb;
				integer iNN=f.neighbors_for_the_internal_node[N_SIDE][0][iN];
				switch (imy_interpol) {
                case LINEAR_INTERPOL:
		                 center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                 center_cord3D(iN, f.nvtx, f.pa, pb,N_SIDE);
		                 f.potent[PAM][iS]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iN], pp.y, pb.y, pp.y-0.5*dy);
				         break;
				case QUAD_INTERPOL:
		                 center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                 center_cord3D(iN, f.nvtx, f.pa, pb,N_SIDE);
						 center_cord3D(iNN, f.nvtx, f.pa, pbb,NN_SIDE);
						 f.potent[PAM][iS]=my_quadratic_interpolation('-', f.potent[PAM][iNN], f.potent[PAM][iN], f.potent[PAM][iP],
	                            pbb.y, pb.y, pp.y, pp.y-0.5*dy);
					break;
				case NO_INTERPOL: f.potent[PAM][iS] = f.potent[PAM][iP];
					break;
				default:
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iN, f.nvtx, f.pa, pb,N_SIDE);
		            f.potent[PAM][iS]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iN], pp.y, pb.y, pp.y-0.5*dy);
					break;
				}
	        }
	        if (bN) {
		        // узел N граничный
				TOCHKA pp,pb,pbb;
				integer iSS=f.neighbors_for_the_internal_node[S_SIDE][0][iS];
				switch (imy_interpol) {
					case LINEAR_INTERPOL:
		                 center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                 center_cord3D(iS, f.nvtx, f.pa, pb,S_SIDE);
		                 f.potent[PAM][iN]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iS], pp.y, pb.y, pp.y+0.5*dy);
				    break;
					case QUAD_INTERPOL:
		                center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                center_cord3D(iS, f.nvtx, f.pa, pb,S_SIDE);
						center_cord3D(iSS, f.nvtx, f.pa, pbb,SS_SIDE);
						f.potent[PAM][iN]=my_quadratic_interpolation('+', f.potent[PAM][iSS], f.potent[PAM][iS], f.potent[PAM][iP],
	                            pbb.y , pb.y, pp.y, pp.y+0.5*dy);
					break;
					case NO_INTERPOL: f.potent[PAM][iN] = f.potent[PAM][iP];
						break;
					default:
						center_cord3D(iP, f.nvtx, f.pa, pp,100);
		                center_cord3D(iS, f.nvtx, f.pa, pb,S_SIDE);
		                f.potent[PAM][iN]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iS], pp.y, pb.y, pp.y+0.5*dy);
					break;
				}
	        }
	       if (bB) {
		        // узел B граничный
			   TOCHKA pp,pb,pbb;
			   integer iTT=f.neighbors_for_the_internal_node[T_SIDE][0][iT];
				switch (imy_interpol) {
		        case LINEAR_INTERPOL:
		               center_cord3D(iP, f.nvtx, f.pa, pp,100);
		               center_cord3D(iT, f.nvtx, f.pa, pb,T_SIDE);
		               f.potent[PAM][iB]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iT], pp.z, pb.z, pp.z-0.5*dz);
				   break;
				case QUAD_INTERPOL:
					center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iT, f.nvtx, f.pa, pb,T_SIDE);
					center_cord3D(iTT, f.nvtx, f.pa, pbb,TT_SIDE);
					f.potent[PAM][iB]=my_quadratic_interpolation('-', f.potent[PAM][iTT], f.potent[PAM][iT], f.potent[PAM][iP],
	                            pbb.z, pb.z, pp.z, pp.z-0.5*dz);

					break;
				case NO_INTERPOL: f.potent[PAM][iB] = f.potent[PAM][iP];
					break;
				default:
					center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iT, f.nvtx, f.pa, pb,T_SIDE);
		            f.potent[PAM][iB]=my_linear_interpolation('-', f.potent[PAM][iP], f.potent[PAM][iT], pp.z, pb.z, pp.z-0.5*dz);
					break;
				}
	       }
	       if (bT) { 
		        // узел T граничный
		        TOCHKA pp,pb, pbb;
				integer iBB=f.neighbors_for_the_internal_node[B_SIDE][0][iB];
				switch (imy_interpol) {
				case LINEAR_INTERPOL:
		              center_cord3D(iP, f.nvtx, f.pa, pp,100);
		              center_cord3D(iB, f.nvtx, f.pa, pb,B_SIDE);
		              f.potent[PAM][iT]=my_linear_interpolation('+', f.potent[PAM][iP], f.potent[PAM][iB], pp.z, pb.z, pp.z+0.5*dz);
					  break;
                case QUAD_INTERPOL:
					  center_cord3D(iP, f.nvtx, f.pa, pp,100);
		              center_cord3D(iB, f.nvtx, f.pa, pb,B_SIDE);
					  center_cord3D(iBB, f.nvtx, f.pa, pbb,BB_SIDE);
					  f.potent[PAM][iT]=my_quadratic_interpolation('+', f.potent[PAM][iBB], f.potent[PAM][iB], f.potent[PAM][iP],
	                            pbb.z , pb.z, pp.z, pp.z+0.5*dz);
					break;
				case NO_INTERPOL: f.potent[PAM][iT] = f.potent[PAM][iP];
					break;
				default:
					  center_cord3D(iP, f.nvtx, f.pa, pp,100);
		              center_cord3D(iB, f.nvtx, f.pa, pb,B_SIDE);
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
	        iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP];
	        iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];
	        

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
		       if ((f.border_neighbor[inumber].MCB<(ls+lw)) && (f.border_neighbor[inumber].MCB>=ls) && (w[f.border_neighbor[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iW],pamout);
			   }
	        }
	        if (bE) {
		        // узел E граничный
		        inumber=iE-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls+lw)) && (f.border_neighbor[inumber].MCB>=ls) && (w[f.border_neighbor[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iE],pamout);
			   }
	        }
	        if (bS) {
		        // узел S граничный
		        inumber=iS-f.maxelm;
                if ((f.border_neighbor[inumber].MCB<(ls+lw)) && (f.border_neighbor[inumber].MCB>=ls) && (w[f.border_neighbor[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iS],pamout);
			   }
	        }
	        if (bN) {
		        // узел N граничный
		        inumber=iN-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls+lw)) && (f.border_neighbor[inumber].MCB>=ls) && (w[f.border_neighbor[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iN],pamout);
			    }
	        }
	       if (bB) {
		        // узел B граничный
		        inumber=iB-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls+lw)) && (f.border_neighbor[inumber].MCB>=ls) && (w[f.border_neighbor[inumber].MCB-ls].bpressure)) {
				   pamout=fmin(f.potent[PAM][iB],pamout);
			   }
	       }
	       if (bT) { 
		        // узел T граничный
		        inumber=iT-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls+lw)) && (f.border_neighbor[inumber].MCB>=ls) && (w[f.border_neighbor[inumber].MCB-ls].bpressure)) {
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
	        iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP];
	        iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];
	        

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
				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure || w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				   // Линейная интерполяция
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iE, f.nvtx, f.pa, pb,E_SIDE);
		           f.potent[iVar][iW]=my_linear_interpolation('-', f.potent[iVar][iP], f.potent[iVar][iE], pp.x, pb.x, pp.x-0.5*dx);
			   }
	        }
	        if (bE) {
		        // узел E граничный
		        inumber=iE-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure || w[f.border_neighbor[inumber].MCB - ls].bopening))) {
					// Линейная интерполяция
				    TOCHKA pp,pb;
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iW, f.nvtx, f.pa, pb,W_SIDE);
		            f.potent[iVar][iE]=my_linear_interpolation('+', f.potent[iVar][iP], f.potent[iVar][iW], pp.x, pb.x, pp.x+0.5*dx);
			   }
	        }
	        if (bS) {
		        // узел S граничный
		        inumber=iS-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure || w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				    TOCHKA pp,pb;
		            center_cord3D(iP, f.nvtx, f.pa, pp,100);
		            center_cord3D(iN, f.nvtx, f.pa, pb,N_SIDE);
		            f.potent[iVar][iS]=my_linear_interpolation('-', f.potent[iVar][iP], f.potent[iVar][iN], pp.y, pb.y, pp.y-0.5*dy);
			   }
	        }
	        if (bN) {
		        // узел N граничный
		        inumber=iN-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure || w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iS, f.nvtx, f.pa, pb,S_SIDE);
		           f.potent[iVar][iN]=my_linear_interpolation('+', f.potent[iVar][iP], f.potent[iVar][iS], pp.y, pb.y, pp.y+0.5*dy);
			    }
	        }
	       if (bB) {
		        // узел B граничный
		        inumber=iB-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure || w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iT, f.nvtx, f.pa, pb,T_SIDE);
		           f.potent[iVar][iB]=my_linear_interpolation('-', f.potent[iVar][iP], f.potent[iVar][iT], pp.z, pb.z, pp.z-0.5*dz);
			   }
	       }
	       if (bT) { 
		        // узел T граничный
		        inumber=iT-f.maxelm;
				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure ||  w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				   TOCHKA pp,pb;
		           center_cord3D(iP, f.nvtx, f.pa, pp,100);
		           center_cord3D(iB, f.nvtx, f.pa, pb,B_SIDE);
		           f.potent[iVar][iT]=my_linear_interpolation('+', f.potent[iVar][iP], f.potent[iVar][iB], pp.z, pb.z, pp.z+0.5*dz);
			   }
	       }
	
		
	}

} // interpolatevel

#endif
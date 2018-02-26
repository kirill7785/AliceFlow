// Файл my_interpolate_v0_07.cpp содержит
// некоторые необходимые для интерполляции функции.
// Перечень функций следующий:
// 1. Формула вычисления первой производной на неравномерной сетке с точностью O(h^2).
// Данная формула заимствована из книги Г.З. Гарбера и обеспечивает второй порядок на неравномерной сетке.
// Формула является трёхточечной и основа на построении параболы по трём точкам и взятии её производной в центральной точке.

// Данные функции объединены в модуль 15 мая 2012 года.

#ifndef MY_INTERPOLATE_V0_07_CPP
#define MY_INTERPOLATE_V0_07_CPP 1

// Производная от искомой величины в точности на грани КО.
// позиция грани xg.
doublereal DFDXg(doublereal x1, doublereal x2, doublereal x3, doublereal x4, doublereal xg,
	      doublereal f1, doublereal f2, doublereal f3, doublereal f4) {

	doublereal DFDXgran=0.0;

	const integer parabola_MNK=0; // производная от параболы построенная по способу наименьших квадратов.
	// Внимание : line_MNK лучше не использовать, т.к. экспериментально установлено, на примере задачи
	// статики tgf2023_01 что сходимости при таком задании производной не наступает, наблюдаются некие автоколебания на протяжени
	// более чем 40 глобальных итераций и более. Возможна даже расходимость.
	const integer line_MNK=1; // тангенс угла наклона прямой, прямая построена по способу МНК по четырём точкам.
	// данный метод также не рекомендуется использовать ! ввиду отсутствия сходимости.
	const integer cubic_parabola=2; // кубическая парабола по 4 точкам и от неё производная.

	integer imetod=parabola_MNK;

	
	if (imetod==parabola_MNK) 
	{

	     // проверено.

	     // Грань одназначно определяется четырьмя соседними Контрольными объёмами.
	     // Построим по 4 точкам параболу используя МНК.

		doublereal* a = NULL;
		a=new doublereal[3];
		doublereal** AM = NULL;
		AM=new doublereal*[3];
	     for (integer i=0; i<3; i++) {
		     AM[i]=new doublereal[3]; // матрица СЛАУ
	     }
		 doublereal* b = NULL;
		 b=new doublereal[3]; // правая часть

	     for (integer i=0; i<3; i++) {
		     for (integer j=0; j<3; j++) {
		    	AM[i][j]=0.0; // инициализация
		     }
		     b[i]=0.0;
		     a[i]=0.0;
	     }

	     AM[0][0]=4.0;
	     AM[0][1]=x1+x2+x3+x4;
	     AM[0][2]=x1*x1+x2*x2+x3*x3+x4*x4;
	     b[0]=f1+f2+f3+f4;
	     AM[1][0]=x1+x2+x3+x4;
	     AM[1][1]=x1*x1+x2*x2+x3*x3+x4*x4;
	     AM[1][2]=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4;
	     b[1]=f1*x1+f2*x2+f3*x3+f4*x4;
	     AM[2][0]=x1*x1+x2*x2+x3*x3+x4*x4;
	     AM[2][1]=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4;
	     AM[2][2]=x1*x1*x1*x1+x2*x2*x2*x2+x3*x3*x3*x3+x4*x4*x4*x4;
	     b[2]=f1*x1*x1+f2*x2*x2+f3*x3*x3+f4*x4*x4;

	     integer nodes=3;
	     eqsolve_simple_gauss(AM, nodes, b, a); // решение СЛАУ методом Гаусса.

	     // Освобождение памяти
		 if (AM != NULL) {
			 for (integer i = 0; i < 3; i++) {
				 if (AM[i] != NULL) {
					 delete[] AM[i];
				 }
			 }
			 delete[] AM;
		 }
		 if (b != NULL) {
			 delete[] b;
		 }

	     //f(x)=a[0]+a[1]*x+a[2]*x^2;

	     //Первая производная в точке xg равна 2.0*a[2]*xg+a[1];
	     DFDXgran=2.0*a[2]*xg+a[1];

		 if (a != NULL) {
			 delete[] a;
		 }

	}

	if (imetod==line_MNK) {  
		// по 4 точкам по способу наименьших квадратов строим линию и находим её производную аналитически.


		doublereal* a = NULL;
		a=new doublereal[2];
		doublereal** AM = NULL;
		AM=new doublereal*[2];
	    for (integer i=0; i<2; i++) {
		    AM[i]=new doublereal[3]; // матрица СЛАУ
	    }
		doublereal* b = NULL;
		b=new doublereal[2]; // правая часть

	    for (integer i=0; i<2; i++) {
		    for (integer j=0; j<2; j++) {
			    AM[i][j]=0.0; // инициализация
		    }
		    b[i]=0.0;
		    a[i]=0.0;
	    }

	    AM[0][0]=4.0;
	    AM[0][1]=x1+x2+x3+x4;
	    b[0]=f1+f2+f3+f4;
	    AM[1][0]=x1+x2+x3+x4;
	    AM[1][1]=x1*x1+x2*x2+x3*x3+x4*x4;
	    b[1]=f1*x1+f2*x2+f3*x3+f4*x4;
	

	    integer nodes=2;
	    eqsolve_simple_gauss(AM, nodes, b, a); // решение СЛАУ методом Гаусса.

	    // Освобождение памяти
		if (AM != NULL) {
			for (integer i = 0; i < 2; i++) {
				if (AM[i] != NULL) {
					delete[] AM[i];
				}
			}
			delete[] AM;
		}
		if (b != NULL) {
			delete[] b;
		}

	    //f(x)=a[0]+a[1]*x;

	    //Первая производная в точке xg равна a[1];
	    DFDXgran=a[1];

		if (a != NULL) {
			delete[] a;
		}

	}

	if (imetod==cubic_parabola) {
		// По четырём точкам строим кубическую параболу и аналитически находим её производную в нужной точке xg.

		doublereal* a = NULL;
		a=new doublereal[4];
		doublereal** AM = NULL;
		AM=new doublereal*[4];
	    for (integer i=0; i<4; i++) {
		    AM[i]=new doublereal[4]; // матрица СЛАУ
	    }
		doublereal* b = NULL;
		b=new doublereal[4]; // правая часть

	    for (integer i=0; i<4; i++) {
		    for (integer j=0; j<4; j++) {
			    AM[i][j]=0.0; // инициализация
		    }
		    b[i]=0.0;
		    a[i]=0.0;
	    }

	     AM[0][0]=1.0;
	     AM[0][1]=x1;
	     AM[0][2]=x1*x1;
	     AM[0][3]=x1*x1*x1;
	     b[0]=f1;
	     AM[1][0]=1.0;
	     AM[1][1]=x2;
	     AM[1][2]=x2*x2;
		 AM[1][3]=x2*x2*x2;
	     b[1]=f2;
		 AM[2][0]=1.0;
	     AM[2][1]=x3;
	     AM[2][2]=x3*x3;
		 AM[2][3]=x3*x3*x3;
	     b[2]=f3;
		 AM[3][0]=1.0;
	     AM[3][1]=x4;
	     AM[3][2]=x4*x4;
		 AM[3][3]=x4*x4*x4;
	     b[3]=f4;
	

	     integer nodes=4;
	     eqsolve_simple_gauss(AM, nodes, b, a); // решение СЛАУ методом Гаусса.

	     // Освобождение памяти
		 if (AM != NULL) {
			 for (integer i = 0; i < 3; i++) {
				 if (AM[i] != NULL) {
					 delete AM[i];
				 }
			 }
			 delete[] AM;
		 }
		 if (b != NULL) {
			 delete[] b;
		 }

	     //f(x)=a[0]+a[1]*x+a[2]*x^2+a[3]*x^3;

	     //Первая производная в точке xg равна 2.0*a[2]*xg+a[1];
	     DFDXgran=3.0*a[3]*xg*xg+2.0*a[2]*xg+a[1];

		 if (a != NULL) {
			 delete[] a;
		 }

	}
	

	return DFDXgran; // производная в точности на грани КО.

} // DFDXg

// Производная от искомой величины в точности на грани КО.
// позиция грани xg. По трём точкам ближайшим к границе.
doublereal DFDXg2(doublereal x1, doublereal x2, doublereal x3, doublereal xg,
	      doublereal f1, doublereal f2, doublereal f3) 
{
	doublereal h1, h2;
	h1=x2-x1; h2=x3-x2;
	doublereal D1, D2, D;
	D=h1*h2*(h1+h2);
	D1=(f3-f2)*h1+(f1-f2)*h2;
	D2=(f3-f2)*h1*h1-(f1-f2)*h2*h2;

	doublereal alpha, beta;
	alpha=D1/D;
	beta=D2/D;

	return 2.0*alpha*xg+beta;

} // DFDXg2

// В книге Г.З. Гарбера на стр 208-210 предложена формула
// для вычисления первой производной на неравномерной сетке с точностью O(h^2).
// По-видимому этот способ будет точнее чем способ предложенный первоначально и 
// его можно будет использовать при коррекции скорости во внутренних контрольных объёмах.
// Математически способ состоит в построении параболы по трём точкам а затем нахождении её 
// первой производной.
doublereal rgradF(doublereal Fback, doublereal Fcenter, doublereal Fforvard, doublereal hback, doublereal hforvard) {
	doublereal r=0.0;

	// Геометрически это выглядит следующим образом:
	//  Fback Fcenter Fforvard
	// -hback    0    hforvard

	// точность данной формулы O(h^2)
	r=((Fforvard-Fcenter)*hback*hback-(Fback-Fcenter)*hforvard*hforvard)/(hback*hforvard*(hback+hforvard));

	return r;
} // rgradF

// линейная интерполяция по двум узлам
doublereal my_linear_interpolation(char chdirect, doublereal PP, doublereal PbG, doublereal posP,
	                         doublereal posbG, doublereal posvirtual) {
	doublereal a=0.0,b=0.0;
	doublereal r=0.0;

	switch (chdirect) {
    	case '+' : a=(PP-PbG)/(posP-posbG);
	    	       b=(PbG*posP-PP*posbG)/(posP-posbG);
		           break;
	    case '-' : a=(PbG-PP)/(posbG-posP);
		    	   b=(PP*posbG-PbG*posP)/(posbG-posP);
		           break;
	}

	r=a*posvirtual+b; // linear interpolation
	return r;

} // my_linear_interpolation

doublereal my_quadratic_interpolation(char chdirect, doublereal PbGG, doublereal PbG, doublereal PP,
	                            doublereal posbGG, doublereal posbG, doublereal posP, doublereal posvirtual) {

    // Для интерполляции давления за границами расчётной области
    doublereal **B3x3=NULL;
	B3x3=new doublereal*[3];
	integer l=0; // счётчик цикла for
	for (l=0; l<3; l++) B3x3[l]=new doublereal[3];

	doublereal posbGG2=0.0, posbG2=0.0, posP2=0.0; // квадраты координат позиций
	posbGG2=posbGG*posbGG; posbG2=posbG*posbG; posP2=posP*posP;

	doublereal r=0.0; // здесь будет формироваться возвращаемое значение

	switch (chdirect) {
	case '+' : B3x3[0][0]=1.0; B3x3[0][1]=posbGG; B3x3[0][2]=posbGG2;
	           B3x3[1][0]=1.0; B3x3[1][1]=posbG; B3x3[1][2]=posbG2;
	           B3x3[2][0]=1.0; B3x3[2][1]=posP; B3x3[2][2]=posP2;

	           inverse_matrix_simple(B3x3, 3, false); // обращает матрицу

			   r=1.0*(B3x3[0][0]*PbGG+B3x3[0][1]*PbG+B3x3[0][2]*PP);
	           r+=posvirtual*(B3x3[1][0]*PbGG+B3x3[1][1]*PbG+B3x3[1][2]*PP);
	           r+=posvirtual*posvirtual*(B3x3[2][0]*PbGG+B3x3[2][1]*PbG+B3x3[2][2]*PP);
		       break;
	case '-' : B3x3[0][0]=1.0; B3x3[0][1]=posP; B3x3[0][2]=posP2;
	           B3x3[1][0]=1.0; B3x3[1][1]=posbG; B3x3[1][2]=posbG2;
	           B3x3[2][0]=1.0; B3x3[2][1]=posbGG; B3x3[2][2]=posbGG2;

	           inverse_matrix_simple(B3x3, 3, false); // обращает матрицу

               r=1.0*(B3x3[0][0]*PP+B3x3[0][1]*PbG+B3x3[0][2]*PbGG);
	           r+=posvirtual*(B3x3[1][0]*PP+B3x3[1][1]*PbG+B3x3[1][2]*PbGG);
	           r+=posvirtual*posvirtual*(B3x3[2][0]*PP+B3x3[2][1]*PbG+B3x3[2][2]*PbGG);
		       break;
	}
		
	if (B3x3 != NULL) {
		for (l = 0; l < 3; l++) {
			if (B3x3[l] != NULL) {
				delete[] B3x3[l];
				B3x3[l] = NULL;
			}
		}
		delete[] B3x3;
		B3x3 = NULL;
	}

	return r;

} // my_quadratic_interpolation

#endif
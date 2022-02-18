// Файл my_approx_convective2.c 
// содержит функции для аппроксимации 
// конвективного и диффузионного члена.

#ifndef MY_APPROX_CONVECTIVE2_C
#define MY_APPROX_CONVECTIVE2_C 1

//#define doublereal double // float

// схемы для аппроксимации конвекции-диффузии
#define CR2 100 // Центральные разности
#define UDS2 101 // Противопоточная первого порядка
#define EXP2 102 // экспоненциальная схема (точная)
#define KUD 103 // схема предложенная в диссертации Кудинова Павла Ивановича

// Higher-Order Schemes
// равномерная сетка
// Linear Higher-Order Schemes
#define QUICK 1000 // схема Леонарда: здесь просто объявляется, а реализуется в модуле my_elmatr_quad_f3D.c
#define LUS 1001 // Linear Upwind scheme
#define CUS 1002 // Cubic upwind scheme
// Non-Linear Higher-Order Schemes
#define SMART 1003
#define H_QUICK 1004
#define UMIST 1005
#define CHARM 1006
#define MUSCL 1007
#define VAN_LEER_HARMONIC 1008
#define OSPRE 1009
#define VAN_ALBADA 1010
#define SUPERBEE 1011
#define MINMOD 1012
#define H_CUS 1013
#define KOREN 1014
#define FROMM 1015
// неравномерная сетка
// неограниченные схемы
#define UNEVENQUICK 1016 // схема Леонарда для неравномерной сетки
// ограниченные схемы
/* // 3 АВГУСТА 2015 в связи с тем что эти схемы стали доступны через GUI.
#define UNEVEN_MUSCL 1017  // van Leer (1977)
#define UNEVEN_SOUCUP 1018 // MINMOD
#define UNEVEN_HLPA 1019
#define UNEVEN_SMART 1020 // Gaskell and Lau (1988)
#define UNEVEN_WACEB 1021
#define UNEVEN_SMARTER 1022
#define UNEVEN_STOIC 1023 // Darwish (1993)
#define UNEVEN_CLAM 1024
#define UNEVEN_OSHER 1025 // Chakravarthy and Osher (1983)
#define UNEVEN_VONOS 1026
#define UNEVEN_LPPA 1027
#define UNEVEN_EXPONENTIAL 1028
#define UNEVEN_SUPER_C 1029
#define UNEVEN_ISNAS 1030
#define UNEVEN_CUBISTA 1031
#define UNEVEN_GAMMA 1032 // схема с параметром beta_m
#define UNEVEN_COPLA 1033 // 1 08 2015
#define UNEVEN_SECBC 1034 // 2 08 2015 Yu et al., (2001b) Сингапур, Малазия.
#define UNEVEN_SGSD 1035 // 3 08 2015 Li and Tao (2002)
*/
// стабильная заглушка (центральные разности для диффузионного члена, UDS схема для конвективного члена).
// Выставляется в крайнем случае, когда все другие варианты не работают.
#define STABILITY 2000

// аппроксимация конвективного члена на существенно неравномерной сетке.
// Реализовано 24 февраля 2012 года.
// по мотивам диссертационной работы Кудинова Павла Ивановича.
doublereal fC(doublereal Pe, integer isheme, bool buneven_grid, doublereal fgplus) {
	doublereal r=1.0;
	const doublereal Pemax=12.0;

	switch (isheme) {
		case CR2: // центральные разности:
			       // Осторожно ! работает только при малых числах Пекле.
			       // есть публикации где центральные разности используются при 
                   // LES моделировании на очень подробных сетках.
			       // Для существенно неравномерных сеток диапазон рабочих 
			       // чисел Пекле не симметричен относительно нуля.
			       if (buneven_grid) {
                       r=fgplus; // для неравномерной сетке. 
				   }
				   else r=0.5; // для равномерной сетки. 
			       break;
		case UDS2: // противопоточная схема первого порядка.
			        if (Pe<0.0) r=1.0; else r=0.0;
			        break;
		case EXP2: // экспоненциальная (точная) схема.
			        if (fabs(Pe) < 0.01) {
						if (buneven_grid) {
                            r=fgplus; // для неравномерной сетке. 
						}
						else r=0.5; // равномерная сетка
					}
					else {
						if (buneven_grid) {
							r=expm1(fgplus*Pe)/expm1(Pe); // для неравномерной сетки.
						}
						else r=expm1(0.5*Pe)/expm1(Pe);   // для равномерной сетки.
					}
					break;
		case KUD: // схема Кудинова Павла Ивановича.
			       if (buneven_grid) {
					   if ((Pe<=2.0*Pemax*fgplus)&&(Pe>=-2.0*Pemax*(1.0-fgplus))) {
						   r=fgplus*(1.0-Pe/(2.0*Pemax*fgplus));
					   }
					   else {
						   r=0.5*(1.0-Pe/fabs(Pe)); // максимальное значение 1.0.
					   }
				   }
				   else {
					   // равномерная сетка.
					   if (fabs(Pe)<=Pemax) {
					       r=0.5*(1.0-Pe/Pemax); 
				       }
				       else{
					       r=0.5*(1.0-Pe/fabs(Pe));
				       }
				   }
			       break;
		case STABILITY: // противопоточная схема первого порядка.
			      if (Pe<0.0) r=1.0; else r=0.0;
			      break;
		default: //  по умолчанию используется
			      // схема Кудинова Павла Ивановича.
			       if (buneven_grid) {
					   if ((Pe<=2.0*Pemax*fgplus)&&(Pe>=-2.0*Pemax*(1.0-fgplus))) {
						   r=fgplus*(1.0-Pe/(2.0*Pemax*fgplus));
					   }
					   else {
						   r=0.5*(1.0-Pe/fabs(Pe)); // максимальное значение 1.0.
					   }
				   }
				   else {
					   // равномерная сетка.
					   if (fabs(Pe)<=Pemax) {
					       r=0.5*(1.0-Pe/Pemax); 
				       }
				       else{
					       r=0.5*(1.0-Pe/fabs(Pe));
				       }
				   } 
			      break;
	}
	return (r);
} // аппроксимация конвекции


// Вычисление экспоненты с большим значением аргумента вызывает появление положительной
// бесконечности, а неопределённость вида бесконечность на бесконечность даёт NaN - Not a Number.
doublereal fDbug(doublereal Pe, integer isheme, bool buneven_grid, doublereal fgplus) {

	// если buneven_grid истина то расчёт ведётся на неравномерной сетке.
	// fgplus - учёт неравномернойсти сетки.

	doublereal r=1.0;
	const doublereal Pemax=12.0;

	// 1. В диссертации Кудинова Павла Ивановича предложено на 
	// на основе аналитического решения одномерной задачи конвекции-диффузии
	// ограничить вклад диффузионного члена для больших чисел Пекле. Это
	// соответствует схеме KUD.
	// 2. В работе И.К. Жарова, Г.В. Кузнецов и др. Исследование взаимодействия
	// импактной струи с поверхностью преграды сложной формы предложено
	// использовать схему Леонарда QUICK для квазиравномерных сеток, причём
	// вклад диффузионной составляющей никак не ограничивается. Это реализовано
	// в схеме QUICK (совпадает для диффузионного члена с центральными разностями CR2).
	// 3. В работе SIMPLE METHOD FOR THE SOLUTION OF INCOMPRESSIBLE FLOW ON NON-STAGGERED GRIDS
	// I. Sezai - Eastern Mediterranean University, January, 2011 также предлагается не ограничивать вклад
	// диффузионного члена. Это реализовано у них для всех схем.
	// Для пунктов 2 и 3 используйте STABILITY. Если согласится с 1 то используйте KUD или default.


	switch (isheme) {
		case CR2: r=1.0; break; // Центральные разности (стабильная схема для диффузии, учитывается полный вклад диффузионного члена).
		case UDS2: // противопоточная схема первого порядка.
			        if (buneven_grid) {
						r=Pemax*exp(Pemax*fgplus)/expm1(Pemax);
					} else r=0.0674; // предел на равномерной сетке при числе Пекле Pe>10.
					break;
		case EXP2: // экспоненциальная (точная) схема.
			        if (fabs(Pe) < 0.001) r=1.0;
					else {
						if (buneven_grid) {
							// неравномерная сетка.
                            r=Pe*exp(Pe*fgplus)/expm1(Pe);
						}
						else r=(Pe*exp(0.5*Pe))/expm1(Pe);
					}
			        break;
		case KUD: // Схема Кудинова Павла Ивановича (аппроксимация экспоненциальной зависимости).
			       // годится только для равномерной сетки. 
			       // для неравномерной сетки непригодна.
			       //-->//if (fabs(Pe)<=Pemax) r=1.0-(Pe*Pe)/(Pemax*Pemax); else r=0.0;
				   // Причина в том, что для неравномерной сетки требуется аппроксимировать семейство кривых (для различных fgplus).
				   // Аппроксимация производится с целью уменьшить число арифметических операций (не вычислять экспоненты). Но
				   // в данном случае удалось найти аппроксимацию кривой только для случая fgplus=0.5 (равномерная сетка). Для остальных fgplus!=0.5
				   // кривые аппроксимировать не удалось по причине их сложного поведения. Возможно для ускорения вычислений здесь поможет табличный
				   // способ аппроксимации.
				   // Вывод: здесь будет автоматически осуществлён переход на экспоненциальную схему: EXP2.
			       if (fabs(Pe) < 0.001) r=1.0;
				   else {
						if (buneven_grid) {
							// неравномерная сетка.
                            r=Pe*exp(Pe*fgplus)/expm1(Pe);
						}
						else r=(Pe*exp(0.5*Pe))/expm1(Pe);
					}
			       break; 
		case STABILITY: r=1.0; break; // диффузионный член учитывается полностью в соответствии с центрально разностной схемой.
		default: // по умолчанию используется экспоненциальная схема на неравномерной сетке, 
			      // т.к. не получилось аппроксимировать семейство кривых (для разных fgplus) 
			      // и добиться уменьшения числа арифметических операций. Кривые в семействе ведут себя слишком по разному.
			      if (fabs(Pe) < 0.001) r=1.0;
				  else {
					 if (buneven_grid) {
						// неравномерная сетка.
                        r=Pe*exp(Pe*fgplus)/expm1(Pe);
					 }
					 else r=(Pe*exp(0.5*Pe))/expm1(Pe); // равномерная сетка.
				  }
			      break; 
	}
	return (r);
} // аппроксимация диффузии

// отладочный диффузионный поток
doublereal fD(doublereal Pe, integer isheme, bool buneven_grid, doublereal fgplus) {

	// если buneven_grid истина то расчёт ведётся на неравномерной сетке.
	// fgplus - учёт неравномерности сетки.

	doublereal r=1.0;
	const doublereal Pemax=12.0;

	// 1. В диссертации Кудинова Павла Ивановича предложено на 
	// на основе аналитического решения одномерной задачи конвекции-диффузии
	// ограничить вклад диффузионного члена для больших чисел Пекле. Это
	// соответствует схеме KUD.
	// 2. В работе И.К. Жарова, Г.В. Кузнецов и др. Исследование взаимодействия
	// импактной струи с поверхностью преграды сложной формы предложено
	// использовать схему Леонарда QUICK для квазиравномерных сеток, причём
	// вклад диффузионной составляющей никак не ограничивается. Это реализовано
	// в схеме QUICK (совпадает для диффузионного члена с центральными разностями CR2).
	// 3. В работе SIMPLE METHOD FOR THE SOLUTION OF INCOMPRESSIBLE FLOW ON NON-STAGGERED GRIDS
	// I. Sezai - Eastern Mediterranean University, January, 2011 также предлагается не ограничивать вклад
	// диффузионного члена. Это реализовано у них для всех схем.
	// Для пунктов 2 и 3 используйте STABILITY. Если согласится с 1 то используйте KUD или default.


	switch (isheme) {
		case CR2: r=1.0; break; // Центральные разности (стабильная схема для диффузии, учитывается полный вклад диффузионного члена).
		case UDS2: // противопоточная схема первого порядка.
			        if (buneven_grid) {
						if (fabs(Pe)>85.0) {
							// неопределённость вида бесконечность делённая на бесконечность.
							// Предельное значение отношения стремится к нулю.
							r=0.0;
						}
						else if ((fabs(Pe)>0.001)&&(fabs(Pe)<=85.0)) {
						    r=Pemax*exp(Pemax*fgplus)/expm1(Pemax);
						}
						else {
							// вычислен предел при малых значениях числа Пекле.
					        r=1.0*exp(Pe*fgplus);
						}
					} else r=0.0674; // предел на равномерной сетке при числе Пекле Pe>10.
					break;
		case EXP2: // экспоненциальная (точная) схема.
			        if (fabs(Pe) < 0.001) {
						// вычислен предел при малых значениях числа Пекле.
					    r=1.0*exp(Pe*fgplus);
					}
					else {
						if (fabs(Pe)>85) {
							// неопределённость вида бесконечность делённая на бесконечность.
							// Предельное значение отношения стремится к нулю.
							r=0.0;
						}
						else {

						     if (buneven_grid) {
							    // неравномерная сетка.
                                r=Pe*exp(Pe*fgplus)/expm1(Pe);
						     }
						     else r=(Pe*exp(0.5*Pe))/expm1(Pe);
						}
					}
			        break;
		case KUD: // Схема Кудинова Павла Ивановича (аппроксимация экспоненциальной зависимости).
			       // годится только для равномерной сетки. 
			       // для неравномерной сетки непригодна.
			       //-->//if (fabs(Pe)<=Pemax) r=1.0-(Pe*Pe)/(Pemax*Pemax); else r=0.0;
				   // Причина в том, что для неравномерной сетки требуется аппроксимировать семейство кривых (для различных fgplus).
				   // Аппроксимация производится с целью уменьшить число арифметических операций (не вычислять экспоненты). Но
				   // в данном случае удалось найти аппроксимацию кривой только для случая fgplus=0.5 (равномерная сетка). Для остальных fgplus!=0.5
				   // кривые аппроксимировать не удалось по причине их сложного поведения. Возможно для ускорения вычислений здесь поможет табличный
				   // способ аппроксимации.
				   // Вывод: здесь будет автоматически осуществлён переход на экспоненциальную схему: EXP2.
			       if (fabs(Pe) < 0.001) {
					   // вычислен предел при малых значениях числа Пекле.
					   r=1.0*exp(Pe*fgplus);
				   }
				   else {

					   if (fabs(Pe)>85) {
							// неопределённость вида бесконечность делённая на бесконечность.
							// Предельное значение отношения стремится к нулю.
							r=0.0;
						}
						else {

						    if (buneven_grid) {
								// неравномерная сетка.
                                r=Pe*exp(Pe*fgplus)/expm1(Pe);
						    }
						    else r=(Pe*exp(0.5*Pe))/expm1(Pe);
						}
					}
			       break; 
		case STABILITY: r=1.0; break; // диффузионный член учитывается полностью в соответствии с центрально разностной схемой.
		default: // по умолчанию используется экспоненциальная схема на неравномерной сетке, 
			      // т.к. не получилось аппроксимировать семейство кривых (для разных fgplus) 
			      // и добиться уменьшения числа арифметических операций. Кривые в семействе ведут себя слишком по разному.
			      if (fabs(Pe) < 0.001) {
					  // вычислен предел при малых значениях числа Пекле.
					  r=1.0*exp(Pe*fgplus);
				  }
				  else {
					  if (fabs(Pe)>85) {
							// неопределённость вида бесконечность делённая на бесконечность.
							// Предельное значение отношения стремится к нулю.
							r=0.0;
						}
						else {

					        if (buneven_grid) {
						        // неравномерная сетка.
                                r=Pe*exp(Pe*fgplus)/expm1(Pe);
					        }
					        else r=(Pe*exp(0.5*Pe))/expm1(Pe); // равномерная сетка.
						}
				  }
			      break; 
	}
	return (r);
} // аппроксимация диффузии

// схема Леонарда [1979]
// См. документ "рабочая схема QUICK" в папке с исходным кодом программы.
// Один и тот же код используется для интерполяции на всех гранях, см. таблицу из вышеприведённого файла.
doublereal workQUICK(doublereal deltaxB, doublereal deltaxC, doublereal xA, doublereal xB, doublereal xC, doublereal xD,
	           doublereal FA, doublereal FB, doublereal FC, doublereal FD, doublereal Fe) {

	/* таблица соответствия:
	*  A	B	C	D	e	+/-
	*  W	P	E	-	e	+
	*  -	    P	E	EE  e   -
	*  WW   W	P	-	w	+
	*  -	    W	P	E	w	-
	*  S	P	N	-	n	+
	*  -	P	N	NN  n	-
	*  SS   S	P	-	s	+
	*  -	S	P	N	s	-
	*  B	P	T	-	t	+
	*  -	P	T	TT  t	-
	*  BB   B	P	-	b	+
	*  -	B	P	T	b	-
	*/

	doublereal r=0.0;
	if (Fe>=0.0) {
		doublereal g1=0.0, g2=0.0; // инициализация
		g1=(deltaxB*(deltaxB+2.0*(xB-xA)))/(4.0*(xC-xB)*(xC-xA));
		g2=((deltaxB*(2.0*(xC-xB)-deltaxB))/(4.0*(xC-xA)*(xB-xA)));
		r=g1*FC-g2*FA+(1.0-g1+g2)*FB;
	}
	else 
	{
		doublereal g3=0.0, g4=0.0; // инициализация
		g3=((deltaxC*(deltaxC+2.0*(xD-xC)))/(4.0*(xC-xB)*(xD-xB)));
		g4=((deltaxC*(2.0*(xC-xB)-deltaxC))/(4.0*(xD-xC)*(xD-xB)));
		r=g3*FB-g4*FD+(1.0-g3+g4)*FC;
	}
	return (r);
} // workQUICK

// минимум из трёх чисел
doublereal fmin3(doublereal fA, doublereal fB, doublereal fC) 
{
	return fmin(fmin(fA, fB), fC);
} // fmin3

// Возвращает минимум из четырёх чисел
doublereal fmin4(doublereal fA, doublereal fB, doublereal fC, doublereal fD)
{
	return fmin(fmin3(fA, fB, fC), fD);
} // fmin4

// минимум из трёх чисел
doublereal fmax3(doublereal fA, doublereal fB, doublereal fC)
{
	return fmax(fmax(fA, fB), fC);
} // fmax3

// см. дискретизация конвективных потоков в уравнениях Навье-Стокса
// на основе разностных схем высокой разрешающей способности. 
// К.Н. Волков. стр. 135.
// Вычислительные методы и программирование. 2004. Т. 5.
doublereal linear_flux_limiter(doublereal kappa, doublereal r)
{
	return (0.5 * ((1.0 + kappa) * r + (1.0 - kappa)));
} // linear_flux_limiter

// по мотивам программы проф. Б. Сполдинга PHOENICS
// К сожалению это многообразие пригодно лишь для равномерной сетки.
doublereal limiter_function(integer ischeme, doublereal r) {
	// Дискретизация конвективных потоков в уравнениях Навье-Стокса
	// на основе разностных схем высокой разрешающей способности К.Н. Волков.
	// Вычислительные методы и программирование. 2004. Т. 5.

	doublereal Br=2.0*(r+fabs(r))/(r+3.0); // H_QUICK - схема по дефолту
	doublereal Konst=0.0;

	switch (ischeme) {
	  case QUICK: Konst=0.5; Br=linear_flux_limiter(Konst, r); break;
	  case LUS: Konst=-1.0; Br=linear_flux_limiter(Konst, r); break;
	  case CUS: Konst=1.0/3.0; Br=linear_flux_limiter(Konst, r); break;
	  case SMART: Br=fmax(0.0,fmin3(2.0*r,0.75*r+0.25,4.0)); break; // ограниченная на основе QUICK
	  case H_QUICK: Br=2.0*(r+fabs(r))/(r+3.0); break; // гладкая на основе QUICK
	  case UMIST:  Br=fmax(0.0,fmin4(2.0*r,0.25+0.75*r,0.75+0.25*r,2.0)); break; // TVD на основе QUICK
	  case CHARM: if (r<=0.0) Br=0.0; else Br=r*(3.0*r+1.0)/((r+1)*(r+1)); break; // ограниченная на основе QUICK
	  case MUSCL: Br=fmax(0.0, fmin3(2.0*r,0.5+0.5*r,2));  break; // TVD на базе Fromm
	  case VAN_LEER_HARMONIC: Br=(r+fabs(r))/(r+1.0); break; // TVD на базе Fromm
	  case OSPRE: Br=3.0*(r*r+r)/(2.0*(r*r+r+1.0)); break; // гладкая на основе Fromm
	  case VAN_ALBADA: Br=(r*r+r)/(r*r+1.0); break; // TVD на основе Fromm
	  case SUPERBEE: Br=fmax3(0.0,fmin(2.0*r,1.0),fmin(r,2.0));  break; // TVD
	  case MINMOD: Br=fmax(0.0,fmin(r,1.0)); break; // (SOUCUP) TVD
	  case H_CUS: Br=1.5*(r+fabs(r))/(r+2.0); break; // (HCUDS) гладкая на основе CUDS
	  case KOREN: Br=fmax(0.0,fmin3(2.0*r,2.0*r/3.0+1.0/3.0,2.0)); break; // TVD на основе CUDS
	  case FROMM: Konst=0.0; Br=0.5*((1.0+Konst)*r+(1.0-Konst)); break;
	  case UNEVENQUICK: Konst=0.5; Br=linear_flux_limiter(Konst, r); break; // здесь дублирует схему QUICK
	  default: Konst=0.5; Br=linear_flux_limiter(Konst, r); break; // QUICK
	}

	return Br;
} // limiter_function

// по мотивам программы проф. Б. Сполдинга PHOENICS
// Возвращает значение искомой величины на границе контрольного объёма.
doublereal cell_face_value_local(integer ischeme, doublereal Fc, doublereal Fd, doublereal Fu) {

	doublereal Ff=Fc; // UDS Br==0 // противопоточная схема.
	bool bcontinue=true;

	if (fabs(Fc-Fu)<1e-30) bcontinue=false;

	if (bcontinue) {

	   doublereal r=(Fd-Fc)/(Fc-Fu);
	

	   switch (ischeme) {
	     case UDS: Ff=Fc; break; // UDS
	     default: Ff=Fc+0.5*limiter_function(ischeme,r)*(Fc-Fu); break; // Higher-Order Scheme
	   }
	}

	return Ff;
} // cell_face_value

// по мотивам программы проф. Сполдинга PHOENICS
// Возвращает значение искомой величины на границе контрольного объёма.
doublereal cell_face_value_global(integer ischeme, doublereal uf, doublereal fa, doublereal fb,
	                        doublereal fc, doublereal fd) {
	// uf - скорость на грани контрольного объёма.
    // fa, fb, fc, fd - значение искомой величины (с предыдущей итерации алгоритма SIMPLE)
	// в центрах контрольного объёма.

    /*  sign(uf)	U	C	f	D 
	*		+		WW	W	f	P
    *
	*	sign(uf)	D	f	C	U
	*		-		W	f	P	E
	*
	*   f - face - поверхность на которой нужно получить искомую функцию (грань контрольного объёма).
	*/

	doublereal ff=0.0;
	doublereal Fc=0.0, Fd=0.0, Fu=0.0;

	if (uf >= 0.0) {
		Fc=fb;
		Fd=fc;
		Fu=fa;
	}
	else {
		Fc=fc;
		Fd=fb;
		Fu=fd;
	}

	ff=cell_face_value_local(ischeme, Fc, Fd, Fu);

	return ff;

} // cell_face_value_global

// Li and Tao(2002)
// 3 08 2015
doublereal UNEVEN_SGSD_SCHEME(doublereal Pe, doublereal f_C, doublereal xQ, doublereal yQ) {
	// Pe - сеточное число Пекле как в схеме С.Патанкара.
	doublereal f_f;
	doublereal beta=2.0/(2.0+fabs(Pe));
	f_f=beta*(((yQ-xQ)/(1.0-xQ))+(((yQ-1.0)*f_C)/(xQ-1.0)))+(1.0-beta)*yQ*f_C/xQ;
	return (f_f);
} // UNEVEN_SGSD

// Yu et. al(2001b)
// 2 08 2015
doublereal UNEVEN_SECBC_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {
	doublereal f_f=0.5*f_C;
	doublereal border=((xQ-yQ)*xQ)/(xQ-yQ+2.0*yQ*yQ-4.0*yQ*xQ);
	if ((f_C>=0.0)&&(f_C<=border)) {
		f_f=(yQ*(1.0-3.0*xQ+2.0*yQ)*f_C)/(xQ*(xQ-1));
	}
	if ((f_C>=border)&&(f_C<=1)) {
		doublereal aw=(yQ-xQ)/(1.0-xQ);
		doublereal bw=(yQ-1.0)/(xQ-1.0);
		f_f=aw+bw*f_C;
	}
	return (f_f);
} // UNEVEN_SECBC

// схема MUSCL на неравномерной сетке. Базовая схема Fromm.
// схема обладает TVD свойством. А поскольку это TVD схема,
// то она также является и ограниченной.
doublereal UNEVEN_MUSCL_SCHEME(doublereal f_C, doublereal xQ,  doublereal yQ) {

	// MUSCL - монотонная противопоточная схема для законов сохранения.
	// Monotonic Upwind Scheme for Conservation Laws.
	// Обладает вторым порядком точности.

	// ссылка: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.
	// Внимание в работе К.Н.Волкова данная схема записана с ошибкой.
	// Здесь ошибки исправлены (диаграмма нормализованных переменных
	// проверена в программе maple на предмет непрерывности).

	doublereal f_f=1.0; // инициализация.
	// проверено в Delphi 29 февраля 2012.

	if ((f_C>=0.0)&&(f_C<0.5*xQ)) {
		f_f=((2.0*yQ-xQ)/(xQ))*f_C;
	}
	else if ((f_C>=0.5*xQ)&&(f_C<1.0+xQ-yQ)) {
		f_f=f_C+(yQ-xQ);
	}
	else if ((f_C>=1.0+xQ-yQ)&&(f_C<=1.0)) {
		f_f=1.0;
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_MUSCL 

// схема COPLA на неравномерной сетке. 
doublereal UNEVEN_COPLA_SCHEME(doublereal f_C, doublereal xQ,  doublereal yQ, doublereal sQ) {

	// COPLA
	// Combination of Piecewise Linear Approximation.

	// ссылка: 

	doublereal f_f=1.0; // инициализация.
	// проверено в Delphi 29 февраля 2012.

	if ((f_C>=0.0)&&(f_C<=0.5*xQ)) {
		doublereal aw=0.0;
		doublereal bw=(2.0*yQ-sQ*xQ)/xQ;
		f_f=aw+bw*f_C;
	}
	else if ((f_C>=0.5*xQ)&&(f_C<=1.5*xQ)) {
		doublereal cw=yQ-sQ*xQ;
		doublereal dw=sQ;
		f_f=cw+dw*f_C;
	}
	else if ((f_C>=1.5*xQ)&&(f_C<=1.0)) {
		doublereal ew=(3.0*xQ-2.0*yQ-sQ*xQ)/(3.0*xQ-2.0);
		doublereal fw=(2.0*yQ+sQ*xQ-2.0)/(3.0*xQ-2.0);
		f_f=ew+fw*f_C;
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_COPLA

// схема SOUCUP (MINMOD) для неравномерной сетки. 
doublereal UNEVEN_SOUCUP_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {
	// Схема SOUCUP, использующая кусочно-линейную интерполяцию искомой 
	// функции, объединяет схемы UDS, CDS-2 и LUDS. Переключение от одной 
	// схемы к другой контролируется критерием конвективной ограниченности.
	// Схема SOUCUP удовлетворяет условию TVD. А поскольку это TVD схема,
    // то она также является и ограниченной.

	// SOUCUP Second-Order Upwind Central difference-first order UPwind
	// Схема с разностями против потока и центральными разностями. 
	// Обеспечивает второй порядок точности.

	// Реализация согласовано совпадает с работой:
	// Normalized Variable and Space Formulation
	// Methodology for High-Resolution Schemes.
	// M.S.Darwish and F.H.Moukalled (год не указан, не ранее 1993 года).

	doublereal f_f=1.0; // инициализация.
	// проверено в Delphi 29 февраля 2012.
	

	if ((f_C>=0.0)&&(f_C<xQ)) {
		doublereal af, bf;
		af=0.0;
	    bf=yQ/xQ;

		// эта часть относится к LUDS
		// LUDS - противопоточная второго порядка.
		f_f=af+bf*f_C;
	}
	else if ((f_C>=xQ)&&(f_C<=1.0)) {
		doublereal cf, df;
		cf=(yQ-xQ)/(1.0-xQ);
	    df=(1.0-yQ)/(1.0-xQ);

		// Эта часть относится к центральным разностям второго порядка CDS-2
		f_f=cf+df*f_C;
	}
	else f_f=f_C; // чистая противопоточная часть 1-ого порядка.

	return (f_f);
} // UNEVEN_SOUCUP_SCHEME

// схема HLPA для неравномерной сетки
doublereal UNEVEN_HLPA_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	// В схеме HLPA используется линейная и параболическая интерполяция искомой функции,
	// при этом она удовлетворяет условию TVD.  А поскольку это TVD схема,
    // то она также является и ограниченной.

	// HLPA гибридная схема с линейно-параболической аппроксимацией.
	// Hybrid Linear/Parabolic Approximation.
	// HLPA схема обеспечивает второй порядок точности.

	// Схема проверена в Delphi 12 марта 2012 года.

	doublereal f_f=1.0; // инициализация.

	if ((f_C>=0.0)&&(f_C<=1.0)) {
		doublereal af, bf, cf;
	    doublereal xQ2=xQ*xQ; // для ускорения вычислений.
	    af=0.0; bf=(yQ-xQ2)/(xQ-xQ2); cf=(xQ-yQ)/(xQ-xQ2);

		f_f=af+bf*f_C+cf*f_C*f_C;
	}
	else f_f=f_C;

    return (f_f);
} // UNEVEN_HLPA_SCHEME

// схема SMART для неравномерной сетки
doublereal UNEVEN_SMART_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	//  При построении схемы SMART в качестве базовой линейной схемы используется схема QUICK. 
	//  Данная схема не обладает TVD свойством но является ограниченной. 
    // 

	// SMART Sharp and Monotonic Algorithm for Realistic Transport
	// Быстрый и монотонный алгоритм для реалистичного переноса.
	// SMART схема обеспечивает третий порядок точности.

	// Исправлена информация в статье К.Н.Волкова. Правильный материал взят
	// из статьи 1: "Normalized Variable and Space Formulation Methodology for High-Resolution Schemes."
	// M.S. Darwish and F.H.Moukalled.

	// ссылка 2: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.

	// Схема проверена в Delphi 12 марта 2012 года.
	// Аналитическая проверка на согласование в точках x1 и x2 
	// также успешно проведена 12 марта 2012 года.

	doublereal f_f=1.0; // инициализация.

	
	doublereal xQ2=xQ*xQ; // для ускорения вычислений.
	doublereal yQ2=yQ*yQ;

	doublereal x1=xQ/3.0;
	doublereal x2=xQ*(1.0-xQ+yQ)/yQ;

	if ((f_C>=0.0)&&(f_C<x1)) {
		doublereal af, bf;
		af=0.0; 
	    bf=(yQ-3.0*xQ*yQ+2.0*yQ2)/(xQ-xQ2);

		f_f=af+bf*f_C;
	}
	else if ((f_C>=x1)&&(f_C<x2)) {
		doublereal cf, df;
		cf=(-xQ*yQ+yQ2)/(1.0-xQ);
	    df=(yQ-yQ2)/(xQ-xQ2);

		f_f=cf+df*f_C;
	}
	else if ((f_C>=x2)&&(f_C<=1.0)) {
		f_f=1.0;
	}
	else f_f=f_C; // UDS

    return (f_f);
} // UNEVEN_SMART_SCHEME

// схема WACEB для неравномерной сетки.
doublereal UNEVEN_WACEB_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	// Схема WACEB представляет собой одну из модификаций схемы SMART
	// отличается от неё лишь в области малых положительных значений f_C.
	// Схема WACEB - удовлетворяет условию TVD.

	// Аналитическая проверка на согласование в точках x1 и x2 
	// успешно проведена 12 марта 2012 года.

	doublereal f_f=1.0; // инициализация.
		
	
	doublereal x1=xQ*yQ*(yQ-xQ)/(2.0*xQ*(1.0-xQ)-yQ*(1-yQ));
	doublereal x2=xQ*(1.0-xQ+yQ)/yQ;

	if ((f_C>=0.0)&&(f_C<x1)) {
		doublereal af, bf;
		af=0.0;
	    bf=2.0;

		f_f=af+bf*f_C;
	}
	else if ((f_C>=x1)&&(f_C<x2)) {
		doublereal xQ2=xQ*xQ; // для ускорения вычислений.
	    doublereal yQ2=yQ*yQ;

		doublereal cf, df;
		cf=(yQ2-xQ*yQ)/(1.0-xQ);
	    df=(yQ-yQ2)/(xQ-xQ2);

		f_f=cf+df*f_C;
	}
	else if ((f_C>=x2)&&(f_C<=1.0)) {
		f_f=1.0;
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_WACEB_SCHEME

// схема SMARTER на неравномерной сетке.
doublereal UNEVEN_SMARTER_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ, doublereal sQ) {
	// Модификация разностной схемы SMART, получившая в некоторых работах 
	// название SMARTER (SMART Efficiently Revised), имеет более хорошие показатели,
	// чем базовая схема. Касательная к кривой, воспроизводящей поведение искомой функции
	// между соседними узлами, имеет в точке Q такой же угол наклона, что и в схеме QUICK.

	// схема SMARTER обладает третьим порядком точности.


	doublereal f_f=1.0; // инициализация.

    doublereal af, bf, cf, df;
	doublereal xQ2=xQ*xQ; // для ускорения вычислений.
	//doublereal yQ2=yQ*yQ;
	doublereal xQ3=xQ2*xQ;
	//doublereal yQ3=yQ2*yQ;
	doublereal zQ=1.0/((xQ-xQ2)*(xQ-xQ2));

	af=0.0;
	bf=(xQ2*xQ2+sQ*(xQ3-xQ2)+yQ*(2.0*xQ-3.0*xQ2))*zQ;
	cf=(-2.0*xQ3+sQ*(xQ-xQ3)+yQ*(3.0*xQ2-1.0))*zQ;
	df=(xQ2+sQ*(xQ2-xQ)+yQ*(1.0-2.0*xQ))*zQ;

	if ((f_C>=0.0)&&(f_C<=1.0)) {
		f_f=af+bf*f_C+cf*f_C*f_C+df*f_C*f_C*f_C;
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_SMARTER_SCHEME

// схема STOIC на неравномерной сетке
doublereal UNEVEN_STOIC_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {
	// ссылка: Normalized Variable and Space Formulation
	// Methodology for High-Resolution Schemes.
	// M.S.Darwish and F.H.Moukalled (год не указан, не ранее 1993 года).

	doublereal f_f=1.0; // инициализация.

	doublereal yQ2=yQ*yQ;
	doublereal x1, x2;
	x1=((xQ-yQ)*xQ)/(xQ+yQ+2.0*yQ2-4.0*yQ*xQ);
	x2=(xQ*(1.0+yQ-xQ))/yQ;

	if ((f_C>0)&&(f_C<x1)) {
		f_f=-(yQ*(1.0-3.0*xQ+2.0*yQ)*f_C)/(xQ*(xQ-1.0));
	}
	else if ((f_C>=x1)&&(f_C<xQ)) {
		f_f=((xQ-yQ)/(xQ-1.0))+((yQ-1.0)/(xQ-1.0))*f_C;
	}
	else if ((f_C>=xQ)&&(f_C<x2)) {
		f_f=((yQ*(yQ-xQ))/(1.0-xQ))+((yQ*(yQ-1.0))/(xQ*(xQ-1.0)))*f_C;
	}
	else if ((f_C>=x2)&&(f_C<1.0)) {
		f_f=1.0;
	} else f_f=f_C;

	return (f_f);
} // UNEVEN_STOIC_SCHEME

// схема CLAM на неравномерной сетке
doublereal UNEVEN_CLAM_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {
	// ссылка 1: Normalized Variable and Space Formulation
	// Methodology for High-Resolution Schemes.
	// M.S.Darwish and F.H.Moukalled (год не указан, не ранее 1993 года).

	// ссылка 2: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.

	doublereal f_f=1.0; // инициализация.
	
	if ((f_C>0.0)&&(f_C<1.0)) {
		doublereal af, bf;
	    doublereal zQ=(xQ*(xQ-1.0));
	    af=(xQ*xQ-yQ)/zQ;
	    bf=(yQ-xQ)/zQ;

		f_f=af*f_C+bf*f_C*f_C;
	} else f_f=f_C;

	return (f_f);
} // UNEVEN_CLAM_SCHEME

// схема OSHER на неравномерной сетке
doublereal UNEVEN_OSHER_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {
	// ссылка 1: Normalized Variable and Space Formulation
	// Methodology for High-Resolution Schemes.
	// M.S.Darwish and F.H.Moukalled (год не указан, не ранее 1993 года).

	// ссылка 2: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.

	doublereal f_f=1.0; // инициализация.

	doublereal x_barrier=xQ/yQ;

	if ((f_C>0.0)&&(f_C<=x_barrier)) {
		f_f=yQ*f_C/xQ;
	} 
	else if ((f_C>x_barrier)&&(f_C<1.0)) {
		f_f=1.0;
	} 
	else f_f=f_C;

    return (f_f);
} // UNEVEN_OSHER_SCHEME

// Схема VONOS на неравномерной сетке.
doublereal UNEVEN_VONOS_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	// Схема VONOS представляет собой одну из модификаций схемы
	// SMART. В качестве базовой схемы используется схема QUICK.

	// VONOS - неосциллирующая схема переменного порядка.
	// Variable-Order Non-Oscillatory Scheme (VONOS)
	// Схема VONOS - обеспечивает третий порядок точности на неравномерной сетке.

	doublereal x1=xQ/3.0;
	doublereal x2=xQ/yQ;
	doublereal xQ2=xQ*xQ;
	doublereal yQ2=yQ*yQ;
	doublereal xyQ=xQ*yQ;
	
	doublereal f_f=1.0; // инициализация.

	if ((f_C>0.0)&&(f_C<x1)) {
		doublereal af=0.0;
	    doublereal bf=(yQ-3.0*xyQ+2.0*yQ2)/(xQ-xQ2);

		f_f=af+bf*f_C;
	}
	else if ((f_C>=x1)&&(f_C<xQ)) {
		doublereal cf=(yQ2-xyQ)/(1.0-xQ);
	    doublereal df=(yQ-yQ2)/(xQ-xQ2);

		f_f=cf+df*f_C;
	}
	else if ((f_C>=xQ)&&(f_C<x2)) {
		doublereal ef=0.0;
	    doublereal ff=yQ/xQ;

		f_f=ef+ff*f_C;
	}
	else if ((f_C>=x2)&&(f_C<1.0)) {
		f_f=1.0;
	}
	else f_f=f_C;

    return (f_f);
} // UNEVEN_VONOS_SCHEME

// Схема LPPA на неравномерной сетке.
doublereal UNEVEN_LPPA_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ, doublereal sQ) {
	// В схеме LPPA используется линейная и кусочно-параболическая интерполяция
	// искомой функции в нормализованных переменных.

	// LPPA - схема с линейной и кусочно-параболической интерполяцией
	// Linear and Piecewise/Parabolic Approximation.
	// Схема имеет третий порядок точности.

	doublereal xQ2=xQ*xQ;
	doublereal xyQ=xQ*yQ;	

    doublereal f_f=1.0; // инициализация.

	if ((f_C>=0.0)&&(f_C<xQ)) {
		doublereal af=0.0;
	    doublereal bf=(-sQ*xQ2+2.0*xyQ)/xQ2;
	    doublereal cf=(sQ*xQ-yQ)/xQ2;

		f_f=af+bf*f_C+cf*f_C*f_C;
	}
	else if ((f_C>=xQ)&&(f_C<=1.0)) {
		doublereal zQ=((1.0-xQ)*(1.0-xQ));

		doublereal df=(xQ2+sQ*(xQ2-xQ)+yQ*(1.0-2.0*xQ))/(zQ);
	    doublereal ef=(-2.0*xQ+sQ*(1.0-xQ2)+2.0*xyQ)/(zQ);
	    doublereal ff=(1.0+sQ*(xQ-1.0)-yQ)/(zQ);

		f_f=df+ef*f_C+ff*f_C*f_C;
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_LPPA_SCHEME

// Экспоненциальная схема на неравномерной сетке.
doublereal UNEVEN_EXPONENTIAL_SCHEME(doublereal f_C) {

	// ссылка: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.

	doublereal f_f=1.0; // инициализация.

	if ((f_C>0.0)&&(f_C<1.0)) {
		f_f=1.125*(1.0-exp(-2.19722*f_C));
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_EXPONENTIAL_SCHEME

// схема SUPER-C на неравномерной сетке.
doublereal UNEVEN_SUPER_C_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	// ссылка 1: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.

	// Данная схема проходит через точку (xQ,yQ).

	doublereal f_f=1.0; // инициализация.

	doublereal x1=2.0*xQ/5.0;
	doublereal x2=xQ*(1.0+yQ-xQ)/yQ;

	// Первоначальный вариант схемы по ссылке 1,
	// приведён в закомментированном виде. Этот вариант 
	// терпит разрыв первого рода в точке x1.
	// Это установлено путём построения графика в системе maple.
	/*
	if ((f_C>0.0)&&(f_C<=x1)) {
		f_f=((yQ*(1.0-3.0*xQ+2.0*yQ))/(xQ*(1.0-xQ)))*f_C;
	}
	else if ((f_C>x1)&&(f_C<xQ)) {
		doublereal af=(1.0-yQ)/(1.0-xQ);
		doublereal bf=(yQ-xQ)/(1.0-xQ);

		f_f=af*f_C+bf;
	}
	else if ((f_C>=xQ)&&(f_C<x2)) {
        doublereal af=(yQ*(1.0-yQ))/(xQ*(1.0-xQ));
		doublereal bf=(yQ*(yQ-xQ))/(1.0-xQ);

		f_f=af*f_C+bf;

	}
	else if ((f_C>=x2)&&(f_C<1.0)) {
		f_f=1.0;
	}
	else f_f=f_C;
	*/

	// Здесь приведён непрерывный вариант он получен
	// изменением функции на втором участке: от x1 до xQ.
	if ((f_C>0.0)&&(f_C<=x1)) {
		f_f=((yQ*(1.0-3.0*xQ+2.0*yQ))/(xQ*(1.0-xQ)))*f_C;
	}
	else if ((f_C>x1)&&(f_C<xQ)) {
		doublereal af=(yQ*(3.0+xQ-4.0*yQ))/(3.0*xQ*(1.0-xQ));
		doublereal bf=(4.0*yQ*(yQ-xQ))/(3.0*(1.0-xQ));

		f_f=af*f_C+bf;
	}
	else if ((f_C>=xQ)&&(f_C<x2)) {
        doublereal af=(yQ*(1.0-yQ))/(xQ*(1.0-xQ));
		doublereal bf=(yQ*(yQ-xQ))/(1.0-xQ);

		f_f=af*f_C+bf;

	}
	else if ((f_C>=x2)&&(f_C<1.0)) {
		f_f=1.0;
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_SUPER_C_SCHEME

// схема ISNAS на неравномерной сетке.
doublereal UNEVEN_ISNAS_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	// ссылка: A High Resolution Pressure Based Algorithm for Fluid Flow at All Speeds.
	// F. Moukalled and M.Darwish 9 июня 1999 года.

	// Interpolation Scheme which is Nonoscillatory for Advected Scalars
	// ISNAS неосциллирующая интерполяционная схема для переноса скалярной величины.

	doublereal f_f=1.0; // инициализация.

	/* 
	// Данная реализация по-видимому содержит ошибку.
	if ((f_C>0.0)&&(f_C<1.0)) {
		doublereal af, bf, cf;
		doublereal zQ=xQ*xQ*(1.0-xQ)*(1.0-xQ);
		af=((yQ-xQ)*(yQ-xQ))/zQ;
		bf=-(2.0*xQ*xQ*xQ+xQ*yQ*yQ-xQ*yQ+yQ*yQ*yQ-3.0*xQ*xQ*yQ)/zQ;
		cf=(xQ*(xQ*xQ*xQ-3.0*xQ*yQ+yQ*yQ+yQ))/zQ;

		f_f=af*f_C*f_C*f_C+bf*f_C*f_C+cf*f_C;
	}
	else f_f=f_C;
	*/

	// Данная схема не прошла проверку в программе maple. График получается лежащим вне области
	// удовлетворяющей критерию конвективной ограниченности. 
	// В интернете оригинал данной схемы найти не удалось.
	printf("ISNAS scheme is incorrect. See the code.\n");
	printf("press any key to exit is programm. \n");
	//system("pause");
	system("pause");
	exit(0); // выход из программы.
	// закомментировал 11.01.2020
	// как недостижимый код.
	//f_f=f_C; // простейшая схема UDS - противопоточная первого порядка.

	return (f_f);
} // UNEVEN_ISNAS_SCHEME

// схема CUBISTA на неравномерной сетке.
doublereal UNEVEN_CUBISTA_SCHEME(doublereal f_C, doublereal xQ, doublereal yQ) {

	// CUBISTA - Convergent and Universally Bounded Interpolation Scheme for the Treatment of Advection.
	// Хорошо сходящаяся и универсально ограниченная схема обращения конвекции.

	// ссылка: A new convergent and universally bounded interpolation scheme for the treatment of advection:
	// application to viscoelastic simulations. M.A. Alves, P.J.Oliveira, F.T.Pinho. (не ранее 2001 года).
	// на cfd wiki данная схема приведена с ошибкой, здесь эта ошибка исправлена.

	doublereal f_f=1.0; // инициализация.

	doublereal x1=3.0*xQ/4.0;
	doublereal x2=xQ*(1.0+2.0*(yQ-xQ))/(2.0*yQ-xQ);

	if ((f_C>0.0)&&(f_C<x1)) {
		f_f=(1.0+(yQ-xQ)/(3.0*(1.0-xQ)))*(yQ/xQ)*f_C;
	}
	else if ((f_C>=x1)&&(f_C<x2)) {
		doublereal af, bf;
		af=(yQ*(1.0-yQ))/(xQ*(1.0-xQ));
		bf=(yQ*(yQ-xQ))/(1.0-xQ);

		f_f=af*f_C+bf;
	}
	else if ((f_C>=x2)&&(f_C<1.0)) {
		f_f=1.0-((1.0-yQ)*(1.0-f_C))/(2.0*(1.0-xQ));
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_CUBISTA_SCHEME


// схема GAMMA на неравномерной сетке. 
// Данная схема содержит один параметр: betam,
// некоторые возможные значения данного параметра 0.45 или 0.5
// 0.5 - значение параметра по умолчанию.
doublereal UNEVEN_GAMMA_SCHEME(doublereal f_C, doublereal xQ,  doublereal yQ, doublereal beta_m=0.5) {

	doublereal f_f=1.0; // инициализация.

	// схема аналитически проверена на непрерывность в точке beta_m
	// 13 марта 2012 года.
	
	if ((f_C>=0.0)&&(f_C<beta_m)) {
		f_f=(1.0+(1.0/beta_m)*((yQ-xQ)/(1.0-xQ))*(1.0-f_C))*f_C;
	}
	else if ((f_C>=beta_m)&&(f_C<1.0)) {
		f_f=((1.0-yQ)/(1.0-xQ))*f_C+((yQ-xQ)/(1.0-xQ));
	}
	else f_f=f_C;

	return (f_f);
} // UNEVEN_GAMMA

// схемы из статьи К.Н. Волкова для неравномерной сетки.
// См. документ "Дискретизация конвективных потоков в уравнениях Навье-Стокса на основе разностных схем
// высокой разрешающей способности."
// Один и тот же код используется для интерполяции на всех гранях, см. таблицу соответствия в этой функции.
doublereal workKN_VOLKOV(doublereal xA, doublereal xB, doublereal xC, doublereal xD,
	           doublereal FA, doublereal FB, doublereal FC, doublereal FD, doublereal Ff, integer isheme) {

	/* таблица соответствия:
	*  A	B	C	D	g	+/-
	*  W	P	E	-	e	+
	*  -	P	E	EE  e   -
	*  WW   W	P	-	w	+
	*  -	W	P	E	w	-
	*  S	P	N	-	n	+
	*  -	P	N	NN  n	-
	*  SS   S	P	-	s	+
	*  -	S	P	N	s	-
	*  B	P	T	-	t	+
	*  -	P	T	TT  t	-
	*  BB   B	P	-	b	+
	*  -	B	P	T	b	-
	*/

	// Ff - скорость на грани.

    doublereal xAB, xBC, xCD;
	xAB=xB-xA;
	xBC=xC-xB;
	xCD=xD-xC;
	doublereal C1, C2, C3;
	C1=xBC/(xBC+xAB);
	C2=xBC/(1.5*xBC+0.5*xCD);
	C3=0.5*(xBC+xCD)/(0.5*xBC+1.5*xCD);

	doublereal denominator = 0.0;

	doublereal r=0.0;
	if (Ff>=0.0) {

		denominator=(FC-FA); // знаменатель

		if (fabs(denominator)<1e-30) {
			// чистый ноль
			// эта ситуация встречается на первой итерации когда скорость всюду равна нулю.

			// нужно включить чистую противопоточную схему 
			r=FB;
		}
		else
		{

		    doublereal xQ, yQ, sQ;
		    xQ=C2/(C1+C2);
		    yQ=C2*(1.0+C1)/(C1+C2);
		    sQ=(1.0+C1)*(1.0-C2);
		    doublereal f_C, f_f;
		    // U - Upwind, C - Center, D - Direct.
		    // f_C=(FC-FU)/(FD-FU) 
		    f_C=(FB-FA)/(FC-FA);
		    switch (isheme) {
				case UNEVEN_MUSCL:  f_f=UNEVEN_MUSCL_SCHEME(f_C, xQ, yQ); break; // TVD second order
		        case UNEVEN_SOUCUP: f_f=UNEVEN_SOUCUP_SCHEME(f_C, xQ, yQ); break; // TVD second order (MINMOD)
		        case UNEVEN_HLPA: f_f=UNEVEN_HLPA_SCHEME(f_C, xQ, yQ); break; // TVD second order
		        case UNEVEN_SMART: f_f=UNEVEN_SMART_SCHEME(f_C,  xQ, yQ); break; // ограниченная 3 порядка.
		        case UNEVEN_WACEB:  f_f=UNEVEN_WACEB_SCHEME(f_C, xQ, yQ); break; // TVD на основе SMART.
		        case UNEVEN_SMARTER: f_f=UNEVEN_SMARTER_SCHEME(f_C, xQ, yQ, sQ); break; // третьего порядка на основе SMART.
		        case UNEVEN_STOIC: f_f=UNEVEN_STOIC_SCHEME( f_C, xQ, yQ); break;
		        case UNEVEN_CLAM: f_f=UNEVEN_CLAM_SCHEME( f_C, xQ, yQ);  break;
		        case UNEVEN_OSHER:  f_f=UNEVEN_OSHER_SCHEME( f_C, xQ, yQ); break;
		        case UNEVEN_VONOS: f_f=UNEVEN_VONOS_SCHEME( f_C, xQ, yQ); break; // неосциллирующая третьего порядка
		        case UNEVEN_LPPA: f_f=UNEVEN_LPPA_SCHEME( f_C, xQ, yQ, sQ); break;
		        case UNEVEN_EXPONENTIAL: f_f=UNEVEN_EXPONENTIAL_SCHEME( f_C); break;
		        case UNEVEN_SUPER_C: f_f=UNEVEN_SUPER_C_SCHEME( f_C, xQ, yQ); break;
		        case UNEVEN_ISNAS: f_f=UNEVEN_ISNAS_SCHEME( f_C, xQ, yQ); break;
		        case UNEVEN_CUBISTA: f_f=UNEVEN_CUBISTA_SCHEME( f_C, xQ, yQ); break;
				case UNEVEN_GAMMA: f_f=UNEVEN_GAMMA_SCHEME( f_C, xQ, yQ); break; // последний параметр по умолчанию.
		        default: f_f=UNEVEN_CUBISTA_SCHEME( f_C, xQ, yQ); break; // по умолчанию хорошо сходящаяся схема  CUBISTA
		    }
		    r=f_f*(FC-FA)+FA;
		}
	}
	else 
	{

		denominator=(FB-FD); // знаменатель

		if (fabs(denominator)<1e-30) {
			// чистый ноль
			// эта ситуация встречается на первой итерации когда скорость всюду равна нулю.

			// нужно включить чистую противопоточную схему 
			r=FC;
		}
		else
		{


		    doublereal xQ, yQ, sQ;
		    xQ=(1.0-C2)/(1.0-C2+C3);
		    yQ=((1.0-C2)*(1.0+C3))/(1.0-C2+C3);
		    sQ=C2*(1.0+C3);
		    doublereal f_C, f_f;
			// U - Upwind, C - Center, D - Direct.
		    // f_C=(FC-FU)/(FD-FU) 
		    f_C=(FC-FD)/(FB-FD);
		    switch (isheme) {
		       case UNEVEN_MUSCL:  f_f=UNEVEN_MUSCL_SCHEME(f_C, xQ, yQ); break; // TVD second order
		       case UNEVEN_SOUCUP: f_f=UNEVEN_SOUCUP_SCHEME(f_C, xQ, yQ); break; // TVD second order (MINMOD)
		       case UNEVEN_HLPA: f_f=UNEVEN_HLPA_SCHEME(f_C, xQ, yQ); break; // TVD second order
               case UNEVEN_SMART: f_f=UNEVEN_SMART_SCHEME(f_C, xQ, yQ); break; // ограниченная 3 порядка.
		       case UNEVEN_WACEB:  f_f=UNEVEN_WACEB_SCHEME(f_C, xQ, yQ); break; // TVD на основе SMART
               case UNEVEN_SMARTER: f_f=UNEVEN_SMARTER_SCHEME(f_C, xQ, yQ, sQ); break; // третьего порядка на основе SMART.
		       case UNEVEN_STOIC: f_f=UNEVEN_STOIC_SCHEME( f_C, xQ, yQ); break;
		       case UNEVEN_CLAM: f_f=UNEVEN_CLAM_SCHEME( f_C, xQ, yQ);  break;
		       case UNEVEN_OSHER: f_f=UNEVEN_OSHER_SCHEME( f_C, xQ, yQ); break;
		       case UNEVEN_VONOS: f_f=UNEVEN_VONOS_SCHEME( f_C, xQ, yQ); break; // неосциллирующая третьего порядка
		       case UNEVEN_LPPA: f_f=UNEVEN_LPPA_SCHEME( f_C, xQ, yQ, sQ); break;
               case UNEVEN_EXPONENTIAL: f_f=UNEVEN_EXPONENTIAL_SCHEME( f_C); break;
		       case UNEVEN_SUPER_C: f_f=UNEVEN_SUPER_C_SCHEME( f_C, xQ, yQ); break;
		       case UNEVEN_ISNAS: f_f=UNEVEN_ISNAS_SCHEME( f_C, xQ, yQ); break;
		       case UNEVEN_CUBISTA: f_f=UNEVEN_CUBISTA_SCHEME( f_C, xQ, yQ); break;
			   case UNEVEN_GAMMA: f_f=UNEVEN_GAMMA_SCHEME( f_C, xQ, yQ); break; // последний параметр по умолчанию.
		       default: f_f=UNEVEN_CUBISTA_SCHEME( f_C, xQ, yQ); break; // по умолчанию хорошо сходящаяся схема CUBISTA
		    }
		    r=f_f*(FB-FD)+FD;
		}
	}

	if (r != r) {

		printf("denominator = %e\n", denominator);
		system("pause");
	}

	// 21.01.2021
	// Для задачи АМПЛИТРОН было выяснено что для температуры получаются
	// значения r гораздо большие 5000.0  и уравнение тепопередачи сильно расходится.
	// Здесь мы искуственно ограничиваем значения r.
	//r = fmin(r, 5000.0); Это не помогает сойтись температуре.
	//if (r > 25000) {

		//std::cout << "convection HO scheme divergence detected\n" << std::endl;
		//std::cout << " r="<< r <<" FA=" << FA << " FB=" << FB << " FC=" << FC << " FD=" << FD << " Ff=" << Ff << " FB -FD=" << FB-FD << " FC -FA=" << FC - FA << std::endl;
		//system("pause");
	//}
	return (r);
} // workKN_VOLKOV

#endif
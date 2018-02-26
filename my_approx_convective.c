// Файл my_approx_convective.c 
// аппроксимация конвективного члена

#ifndef MY_APPROX_CONVECTIVE_C
#define MY_APPROX_CONVECTIVE_C 1

//#include <math.h>

//#define doublereal double

// схемы для аппроксимации конвекции-диффузии
#define CR 1 // Центрально-разностная
#define UDS 2 // Противопоточная первого порядка
#define COMB 3 // Комбинированная 
#define POLY 4 // Полиномиальная C. Патанкара
#define EXP 5 // экспоненциальная схема
#define BULG 6 // схема В.К. Булгакова (23) из статьи
#define POW 7 // показательная

// возвращает максимальное из двух 
// вещественных чисел.
/*
doublereal fmax(doublereal fA, doublereal fB) {
	doublereal r=fB;
	if (fA > fB) r=fA;
	return r;
} // fmax
*/

// Функция A(|P|) 
// для различных схем:
// fabsPe - модуль числа Пекле.
// ishconvection - номер схемы:
// CR==1 - центрально-разностная схема,
// UDS==2 - с разностями против потока,
// COMB==3 - комбинированная,
// POLY==4 - со степенным законом,
// EXP==5 - экспоненциальная (точная),
// BULG==6 - схема В.К. Булгакова (23) из статьи,
// POW==7 - показательная зависимость.
doublereal ApproxConvective(doublereal fabsPe, integer ishconvection) {
	doublereal r=1; // по умолчанию с разностью против потока
	doublereal rCR, rpoly, rpoly2, rpoly5;
	if (ishconvection<6) {
 	   rCR=1.0-0.5*fabsPe;
	   rpoly=1.0-0.1*fabsPe;
	   rpoly2=rpoly*rpoly;
	   rpoly5=rpoly2*rpoly2*rpoly; // со степенным законом
	}

	switch (ishconvection) {
		case CR : // Центрально разностная схема пригодна только 
			    // для небольших чисел Пекле меньших по модулю 1.0.
			    if (fabsPe<1.0) r=rCR; 
				// поэтому при числах Пекле по модулю больших 1.0 но
				// меньших 10.0 осуществляется переход на рекомендуемую
				// Патанкаром схему со степенным законом.
			    if ((fabsPe>=1.0)&&(fabsPe<10.0)) r=rpoly5;
				else r=0.0; // при числах Пекле больших 10.0 аппроксимация нулём.
			    break;
		case UDS : // с разностями против потока
			    //if (fabsPe<1.0) r=1.0; else r=0.0;
				r=1.0;
				//getchar();
			    break;
		case COMB: // комбинированная
			    r=fmax(0.0, rCR);
			    break;
		case POLY : // со степенным законом
			    if (fabsPe<10.0) r=rpoly5; else r=0.0;
			    break;
		case EXP : // экспоненциальная (точная)
			    if (fabsPe<10.0) {
			    	if (fabsPe<0.01) r=rCR; else r=fabsPe/(exp(fabsPe)-1.0);
			    }
			    else r=0.0;
			    break;
		case BULG : // В.К. Булгаков, Н.В. Булгаков
			    // Хабаровкий Государственный Технический Университет
			    // зависимость (23) из статьи.
			    r=1.0/(1.0+0.6712*fabsPe*fabsPe);
			    break;
		case POW : // показательная зависимость
			    r=pow(0.553,fabsPe);
			    break;
		default: r=1.0/(1.0+0.6712*fabsPe*fabsPe); break;
	}
	return r;
} // ApproxConvective

#endif
// 24 ноября 2011 вода и воздух
// 25 ноября 2011 года 15 твёрдых тел.
// 26-27 ноября 2011 года написание интерфейса AliceMeshv0_08.
// 28 ноября 2011 года написание update_temp_properties только
// для внутренних контрольных объёмов.
// В модуле my_material_properties.c содержится
// внутрипрограммная библиотека нелинейных свойств материалов.

#ifndef  MY_MATERIAL_PROPERTIES_CPP
#define  MY_MATERIAL_PROPERTIES_CPP 1

// Объявление модели вещественной арифметики содержится в самом начале программы в главном модуле 
// AliceFlow_v0_48
//#define doublereal double
#include <math.h>
#include <time.h>

// Осуществляет выход из программы в случае превышения 
// критической температуры определённой в паспортных данных
// на мощный 110 ватный транзистор TGF2023-20: TEMPERATURE_FAILURE_DC.
void diagnostic_critical_temperature(doublereal TiP,  FLOW* &fglobal, TEMPER &t, 
	BLOCK* b, int lb) 
{
	 if (TiP > TEMPERATURE_FAILURE_DC) {
		   printf("Attantion! unit burned...\n");
		   printf("temperature exceeds a critical TEMPERATURE_FAILURE_DC==%3.2f",TEMPERATURE_FAILURE_DC);



		   if (1) {
			   if (!b_on_adaptive_local_refinement_mesh) {

				   bool bextendedprint = false;
				   integer flow_interior = 0;
				   // экспорт результата вычисления в программу tecplot360:
				   exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, bextendedprint, 0, b, lb);
			   }
			   else {
				   
				   // Экспорт в программу tecplot температуры.
				   //С АЛИС сетки.
				   ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t, fglobal, 0, b, lb);
			   }
		   }


		   unsigned int calculation_main_end_time = clock();
		   unsigned int calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;
		   
		   // Общее время вычисления до возникновения критического сообщения.
		   int im = 0, is = 0, ims = 0;
		   im = (int)(calculation_main_seach_time / 60000); // минуты
		   is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // секунды
		   ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

		   printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);
		   //  getchar();
		   system("pause");
		   char ch;
		   printf("please enter c)ontinue h)alt.\n");
		   ch=getchar();
		   if ((ch=='h')||(ch=='H')) exit(0);
	   }
} // diagnostic_critical_temperature

// значение 0 для материала заданного пользователем.
#define MY_AIR 1
#define MY_WATER 2

// В основном в качестве жидких рабочих веществ используются вода или воздух.
// Ну и совсем редко (очень и очень редко) что-либо другое.
// Поэтому здесь подобраны максимально приближенные к реальности свойства воды и воздуха.

// Воздух используется по умолчанию во всех расчётах,
// если пользователь не захотел иначе.
void my_air_properties(doublereal TiP, doublereal PiP, 
	                   float &rho, float &cp,
					   float &lam, float &mu,
					   float &beta) {
	// возвращает реальные свойства воздуха !
	// Входные параметры: 
    // TiP - температура в градусах цельсия в центре контрольного объёма.
	// PiP - давление в избыточных (отсчитываемое от нуля, 0 - опорное значение). 

	// Данная приближённая аппроксимация справедлива в диапазоне
	// температур 273-473K.
	// Ссылка: HighExpert.Ru Физические свойства воздуха
	// www.engineeringtoolbox.com.

    doublereal TK=TiP+273.15; // K
	//doublereal Pair=101325+PiP; // Па

	// скорость звука в воздухе:
	/* temp grad C  Vel m/s  Vel km/h
	*  -150			216.7	  780.1
	*  -100			263.7	  949.2
	*  -50			299.3	  1077.6
	*  -20			318.8	  1147.8
	*  -10			325.1	  1170.3
	*  0			331.5	  1193.4
	*  10			337.3	  1214.1
	*  20			343.1	  1235.2
	*  30			348.9	  1256.2
	*  50			360.3	  1296.9
	*  100			387.1     1393.7
	*  200			436.0     1569.5
	*  300			479.8     1727.4
	*  400			520.0     1872.1
	*  500			557.3     2006.4
	*  1000			715.2     2574.8
	*/
	// уравнения справедливы при числах Маха < 0.3.
	// Пороговая скорость 102 м/с.
	//rho=Pair/(287.4*TK); // kg/m^3
	// 27. 07. 2016.
	// В приближении Буссинеска опорное значение плотности константа,
	// зависимость плотности от температуры учитывается через коэффициент линейного температурного расширения.
	rho = 1.1614f;
	cp=(float)(1000.0*(1.0005+1.1904e-4*(TiP))); // J/(kg*GC)
	lam= (float)(2.44e-2*exp(0.82*log(TK/273.15))); // W/(m*K)
	mu= (float)(1.717e-5*exp(0.683*log(TK/273.15))); // Pa*s
	//printf("mu=%e, TK=%e\n",mu,TK); // debug Ok
	//getchar();
	// в maple по способу наименьших квадратов получена следующая зависимость:
    // Исправлено 27.07.2016
	beta = (float)(0.001*(-0.06810259 + 0.00013411*TiP - 8.287412214e-8*TiP*TiP + 1022.18 / (273.15 + TiP)));

	// возвращаем постоянные свойства
	//rho=1.1614; // несжимаемость !!!
	//mu=1.84e-5; 
	//beta=0.003331;
	//cp=1005;
	//lam=0.025;
	
} // my_air_properties

// Вода - это второе после воздуха жидкое рабочее тело,
// для расчётов. Его свойства также важны.
void my_water_properties(doublereal TiP,  
	                   float &rho, float &cp,
					   float &lam, float &mu,
					   float &beta) {
	// возвращает реальные свойства воды !
	// Входные параметры: 
    // TiP - температура в градусах цельсия в центре контрольного объёма.

	// Данная приближённая аппроксимация справедлива в диапазоне
	// температур 283-373K.
	// Ссылка: HighExpert.Ru Физические свойства воды
	// www.engineeringtoolbox.com.

    doublereal TK=TiP+273.15;

	// Скорость звука в воде:
	/* temp grad C  Vel Mag m/s
	*		 0		1.403
	*		 5		1.427
	*		10		1.447
	*		20		1.481
	*		30		1.507
	*		40		1.526
	*		50		1.541
	*		60		1.552
	*		70		1.555
	*		80		1.555
	*		90		1.550
	*	   100		1.543
	*/
	// Уравнения справедливы при числах Маха < 0.3.

	//rho=995.7/(0.984+0.483e-3*TiP); // kg/m^3
	// 27. 07. 2016.
	// В приближении Буссинеска опорное значение плотности константа,
	// зависимость плотности от температуры учитывается через коэффициент линейного температурного расширения.
	rho = (float)(995.7 / (0.984 + 0.483e-3*(20.0))); // kg/m!3
	cp= (float)((4194-1.15*TiP+1.5e-2*TiP*TiP)); // J/(kg*GC)
	lam= (float)(0.553*(1.0+0.003*TiP)); // W/(m*K)
	mu= (float)((1.78e-6/(1.0+0.0337*TiP+0.000221*TiP*TiP))*rho); // Pa*s
	// в maple по способу наименьших квадратов получена следующая зависимость:
	beta= (float)(0.02289718770+0.142992772827544502e-6*TiP*TiP-0.67623235e-4*TiP-(6.27314791869)/(TK));

	// возвращаем постоянные свойства при 20 град С.
	//rho=998.2; // несжимаемость !!!
	//mu=1e-3;
	//cp=4177;
	//lam=0.58618;
	//beta=0.192166e-3;
} // my_water_properties

// SOLID твёрдые тела:
// Главный принцип: приводить только проверенные данные.

// Значение 100 для материала заданного пользователем.
//#define MY_AUSN 101    // Припой золото 80% станум 20%
#define MY_ALUMINA 101   // поликор Alumina, Al2O3 -ref. MatWeb[2] 99.9% Al2O3, polycrystalline aluminum oxide.
#define MY_SI 102        // Кремний
#define MY_GAAS 103      // Арсенид Галия
#define MY_GAN 104       // Нитрид Галия
#define MY_SIC4H 105     // Карбид кремния SiC4H
#define MY_SAPPHIRE 106  // Сапфир
#define MY_DIAMOND 107   // Алмаз
#define MY_MD40 108      // Псевдосплав МД40
#define MY_AU 109        // Золото
#define MY_SIO2 110      // оксид кремния
#define MY_COPPER 111    // Медь
#define MY_KOVAR 112     // Ковар
#define MY_BRASS 113     // Латунь ЛС-59-1-Л
#define MY_DURALUMIN 114 // Дюралюминий Д16
#define MY_ALUMINIUM_NITRIDE 115 // Нитрид аллюминия (кристаллический).
#define MY_GLUE_ECHES 116 // Клей ЭЧЭС

// Припой золото станум (олово) 
// 21 марта 2016 года сделаны температурные зависимости на основе
// известных свойств компонентов двухкомпонентного соединения.
void my_ausn_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // AuSn 80 золота/20 олова
       //doublereal TK=TiP+273.15;
	   rho=14510; // kg/m^3
	   // cp=150 по данным программы SYMMIC
	   cp=143; // J/(kg*K) по данным TDIM - Master
	  
	   // Теплопроводность в литературе известна только 
	   // для двух точек: 68 при 300К и 57 при 85 град. С.
	   lam=57; // т.к. нигде не нашёл достоверное значение показателя степени.
	   // формула на основе показателя степени верна по-видимости в основном для 
	   // полупроводников т.е. она не универсальна. Есть металлы и сплавы теплопроводность
	   // которых меняется по закону близкому к линейному (полиномиальному) есть даже вещества у 
	   // которых теплопроводность растёт с ростом температуры.
	   //lam=68*exp(-0.99678*log(TK/300)); // 57 W/(m*K) неверный подход, читай ниже

	   // Теплопроводность и теплоёмкость золото-оловянного припоя нельзя 
	   // вычислять на основе свойств его компонентов в отдельности зная 
	   // процентный состав сплава, т.к. структура соединения сильно не упорядочена
	   // имеются большие случайные вкрапления олова в золото,  поэтому теплопроводность и
	   // теплоёмкость данного соединения хуже чем каждая из  компонент сплава в отдельности.
	   // Данные величины являются экспериментальными данными, теоретически рассчитать которые невозможно.

} // AuSn

// Кремний Si Silicon 
void my_si_properties(doublereal TiP,
	                    float &rho, float &cp,
						float &lam) {
       // Silicon
       float TK=(float)(TiP+273.15);
	   rho=2330; // kg/m^3
	   //cp=711; // J/(kg*K)
	   // Формула для теплоёмкости справедлива в диапазоне температур: 77K-773K.
	   // TK  77 173 273 373 573 773
	   // cp 180 490 680 770 850 880
	   if (TK < 800.0) {
		   cp= (float)(83.23167276-11801.68987/TK+0.3123808545e-5*TK*TK*TK+3.672005743*TK-0.5805489801e-2*TK*TK);
	   }
	   else {
		   // ограничитель сверху на ту область где аппроксимация по способу наименьших квадратов не работает.
		   cp=880;
	   }
	   // Функция изменения теплопроводности от температуры взята из
	   // книги Р. Квэя Электроника нитрида галия.
	   //lam=120; // данные В.Д. Красильникова
	   lam= (float)(148*exp(-1.35*log(TK/300))); // 120 W/(m*K)
} // Si

// 18 october 2016.
// Поликор понадобился при расчётах резистора на поликоровой подложке.
// Alumina, Al2O3-ref. Mat Web[2] 99.9% Al2O3, polycrystalline aluminum oxide.
void my_alumina_properties(doublereal TiP,
	float &rho, float &cp,
	float &lam) {

	float TK = (float)(TiP + 273.15);
	rho = 2330; // kg/m^3

	//cp=753; // J/(kg*K)
	// Формула для теплоёмкости справедлива в диапазоне температур: 200K-800K.
	// TK  200 298 400 500 600 800
	// cp 502 753 920 1046 1088 1172
	if (TK < 800.0) {
		cp = (float)(0.00000505903*TK*TK*TK - 0.01042874335*TK*TK + 7.612830818*TK + (50982.822 / TK) - 898.2525);
	}
	else {
		// ограничитель сверху на ту область где аппроксимация по способу наименьших квадратов не работает.
		cp = 1172;
	}

	// Теплопроводность по данным программы SYMMIC.
	// TK lam
	// 200 82
	// 298 46
	// 400 32.3
	// 500 24.2
	// 600 18.9
	// 800 13.0
	lam = (float)(47.84832382*exp(-1.328556144*log(TK / 300.0))); // 46

} // alumina

// Арсенид Галия GaAs 
void my_gaas_properties(doublereal TiP,
	                    float &rho, float &cp,
						float &lam) {
       // GaAs
       float TK=(float)(TiP+273.15);
	   rho=5300.0f; // kg/m^3
	   cp=322.0f; // J/(kg*K)
	   if ((TK > 173.0) && (TK < 1514.0))
	   {
		   cp = (float)(4.059e-8*TK*TK*TK - 1.296e-4*TK*TK + 1.854e-1*TK + 279.0);
	   }
	   else if (TK <= 173.0) cp = 307.4f;
	   else cp = 402.9f;
	   // Функция изменения теплопроводности от температуры взята из
	   // книги Р. Квэя Электроника нитрида галия.
	   //lam=47; // данные В.Д. Красильникова
	   lam= (float)(54*exp(-1.25*log(TK/300))); // 54 W/(m*K)

	   
} // GaAs

// Нитрид Галия GaN
void my_gan_properties(doublereal TiP,
	                    float &rho, float &cp,
					     float &lam) {
       // GaN
       float TK=(float)(TiP+273.15);
	   rho=6150; // kg/m^3
	   cp=491; // J/(kg*K)
	   // Формула для теплоёмкости справедлива в диапазоне температур: 200K-1300K.
	   // TK  200 300 400 500 600 700 800 900 1000 1100 1200 1300
	   // cp 322.3 431.3 501.2 543.8 572.17 592.8 608.9 622.2 633.7 643.99 653.4 662.22
	   if ((TK > 200) && (TK < 1300)) {
		   cp = (float)(591.807 + 0.064971*TK - 2.6e7 / (TK*TK) + 2.94e9 / (TK*TK*TK));
	   }
	   // Функция изменения теплопроводности от температуры взята из
	   // книги Р. Квэя Электроника нитрида галия.
	   //lam=130; // данные В.Д. Красильникова
	   lam= (float)(130*exp(-0.43*log(TK/300))); // 130 W/(m*K)
} // GaN

// Карбид Кремния SiC4H
void my_sic4h_properties(doublereal TiP,
	                    float &rho, float &cp,
						float &lam) {
       // SiC4H
       float TK=(float)(TiP+273.15);
	   rho=3210; // kg/m^3
	   //cp=690; // J/(kg*K)
	   // Зависимость теплоёмкости карбида от температуры
	   // построена по следующим опорным точкам:
	   // TK 77 300 373 573 773
	   // Cp 50 690 820 1010 1120
	   if ((TK>=77.0)&&(TK<=773.2)) {
		   cp= (float)(-600.9595790625+6.41266663224237*TK+16052.9931969468/TK-0.00900588312450910*TK*TK+0.459922567674634e-5*TK*TK*TK);
	   } else if (TK<77.0) cp=50;
	   else cp=1120;

	   // Функции изменения теплопроводности от температуры 
	   // в книге Р. Квэя Электроника нитрида галия. НЕТ!!!
	   //lam=370; // V-легир. SiC
	   lam= (float)(370*exp(-1.411*log(TK/300))); // 370 W/(m*K)
} // SiC4H

// Сапфир Sapphire
void my_sapphire_properties(doublereal TiP,
	                    float &rho, float &cp,
						float &lam) {
       // Sapphire Al2O3
       float TK=(float)(TiP+273.15);
	   rho=3980; // kg/m^3
	   //cp=750; // J/(kg*K)
	   // Формула для теплоёмкости справедлива в диапазоне температур 77K - 773K
	   // TK 77 173 273 373 573 773
	   // Cp 60 403 718 907 1089 1168
	   if (TK<774) {
		   cp= (float)(-703.3920492+19101.78790/TK+0.5623375490e-5*TK*TK*TK+7.499240665*TK-0.01095728614*TK*TK);
	   }
	   else cp=1168;

	   // Функции изменения теплопроводности от температуры 
	   // в  книге Р. Квэя Электроника нитрида галия. НЕТ!!!
	  // lam=23; // Данные В.Д. Красильникова
	   // По данным программы SYMMIC.
	   lam= (float)(42*exp((-1.134821)*log(TK/298))); // 42 W/(m*K) при 298К

} // Sapphire

// Алмаз Diamond
void my_diamond_properties(doublereal TiP,
	                    float &rho, float &cp,
						float &lam) {
       // Diamond
       float TK=(float)(TiP+273.15);
	   rho=3500; // kg/m^3
	   //cp=520; // при 300К J/(kg*K)
	   // Данные по теплоёмкости справедливы в диапазоне температур:
	   // 77K - 773K
	   // TK 77  173 273 373 573 773
	   // Cp 8  140  420 770 1300 1590 и продолжает расти, но данных нет.
	   if (TK < 774) {
		   cp= (float)(-696.4289+32660.38060/TK-0.2774926361e-5*TK*TK*TK+3.566362474*TK+0.1285766705e-2*TK*TK);
	   }
	   else cp=1590;
	   // cp=1260; // По данным В.Д. Красильникова
	   // Функция изменения теплопроводности от температуры взята из
	   // книги Р. Квэя Электроника нитрида галия.
	   //lam=800; // данные В.Д. Красильникова
	   lam= (float)(2300*exp(-1.85*log(TK/300))); // 130 W/(m*K)
	   // По данным Р. Квэя теплопроводность алмаза при 300К лежит между 2000 и 2500.
} // Diamond

// Псевдосплав CuMo с 40% содержанием меди: МД40
void my_md40_properties(doublereal TiP,
	                    float &rho, float &cp,
						float &lam) {
       // Из материала этого псевдосплава сделан корпус прибора:
	   // Основание, боковые стенки и крышка. Крышка приклеивается.
       // MD40
       float TK=(float)(TiP+273.15);
	   rho=9660.0f; // kg/m^3
	   // cp=540 по данным программы Захарова, Асвадуровой.
	   // Почему в программе Захарова и Асвадуровой такая большая теплоёмкость непонятно.
	   cp=318.2f; // J/(kg*K)
	   // heat capacity cp
	   // TK 173 273 373 573 
	   // cp J/(kg-K) 268.8 302.5 318.17 335.1
	   if ((TK > 173) && (TK < 573))
	   {
		   cp = (float)(1.66e-6*TK*TK*TK - 2.263e-3*TK*TK + 1.095*TK + 138.5);
	   }
	   else if (TK <= 173) cp = 268.77f;
	   else cp = 335.111f;
	   lam=210.0f; // W/(m*K)
	   if ((TK >= 273) && (TK < 1600)) {
		   lam = (float)(7.15e-9*TK*TK*TK - 2.18e-6*TK*TK - 6.03e-2*TK + 271.0);
	   }
	   else if (TK < 273) lam = 254.52f;
	   else lam = 198.0f;
} // MD40

// Золото Au
void my_au_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // Au
       float TK=(float)(TiP+273.15);
	   rho=19300; // kg/m^3
	   //cp=126; // J/(kg*K)

	   // Функциональная зависимость получена по следующим точкам
	   // с помощью МНК:
	   // TK 77 173 273 373 573 773
	   // Cp 97 121 128 131 135 140
	   if ((TK>77.0)&&(TK<773.0)) {
		   cp= (float)(140.2979168-3345.749125/TK+0.3399388027e-7*TK*TK*TK+0.3605589168e-2*TK-0.2418961279e-4*TK*TK);
	   }
	   else if (TK<=77.0) cp=97;
	   else cp=140;

	   //lam=293; // W/(m*K)
	   // Вид зависимости, предложенный в книге Р. Квэя не годится для
	   // теплопроводности золота.
	   // Поэтому была получена своя оригинальная формула с помощью МНК.
	   if (TK<974) {
		   lam= (float)(327.5268055+146.5916497/TK+0.5614880962e-7*TK*TK*TK-0.008865308782*TK-0.1043207287e-3*TK*TK);
	   }
	   else lam=272;

} // Au

// SiO2
// Термически выращенный из газовой фазы.
// Достоверность данных свойств гарантируется программой SYMMIC.
void my_sio2_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // SiO2
       float TK=(float)(TiP+273.15);
	   rho=2100.0f; // kg/m^3
	   cp=741.0f; // J/(kg*K)
	   if ((TK >= 273.0) && (TK <= 773))
	   {
		   cp = (float)(-1.6667e-7*TK*TK*TK - 9.635e-4*TK*TK + 1.975*TK + 236.02);
	   }
	   else if (TK < 273.0) cp = 700.0f;
	   else cp = 1110.0f;

	   //lam=7; // 7-11.5 W/(m*K) quartz crystal
	   lam = 1.27f; // термически выращенный из газовой фазы.

} // SiO2

// Медь Copper
void my_copper_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // Copper
       float TK=(float)(TiP+273.15);
	   rho=8930; // kg/m^3
	   //cp=390; // J/(kg*K)
	   // Функциональная зависимость теплоёмкости от температуры
	   // получена по следующим точкам с помощью МНК.
	   // TK    77  173 273 373 573 773
	   // Cp    195 341 379 397 419 430 
	   if ((TK>77.0)&&(TK<773.0)) {
		   cp= (float)(506.2148755-22501.73085/TK-0.3085157062e-6*TK*TK*TK-0.2854473337*TK+0.5289207596e-3*TK*TK);
	   }
	   else if (TK<77.0) cp=195;
	   else cp=430;

	   //lam=390; // W/(m*K)
	   // Функциональная зависимость теплопроводности от температуры
	   // получена по следующим точкам с помощью МНК.
	   // TK      173.2 273.2 373.2 573.2 973.2
	   // Lambda  420.0 403.0 395.0 381   354
	   if ((TK>173.2) && (TK<973.2)) {
		   lam= (float)(318.4875333+12313.69842/TK+0.1748333335e-6*TK*TK*TK+0.2380247005*TK-0.3905912003e-3*TK*TK);
	   }
	   else if (TK<=173.2) lam=420.0;
	   else lam=354.0;
} // Copper

// Ковар Kovar
void my_kovar_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // Kovar
	   // Ковар это сплав четырёх веществ: C, Fe 54%, Mn, Ni 29%
       float TK=(float)(TiP+273.15);
	   rho=8300; // kg/m^3
	   cp=669; // J/(kg*K)
	   // heat capacity kovar
	   // TK 273 703
	   // cp 439.614 648.954
	   if ((TK>273) && (TK < 703))
	   {
		   cp = (float)(0.486*TK + 306.6);
	   }
	   else if (TK < 273) cp = 439.6;
	   else cp = 648.954;
	   //lam=19; // W/(m*K) 
	   // Функциональная зависимость теплопроводности от температуры
	   // получена по следующим точкам с помощью МНК.
	   // TK      273  373  573  773  973
	   // Lambda  14.1 14.7 15.6 17.5 19.3
	   if ((TK>273.0) && (TK<973.0)) {
		   lam= (float)(44.13357753-3650.246187/TK-0.4572641192e-7*TK*TK*TK-0.8854732756e-1*TK+0.1132267331e-3*TK*TK);
	   } else if (TK<273.0) lam=14.1;
	   else lam=19.3;
} // Kovar

// Латунь ЛС-59-1-Л Brass
void my_brass_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // Brass
	   // Латунь Cu 70%, Zn 30%
       float TK=(float)(TiP+273.15);
	   rho=8440; // kg/m^3
	   cp=377; // J/(kg*K)
	   //lam=108.784; //  W/(m*K)
	   // Данные о температурной зависимости теплопроводности Латуни:
	   // T GC   -100 0 100 200 300 400
	   // TK      173.2 273.2 373.2 473.2 573.2 673.2
	   // Lambda  90 106 131 143 145 148
	   if ((TK>=173.2)&&(TK<=673.2)) {
		   lam= (float)(-591.598556543738+2.90903837477216*TK+52133.0458047673/TK-0.00454066869215317*TK*TK+0.249629141483149e-5*TK*TK*TK);
	   } else if (TK<173.2) lam=90.0;
	   else lam=148.0;

} // Brass

// Дюралюминий Duralumin Д16
void my_duralumin_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // Duralumin
       float TK=(float)(TiP+273.15);
	   rho=2780.0f; // kg/m^3
	   cp=922.0f; // J/(kg*K)
	   // heat capacity duralumin
	   // TK  273 373 473 573 673
	   // cp J/(kg-K) 750.13 921 1047 1130 1172
	   if ((TK>273.0) && (TK < 673.0)) {
		   cp = (float)(3.33e-7*TK*TK*TK - 2.62e-3*TK*TK + 3.3*TK + 38.0);
	   }
	   else if (TK <= 273.0) cp = 750.13f;
	   else cp = 1172.0f;
	   lam=130.0f; // W/(m*K)
	   if ((TK>150.0) && (TK < 573)) {
		   lam = (float)(0.171*TK + 66.1);
	   }
	   else if (TK <= 150.0) lam = 90.0f;
	   else lam = 164.0f;
} // Duralumin

// Клей ЭЧЭС
void my_glueeches_properties(doublereal TiP,
	                  float &rho, float &cp,
					  float &lam) {
       // glue ECHES
	   // По данным В.Д. Красильникова.
       //doublereal TK=TiP+273.15;
	   rho=1000; // kg/m^3
	   cp=100; // J/(kg*K)
	   lam=4; // W/(m*K)
} // glue ECHES

// Нитрид Алюминия (кристаллический)
void my_aluminium_nitride_properties(doublereal TiP,
	float &rho, float &cp,
	float &lam) {
	// aluminium nitride
	// По данным программы SYMMIC.
	float TK = (float)(TiP + 273.15);
	rho = 3255; // kg/m^3
	cp = 600; // J/(kg*K) // SYMMIC
	//cp=748; // J/(kg*K) // Р. Куэй.
	//lam = 319; // W/(m*K)
	// lam=285; // W/(m*K) Р. Куэй

	// TSIDE,K lam W/(mxK)
	// 200 780
	// 300 319
	// 400 195
	// 600 100
	// 1000 49

	// Функции изменения теплопроводности от температуры 
	// в  книге Р. Куэя Электроника нитрида галия. присутствует, как для кристаллического так и для керамики.
	// lam=319; // Данные программы SYMMIC
	// По данным программы SYMMIC.
	lam = (float)(319 * exp((-1.57)*log(TK / 300))); // 319 W/(m*K) при 300К
} // aluminium nitride


// внутрипрограммная библиотека жидких материалов.
void my_fluid_properties(doublereal TiP, doublereal PiP, 
	                   float &rho, float &cp,
					   float &lam, float &mu,
					   float &beta, integer ilibident) {
	switch (ilibident) {
	  case MY_AIR: my_air_properties(TiP, PiP, rho, cp, lam, mu, beta); break;
	  case MY_WATER: my_water_properties(TiP, rho, cp,lam, mu, beta); break;
	  default: my_air_properties(TiP, PiP, rho, cp, lam, mu, beta);  break;
	} // end switch
} // my_fluid_properties

doublereal dsic=1.0e30, dgan=1.0e30, dcu=1.0e30;

 // внутрипрограммная библиотека твёрдых материалов.
void my_solid_properties(doublereal TiP, float &rho, float &cp,
	                     float &lam, integer ilibident) {
	switch (ilibident) {
	  //case MY_AUSN: my_ausn_properties(TiP,rho,cp,lam); break; // припой золото станум
	  case MY_ALUMINA:  my_alumina_properties(TiP, rho, cp, lam); break; // поликор Al2O3.
	  case MY_SI: my_si_properties(TiP,rho,cp,lam); break; // кремний
	  case MY_GAAS: my_gaas_properties(TiP,rho,cp,lam); break; // арсенид галия
	  case MY_GAN: my_gan_properties(TiP, rho, cp, lam); if (lam < dgan) { dgan = lam; }  break; // нитрид галия
	  case MY_SIC4H: my_sic4h_properties(TiP, rho, cp, lam); if (lam < dsic) { dsic = lam; } break; // карбид кремния
	  case MY_SAPPHIRE: my_sapphire_properties(TiP,rho,cp,lam); break; // сапфир
	  case MY_DIAMOND: my_diamond_properties(TiP,rho,cp,lam); break; // алмаз
	  case MY_MD40: my_md40_properties(TiP,rho,cp,lam); break; // Псевдосплав МД40
	  case MY_AU: my_au_properties(TiP,rho,cp,lam); break; // Золото
	  case MY_SIO2: my_sio2_properties(TiP,rho,cp,lam); break; // SiO2
	  case MY_COPPER: my_copper_properties(TiP, rho, cp, lam); if (lam < dcu) { dcu = lam; }  break; // Медь
	  case MY_KOVAR: my_kovar_properties(TiP,rho,cp,lam); break; // Ковар
	  case MY_BRASS: my_brass_properties(TiP,rho,cp,lam); break; // Латунь ЛС-59-1-Л Brass
	  case MY_DURALUMIN: my_duralumin_properties(TiP,rho,cp,lam); break; // дюралюминий
	  case MY_ALUMINIUM_NITRIDE: my_aluminium_nitride_properties(TiP, rho, cp, lam); break; // aluminium nitride
	  case MY_GLUE_ECHES: my_glueeches_properties(TiP, rho, cp, lam); break; // клей ЭЧЭС
	  default: my_duralumin_properties(TiP,rho,cp,lam); break; // дюралюминий - библиотечный материал по умолчанию
	} // end switch
} // my_solid_properties

// Меняет тепловые свойства материала в граничном КО.
void gran_prop(TEMPER &t, FLOW* &f, BLOCK* b, int lb, integer iP, integer G, int ib, TPROP* matlist) {

	// t.ptr[1][iP] - идентификатор Fluid domain или -1 если это чистый Solid.

	float rho, cp, lam;
	integer iG; // номер соседа который может оказаться граничным узлом.
	iG=t.neighbors_for_the_internal_node[G][0][iP];
	if (iG>=t.maxelm) {
		// это граничный узел
        rho=1.1614f; cp=1005.0f; lam=0.025f; // инициализация default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat==1) {
			if (b[ib].itype== PHYSICS_TYPE_IN_BODY::SOLID) {
		        my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // подставляется температура в граничном узле
				// проверка на допустимость температур.
				diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
		    } // SOLID
		    else if (b[ib].itype== PHYSICS_TYPE_IN_BODY::FLUID) {
		       float mu, beta_t; // значения не используются но требуются.
		       doublereal pressure;
		       if (t.ptr[1][iP]==-1) {
			       pressure=0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
		       }
			   else {
				   if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][0][t.ptr[0][iP]] > -1) {
					   pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][0][t.ptr[0][iP]]];
				   }
				   else {
					   // STUB 15.09.2018
					   pressure = 0.0;  // В твёрдом теле.

#if doubleintprecision == 1
					  // printf("error in gran_prop in my_material_properties.c NODE1 G=%lld\n", G);
#else
					   //printf("error in gran_prop in my_material_properties.c NODE1 G=%d\n", G);
#endif

					   
					   //getchar();
					   //system("PAUSE");
					   //exit(1);
				   }
			   }
			   my_fluid_properties(t.potent[iG], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
		    } // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat==0) {
            // материал определённый пользователем:
			// постоянные свойства.
			rho=matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

			// переносим постоянные свойства материала из внутренней точки области
			// в ближайшую граничную точку.
		}
		// Свойства для внутреннего контрольного объёма.
		t.prop_b[RHO][iG-t.maxelm]=rho;
		t.prop_b[HEAT_CAPACITY][iG-t.maxelm]=cp;
		t.prop_b[LAM][iG-t.maxelm]=lam;
	} // G Side

	if (t.neighbors_for_the_internal_node[G][1] != nullptr) {
		iG = t.neighbors_for_the_internal_node[G][1][iP];
		if (iG >= t.maxelm) {
			// это граничный узел
			rho = 1.1614f; cp = 1005.0f; lam = 0.025f; // инициализация default  dry air 300K 1atm properties
			if (matlist[b[ib].imatid].blibmat == 1) {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // подставляется температура в граничном узле
					// проверка на допустимость температур.
					diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
				} // SOLID
				else if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					float mu, beta_t; // значения не используются но требуются.
					doublereal pressure;
					if (t.ptr[1][iP] == -1) {
						pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
					}
					else {
						if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][1][t.ptr[0][iP]] > -1) {
							pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][1][t.ptr[0][iP]]];
						}
						else {
							// STUB 15.09.2018
							pressure = 0.0; // В твёрдом теле.

#if doubleintprecision == 1
						//printf("error in gran_prop in my_material_properties.c NODE2 G=%lld\n", G);
#else
						//printf("error in gran_prop in my_material_properties.c NODE2 G=%d\n", G);
#endif

						//getchar();
						//system("PAUSE");
						//exit(1);
						}
					}
					my_fluid_properties(t.potent[iG], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				} // FLUID
			}
			else if (matlist[b[ib].imatid].blibmat == 0) {
				// материал определённый пользователем:
				// постоянные свойства.
				rho = matlist[b[ib].imatid].rho;
				//cp = matlist[b[ib].imatid].cp;
				//lam = matlist[b[ib].imatid].lam;
				cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
				lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

				// переносим постоянные свойства материала из внутренней точки области
				// в ближайшую граничную точку.
			}
			// Свойства для внутреннего контрольного объёма.
			t.prop_b[RHO][iG - t.maxelm] = rho;
			t.prop_b[HEAT_CAPACITY][iG - t.maxelm] = cp;
			t.prop_b[LAM][iG - t.maxelm] = lam;
		} // G Side
	}


	if (t.neighbors_for_the_internal_node[G][2] != nullptr) {
		iG = t.neighbors_for_the_internal_node[G][2][iP];
		if (iG >= t.maxelm) {
			// это граничный узел
			rho = 1.1614f; cp = 1005.0f; lam = 0.025f; // инициализация default  dry air 300K 1atm properties
			if (matlist[b[ib].imatid].blibmat == 1) {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // подставляется температура в граничном узле
					// проверка на допустимость температур.
					diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
				} // SOLID
				else if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					float mu, beta_t; // значения не используются но требуются.
					doublereal pressure;
					if (t.ptr[1][iP] == -1) {
						pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
					}
					else {
						if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][2][t.ptr[0][iP]] > -1) {
							pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][2][t.ptr[0][iP]]];
						}
						else {
							// STUB 15.09.2018
							pressure = 0.0; // В твёрдом теле.

#if doubleintprecision == 1
						//printf("error in gran_prop in my_material_properties.c NODE3 G=%lld\n", G);
#else
						//printf("error in gran_prop in my_material_properties.c NODE3 G=%d\n", G);
#endif

						//getchar();
						//system("PAUSE");
						//exit(1);
						}
					}
					my_fluid_properties(t.potent[iG], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				} // FLUID
			}
			else if (matlist[b[ib].imatid].blibmat == 0) {
				// материал определённый пользователем:
				// постоянные свойства.
				rho = matlist[b[ib].imatid].rho;
				//cp = matlist[b[ib].imatid].cp;
				//lam = matlist[b[ib].imatid].lam;
				cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
				lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

				// переносим постоянные свойства материала из внутренней точки области
				// в ближайшую граничную точку.
			}
			// Свойства для внутреннего контрольного объёма.
			t.prop_b[RHO][iG - t.maxelm] = rho;
			t.prop_b[HEAT_CAPACITY][iG - t.maxelm] = cp;
			t.prop_b[LAM][iG - t.maxelm] = lam;
		} // G Side
	}

	if (t.neighbors_for_the_internal_node[G][3] != nullptr) {
		iG = t.neighbors_for_the_internal_node[G][3][iP];
		if (iG >= t.maxelm) {
			// это граничный узел
			rho = 1.1614f; cp = 1005.0f; lam = 0.025f; // инициализация default  dry air 300K 1atm properties
			if (matlist[b[ib].imatid].blibmat == 1) {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // подставляется температура в граничном узле
					// проверка на допустимость температур.
					diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
				} // SOLID
				else if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					float mu, beta_t; // значения не используются но требуются.
					doublereal pressure;
					if (t.ptr[1][iP] == -1) {
						pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
					}
					else {
						if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][3][t.ptr[0][iP]] > -1) {
							pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][3][t.ptr[0][iP]]];
						}
						else {
							// STUB 15.09.2018
							pressure = 0.0; // В твёрдом теле.

#if doubleintprecision == 1
						//printf("error in gran_prop in my_material_properties.c NODE4 G=%lld\n", G);
#else
						//printf("error in gran_prop in my_material_properties.c NODE4 G=%d\n", G);
#endif

						//getchar();
						//system("PAUSE");
						//exit(1);
						}
					}
					my_fluid_properties(t.potent[iG], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				} // FLUID
			}
			else if (matlist[b[ib].imatid].blibmat == 0) {
				// материал определённый пользователем:
				// постоянные свойства.
				rho = matlist[b[ib].imatid].rho;
				//cp = matlist[b[ib].imatid].cp;
				//lam = matlist[b[ib].imatid].lam;
				cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
				lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

				// переносим постоянные свойства материала из внутренней точки области
				// в ближайшую граничную точку.
			}
			// Свойства для внутреннего контрольного объёма.
			t.prop_b[RHO][iG - t.maxelm] = rho;
			t.prop_b[HEAT_CAPACITY][iG - t.maxelm] = cp;
			t.prop_b[LAM][iG - t.maxelm] = lam;
		} // G Side
	}
} // gran_prop

bool bswitch_print_message = true;

// обновление свойств материалов
// для уравнения теплопроводности
void update_temp_properties(TEMPER &t, FLOW* &f, BLOCK* b, int lb, TPROP* matlist) {
	// Свойства материалов, такие как плотность, теплоёмкость и теплопроводность
	// значительно зависят от температуры. В данной функции производится обновление
	// свойств материалов с учётом их зависимости от текущего значения температуры.
	// Если теплопроводность материала падает с ростом температуры, то мы имеем систему
	// с положительной обратной связью в которой расчётная максимальная температура 
	// может быть намного выше чем соответствующая расчётная температура в линейной 
	// системе с постоянными свойствами.

	//TOCHKA p; // точка - центр рассматриваемого КО.	
	
	const integer ISIZE = t.maxelm;

#pragma omp parallel for
	for (integer iP=0; iP<ISIZE; ++iP) {
		// проход по всем внутренним контрольным объёмам расчётной области.

		doublereal dmin = 1.0e30;
		doublereal dmax = -1.0e30;
		dgan = 1.0e30;
		dsic = 1.0e30;
		dcu = 1.0e30;

		//center_cord3D(iP, t.nvtx, t.pa, p); // вычисление координат центра КО.
		//in_model_temp(p,ib,b,lb); // возвращает номер блока ib которому принадлежит контрольный объём с номером iP.
		// 8 января 2016 Ускорение вычисления ib
		

		int ib; // номер блока которому принадлежит контрольный объём.
		float rho, cp, lam;

		ib = t.whot_is_block[iP];
		rho=1.1614f; cp=1005.0f; lam=0.025f; // инициализация default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat==1) {
			// библиотечный, находящийся внутри программы AliceFlow материал.
			if (b[ib].itype== PHYSICS_TYPE_IN_BODY::SOLID) {
			    my_solid_properties(t.potent[iP], rho, cp, lam, matlist[b[ib].imatid].ilibident);
				// проверка на допустимость температур.
				diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
		    } // SOLID
		    if (b[ib].itype== PHYSICS_TYPE_IN_BODY::FLUID) {
			   float mu, beta_t; // значения не используются но требуются.
		       doublereal pressure;
			   if ((t.ptr==nullptr)||(t.ptr[1][iP]==-1)) {
				   pressure=0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
			   }
			   else pressure=f[t.ptr[1][iP]].potent[PRESS][t.ptr[0][iP]];
			   my_fluid_properties(t.potent[iP], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
		    } // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat==0) {
			// материал определённый пользователем:
			// постоянные свойства.
			rho=matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iP]);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iP]);

			// Механические свойства материала, зависящие от температуры.
			t.prop[BETA_T_MECHANICAL][iP]= get_beta_t_solid(matlist[b[ib].imatid].n_beta_t_solid, matlist[b[ib].imatid].temp_beta_t_solid, matlist[b[ib].imatid].arr_beta_t_solid, t.potent[iP]);
			t.prop[YOUNG_MODULE][iP] = get_Young_Module(matlist[b[ib].imatid].n_YoungModule, matlist[b[ib].imatid].temp_Young_Module, matlist[b[ib].imatid].arr_Young_Module, t.potent[iP]);
		    t.prop[POISSON_RATIO][iP]= get_Poisson_ratio(matlist[b[ib].imatid].n_Poisson_ratio, matlist[b[ib].imatid].temp_Poisson_ratio, matlist[b[ib].imatid].arr_Poisson_ratio, t.potent[iP]);
		}
		// Свойства для внутреннего контрольного объёма.
		if (t.ptr != NULL) {
			if (f[0].maxelm == 0) {
				if (lam > dmax) dmax = lam;
				if (lam < dmin) dmin = lam;
			}
			else {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					if (lam > dmax) dmax = lam;
					if (lam < dmin) dmin = lam;
				}
			}
		}
		else {
			if (lam > dmax) dmax = lam;
			if (lam < dmin) dmin = lam;
		}
		t.prop[RHO][iP]=rho;
		t.prop[HEAT_CAPACITY][iP]=cp;
		t.prop[LAM][iP]=lam;

		// Теперь требуется обработать граничные контрольные объёмы,
		// которые соседствуют с данным внутренним КО.
		//
		// 25 сентября 2016 gran_prop теперь работает на АЛИС сетке. 
		gran_prop(t, f, b, lb, iP, E_SIDE, ib, matlist); // East Side
		gran_prop(t, f, b, lb, iP, W_SIDE, ib, matlist); // West Side
		gran_prop(t, f, b, lb, iP, N_SIDE, ib, matlist); // North Side
		gran_prop(t, f, b, lb, iP, S_SIDE, ib, matlist); // South Side
		gran_prop(t, f, b, lb, iP, T_SIDE, ib, matlist); // Top Side
		gran_prop(t, f, b, lb, iP, B_SIDE, ib, matlist); // Bottom Side

		//if (bswitch_print_message) {
			//printf("lam_min=%e lam_max=%e \n", dmin, dmax);
			//if (fabs(dmin - dmax) > 1.0e-10) {
				//std::cout << "thermal conductivity minimum=" << dmin << " thermal conductivity maximum=" << dmax << std::endl;
			//}
		//}
		//if (dgan < 1.0e29) {
			//printf("GaN nonlinear programm library Ok. %e\n",dgan);
			//std::cout << "GaN nonlinear programm library Ok. " << dgan << std::endl;
		//}
		//if (dsic < 1.0e29) {
			//printf("SiC4H nonlinear programm library Ok. %e\n",dsic);
			//std::cout << "SiC4H nonlinear programm library Ok.  " << dsic << std::endl;
		//}
		//if (dcu < 1.0e29) {
			//printf("Cu nonlinear programm library Ok.%e\n",dcu);
			//std::cout << "Cu nonlinear programm library Ok." << dcu << std::endl;
		//}

	}
	
	//getchar();
	bswitch_print_message = !bswitch_print_message;
} // update_temp_properties 

  // обновление свойств материалов
  // для уравнения теплопроводности
void update_temp_properties1(TEMPER &t, FLOW* &f, BLOCK* b, int lb, 
	TPROP* matlist, doublereal* &temperature, integer iadd, integer iu,
	doublereal* &lam_export, int** &nvtx_global) {
	// Свойства материалов, такие как плотность, теплоёмкость и теплопроводность
	// значительно зависят от температуры. В данной функции производится обновление
	// свойств материалов с учётом их зависимости от текущего значения температуры.
	// Если теплопроводность материала падает с ростом температуры, то мы имеем систему
	// с положительной обратной связью в которой расчётная максимальная температура 
	// может быть намного выше чем соответствующая расчётная температура в линейной 
	// системе с постоянными свойствами.
	 
	//TOCHKA p; // точка - центр рассматриваемого КО.
	
	doublereal dmin = 1.0e30;
	doublereal dmax = -1.0e30;
	for (integer iP = iadd; iP<iadd+t.maxelm; iP++) {
		// проход по всем внутренним контрольным объёмам расчётной области.

		
		doublereal Temperature_in_cell = 0.125*(temperature[nvtx_global[0][iP]-1]+ temperature[nvtx_global[1][iP] - 1] + temperature[nvtx_global[2][iP] - 1] + temperature[nvtx_global[3][iP] - 1] + temperature[nvtx_global[4][iP] - 1] + temperature[nvtx_global[5][iP] - 1] + temperature[nvtx_global[6][iP] - 1] + temperature[nvtx_global[7][iP] - 1]);

		int ib; // номер блока которому принадлежит контрольный объём.
		float rho, cp, lam;

		if (iu == -1) {
			//center_cord3D(iP, t.nvtx, t.pa, p); // вычисление координат центра КО.
			//in_model_temp(p,ib,b,lb); // возвращает номер блока ib которому принадлежит контрольный объём с номером iP.
			// 8 января 2016 Ускорение вычисления ib

			ib = t.whot_is_block[iP];

		}
		else {
			ib = t.whot_is_block[iP - iadd];
		}
		rho = 1.1614f; cp = 1005.0f; lam = 0.025f; // инициализация default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat == 1) {
			// библиотечный, находящийся внутри программы AliceFlow материал.
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				my_solid_properties(Temperature_in_cell, rho, cp, lam, matlist[b[ib].imatid].ilibident);
				// проверка на допустимость температур.
				diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
			} // SOLID
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				float mu, beta_t; // значения не используются но требуются.
				doublereal pressure;
				if ((t.ptr==nullptr)||(t.ptr[1][iP-iadd] == -1)) {
					pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
				}
				else pressure = f[t.ptr[1][iP-iadd]].potent[PRESS][t.ptr[0][iP-iadd]];
				my_fluid_properties(Temperature_in_cell, pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
			} // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat == 0) {
			// материал определённый пользователем:
			// постоянные свойства.
			rho = matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, Temperature_in_cell);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, Temperature_in_cell);

			// Загружает механические свойства материала зависящие от температуры.
			t.prop[POISSON_RATIO][iP - iadd] = get_Poisson_ratio(matlist[b[ib].imatid].n_Poisson_ratio, matlist[b[ib].imatid].temp_Poisson_ratio, matlist[b[ib].imatid].arr_Poisson_ratio, Temperature_in_cell);
			t.prop[YOUNG_MODULE][iP - iadd] = get_Young_Module(matlist[b[ib].imatid].n_YoungModule, matlist[b[ib].imatid].temp_Young_Module, matlist[b[ib].imatid].arr_Young_Module, Temperature_in_cell);
			t.prop[BETA_T_MECHANICAL][iP - iadd] = get_beta_t_solid(matlist[b[ib].imatid].n_beta_t_solid, matlist[b[ib].imatid].temp_beta_t_solid, matlist[b[ib].imatid].arr_beta_t_solid, Temperature_in_cell);

		}
		// Свойства для внутреннего контрольного объёма.
		if (t.ptr != NULL) {
			if (f[0].maxelm == 0) {
				if (lam > dmax) dmax = lam;
				if (lam < dmin) dmin = lam;
			}
			else {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					if (lam > dmax) dmax = lam;
					if (lam < dmin) dmin = lam;
				}
			}
		}
		else {
			if (lam > dmax) dmax = lam;
			if (lam < dmin) dmin = lam;
		}
		t.prop[RHO][iP - iadd] = rho;
		t.prop[HEAT_CAPACITY][iP - iadd] = cp;
		t.prop[LAM][iP - iadd] = lam;
		lam_export[iP] = lam;

		// Теперь требуется обработать граничные контрольные объёмы,
		// которые соседствуют с данным внутренним КО.
		//
		// 25 сентября 2016 gran_prop теперь работает на АЛИС сетке. 
		gran_prop(t, f, b, lb, iP - iadd, E_SIDE, ib, matlist); // East Side
		gran_prop(t, f, b, lb, iP - iadd, W_SIDE, ib, matlist); // West Side
		gran_prop(t, f, b, lb, iP - iadd, N_SIDE, ib, matlist); // North Side
		gran_prop(t, f, b, lb, iP - iadd, S_SIDE, ib, matlist); // South Side
		gran_prop(t, f, b, lb, iP - iadd, T_SIDE, ib, matlist); // Top Side
		gran_prop(t, f, b, lb, iP - iadd, B_SIDE, ib, matlist); // Bottom Side

	}
	if (bswitch_print_message) {
		//printf("lam_min=%e lam_max=%e \n", dmin, dmax);
		if (fabs(dmin - dmax) > 1.0e-10) {
			std::cout << "thermal conductivity minimum=" << dmin << " thermal conductivity maximum=" << dmax << std::endl;
		}
	}
	//bswitch_print_message = !bswitch_print_message;
} // update_temp_properties1 

// возвращает коэффициент динамической вязкости заданный пользователем
// 26 января 2012 реализация не Ньютоновских жидкостей.
float return_dynamic_viscosity(TPROP* matlist, integer imatid, 
	                          bool bfirst_start, doublereal SInvariantStrainRateTensor) {
	float mu;

	// Вычисление динамической вязкости:
	switch (matlist[imatid].ilawmu) {
	case 0: // постоянное значение
		mu=matlist[imatid].mu;
		break;
	case 1: // power-law fluid
		// Закон Оствальда-де Вела (Ostwald de Vel)
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
		   mu=(float)(fmax(matlist[imatid].mumin,fmin(2.0*matlist[imatid].Amu*exp((matlist[imatid].degreennmu-1.0)*log(SInvariantStrainRateTensor)),matlist[imatid].mumax)));
		}
		else {
			if (matlist[imatid].degreennmu>1.0) {
				mu=matlist[imatid].mumin;
			}
			else if ((matlist[imatid].degreennmu<1.0)&&(matlist[imatid].degreennmu>0.0)) {
				mu=matlist[imatid].mumax;
			} else mu=matlist[imatid].mu;
		}
		break;
	case 2: // закон Кессона Caisson
		// Согласно закону Кессона вязкость убывает с возрастанием напряжения сдвига.
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
			if (matlist[imatid].Amu>=0.0) {
				mu= (float)(fmax(matlist[imatid].mumin,fmin(((sqrt(matlist[imatid].Amu/SInvariantStrainRateTensor))+matlist[imatid].Bmu)*((sqrt(matlist[imatid].Amu/SInvariantStrainRateTensor))+matlist[imatid].Bmu),matlist[imatid].mumax)));
			}
			else {
				printf("A constant is given incorrectly on the law of the caisson\n");
					printf("for the non-Newtonian fluid");
					printf("please, press any key to exit...\n");
					exit(0);
			}
		}
		else mu=matlist[imatid].mumax;
		break;
	case 3: // закон Прандтля
		// согласно закону Прандтля при B>0 вязкость возрастает с увеличением напряжения сдвига,
		// а при B<0 динамическая вязкость убывает с ростом напряжения сдвига.
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
			if ( (fabs(matlist[imatid].Amu)>1e-30) && (fabs(matlist[imatid].Bmu)>1e-30) && ((SInvariantStrainRateTensor/matlist[imatid].Bmu>-1.0)&&(SInvariantStrainRateTensor/matlist[imatid].Bmu<1.0))) {
				mu= (float)(fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,(matlist[imatid].Amu/SInvariantStrainRateTensor)*asin(SInvariantStrainRateTensor/matlist[imatid].Bmu))));
			}
			else {
				printf("A or B constant is given incorrectly on the law of the Prandtl\n");
				printf("for the non-Newtonian fluid");
				printf("please, press any key to exit...\n");
				exit(0);
			}
		}
		else {
			if (matlist[imatid].Bmu>0.0) {
	    		mu=matlist[imatid].mumin;
			}
			else mu=matlist[imatid].mumax;
		}
		break;
	case 4: // Carreau fluid
		mu= (float)(fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,matlist[imatid].Amu+(matlist[imatid].Bmu-matlist[imatid].Amu)*exp(0.5*(matlist[imatid].degreennmu-1.0)*log(1.0+(matlist[imatid].Cmu*SInvariantStrainRateTensor)*(matlist[imatid].Cmu*SInvariantStrainRateTensor))))));
		break;
	case 5: // Пауэлл-Эйринг
		// Если B*C>0.0 то динамическая вязкость убывает с ростом напряжения сдвига,
		// если B*C<0.0 то динамическая вязкость возрастает с ростом напряжения сдвига.
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
			// Ареа-синус Arsh(x)=log(x+sqrt(x*x+1.0));
			mu= (float)(fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,matlist[imatid].Amu+(matlist[imatid].Bmu/SInvariantStrainRateTensor)*log(matlist[imatid].Cmu*SInvariantStrainRateTensor+sqrt(1.0+(matlist[imatid].Cmu*SInvariantStrainRateTensor)*(matlist[imatid].Cmu*SInvariantStrainRateTensor))))));
		} else 
		{
			if (matlist[imatid].Bmu*matlist[imatid].Cmu>0.0) {
				mu=matlist[imatid].mumax;
			} else mu=matlist[imatid].mumin;
		}
		break;
	case 6: // Williamson non-Newtonian fluid
		// Здесь априори предполагается что параметры формулы заданы правильно (ответственность целиком и полностью ложится на пользователя).
		mu= (float)(fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,matlist[imatid].Amu/(matlist[imatid].Bmu+SInvariantStrainRateTensor)+matlist[imatid].Cmu)));
		break;
	default: // постоянное значение
		mu=matlist[imatid].mu;
		break;
	}

	return mu;
} // return_dynamic_viscosity


// Граничные свойства материалов в задаче гидродинамики
// Меняет свойства материалов в граничном КО.
void gran_prop_flow(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, integer iflow,
	                integer iP, integer G, integer ib, TPROP* matlist, bool bfirst_start) {

       float rho=0.0, mu=0.0, beta_t=0.0;
	   integer iG=0; // номер соседа который может оказаться граничным узлом.
	   iG=f[iflow].neighbors_for_the_internal_node[G][0][iP];
	   if (iG>=f[iflow].maxelm) {
		   // это граничный узел
           rho=1.1614f; mu=1.84e-5f; beta_t=0.003331f; // инициализация default dry air 300K 1atm properties
		   if (matlist[b[ib].imatid].blibmat==1) {
				float cp, lam;
				cp=1005.0f; lam=0.025f;
				// библиотечный находящийся внутри программы AliceFlow материал
				doublereal temperature=20.0;
				if (t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]]>=t.maxelm) {
					// граничная для температуры точка
                    temperature=t.potent[t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]]];
				} else {
					if (t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]] > -1) {
						// температура на грани восстанавливается линейной интерполяцией:
						TOCHKA p1, p2;
						center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1,100); // вычисление координат центра КО.
						center_cord3D(t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]], t.nvtx, t.pa, p2, G);
						doublereal fgplus = 0.0;
						doublereal dx = 0.0, dy = 0.0, dz = 0.0;
						volume3D(f[iflow].ptr[iP], t.nvtx, t.pa, dx, dy, dz);
						switch (G) {
						case E_SIDE: case W_SIDE: fgplus = dx / fabs(p1.x - p2.x); break;
						case N_SIDE: case S_SIDE: fgplus = dy / fabs(p1.y - p2.y); break;
						case T_SIDE: case B_SIDE: fgplus = dz / fabs(p1.z - p2.z); break;
						}

						temperature = (1.0 - fgplus)*t.potent[f[iflow].ptr[iP]] + fgplus*t.potent[t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]]];
					}
					else {
						temperature = t.potent[iP];
#if doubleintprecision == 1
						//printf("error in gran_prop_flow in my_material_properties.c NODE1 G=%lld\n", G);
#else
						//printf("error in gran_prop_flow in my_material_properties.c NODE1 G=%d\n", G);
#endif
						
						//getchar();
						//system("PAUSE");
						//exit(1);
					}
				}
                // библиотечный находящийся внутри программы AliceFlow материал
				if (b[ib].itype== PHYSICS_TYPE_IN_BODY::FLUID) {
					my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				}
		   }
		   else if (matlist[b[ib].imatid].blibmat==0) {
				// материал определённый пользователем:
				// постоянные свойства.

			    // плотность
				rho=matlist[b[ib].imatid].rho;
				// динамическая вязкость
                mu=return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				// коэффициент линейного температурного расширения.
				beta_t=matlist[b[ib].imatid].beta_t;
				// переносим постоянные свойства материала из внутренней точки области
				// в ближайшую граничную точку.
			}
		    // Свойства граничного контрольного объёма:
		    f[iflow].prop_b[RHO][iG-f[iflow].maxelm]=rho;
		    f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG-f[iflow].maxelm]=mu;
		    f[iflow].prop_b[BETA_T][iG-f[iflow].maxelm]=beta_t;

	   } // G Side

	   if (f[iflow].neighbors_for_the_internal_node[G][1] != nullptr) {
		   iG = f[iflow].neighbors_for_the_internal_node[G][1][iP];
		   if (iG >= f[iflow].maxelm) {
			   // это граничный узел
			   rho = 1.1614f; mu = 1.84e-5f; beta_t = 0.003331f; // инициализация default dry air 300K 1atm properties
			   if (matlist[b[ib].imatid].blibmat == 1) {
				    float  cp, lam;
				   cp = 1005; lam = 0.025;
				   // библиотечный находящийся внутри программы AliceFlow материал
				   doublereal temperature = 20.0;
				   if (t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]] >= t.maxelm) {
					   // граничная для температуры точка
					   temperature = t.potent[t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]]];
				   }
				   else {
					   if (t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]] > -1) {
						   // температура на грани восстанавливается линейной интерполяцией:
						   TOCHKA p1, p2;
						   center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1, 100); // вычисление координат центра КО.
						   center_cord3D(t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]], t.nvtx, t.pa, p2, G);
						   doublereal fgplus = 0.0;
						   doublereal dx = 0.0, dy = 0.0, dz = 0.0;
						   volume3D(f[iflow].ptr[iP], t.nvtx, t.pa, dx, dy, dz);
						   switch (G) {
						   case E_SIDE: case W_SIDE: fgplus = dx / fabs(p1.x - p2.x); break;
						   case N_SIDE: case S_SIDE: fgplus = dy / fabs(p1.y - p2.y); break;
						   case T_SIDE: case B_SIDE: fgplus = dz / fabs(p1.z - p2.z); break;
						   }

						   temperature = (1.0 - fgplus) * t.potent[f[iflow].ptr[iP]] + fgplus * t.potent[t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]]];
					   }
					   else {
						   temperature = t.potent[iP];
#if doubleintprecision == 1
						   // printf("error in gran_prop_flow in my_material_properties.c NODE2 G=%lld\n", G);
#else
						   //printf("error in gran_prop_flow in my_material_properties.c NODE2 G=%d\n", G);
#endif

					   //getchar();
					   //system("PAUSE");
					   //exit(1);
					   }
				   }
				   // библиотечный находящийся внутри программы AliceFlow материал
				   if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					   my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				   }
			   }
			   else if (matlist[b[ib].imatid].blibmat == 0) {
				   // материал определённый пользователем:
				   // постоянные свойства.

				   // плотность
				   rho = matlist[b[ib].imatid].rho;
				   // динамическая вязкость
				   mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				   // коэффициент линейного температурного расширения.
				   beta_t = matlist[b[ib].imatid].beta_t;
				   // переносим постоянные свойства материала из внутренней точки области
				   // в ближайшую граничную точку.
			   }
			   // Свойства граничного контрольного объёма:
			   f[iflow].prop_b[RHO][iG - f[iflow].maxelm] = rho;
			   f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG - f[iflow].maxelm] = mu;
			   f[iflow].prop_b[BETA_T][iG - f[iflow].maxelm] = beta_t;

		   } // G Side
	   }

	   if (f[iflow].neighbors_for_the_internal_node[G][2] != nullptr) {
		   iG = f[iflow].neighbors_for_the_internal_node[G][2][iP];
		   if (iG >= f[iflow].maxelm) {
			   // это граничный узел
			   rho = 1.1614f; mu = 1.84e-5f; beta_t = 0.003331f; // инициализация default dry air 300K 1atm properties
			   if (matlist[b[ib].imatid].blibmat == 1) {
				   float cp, lam;
				   cp = 1005.0f; lam = 0.025f;
				   // библиотечный находящийся внутри программы AliceFlow материал
				   doublereal temperature = 20.0;
				   if (t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]] >= t.maxelm) {
					   // граничная для температуры точка
					   temperature = t.potent[t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]]];
				   }
				   else {
					   if (t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]] > -1) {
						   // температура на грани восстанавливается линейной интерполяцией:
						   TOCHKA p1, p2;
						   center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1, 100); // вычисление координат центра КО.
						   center_cord3D(t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]], t.nvtx, t.pa, p2, G);
						   doublereal fgplus = 0.0;
						   doublereal dx = 0.0, dy = 0.0, dz = 0.0;
						   volume3D(f[iflow].ptr[iP], t.nvtx, t.pa, dx, dy, dz);
						   switch (G) {
						   case E_SIDE: case W_SIDE: fgplus = dx / fabs(p1.x - p2.x); break;
						   case N_SIDE: case S_SIDE: fgplus = dy / fabs(p1.y - p2.y); break;
						   case T_SIDE: case B_SIDE: fgplus = dz / fabs(p1.z - p2.z); break;
						   }

						   temperature = (1.0 - fgplus) * t.potent[f[iflow].ptr[iP]] + fgplus * t.potent[t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]]];
					   }
					   else {
						   temperature = t.potent[iP];

#if doubleintprecision == 1
						   // printf("error in gran_prop_flow in my_material_properties.c NODE3 G=%lld\n", G);
#else
						   //printf("error in gran_prop_flow in my_material_properties.c NODE3 G=%d\n", G);
#endif

					   //getchar();
					   //system("PAUSE");
					   //exit(1);
					   }
				   }
				   // библиотечный находящийся внутри программы AliceFlow материал
				   if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					   my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				   }
			   }
			   else if (matlist[b[ib].imatid].blibmat == 0) {
				   // материал определённый пользователем:
				   // постоянные свойства.

				   // плотность
				   rho = matlist[b[ib].imatid].rho;
				   // динамическая вязкость
				   mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				   // коэффициент линейного температурного расширения.
				   beta_t = matlist[b[ib].imatid].beta_t;
				   // переносим постоянные свойства материала из внутренней точки области
				   // в ближайшую граничную точку.
			   }
			   // Свойства граничного контрольного объёма:
			   f[iflow].prop_b[RHO][iG - f[iflow].maxelm] = rho;
			   f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG - f[iflow].maxelm] = mu;
			   f[iflow].prop_b[BETA_T][iG - f[iflow].maxelm] = beta_t;

		   } // G Side
	   }

	   if (f[iflow].neighbors_for_the_internal_node[G][3] != nullptr) {
		   iG = f[iflow].neighbors_for_the_internal_node[G][3][iP];
		   if (iG >= f[iflow].maxelm) {
			   // это граничный узел
			   rho = 1.1614f; mu = 1.84e-5f; beta_t = 0.003331f; // инициализация default dry air 300K 1atm properties
			   if (matlist[b[ib].imatid].blibmat == 1) {
				   float cp, lam;
				   cp = 1005; lam = 0.025;
				   // библиотечный находящийся внутри программы AliceFlow материал
				   doublereal temperature = 20.0;
				   if (t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]] >= t.maxelm) {
					   // граничная для температуры точка
					   temperature = t.potent[t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]]];
				   }
				   else {
					   if (t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]] > -1) {
						   // температура на грани восстанавливается линейной интерполяцией:
						   TOCHKA p1, p2;
						   center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1, 100); // вычисление координат центра КО.
						   center_cord3D(t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]], t.nvtx, t.pa, p2, G);
						   doublereal fgplus = 0.0;
						   doublereal dx = 0.0, dy = 0.0, dz = 0.0;
						   volume3D(f[iflow].ptr[iP], t.nvtx, t.pa, dx, dy, dz);
						   switch (G) {
						   case E_SIDE: case W_SIDE: fgplus = dx / fabs(p1.x - p2.x); break;
						   case N_SIDE: case S_SIDE: fgplus = dy / fabs(p1.y - p2.y); break;
						   case T_SIDE: case B_SIDE: fgplus = dz / fabs(p1.z - p2.z); break;
						   }

						   temperature = (1.0 - fgplus) * t.potent[f[iflow].ptr[iP]] + fgplus * t.potent[t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]]];
					   }
					   else {
						   temperature = t.potent[iP];
#if doubleintprecision == 1
						   //printf("error in gran_prop_flow in my_material_properties.c NODE4 G=%lld\n", G);
#else
						   //printf("error in gran_prop_flow in my_material_properties.c NODE4 G=%d\n", G);
#endif

					   //getchar();
					   //system("PAUSE");
					   //exit(1);
					   }
				   }
				   // библиотечный находящийся внутри программы AliceFlow материал
				   if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					   my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				   }
			   }
			   else if (matlist[b[ib].imatid].blibmat == 0) {
				   // материал определённый пользователем:
				   // постоянные свойства.

				   // плотность
				   rho = matlist[b[ib].imatid].rho;
				   // динамическая вязкость
				   mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				   // коэффициент линейного температурного расширения.
				   beta_t = matlist[b[ib].imatid].beta_t;
				   // переносим постоянные свойства материала из внутренней точки области
				   // в ближайшую граничную точку.
			   }
			   // Свойства граничного контрольного объёма:
			   f[iflow].prop_b[RHO][iG - f[iflow].maxelm] = rho;
			   f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG - f[iflow].maxelm] = mu;
			   f[iflow].prop_b[BETA_T][iG - f[iflow].maxelm] = beta_t;

		   } // G Side
	   }
} // gran_prop_flow

// Обновление свойств материалов для гидродинамики
void update_flow_properties(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, integer flow_interior, TPROP* matlist, bool bfirst_start) {

	// если жидких зон нет то никакого обновления свойств материалов происходить не будет.
	for (integer ifi=0; ifi<flow_interior; ifi++) {
		// проход по всем жидким зонам.
		
		//TOCHKA p; // точка центр рассматриваемого КО.
		
		doublereal dmin = 1.0e30;
		doublereal dmax = -1.0e30;

		for (integer iP=0; iP<f[ifi].maxelm; iP++) {


			integer ib; // номер блока которому принадлежит рассматриваемый контрольный объём.
		    float rho, mu, beta_t;

			// проход по всем внутренним контрольным объёмам
			//center_cord3D(iP, f[ifi].nvtx, f[ifi].pa, p,100); // вычисление координат центра КО.
			//in_model_flow(p, ib, b, lb); // возвращает номер блока ib которому принадлежит контрольный объём с номером iP.
			// Более быстрая операция - индексирование в hash таблице.
			ib=t.whot_is_block[f[ifi].ptr[iP]];

			rho=1.1614f; mu=1.84e-5f; beta_t=0.003331f; // инициализация default dry air 300K 1atm properties

			if (ib > -1) {
				if (matlist[b[ib].imatid].blibmat == 1) {
					float cp, lam;
					cp = 1005; lam = 0.025;
					// библиотечный находящийся внутри программы AliceFlow материал
					if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						//printf("rho=%e, mu=%e",rho,mu); // debug
						//getchar();
						my_fluid_properties(t.potent[f[ifi].ptr[iP]], f[ifi].potent[PRESS][iP], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
						//printf("rho=%e, mu=%e",rho,mu); // debug
						//getchar();
					}
				}
				else if (matlist[b[ib].imatid].blibmat == 0) {
					// материал определённый пользователем:
					// постоянные свойства.
					rho = matlist[b[ib].imatid].rho;
					// Вычисление динамической вязкости:
					// в том числе и для неньютоновских жидкостей.
					mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[ifi].SInvariantStrainRateTensor[iP]);
					// коэффициент линейного температурного расширения
					beta_t = matlist[b[ib].imatid].beta_t;
				}
			}
			// Свойства внутреннего контрольного объёма:
			if (mu > dmax) dmax = mu;
			if (mu < dmin) dmin = mu;
			f[ifi].prop[RHO][iP]=rho;
			f[ifi].prop[MU_DYNAMIC_VISCOSITY][iP]=mu;
			f[ifi].prop[BETA_T][iP]=beta_t;
			//printf("rho=%e, mu=%e",f[ifi].prop[RHO][iP],f[ifi].prop[MU][iP]); // debug
			//getchar();

			// Теперь требуется обработать граничные контрольные объёмы,
			// которые соседствуют с данным внутренним КО.
			//
			gran_prop_flow(t, f, b, lb, ifi, iP, E_SIDE, ib, matlist, bfirst_start); // East Side
            gran_prop_flow(t, f, b, lb, ifi, iP, W_SIDE, ib, matlist, bfirst_start); // West Side
			gran_prop_flow(t, f, b, lb, ifi, iP, N_SIDE, ib, matlist, bfirst_start); // North Side
            gran_prop_flow(t, f, b, lb, ifi, iP, S_SIDE, ib, matlist, bfirst_start); // South Side
			gran_prop_flow(t, f, b, lb, ifi, iP, T_SIDE, ib, matlist, bfirst_start); // Top Side
            gran_prop_flow(t, f, b, lb, ifi, iP, B_SIDE, ib, matlist, bfirst_start); // Bottom Side
		}

		if (!b_on_adaptive_local_refinement_mesh) {
			//printf("\nmu_min=%e mu_max=%e \n", dmin, dmax);
			if (fabs(dmin - dmax) > 1.0e-10) {
				std::cout << std::endl << "dynamic viscosity minimum=" << dmin << " dynamic viscosity maximum=" << dmax << std::endl;
			}
		}
	}
} // update_flow_properties

// Экспортирует внешнюю поверхность геометрии пользователя в .stl формате.
// 09.09.2019.
void export_User_Geom_in_STL_format(TEMPER& t) {

	FILE* fp = NULL;
	

#ifdef MINGW_COMPILLER
	int  err = 0;
	fp = fopen64("user_geom.stl", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "user_geom.stl", "w");
#endif

	if (err != 0) {
		printf("Error open file user_geom.stl\n");
		system("pause");
		exit(0);
	}

	if (fp != NULL) {
		fprintf(fp,"solid user_geom.stl\n");

		for (integer iP = 0; iP < t.maxelm; iP++) {
			if (t.neighbors_for_the_internal_node[T_SIDE][0][iP] >= t.maxelm) {
				// По часовой стрелке

				// Пишем две треугольных грани
				integer ind;
				fprintf(fp, "facet normal %e %e %e\n",0.0,0.0,1.0);
				fprintf(fp, "outer loop\n");
				ind = 4;
				fprintf(fp, "vertex %e %e %e\n",t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 5;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 7;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");


				fprintf(fp, "facet normal %e %e %e\n", 0.0, 0.0, 1.0);
				fprintf(fp, "outer loop\n");
				ind = 4;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 7;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 6;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");
			}

			if (t.neighbors_for_the_internal_node[B_SIDE][0][iP] >= t.maxelm) {
				// Против часовой стрелки.

				// Пишем две треугольных грани
				integer ind;
				fprintf(fp, "facet normal %e %e %e\n", 0.0, 0.0, -1.0);
				fprintf(fp, "outer loop\n");
				ind = 0;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 3;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 1;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");


				fprintf(fp, "facet normal %e %e %e\n", 0.0, 0.0, -1.0);
				fprintf(fp, "outer loop\n");
				ind = 2;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 3;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 0;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");
			}

			if (t.neighbors_for_the_internal_node[E_SIDE][0][iP] >= t.maxelm) {
				// Вращение по часовой стрелке.

				// Пишем две треугольных грани
				integer ind;
				fprintf(fp, "facet normal %e %e %e\n", 1.0, 0.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 3;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 5;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 1;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");


				fprintf(fp, "facet normal %e %e %e\n", 1.0, 0.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 3;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 7;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 5;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");
			}

			if (t.neighbors_for_the_internal_node[W_SIDE][0][iP] >= t.maxelm) {
				// Вращение против часовой стрелки.

				// Пишем две треугольных грани
				integer ind;
				fprintf(fp, "facet normal %e %e %e\n", -1.0, 0.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 0;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 4;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 6;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");


				fprintf(fp, "facet normal %e %e %e\n", -1.0, 0.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 6;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 2;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 0;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");
			}


			if (t.neighbors_for_the_internal_node[N_SIDE][0][iP] >= t.maxelm) {
				// Вращение по часовой стрелке.

				// Пишем две треугольных грани
				integer ind;
				fprintf(fp, "facet normal %e %e %e\n", 0.0, 1.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 6;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 3;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 2;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");


				fprintf(fp, "facet normal %e %e %e\n", 0.0, 1.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 6;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 7;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 3;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");
			}

			if (t.neighbors_for_the_internal_node[S_SIDE][0][iP] >= t.maxelm) {
				// Вращение против часовой стрелки.

				// Пишем две треугольных грани
				integer ind;
				fprintf(fp, "facet normal %e %e %e\n", 0.0, -1.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 0;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 1;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 5;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");


				fprintf(fp, "facet normal %e %e %e\n", 0.0, -1.0, 0.0);
				fprintf(fp, "outer loop\n");
				ind = 5;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 4;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				ind = 0;
				fprintf(fp, "vertex %e %e %e\n", t.pa[t.nvtx[ind][iP] - 1].x, t.pa[t.nvtx[ind][iP] - 1].y, t.pa[t.nvtx[ind][iP] - 1].z);
				fprintf(fp, "endloop\n");
				fprintf(fp, "endfacet\n");
			}
		}

		fprintf(fp, "endsolid user_geom.stl");

		fclose(fp);
	}
} // export_User_Geom_in_STL_format

// Вычисляет массу рассчитываемой модели.
// 12.03.2017
// Вычисляет объем рассчитываемой модели.
// 24.06.2021
doublereal massa_cabinet(TEMPER &t, FLOW* &f, 
	int &flow_interior,
	BLOCK* b, int lb, doublereal temp_ref,
	TPROP* matlist) {

	doublereal massa = 0.0;
	doublereal volume_model = 0.0;

	for (integer iP = 0; iP < t.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		//volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
		volume3D_q(iP, t.nvtx, t.pa, dx, dy, dz);
		int ib = t.whot_is_block[iP];
	    float rho, cp, lam;
		rho = 1.1614f; cp = 1005.0f; lam = 0.025f; // инициализация default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat == 1) {
			// библиотечный, находящийся внутри программы AliceFlow материал.
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				my_solid_properties(t.potent[iP], rho, cp, lam, matlist[b[ib].imatid].ilibident);
				// проверка на допустимость температур.
				diagnostic_critical_temperature(t.potent[iP],f,t,b,lb);
			} // SOLID
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				float mu, beta_t; // значения не используются но требуются.
				doublereal pressure;
				if (t.ptr[1][iP] == -1) {
					pressure = 0.0; // давление внутри твёрдого тела (этого не может быть, т.к. здесь обязательно жидкость).
				}
				else pressure = f[t.ptr[1][iP]].potent[PRESS][t.ptr[0][iP]];
				my_fluid_properties(t.potent[iP], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
			} // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat == 0) {
			// материал определённый пользователем:
			// постоянные свойства.
			rho = matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iP]);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iP]);

		}
		// Свойства для внутреннего контрольного объёма.
		massa += rho*dx*dy*dz;
		volume_model += dx * dy * dz; // объём модели в кубических метрах.
	}

	//printf("massa=%1.3f kg\n", massa);
	std::cout << "massa=" << massa << " kg" << std::endl;
	std::cout << "volume model=" << volume_model << " m^3" << std::endl;

	return massa;
}


// Печатает расход через выходную границу потока.
// Печатает сколько Вт тепла покидает расчётную область через выходную границу потока.
// 28.10.2019
void report_out_boundary(FLOW &f, TEMPER &t, int ls, int lw, WALL* &w, BLOCK* &b, int lb,
	TPROP* &matlist, doublereal Tamb) 
{

	integer* number_control_volume_on_wall = new integer[lw];
	doublereal* wall_power = new doublereal[lw];
	for (int iwall_scan = 0; iwall_scan < lw; iwall_scan++) {
		wall_power[iwall_scan] = 0.0;
		number_control_volume_on_wall[iwall_scan] = 0;
		for (integer j = 0; j < t.maxbound; j++) {

			if (t.border_neighbor[j].MCB == (ls + iwall_scan)) {

				number_control_volume_on_wall[iwall_scan]++;

				integer iP = t.border_neighbor[j].iI;// f.maxelm + j;

				TOCHKA p;
				center_cord3D(iP, t.nvtx, t.pa, p, 100);
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;
				volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
				dx = fabs(dx);
				dy = fabs(dy);
				dz = fabs(dz);
				doublereal dx1 = 0.0, dy1 = 0.0, dz1 = 0.0;
				if (t.border_neighbor[j].iII > -1) {
					volume3D(t.border_neighbor[j].iII, t.nvtx, t.pa, dx1, dy1, dz1);
					dx1 = fabs(dx1);
					dy1 = fabs(dy1);
					dz1 = fabs(dz1);
				}
				int ib; // номер искомого блока
				in_model_temp(p, ib, b, lb);

				doublereal lam= t.prop[LAM][iP]; // значения не используются но требуются
				doublereal temperature_i = t.potent[iP]; // но на самом деле давление требуется с предыдущего временного слоя.
				doublereal temperature_ii = temperature_i;
				if (t.border_neighbor[j].iII > -1) {
					temperature_ii = t.potent[t.border_neighbor[j].iII];
				}
				doublereal temperature_w=t.potent[t.border_neighbor[j].iB];

				switch (w[iwall_scan].iPlane) {
				case XY_PLANE: if (t.border_neighbor[j].Norm == T_SIDE) {// Низ, внутренняя номаль.
					//+ втекает
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dy) / (0.5 * (dz + dz1))) >
						fabs((lam * (temperature_i-temperature_w) * dx * dy) / (0.5 * (dz)))) {
					wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dy) / (0.5 * (dz + dz1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dy) / (0.5 * (dz));
					}
				}
				if (t.border_neighbor[j].Norm == B_SIDE) {// Верх, внутренняя номаль.
					//+ втекает
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dy) / (0.5 * (dz + dz1)))>
						fabs((lam * (temperature_i-temperature_w) * dx * dy) / (0.5 * (dz)))) {
					    wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dy) / (0.5 * (dz + dz1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dy) / (0.5 * (dz));
					}
				}
				break;
				case XZ_PLANE: 
				if (t.border_neighbor[j].Norm == N_SIDE) {// Юг, внутренняя номаль.
					//+ втекает
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1)))>
						fabs((lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy)))) {
					    wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy));
					}
				}
				if (t.border_neighbor[j].Norm == S_SIDE) {// Север, внутренняя номаль.
					//+ втекает
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1))) >
						fabs((lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy)))) {
					wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy));
					}
				}
				break;
				case YZ_PLANE:  if (t.border_neighbor[j].Norm == E_SIDE) {// запад, внутренняя номаль.
					//+ втекает
					if (fabs((lam * (temperature_ii - temperature_i) * dy * dz) / (0.5 * (dx + dx1)))>
						fabs((lam * (temperature_i-temperature_w) * dy * dz) / (0.5 * (dx)))) {
						wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dy * dz) / (0.5 * (dx + dx1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dy * dz) / (0.5 * (dx));
					}
				}
				if (t.border_neighbor[j].Norm == W_SIDE) {// Восток, внутренняя номаль.
					//+ втекает
					if (fabs((lam * (temperature_ii - temperature_i) * dy * dz) / (0.5 * (dx + dx1))) >
						fabs((lam * (temperature_i-temperature_w) * dy * dz) / (0.5 * (dx)))) {
						wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dy * dz) / (0.5 * (dx + dx1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dy * dz) / (0.5 * (dx));
					}
				}
				break;
				}
			}
		}
	}

	printf("\n");
	for (int iwall_scan = 0; iwall_scan < lw; iwall_scan++) {
		//printf("wall[%d].name = %s power is %e W. Number control volume in wall=%lld\n", iwall_scan, w[iwall_scan].name,  wall_power[iwall_scan], number_control_volume_on_wall[iwall_scan]);
		std::cout << "wall[" << iwall_scan << "].name = " << w[iwall_scan].name << " power is " << wall_power[iwall_scan] << " W. Number control volume in wall=" << number_control_volume_on_wall[iwall_scan] << std::endl;
	}

	delete[] number_control_volume_on_wall;
	delete[] wall_power;

	int idlw = 0;
	if (lw > 0) {
		for (int iwall_scan = 0; iwall_scan < lw; iwall_scan++) {
			if (w[iwall_scan].bpressure) idlw = 1;
			if (w[iwall_scan].bopening) idlw = 1;
		}
	}

	if (idlw > 0) {
		printf("\n");
		printf("\n        rashod m!3/s;  rashod kg/s;  Power out, W;  type wall out; delta Tavg_wall, oC; Tmax_wall, oC.\n");
	}

	doublereal Tamb0 = 1.0e30;// Tamb; // Tamb
	for (int iwall_scan = 0; iwall_scan < lw; iwall_scan++) {
		// Определяем минимальную заданную температуру.
		if ((!w[iwall_scan].bpressure) && (!w[iwall_scan].bsymmetry)) {
			if (w[iwall_scan].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
				Tamb0 = fmin(Tamb0, w[iwall_scan].Tamb);
			}
		}
	}
	if (Tamb0 > 1.0e10) Tamb0 = Tamb;

	float cp=0.0;

	for (int iwall_scan = 0; iwall_scan < lw; iwall_scan++) {

		doublereal Tmax_wall = -1.0e30;

		//if (iwall_scan == 0) printf("\n");

		if (w[iwall_scan].bpressure) {
			doublereal rashod = 0.0; // m!3/s
			doublereal rashod2 = 0.0; // kg/s
									  // Мощность в Вт которая уносится через выходную границу потока.
									  // Формулу нашел Ионов В.Е. Позволяет найти мощность теплосъёма
									  // если внутри области заданы только условия Дирихле по температуре.
			doublereal Qout = 0.0; // Вт
			for (integer j = 0; j < f.maxbound; j++) {

				if (f.border_neighbor[j].MCB == (ls + iwall_scan)) {

				integer iP = f.border_neighbor[j].iI;// f.maxelm + j;

				TOCHKA p;
				center_cord3D(iP, f.nvtx, f.pa, p, 100);
				int ib; // номер искомого блока
				in_model_flow(p, ib, b, lb);

				float rho, mu, beta_t, lam; // значения не используются но требуются
				doublereal pressure = f.potent[PRESS][iP]; // но на самом деле давление требуется с предыдущего временного слоя.
				my_fluid_properties(t.potent[f.ptr[iP]], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				Tmax_wall = fmax(Tmax_wall, t.potent[f.ptr[iP]]);

				// Внимание!!! f.prop_b[HEAT_CAPACITY][j] использовать нельзя, т.к. для fluid это не определено.

				
					switch (w[iwall_scan].iPlane) {
					case XY_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == T_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);

						}
						else {
							rashod += f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					case XZ_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == N_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						else {
							rashod += f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					case YZ_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == E_SIDE) {
							rashod += -f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);

						}
						else {
							rashod += f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					}
				}
			}

			//printf("wall[%lld] out of boundary (bpressure==true).\n", iwall_scan);
			//printf("rashod = %e m!3/s; rashod = %e kg/s; Power out = %e W, Tavg_wall= %e, Tmax_wall= %e. \n", rashod, rashod2, Qout, Qout / (cp*rashod2), Tmax_wall);
			if (rashod > 0) {
				printf("wall[%d]  %e  %e %e bpressure %e %e\n", iwall_scan, rashod, rashod2, Qout, fabs(Qout/(cp*rashod2)), Tmax_wall);
			}
			else {
				printf("wall[%d] %e %e  %e bpressure %e %e\n", iwall_scan, rashod, rashod2, Qout, fabs(Qout / (cp*rashod2)), Tmax_wall);
			}
		}
		if (w[iwall_scan].bopening) {
			doublereal rashod = 0.0; // m!3/s
			doublereal rashod2 = 0.0; // kg/s
									  // Мощность в Вт которая уносится через выходную границу потока.
									  // Формулу нашел Ионов В.Е. Позволяет найти мощность теплосъёма
									  // если внутри области заданы только условия Дирихле по температуре.
			doublereal Qout = 0.0; // Вт
			for (integer j = 0; j < f.maxbound; j++) {

				if (f.border_neighbor[j].MCB == (ls + iwall_scan)) {

				integer iP = f.border_neighbor[j].iI;// f.maxelm + j;

				

				TOCHKA p;
				center_cord3D(iP, f.nvtx, f.pa, p, 100);
				int ib; // номер искомого блока
				in_model_flow(p, ib, b, lb);

				float rho, mu, beta_t, lam; // значения не используются но требуются
				doublereal pressure = f.potent[PRESS][iP]; // но на самом деле давление требуется с предыдущего временного слоя.
				my_fluid_properties(t.potent[f.ptr[iP]], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				Tmax_wall = fmax(Tmax_wall, t.potent[f.ptr[iP]]);

				// 1 - minx, 2-maxx; 3 - miny, 4-maxy; 5 - minz, 6-maxz;
				//if ((iwall_scan == 2)||(iwall_scan == 5)) {
					//printf("iP=%ld rho=%e cp=%e ptr=%d T=%e iPlane=%d Norm=%d\n", iP,rho,cp, f.ptr[iP], t.potent[f.ptr[iP]], w[iwall_scan].iPlane, f.border_neighbor[j].Norm);
					//printf("xE=%e xS=%e\n", w[iwall_scan].g.xE, w[iwall_scan].g.xS);
					//getchar();
				//}

				// Внимание!!! f.prop_b[HEAT_CAPACITY][j] использовать нельзя, т.к. для fluid это не определено.

				//printf("HEAT_CAPACITY=%e\n", cp); // debug
				//system("pause"); // debug
				
					switch (w[iwall_scan].iPlane) {
					case XY_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == T_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);

						}
						else {
							//printf("VZ=%e dS=%e BSIDE internal normal TOP boundary\n", f.potent[VZCOR][f.maxelm + j], f.border_neighbor[j].dS);
							//if (j % 10 == 0) system("pause");
							rashod += f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VZCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					case XZ_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == N_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						else {
							rashod += f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VYCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					case YZ_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == E_SIDE) {
							rashod += -f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);

						}
						else {
							rashod += f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS*f.potent[VXCOR][f.maxelm + j] *
								cp * ( t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					}
				}
			}

			//printf("wall[%lld] out of boundary (bopening==true).\n", iwall_scan);
			//printf("rashod = %e m!3/s; rashod = %e kg/s; Power out = %e W, Tavg_wall= %e, Tmax_wall= %e. \n", rashod, rashod2, Qout, Qout / (cp*rashod2), Tmax_wall);
			if (rashod > 0) {
				printf("wall[%d]  %e  %e %e bopening  %e %e\n", iwall_scan, rashod, rashod2, Qout, fabs(Qout / (cp*rashod2)), Tmax_wall);
			}
			else {
				printf("wall[%d] %e %e  %e bopening  %e %e\n", iwall_scan, rashod, rashod2, Qout, fabs(Qout / (cp*rashod2)), Tmax_wall);
			}
		}
	
		if ((fabs(w[iwall_scan].Vx)>1.0e-20)||(fabs(w[iwall_scan].Vy) > 1.0e-20)||(fabs(w[iwall_scan].Vz) > 1.0e-20)) {
			// Входная граница потока.

			doublereal rashod = 0.0; // m!3/s
			doublereal rashod2 = 0.0; // kg/s
									  // Мощность в Вт которая уносится через выходную границу потока.
									  // Формулу нашел Ионов В.Е. Позволяет найти мощность теплосъёма
									  // если внутри области заданы только условия Дирихле по температуре.
			doublereal Qout = 0.0; // Вт
			for (integer j = 0; j < f.maxbound; j++) {

				if (f.border_neighbor[j].MCB == (ls + iwall_scan)) {

					integer iP = f.border_neighbor[j].iI;// f.maxelm + j;

					TOCHKA p;
					center_cord3D(iP, f.nvtx, f.pa, p, 100);
					int ib; // номер искомого блока
					in_model_flow(p, ib, b, lb);

					float rho, mu, beta_t, lam; // значения не используются но требуются
					doublereal pressure = f.potent[PRESS][iP]; // но на самом деле давление требуется с предыдущего временного слоя.
					my_fluid_properties(t.potent[f.ptr[iP]], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
					Tmax_wall = fmax(Tmax_wall, t.potent[f.ptr[iP]]);

					// Внимание!!! f.prop_b[HEAT_CAPACITY][j] использовать нельзя, т.к. для fluid это не определено.


					switch (w[iwall_scan].iPlane) {
					case XY_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == T_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);

						}
						else {
							rashod += f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VZCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					case XZ_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == N_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						else {
							rashod += f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VYCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					case YZ_PLANE:
						// Нормаль внутренняя
						// Расход: то что вытекает то с плюсом.
						if (f.border_neighbor[j].Norm == E_SIDE) {
							rashod += -f.border_neighbor[j].dS * f.potent[VXCOR][f.maxelm + j];
							rashod2 += -f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VXCOR][f.maxelm + j];
							Qout += f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VXCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);

						}
						else {
							rashod += f.border_neighbor[j].dS * f.potent[VXCOR][f.maxelm + j];
							rashod2 += f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VXCOR][f.maxelm + j];
							Qout += -f.prop_b[RHO][j] * f.border_neighbor[j].dS * f.potent[VXCOR][f.maxelm + j] *
								cp * (t.potent[f.ptr[f.border_neighbor[j].iI]] - Tamb0);
						}
						break;
					}
				}
			}


			//printf("wall[%lld] out of boundary (bpressure==true).\n", iwall_scan);
			//printf("rashod = %e m!3/s; rashod = %e kg/s; Power out = %e W, Tavg_wall= %e, Tmax_wall= %e. \n", rashod, rashod2, Qout, Qout / (cp*rashod2), Tmax_wall);
			if (rashod > 0) {
				printf("wall[%d]  %e  %e %e inflow %e %e\n", iwall_scan, rashod, rashod2, Qout, fabs(Qout / (cp * rashod2)), Tmax_wall);
			}
			else {
				printf("wall[%d] %e %e  %e inflow %e %e\n", iwall_scan, rashod, rashod2, Qout, fabs(Qout / (cp * rashod2)), Tmax_wall);
			}

		}

    }
	if (idlw > 0) {
		printf("\n");
	}
} // report_out_boundary

#endif
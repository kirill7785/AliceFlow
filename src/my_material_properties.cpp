// 24 ������ 2011 ���� � ������
// 25 ������ 2011 ���� 15 ������ ���.
// 26-27 ������ 2011 ���� ��������� ���������� AliceMeshv0_08.
// 28 ������ 2011 ���� ��������� update_temp_properties ������
// ��� ���������� ����������� �������.
// � ������ my_material_properties.c ����������
// ����������������� ���������� ���������� ������� ����������.

#ifndef  MY_MATERIAL_PROPERTIES_CPP
#define  MY_MATERIAL_PROPERTIES_CPP 1

// ���������� ������ ������������ ���������� ���������� � ����� ������ ��������� � ������� ������ 
// AliceFlow_v0_48
//#define doublereal double
#include <math.h>
#include <time.h>

// ������������ ����� �� ��������� � ������ ���������� 
// ����������� ����������� ����������� � ���������� ������
// �� ������ 110 ������ ���������� TGF2023-20: TEMPERATURE_FAILURE_DC.
void diagnostic_critical_temperature(doublereal TiP,  FLOW* &fglobal, TEMPER &t, 
	BLOCK* b, integer lb) 
{
	 if (TiP > TEMPERATURE_FAILURE_DC) {
		   printf("Attantion! unit burned...\n");
		   printf("temperature exceeds a critical TEMPERATURE_FAILURE_DC==%3.2f",TEMPERATURE_FAILURE_DC);



		   if (1) {
			   if (!b_on_adaptive_local_refinement_mesh) {

				   bool bextendedprint = false;
				   integer flow_interior = 0;
				   // ������� ���������� ���������� � ��������� tecplot360:
				   exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, bextendedprint, 0, b, lb);
			   }
			   else {
				   
				   // ������� � ��������� tecplot �����������.
				   //� ���� �����.
				   ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t, fglobal, 0, b, lb);
			   }
		   }


		   unsigned int calculation_main_end_time = clock();
		   unsigned int calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;
		   
		   // ����� ����� ���������� �� ������������� ������������ ���������.
		   int im = 0, is = 0, ims = 0;
		   im = (int)(calculation_main_seach_time / 60000); // ������
		   is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // �������
		   ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

		   printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);
		   //  getchar();
		   system("pause");
		   char ch;
		   printf("please enter c)ontinue h)alt.\n");
		   ch=getchar();
		   if ((ch=='h')||(ch=='H')) exit(0);
	   }
} // diagnostic_critical_temperature

// �������� 0 ��� ��������� ��������� �������������.
#define MY_AIR 1
#define MY_WATER 2

// � �������� � �������� ������ ������� ������� ������������ ���� ��� ������.
// �� � ������ ����� (����� � ����� �����) ���-���� ������.
// ������� ����� ��������� ����������� ������������ � ���������� �������� ���� � �������.

// ������ ������������ �� ��������� �� ���� ��������,
// ���� ������������ �� ������� �����.
void my_air_properties(doublereal TiP, doublereal PiP, 
	                   doublereal &rho, doublereal &cp,
					   doublereal &lam, doublereal &mu,
					   doublereal &beta) {
	// ���������� �������� �������� ������� !
	// ������� ���������: 
    // TiP - ����������� � �������� ������� � ������ ������������ ������.
	// PiP - �������� � ���������� (������������� �� ����, 0 - ������� ��������). 

	// ������ ����������� ������������� ����������� � ���������
	// ���������� 273-473K.
	// ������: HighExpert.Ru ���������� �������� �������
	// www.engineeringtoolbox.com.

    doublereal TK=TiP+273.15; // K
	//doublereal Pair=101325+PiP; // ��

	// �������� ����� � �������:
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
	// ��������� ����������� ��� ������ ���� < 0.3.
	// ��������� �������� 102 �/�.
	//rho=Pair/(287.4*TK); // kg/m^3
	// 27. 07. 2016.
	// � ����������� ���������� ������� �������� ��������� ���������,
	// ����������� ��������� �� ����������� ����������� ����� ����������� ��������� �������������� ����������.
	rho = 1.1614;
	cp=1000.0*(1.0005+1.1904e-4*(TiP)); // J/(kg*GC)
	lam=2.44e-2*exp(0.82*log(TK/273.15)); // W/(m*K)
	mu=1.717e-5*exp(0.683*log(TK/273.15)); // Pa*s
	//printf("mu=%e, TK=%e\n",mu,TK); // debug Ok
	//getchar();
	// � maple �� ������� ���������� ��������� �������� ��������� �����������:
    // ���������� 27.07.2016
	beta = 0.001*(-0.06810259 + 0.00013411*TiP - 8.287412214e-8*TiP*TiP + 1022.18 / (273.15 + TiP));

	// ���������� ���������� ��������
	//rho=1.1614; // ������������� !!!
	//mu=1.84e-5; 
	//beta=0.003331;
	//cp=1005;
	//lam=0.025;
	
} // my_air_properties

// ���� - ��� ������ ����� ������� ������ ������� ����,
// ��� ��������. ��� �������� ����� �����.
void my_water_properties(doublereal TiP,  
	                   doublereal &rho, doublereal &cp,
					   doublereal &lam, doublereal &mu,
					   doublereal &beta) {
	// ���������� �������� �������� ���� !
	// ������� ���������: 
    // TiP - ����������� � �������� ������� � ������ ������������ ������.

	// ������ ����������� ������������� ����������� � ���������
	// ���������� 283-373K.
	// ������: HighExpert.Ru ���������� �������� ����
	// www.engineeringtoolbox.com.

    doublereal TK=TiP+273.15;

	// �������� ����� � ����:
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
	// ��������� ����������� ��� ������ ���� < 0.3.

	//rho=995.7/(0.984+0.483e-3*TiP); // kg/m^3
	// 27. 07. 2016.
	// � ����������� ���������� ������� �������� ��������� ���������,
	// ����������� ��������� �� ����������� ����������� ����� ����������� ��������� �������������� ����������.
	rho = 995.7 / (0.984 + 0.483e-3*(20.0)); // kg/m!3
	cp=(4194-1.15*TiP+1.5e-2*TiP*TiP); // J/(kg*GC)
	lam=0.553*(1.0+0.003*TiP); // W/(m*K)
	mu=(1.78e-6/(1.0+0.0337*TiP+0.000221*TiP*TiP))*rho; // Pa*s
	// � maple �� ������� ���������� ��������� �������� ��������� �����������:
	beta=0.02289718770+0.142992772827544502e-6*TiP*TiP-0.67623235e-4*TiP-(6.27314791869)/(TK);

	// ���������� ���������� �������� ��� 20 ���� �.
	//rho=998.2; // ������������� !!!
	//mu=1e-3;
	//cp=4177;
	//lam=0.58618;
	//beta=0.192166e-3;
} // my_water_properties

// SOLID ������ ����:
// ������� �������: ��������� ������ ����������� ������.

// �������� 100 ��� ��������� ��������� �������������.
//#define MY_AUSN 101    // ������ ������ 80% ������ 20%
#define MY_ALUMINA 101   // ������� Alumina, Al2O3 -ref. MatWeb[2] 99.9% Al2O3, polycrystalline aluminum oxide.
#define MY_SI 102        // �������
#define MY_GAAS 103      // ������� �����
#define MY_GAN 104       // ������ �����
#define MY_SIC4H 105     // ������ ������� SiC4H
#define MY_SAPPHIRE 106  // ������
#define MY_DIAMOND 107   // �����
#define MY_MD40 108      // ����������� ��40
#define MY_AU 109        // ������
#define MY_SIO2 110      // ����� �������
#define MY_COPPER 111    // ����
#define MY_KOVAR 112     // �����
#define MY_BRASS 113     // ������ ��-59-1-�
#define MY_DURALUMIN 114 // ����������� �16
#define MY_ALUMINIUM_NITRIDE 115 // ������ ��������� (���������������).
#define MY_GLUE_ECHES 116 // ���� ����

// ������ ������ ������ (�����) 
// 21 ����� 2016 ���� ������� ������������� ����������� �� ������
// ��������� ������� ����������� ����������������� ����������.
void my_ausn_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // AuSn 80 ������/20 �����
       //doublereal TK=TiP+273.15;
	   rho=14510; // kg/m^3
	   // cp=150 �� ������ ��������� SYMMIC
	   cp=143; // J/(kg*K) �� ������ TDIM - Master
	  
	   // ���������������� � ���������� �������� ������ 
	   // ��� ���� �����: 68 ��� 300� � 57 ��� 85 ����. �.
	   lam=57; // �.�. ����� �� ����� ����������� �������� ���������� �������.
	   // ������� �� ������ ���������� ������� ����� ��-��������� � �������� ��� 
	   // ��������������� �.�. ��� �� ������������. ���� ������� � ������ ����������������
	   // ������� �������� �� ������ �������� � ��������� (���������������) ���� ���� �������� � 
	   // ������� ���������������� ����� � ������ �����������.
	   //lam=68*exp(-0.99678*log(TK/300)); // 57 W/(m*K) �������� ������, ����� ����

	   // ���������������� � ����������� ������-���������� ������ ������ 
	   // ��������� �� ������ ������� ��� ����������� � ����������� ���� 
	   // ���������� ������ ������, �.�. ��������� ���������� ������ �� �����������
	   // ������� ������� ��������� ���������� ����� � ������,  ������� ���������������� �
	   // ����������� ������� ���������� ���� ��� ������ ��  ��������� ������ � �����������.
	   // ������ �������� �������� ������������������ �������, ������������ ���������� ������� ����������.

} // AuSn

// ������� Si Silicon 
void my_si_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // Silicon
       doublereal TK=TiP+273.15;
	   rho=2330; // kg/m^3
	   //cp=711; // J/(kg*K)
	   // ������� ��� ����������� ����������� � ��������� ����������: 77K-773K.
	   // TK  77 173 273 373 573 773
	   // cp 180 490 680 770 850 880
	   if (TK < 800.0) {
		   cp=83.23167276-11801.68987/TK+0.3123808545e-5*TK*TK*TK+3.672005743*TK-0.5805489801e-2*TK*TK;
	   }
	   else {
		   // ������������ ������ �� �� ������� ��� ������������� �� ������� ���������� ��������� �� ��������.
		   cp=880;
	   }
	   // ������� ��������� ���������������� �� ����������� ����� ��
	   // ����� �. ���� ����������� ������� �����.
	   //lam=120; // ������ �.�. �������������
	   lam=148*exp(-1.35*log(TK/300)); // 120 W/(m*K)
} // Si

// 18 october 2016.
// ������� ����������� ��� �������� ��������� �� ����������� ��������.
// Alumina, Al2O3-ref. Mat Web[2] 99.9% Al2O3, polycrystalline aluminum oxide.
void my_alumina_properties(doublereal TiP,
	doublereal &rho, doublereal &cp,
	doublereal &lam) {

	doublereal TK = TiP + 273.15;
	rho = 2330; // kg/m^3

	//cp=753; // J/(kg*K)
	// ������� ��� ����������� ����������� � ��������� ����������: 200K-800K.
	// TK  200 298 400 500 600 800
	// cp 502 753 920 1046 1088 1172
	if (TK < 800.0) {
		cp = 0.00000505903*TK*TK*TK - 0.01042874335*TK*TK + 7.612830818*TK + (50982.822 / TK) - 898.2525;
	}
	else {
		// ������������ ������ �� �� ������� ��� ������������� �� ������� ���������� ��������� �� ��������.
		cp = 1172;
	}

	// ���������������� �� ������ ��������� SYMMIC.
	// TK lam
	// 200 82
	// 298 46
	// 400 32.3
	// 500 24.2
	// 600 18.9
	// 800 13.0
	lam = 47.84832382*exp(-1.328556144*log(TK / 300.0)); // 46

} // alumina

// ������� ����� GaAs 
void my_gaas_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // GaAs
       doublereal TK=TiP+273.15;
	   rho=5300; // kg/m^3
	   cp=322; // J/(kg*K)
	   if ((TK > 173.0) && (TK < 1514.0))
	   {
		   cp = 4.059e-8*TK*TK*TK - 1.296e-4*TK*TK + 1.854e-1*TK + 279.0;
	   }
	   else if (TK <= 173.0) cp = 307.4;
	   else cp = 402.9;
	   // ������� ��������� ���������������� �� ����������� ����� ��
	   // ����� �. ���� ����������� ������� �����.
	   //lam=47; // ������ �.�. �������������
	   lam=54*exp(-1.25*log(TK/300)); // 54 W/(m*K)

	   
} // GaAs

// ������ ����� GaN
void my_gan_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // GaN
       doublereal TK=TiP+273.15;
	   rho=6150; // kg/m^3
	   cp=491; // J/(kg*K)
	   // ������� ��� ����������� ����������� � ��������� ����������: 200K-1300K.
	   // TK  200 300 400 500 600 700 800 900 1000 1100 1200 1300
	   // cp 322.3 431.3 501.2 543.8 572.17 592.8 608.9 622.2 633.7 643.99 653.4 662.22
	   if ((TK > 200) && (TK < 1300)) {
		   cp = 591.807 + 0.064971*TK - 2.6e7 / (TK*TK) + 2.94e9 / (TK*TK*TK);
	   }
	   // ������� ��������� ���������������� �� ����������� ����� ��
	   // ����� �. ���� ����������� ������� �����.
	   //lam=130; // ������ �.�. �������������
	   lam=130*exp(-0.43*log(TK/300)); // 130 W/(m*K)
} // GaN

// ������ ������� SiC4H
void my_sic4h_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // SiC4H
       doublereal TK=TiP+273.15;
	   rho=3210; // kg/m^3
	   //cp=690; // J/(kg*K)
	   // ����������� ����������� ������� �� �����������
	   // ��������� �� ��������� ������� ������:
	   // TK 77 300 373 573 773
	   // Cp 50 690 820 1010 1120
	   if ((TK>=77.0)&&(TK<=773.2)) {
		   cp=-600.9595790625+6.41266663224237*TK+16052.9931969468/TK-0.00900588312450910*TK*TK+0.459922567674634e-5*TK*TK*TK;
	   } else if (TK<77.0) cp=50;
	   else cp=1120;

	   // ������� ��������� ���������������� �� ����������� 
	   // � ����� �. ���� ����������� ������� �����. ���!!!
	   //lam=370; // V-�����. SiC
	   lam=370*exp(-1.411*log(TK/300)); // 370 W/(m*K)
} // SiC4H

// ������ Sapphire
void my_sapphire_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // Sapphire Al2O3
       doublereal TK=TiP+273.15;
	   rho=3980; // kg/m^3
	   //cp=750; // J/(kg*K)
	   // ������� ��� ����������� ����������� � ��������� ���������� 77K - 773K
	   // TK 77 173 273 373 573 773
	   // Cp 60 403 718 907 1089 1168
	   if (TK<774) {
		   cp=-703.3920492+19101.78790/TK+0.5623375490e-5*TK*TK*TK+7.499240665*TK-0.01095728614*TK*TK;
	   }
	   else cp=1168;

	   // ������� ��������� ���������������� �� ����������� 
	   // �  ����� �. ���� ����������� ������� �����. ���!!!
	  // lam=23; // ������ �.�. �������������
	   // �� ������ ��������� SYMMIC.
	   lam=42*exp((-1.134821)*log(TK/298)); // 42 W/(m*K) ��� 298�

} // Sapphire

// ����� Diamond
void my_diamond_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // Diamond
       doublereal TK=TiP+273.15;
	   rho=3500; // kg/m^3
	   //cp=520; // ��� 300� J/(kg*K)
	   // ������ �� ����������� ����������� � ��������� ����������:
	   // 77K - 773K
	   // TK 77  173 273 373 573 773
	   // Cp 8  140  420 770 1300 1590 � ���������� �����, �� ������ ���.
	   if (TK < 774) {
		   cp=-696.4289+32660.38060/TK-0.2774926361e-5*TK*TK*TK+3.566362474*TK+0.1285766705e-2*TK*TK;
	   }
	   else cp=1590;
	   // cp=1260; // �� ������ �.�. �������������
	   // ������� ��������� ���������������� �� ����������� ����� ��
	   // ����� �. ���� ����������� ������� �����.
	   //lam=800; // ������ �.�. �������������
	   lam=2300*exp(-1.85*log(TK/300)); // 130 W/(m*K)
	   // �� ������ �. ���� ���������������� ������ ��� 300� ����� ����� 2000 � 2500.
} // Diamond

// ����������� CuMo � 40% ����������� ����: ��40
void my_md40_properties(doublereal TiP,
	                    doublereal &rho, doublereal &cp,
						doublereal &lam) {
       // �� ��������� ����� ������������ ������ ������ �������:
	   // ���������, ������� ������ � ������. ������ �������������.
       // MD40
       doublereal TK=TiP+273.15;
	   rho=9660; // kg/m^3
	   // cp=540 �� ������ ��������� ��������, �����������.
	   // ������ � ��������� �������� � ����������� ����� ������� ����������� ���������.
	   cp=318.2; // J/(kg*K)
	   // heat capacity cp
	   // TK 173 273 373 573 
	   // cp J/(kg-K) 268.8 302.5 318.17 335.1
	   if ((TK > 173) && (TK < 573))
	   {
		   cp = 1.66e-6*TK*TK*TK - 2.263e-3*TK*TK + 1.095*TK + 138.5;
	   }
	   else if (TK <= 173) cp = 268.77;
	   else cp = 335.111;
	   lam=210; // W/(m*K)
	   if ((TK >= 273) && (TK < 1600)) {
		   lam = 7.15e-9*TK*TK*TK - 2.18e-6*TK*TK - 6.03e-2*TK + 271.0;
	   }
	   else if (TK < 273) lam = 254.52;
	   else lam = 198.0;
} // MD40

// ������ Au
void my_au_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // Au
       doublereal TK=TiP+273.15;
	   rho=19300; // kg/m^3
	   //cp=126; // J/(kg*K)

	   // �������������� ����������� �������� �� ��������� ������
	   // � ������� ���:
	   // TK 77 173 273 373 573 773
	   // Cp 97 121 128 131 135 140
	   if ((TK>77.0)&&(TK<773.0)) {
		   cp=140.2979168-3345.749125/TK+0.3399388027e-7*TK*TK*TK+0.3605589168e-2*TK-0.2418961279e-4*TK*TK;
	   }
	   else if (TK<=77.0) cp=97;
	   else cp=140;

	   //lam=293; // W/(m*K)
	   // ��� �����������, ������������ � ����� �. ���� �� ������� ���
	   // ���������������� ������.
	   // ������� ���� �������� ���� ������������ ������� � ������� ���.
	   if (TK<974) {
		   lam=327.5268055+146.5916497/TK+0.5614880962e-7*TK*TK*TK-0.008865308782*TK-0.1043207287e-3*TK*TK;
	   }
	   else lam=272;

} // Au

// SiO2
// ���������� ���������� �� ������� ����.
// ������������� ������ ������� ������������� ���������� SYMMIC.
void my_sio2_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // SiO2
       doublereal TK=TiP+273.15;
	   rho=2100; // kg/m^3
	   cp=741; // J/(kg*K)
	   if ((TK >= 273.0) && (TK <= 773))
	   {
		   cp = -1.6667e-7*TK*TK*TK - 9.635e-4*TK*TK + 1.975*TK + 236.02;
	   }
	   else if (TK < 273.0) cp = 700.0;
	   else cp = 1110.0;

	   //lam=7; // 7-11.5 W/(m*K) quartz crystal
	   lam = 1.27; // ���������� ���������� �� ������� ����.

} // SiO2

// ���� Copper
void my_copper_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // Copper
       doublereal TK=TiP+273.15;
	   rho=8930; // kg/m^3
	   //cp=390; // J/(kg*K)
	   // �������������� ����������� ����������� �� �����������
	   // �������� �� ��������� ������ � ������� ���.
	   // TK    77  173 273 373 573 773
	   // Cp    195 341 379 397 419 430 
	   if ((TK>77.0)&&(TK<773.0)) {
		   cp=506.2148755-22501.73085/TK-0.3085157062e-6*TK*TK*TK-0.2854473337*TK+0.5289207596e-3*TK*TK;
	   }
	   else if (TK<77.0) cp=195;
	   else cp=430;

	   //lam=390; // W/(m*K)
	   // �������������� ����������� ���������������� �� �����������
	   // �������� �� ��������� ������ � ������� ���.
	   // TK      173.2 273.2 373.2 573.2 973.2
	   // Lambda  420.0 403.0 395.0 381   354
	   if ((TK>173.2) && (TK<973.2)) {
		   lam=318.4875333+12313.69842/TK+0.1748333335e-6*TK*TK*TK+0.2380247005*TK-0.3905912003e-3*TK*TK;
	   }
	   else if (TK<=173.2) lam=420.0;
	   else lam=354.0;
} // Copper

// ����� Kovar
void my_kovar_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // Kovar
	   // ����� ��� ����� ������ �������: C, Fe 54%, Mn, Ni 29%
       doublereal TK=TiP+273.15;
	   rho=8300; // kg/m^3
	   cp=669; // J/(kg*K)
	   // heat capacity kovar
	   // TK 273 703
	   // cp 439.614 648.954
	   if ((TK>273) && (TK < 703))
	   {
		   cp = 0.486*TK + 306.6;
	   }
	   else if (TK < 273) cp = 439.6;
	   else cp = 648.954;
	   //lam=19; // W/(m*K) 
	   // �������������� ����������� ���������������� �� �����������
	   // �������� �� ��������� ������ � ������� ���.
	   // TK      273  373  573  773  973
	   // Lambda  14.1 14.7 15.6 17.5 19.3
	   if ((TK>273.0) && (TK<973.0)) {
		   lam=44.13357753-3650.246187/TK-0.4572641192e-7*TK*TK*TK-0.8854732756e-1*TK+0.1132267331e-3*TK*TK;
	   } else if (TK<273.0) lam=14.1;
	   else lam=19.3;
} // Kovar

// ������ ��-59-1-� Brass
void my_brass_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // Brass
	   // ������ Cu 70%, Zn 30%
       doublereal TK=TiP+273.15;
	   rho=8440; // kg/m^3
	   cp=377; // J/(kg*K)
	   //lam=108.784; //  W/(m*K)
	   // ������ � ������������� ����������� ���������������� ������:
	   // T GC   -100 0 100 200 300 400
	   // TK      173.2 273.2 373.2 473.2 573.2 673.2
	   // Lambda  90 106 131 143 145 148
	   if ((TK>=173.2)&&(TK<=673.2)) {
		   lam=-591.598556543738+2.90903837477216*TK+52133.0458047673/TK-0.00454066869215317*TK*TK+0.249629141483149e-5*TK*TK*TK;
	   } else if (TK<173.2) lam=90.0;
	   else lam=148.0;

} // Brass

// ����������� Duralumin �16
void my_duralumin_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // Duralumin
       doublereal TK=TiP+273.15;
	   rho=2780; // kg/m^3
	   cp=922; // J/(kg*K)
	   // heat capacity duralumin
	   // TK  273 373 473 573 673
	   // cp J/(kg-K) 750.13 921 1047 1130 1172
	   if ((TK>273.0) && (TK < 673.0)) {
		   cp = 3.33e-7*TK*TK*TK - 2.62e-3*TK*TK + 3.3*TK + 38.0;
	   }
	   else if (TK <= 273.0) cp = 750.13;
	   else cp = 1172.0;
	   lam=130; // W/(m*K)
	   if ((TK>150.0) && (TK < 573)) {
		   lam = 0.171*TK + 66.1;
	   }
	   else if (TK <= 150.0) lam = 90.0;
	   else lam = 164.0;
} // Duralumin

// ���� ����
void my_glueeches_properties(doublereal TiP,
	                  doublereal &rho, doublereal &cp,
					  doublereal &lam) {
       // glue ECHES
	   // �� ������ �.�. �������������.
       //doublereal TK=TiP+273.15;
	   rho=1000; // kg/m^3
	   cp=100; // J/(kg*K)
	   lam=4; // W/(m*K)
} // glue ECHES

// ������ �������� (���������������)
void my_aluminium_nitride_properties(doublereal TiP,
	doublereal &rho, doublereal &cp,
	doublereal &lam) {
	// aluminium nitride
	// �� ������ ��������� SYMMIC.
	doublereal TK = TiP + 273.15;
	rho = 3255; // kg/m^3
	cp = 600; // J/(kg*K) // SYMMIC
	//cp=748; // J/(kg*K) // �. ����.
	//lam = 319; // W/(m*K)
	// lam=285; // W/(m*K) �. ����

	// TSIDE,K lam W/(mxK)
	// 200 780
	// 300 319
	// 400 195
	// 600 100
	// 1000 49

	// ������� ��������� ���������������� �� ����������� 
	// �  ����� �. ���� ����������� ������� �����. ������������, ��� ��� ���������������� ��� � ��� ��������.
	// lam=319; // ������ ��������� SYMMIC
	// �� ������ ��������� SYMMIC.
	lam = 319 * exp((-1.57)*log(TK / 300)); // 319 W/(m*K) ��� 300�
} // aluminium nitride


// ����������������� ���������� ������ ����������.
void my_fluid_properties(doublereal TiP, doublereal PiP, 
	                   doublereal &rho, doublereal &cp,
					   doublereal &lam, doublereal &mu,
					   doublereal &beta, integer ilibident) {
	switch (ilibident) {
	  case MY_AIR: my_air_properties(TiP, PiP, rho, cp, lam, mu, beta); break;
	  case MY_WATER: my_water_properties(TiP, rho, cp,lam, mu, beta); break;
	  default: my_air_properties(TiP, PiP, rho, cp, lam, mu, beta);  break;
	} // end switch
} // my_fluid_properties

doublereal dsic=1.0e30, dgan=1.0e30, dcu=1.0e30;

 // ����������������� ���������� ������ ����������.
void my_solid_properties(doublereal TiP, doublereal &rho, doublereal &cp,
	                     doublereal &lam, integer ilibident) {
	switch (ilibident) {
	  //case MY_AUSN: my_ausn_properties(TiP,rho,cp,lam); break; // ������ ������ ������
	  case MY_ALUMINA:  my_alumina_properties(TiP, rho, cp, lam); break; // ������� Al2O3.
	  case MY_SI: my_si_properties(TiP,rho,cp,lam); break; // �������
	  case MY_GAAS: my_gaas_properties(TiP,rho,cp,lam); break; // ������� �����
	  case MY_GAN: my_gan_properties(TiP, rho, cp, lam); if (lam < dgan) { dgan = lam; }  break; // ������ �����
	  case MY_SIC4H: my_sic4h_properties(TiP, rho, cp, lam); if (lam < dsic) { dsic = lam; } break; // ������ �������
	  case MY_SAPPHIRE: my_sapphire_properties(TiP,rho,cp,lam); break; // ������
	  case MY_DIAMOND: my_diamond_properties(TiP,rho,cp,lam); break; // �����
	  case MY_MD40: my_md40_properties(TiP,rho,cp,lam); break; // ����������� ��40
	  case MY_AU: my_au_properties(TiP,rho,cp,lam); break; // ������
	  case MY_SIO2: my_sio2_properties(TiP,rho,cp,lam); break; // SiO2
	  case MY_COPPER: my_copper_properties(TiP, rho, cp, lam); if (lam < dcu) { dcu = lam; }  break; // ����
	  case MY_KOVAR: my_kovar_properties(TiP,rho,cp,lam); break; // �����
	  case MY_BRASS: my_brass_properties(TiP,rho,cp,lam); break; // ������ ��-59-1-� Brass
	  case MY_DURALUMIN: my_duralumin_properties(TiP,rho,cp,lam); break; // �����������
	  case MY_ALUMINIUM_NITRIDE: my_aluminium_nitride_properties(TiP, rho, cp, lam); break; // aluminium nitride
	  case MY_GLUE_ECHES: my_glueeches_properties(TiP, rho, cp, lam); break; // ���� ����
	  default: my_duralumin_properties(TiP,rho,cp,lam); break; // ����������� - ������������ �������� �� ���������
	} // end switch
} // my_solid_properties

// ������ �������� �������� ��������� � ��������� ��.
void gran_prop(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, integer iP, integer G, integer ib, TPROP* matlist) {

	// t.ptr[1][iP] - ������������� Fluid domain ��� -1 ���� ��� ������ Solid.

	doublereal rho, cp, lam;
	integer iG; // ����� ������ ������� ����� ��������� ��������� �����.
	iG=t.neighbors_for_the_internal_node[G][0][iP];
	if (iG>=t.maxelm) {
		// ��� ��������� ����
        rho=1.1614; cp=1005; lam=0.025; // ������������� default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat==1) {
			if (b[ib].itype== PHYSICS_TYPE_IN_BODY::SOLID) {
		        my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // ������������� ����������� � ��������� ����
				// �������� �� ������������ ����������.
				diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
		    } // SOLID
		    else if (b[ib].itype== PHYSICS_TYPE_IN_BODY::FLUID) {
		       doublereal mu, beta_t; // �������� �� ������������ �� ���������.
		       doublereal pressure;
		       if (t.ptr[1][iP]==-1) {
			       pressure=0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
		       }
			   else {
				   if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][0][t.ptr[0][iP]] > -1) {
					   pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][0][t.ptr[0][iP]]];
				   }
				   else {
					   // STUB 15.09.2018
					   pressure = 0.0;  // � ������ ����.

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
            // �������� ����������� �������������:
			// ���������� ��������.
			rho=matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

			// ��������� ���������� �������� ��������� �� ���������� ����� �������
			// � ��������� ��������� �����.
		}
		// �������� ��� ����������� ������������ ������.
		t.prop_b[RHO][iG-t.maxelm]=rho;
		t.prop_b[HEAT_CAPACITY][iG-t.maxelm]=cp;
		t.prop_b[LAM][iG-t.maxelm]=lam;
	} // G Side

	if (t.neighbors_for_the_internal_node[G][1] != nullptr) {
		iG = t.neighbors_for_the_internal_node[G][1][iP];
		if (iG >= t.maxelm) {
			// ��� ��������� ����
			rho = 1.1614; cp = 1005; lam = 0.025; // ������������� default  dry air 300K 1atm properties
			if (matlist[b[ib].imatid].blibmat == 1) {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // ������������� ����������� � ��������� ����
					// �������� �� ������������ ����������.
					diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
				} // SOLID
				else if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					doublereal mu, beta_t; // �������� �� ������������ �� ���������.
					doublereal pressure;
					if (t.ptr[1][iP] == -1) {
						pressure = 0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
					}
					else {
						if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][1][t.ptr[0][iP]] > -1) {
							pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][1][t.ptr[0][iP]]];
						}
						else {
							// STUB 15.09.2018
							pressure = 0.0; // � ������ ����.

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
				// �������� ����������� �������������:
				// ���������� ��������.
				rho = matlist[b[ib].imatid].rho;
				//cp = matlist[b[ib].imatid].cp;
				//lam = matlist[b[ib].imatid].lam;
				cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
				lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

				// ��������� ���������� �������� ��������� �� ���������� ����� �������
				// � ��������� ��������� �����.
			}
			// �������� ��� ����������� ������������ ������.
			t.prop_b[RHO][iG - t.maxelm] = rho;
			t.prop_b[HEAT_CAPACITY][iG - t.maxelm] = cp;
			t.prop_b[LAM][iG - t.maxelm] = lam;
		} // G Side
	}


	if (t.neighbors_for_the_internal_node[G][2] != nullptr) {
		iG = t.neighbors_for_the_internal_node[G][2][iP];
		if (iG >= t.maxelm) {
			// ��� ��������� ����
			rho = 1.1614; cp = 1005; lam = 0.025; // ������������� default  dry air 300K 1atm properties
			if (matlist[b[ib].imatid].blibmat == 1) {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // ������������� ����������� � ��������� ����
					// �������� �� ������������ ����������.
					diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
				} // SOLID
				else if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					doublereal mu, beta_t; // �������� �� ������������ �� ���������.
					doublereal pressure;
					if (t.ptr[1][iP] == -1) {
						pressure = 0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
					}
					else {
						if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][2][t.ptr[0][iP]] > -1) {
							pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][2][t.ptr[0][iP]]];
						}
						else {
							// STUB 15.09.2018
							pressure = 0.0; // � ������ ����.

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
				// �������� ����������� �������������:
				// ���������� ��������.
				rho = matlist[b[ib].imatid].rho;
				//cp = matlist[b[ib].imatid].cp;
				//lam = matlist[b[ib].imatid].lam;
				cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
				lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

				// ��������� ���������� �������� ��������� �� ���������� ����� �������
				// � ��������� ��������� �����.
			}
			// �������� ��� ����������� ������������ ������.
			t.prop_b[RHO][iG - t.maxelm] = rho;
			t.prop_b[HEAT_CAPACITY][iG - t.maxelm] = cp;
			t.prop_b[LAM][iG - t.maxelm] = lam;
		} // G Side
	}

	if (t.neighbors_for_the_internal_node[G][3] != nullptr) {
		iG = t.neighbors_for_the_internal_node[G][3][iP];
		if (iG >= t.maxelm) {
			// ��� ��������� ����
			rho = 1.1614; cp = 1005; lam = 0.025; // ������������� default  dry air 300K 1atm properties
			if (matlist[b[ib].imatid].blibmat == 1) {
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					my_solid_properties(t.potent[iG], rho, cp, lam, matlist[b[ib].imatid].ilibident); // ������������� ����������� � ��������� ����
					// �������� �� ������������ ����������.
					diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
				} // SOLID
				else if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					doublereal mu, beta_t; // �������� �� ������������ �� ���������.
					doublereal pressure;
					if (t.ptr[1][iP] == -1) {
						pressure = 0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
					}
					else {
						if (f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][3][t.ptr[0][iP]] > -1) {
							pressure = f[t.ptr[1][iP]].potent[PRESS][f[t.ptr[1][iP]].neighbors_for_the_internal_node[G][3][t.ptr[0][iP]]];
						}
						else {
							// STUB 15.09.2018
							pressure = 0.0; // � ������ ����.

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
				// �������� ����������� �������������:
				// ���������� ��������.
				rho = matlist[b[ib].imatid].rho;
				//cp = matlist[b[ib].imatid].cp;
				//lam = matlist[b[ib].imatid].lam;
				cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iG]);
				lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iG]);

				// ��������� ���������� �������� ��������� �� ���������� ����� �������
				// � ��������� ��������� �����.
			}
			// �������� ��� ����������� ������������ ������.
			t.prop_b[RHO][iG - t.maxelm] = rho;
			t.prop_b[HEAT_CAPACITY][iG - t.maxelm] = cp;
			t.prop_b[LAM][iG - t.maxelm] = lam;
		} // G Side
	}
} // gran_prop

bool bswitch_print_message = true;

// ���������� ������� ����������
// ��� ��������� ����������������
void update_temp_properties(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, TPROP* matlist) {
	// �������� ����������, ����� ��� ���������, ����������� � ����������������
	// ����������� ������� �� �����������. � ������ ������� ������������ ����������
	// ������� ���������� � ������ �� ����������� �� �������� �������� �����������.
	// ���� ���������������� ��������� ������ � ������ �����������, �� �� ����� �������
	// � ������������� �������� ������ � ������� ��������� ������������ ����������� 
	// ����� ���� ������� ���� ��� ��������������� ��������� ����������� � �������� 
	// ������� � ����������� ����������.

	//TOCHKA p; // ����� - ����� ���������������� ��.	
	
	const integer ISIZE = t.maxelm;

#pragma omp parallel for
	for (integer iP=0; iP<ISIZE; iP++) {
		// ������ �� ���� ���������� ����������� ������� ��������� �������.

		doublereal dmin = 1.0e30;
		doublereal dmax = -1.0e30;
		dgan = 1.0e30;
		dsic = 1.0e30;
		dcu = 1.0e30;

		//center_cord3D(iP, t.nvtx, t.pa, p); // ���������� ��������� ������ ��.
		//in_model_temp(p,ib,b,lb); // ���������� ����� ����� ib �������� ����������� ����������� ����� � ������� iP.
		// 8 ������ 2016 ��������� ���������� ib
		

		integer ib; // ����� ����� �������� ����������� ����������� �����.
		doublereal rho, cp, lam;

		ib = t.whot_is_block[iP];
		rho=1.1614; cp=1005; lam=0.025; // ������������� default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat==1) {
			// ������������, ����������� ������ ��������� AliceFlow ��������.
			if (b[ib].itype== PHYSICS_TYPE_IN_BODY::SOLID) {
			    my_solid_properties(t.potent[iP], rho, cp, lam, matlist[b[ib].imatid].ilibident);
				// �������� �� ������������ ����������.
				diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
		    } // SOLID
		    if (b[ib].itype== PHYSICS_TYPE_IN_BODY::FLUID) {
			   doublereal mu, beta_t; // �������� �� ������������ �� ���������.
		       doublereal pressure;
			   if ((t.ptr==nullptr)||(t.ptr[1][iP]==-1)) {
				   pressure=0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
			   }
			   else pressure=f[t.ptr[1][iP]].potent[PRESS][t.ptr[0][iP]];
			   my_fluid_properties(t.potent[iP], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
		    } // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat==0) {
			// �������� ����������� �������������:
			// ���������� ��������.
			rho=matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iP]);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iP]);

			// ������������ �������� ���������, ��������� �� �����������.
			t.prop[BETA_T_MECHANICAL][iP]= get_beta_t_solid(matlist[b[ib].imatid].n_beta_t_solid, matlist[b[ib].imatid].temp_beta_t_solid, matlist[b[ib].imatid].arr_beta_t_solid, t.potent[iP]);
			t.prop[YOUNG_MODULE][iP] = get_Young_Module(matlist[b[ib].imatid].n_YoungModule, matlist[b[ib].imatid].temp_Young_Module, matlist[b[ib].imatid].arr_Young_Module, t.potent[iP]);
		    t.prop[POISSON_RATIO][iP]= get_Poisson_ratio(matlist[b[ib].imatid].n_Poisson_ratio, matlist[b[ib].imatid].temp_Poisson_ratio, matlist[b[ib].imatid].arr_Poisson_ratio, t.potent[iP]);
		}
		// �������� ��� ����������� ������������ ������.
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

		// ������ ��������� ���������� ��������� ����������� ������,
		// ������� ����������� � ������ ���������� ��.
		//
		// 25 �������� 2016 gran_prop ������ �������� �� ���� �����. 
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

  // ���������� ������� ����������
  // ��� ��������� ����������������
void update_temp_properties1(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, 
	TPROP* matlist, doublereal* &temperature, integer iadd, integer iu,
	doublereal* &lam_export, int** &nvtx_global) {
	// �������� ����������, ����� ��� ���������, ����������� � ����������������
	// ����������� ������� �� �����������. � ������ ������� ������������ ����������
	// ������� ���������� � ������ �� ����������� �� �������� �������� �����������.
	// ���� ���������������� ��������� ������ � ������ �����������, �� �� ����� �������
	// � ������������� �������� ������ � ������� ��������� ������������ ����������� 
	// ����� ���� ������� ���� ��� ��������������� ��������� ����������� � �������� 
	// ������� � ����������� ����������.
	 
	//TOCHKA p; // ����� - ����� ���������������� ��.
	
	doublereal dmin = 1.0e30;
	doublereal dmax = -1.0e30;
	for (integer iP = iadd; iP<iadd+t.maxelm; iP++) {
		// ������ �� ���� ���������� ����������� ������� ��������� �������.

		
		doublereal Temperature_in_cell = 0.125*(temperature[nvtx_global[0][iP]-1]+ temperature[nvtx_global[1][iP] - 1] + temperature[nvtx_global[2][iP] - 1] + temperature[nvtx_global[3][iP] - 1] + temperature[nvtx_global[4][iP] - 1] + temperature[nvtx_global[5][iP] - 1] + temperature[nvtx_global[6][iP] - 1] + temperature[nvtx_global[7][iP] - 1]);

		integer ib; // ����� ����� �������� ����������� ����������� �����.
		doublereal rho, cp, lam;

		if (iu == -1) {
			//center_cord3D(iP, t.nvtx, t.pa, p); // ���������� ��������� ������ ��.
			//in_model_temp(p,ib,b,lb); // ���������� ����� ����� ib �������� ����������� ����������� ����� � ������� iP.
			// 8 ������ 2016 ��������� ���������� ib

			ib = t.whot_is_block[iP];

		}
		else {
			ib = t.whot_is_block[iP - iadd];
		}
		rho = 1.1614; cp = 1005; lam = 0.025; // ������������� default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat == 1) {
			// ������������, ����������� ������ ��������� AliceFlow ��������.
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				my_solid_properties(Temperature_in_cell, rho, cp, lam, matlist[b[ib].imatid].ilibident);
				// �������� �� ������������ ����������.
				diagnostic_critical_temperature(t.potent[iP], f, t, b, lb);
			} // SOLID
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				doublereal mu, beta_t; // �������� �� ������������ �� ���������.
				doublereal pressure;
				if ((t.ptr==nullptr)||(t.ptr[1][iP-iadd] == -1)) {
					pressure = 0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
				}
				else pressure = f[t.ptr[1][iP-iadd]].potent[PRESS][t.ptr[0][iP-iadd]];
				my_fluid_properties(Temperature_in_cell, pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
			} // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat == 0) {
			// �������� ����������� �������������:
			// ���������� ��������.
			rho = matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, Temperature_in_cell);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, Temperature_in_cell);

			// ��������� ������������ �������� ��������� ��������� �� �����������.
			t.prop[POISSON_RATIO][iP - iadd] = get_Poisson_ratio(matlist[b[ib].imatid].n_Poisson_ratio, matlist[b[ib].imatid].temp_Poisson_ratio, matlist[b[ib].imatid].arr_Poisson_ratio, Temperature_in_cell);
			t.prop[YOUNG_MODULE][iP - iadd] = get_Young_Module(matlist[b[ib].imatid].n_YoungModule, matlist[b[ib].imatid].temp_Young_Module, matlist[b[ib].imatid].arr_Young_Module, Temperature_in_cell);
			t.prop[BETA_T_MECHANICAL][iP - iadd] = get_beta_t_solid(matlist[b[ib].imatid].n_beta_t_solid, matlist[b[ib].imatid].temp_beta_t_solid, matlist[b[ib].imatid].arr_beta_t_solid, Temperature_in_cell);

		}
		// �������� ��� ����������� ������������ ������.
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

		// ������ ��������� ���������� ��������� ����������� ������,
		// ������� ����������� � ������ ���������� ��.
		//
		// 25 �������� 2016 gran_prop ������ �������� �� ���� �����. 
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

// ���������� ����������� ������������ �������� �������� �������������
// 26 ������ 2012 ���������� �� ������������ ���������.
doublereal return_dynamic_viscosity(TPROP* matlist, integer imatid, 
	                          bool bfirst_start, doublereal SInvariantStrainRateTensor) {
	doublereal mu;

	// ���������� ������������ ��������:
	switch (matlist[imatid].ilawmu) {
	case 0: // ���������� ��������
		mu=matlist[imatid].mu;
		break;
	case 1: // power-law fluid
		// ����� ���������-�� ���� (Ostwald de Vel)
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
		   mu=fmax(matlist[imatid].mumin,fmin(2.0*matlist[imatid].Amu*exp((matlist[imatid].degreennmu-1.0)*log(SInvariantStrainRateTensor)),matlist[imatid].mumax));
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
	case 2: // ����� ������� Caisson
		// �������� ������ ������� �������� ������� � ������������ ���������� ������.
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
			if (matlist[imatid].Amu>=0.0) {
				mu=fmax(matlist[imatid].mumin,fmin(((sqrt(matlist[imatid].Amu/SInvariantStrainRateTensor))+matlist[imatid].Bmu)*((sqrt(matlist[imatid].Amu/SInvariantStrainRateTensor))+matlist[imatid].Bmu),matlist[imatid].mumax));
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
	case 3: // ����� ��������
		// �������� ������ �������� ��� B>0 �������� ���������� � ����������� ���������� ������,
		// � ��� B<0 ������������ �������� ������� � ������ ���������� ������.
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
			if ( (fabs(matlist[imatid].Amu)>1e-30) && (fabs(matlist[imatid].Bmu)>1e-30) && ((SInvariantStrainRateTensor/matlist[imatid].Bmu>-1.0)&&(SInvariantStrainRateTensor/matlist[imatid].Bmu<1.0))) {
				mu=fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,(matlist[imatid].Amu/SInvariantStrainRateTensor)*asin(SInvariantStrainRateTensor/matlist[imatid].Bmu)));
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
		mu=fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,matlist[imatid].Amu+(matlist[imatid].Bmu-matlist[imatid].Amu)*exp(0.5*(matlist[imatid].degreennmu-1.0)*log(1.0+(matlist[imatid].Cmu*SInvariantStrainRateTensor)*(matlist[imatid].Cmu*SInvariantStrainRateTensor)))));
		break;
	case 5: // ������-������
		// ���� B*C>0.0 �� ������������ �������� ������� � ������ ���������� ������,
		// ���� B*C<0.0 �� ������������ �������� ���������� � ������ ���������� ������.
		if ((!bfirst_start)&&(SInvariantStrainRateTensor>0.0)) {
			// ����-����� Arsh(x)=log(x+sqrt(x*x+1.0));
			mu=fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,matlist[imatid].Amu+(matlist[imatid].Bmu/SInvariantStrainRateTensor)*log(matlist[imatid].Cmu*SInvariantStrainRateTensor+sqrt(1.0+(matlist[imatid].Cmu*SInvariantStrainRateTensor)*(matlist[imatid].Cmu*SInvariantStrainRateTensor)))));
		} else 
		{
			if (matlist[imatid].Bmu*matlist[imatid].Cmu>0.0) {
				mu=matlist[imatid].mumax;
			} else mu=matlist[imatid].mumin;
		}
		break;
	case 6: // Williamson non-Newtonian fluid
		// ����� ������� �������������� ��� ��������� ������� ������ ��������� (��������������� ������� � ��������� ������� �� ������������).
		mu=fmax(matlist[imatid].mumin,fmin(matlist[imatid].mumax,matlist[imatid].Amu/(matlist[imatid].Bmu+SInvariantStrainRateTensor)+matlist[imatid].Cmu));
		break;
	default: // ���������� ��������
		mu=matlist[imatid].mu;
		break;
	}

	return mu;
} // return_dynamic_viscosity


// ��������� �������� ���������� � ������ �������������
// ������ �������� ���������� � ��������� ��.
void gran_prop_flow(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, integer iflow,
	                integer iP, integer G, integer ib, TPROP* matlist, bool bfirst_start) {

       doublereal rho=0.0, mu=0.0, beta_t=0.0;
	   integer iG=0; // ����� ������ ������� ����� ��������� ��������� �����.
	   iG=f[iflow].neighbors_for_the_internal_node[G][0][iP];
	   if (iG>=f[iflow].maxelm) {
		   // ��� ��������� ����
           rho=1.1614; mu=1.84e-5; beta_t=0.003331; // ������������� default dry air 300K 1atm properties
		   if (matlist[b[ib].imatid].blibmat==1) {
				doublereal cp, lam;
				cp=1005; lam=0.025;
				// ������������ ����������� ������ ��������� AliceFlow ��������
				doublereal temperature=20.0;
				if (t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]]>=t.maxelm) {
					// ��������� ��� ����������� �����
                    temperature=t.potent[t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]]];
				} else {
					if (t.neighbors_for_the_internal_node[G][0][f[iflow].ptr[iP]] > -1) {
						// ����������� �� ����� ����������������� �������� �������������:
						TOCHKA p1, p2;
						center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1,100); // ���������� ��������� ������ ��.
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
                // ������������ ����������� ������ ��������� AliceFlow ��������
				if (b[ib].itype== PHYSICS_TYPE_IN_BODY::FLUID) {
					my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				}
		   }
		   else if (matlist[b[ib].imatid].blibmat==0) {
				// �������� ����������� �������������:
				// ���������� ��������.

			    // ���������
				rho=matlist[b[ib].imatid].rho;
				// ������������ ��������
                mu=return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				// ����������� ��������� �������������� ����������.
				beta_t=matlist[b[ib].imatid].beta_t;
				// ��������� ���������� �������� ��������� �� ���������� ����� �������
				// � ��������� ��������� �����.
			}
		    // �������� ���������� ������������ ������:
		    f[iflow].prop_b[RHO][iG-f[iflow].maxelm]=rho;
		    f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG-f[iflow].maxelm]=mu;
		    f[iflow].prop_b[BETA_T][iG-f[iflow].maxelm]=beta_t;

	   } // G Side

	   if (f[iflow].neighbors_for_the_internal_node[G][1] != nullptr) {
		   iG = f[iflow].neighbors_for_the_internal_node[G][1][iP];
		   if (iG >= f[iflow].maxelm) {
			   // ��� ��������� ����
			   rho = 1.1614; mu = 1.84e-5; beta_t = 0.003331; // ������������� default dry air 300K 1atm properties
			   if (matlist[b[ib].imatid].blibmat == 1) {
				   doublereal cp, lam;
				   cp = 1005; lam = 0.025;
				   // ������������ ����������� ������ ��������� AliceFlow ��������
				   doublereal temperature = 20.0;
				   if (t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]] >= t.maxelm) {
					   // ��������� ��� ����������� �����
					   temperature = t.potent[t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]]];
				   }
				   else {
					   if (t.neighbors_for_the_internal_node[G][1][f[iflow].ptr[iP]] > -1) {
						   // ����������� �� ����� ����������������� �������� �������������:
						   TOCHKA p1, p2;
						   center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1, 100); // ���������� ��������� ������ ��.
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
				   // ������������ ����������� ������ ��������� AliceFlow ��������
				   if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					   my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				   }
			   }
			   else if (matlist[b[ib].imatid].blibmat == 0) {
				   // �������� ����������� �������������:
				   // ���������� ��������.

				   // ���������
				   rho = matlist[b[ib].imatid].rho;
				   // ������������ ��������
				   mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				   // ����������� ��������� �������������� ����������.
				   beta_t = matlist[b[ib].imatid].beta_t;
				   // ��������� ���������� �������� ��������� �� ���������� ����� �������
				   // � ��������� ��������� �����.
			   }
			   // �������� ���������� ������������ ������:
			   f[iflow].prop_b[RHO][iG - f[iflow].maxelm] = rho;
			   f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG - f[iflow].maxelm] = mu;
			   f[iflow].prop_b[BETA_T][iG - f[iflow].maxelm] = beta_t;

		   } // G Side
	   }

	   if (f[iflow].neighbors_for_the_internal_node[G][2] != nullptr) {
		   iG = f[iflow].neighbors_for_the_internal_node[G][2][iP];
		   if (iG >= f[iflow].maxelm) {
			   // ��� ��������� ����
			   rho = 1.1614; mu = 1.84e-5; beta_t = 0.003331; // ������������� default dry air 300K 1atm properties
			   if (matlist[b[ib].imatid].blibmat == 1) {
				   doublereal cp, lam;
				   cp = 1005; lam = 0.025;
				   // ������������ ����������� ������ ��������� AliceFlow ��������
				   doublereal temperature = 20.0;
				   if (t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]] >= t.maxelm) {
					   // ��������� ��� ����������� �����
					   temperature = t.potent[t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]]];
				   }
				   else {
					   if (t.neighbors_for_the_internal_node[G][2][f[iflow].ptr[iP]] > -1) {
						   // ����������� �� ����� ����������������� �������� �������������:
						   TOCHKA p1, p2;
						   center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1, 100); // ���������� ��������� ������ ��.
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
				   // ������������ ����������� ������ ��������� AliceFlow ��������
				   if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					   my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				   }
			   }
			   else if (matlist[b[ib].imatid].blibmat == 0) {
				   // �������� ����������� �������������:
				   // ���������� ��������.

				   // ���������
				   rho = matlist[b[ib].imatid].rho;
				   // ������������ ��������
				   mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				   // ����������� ��������� �������������� ����������.
				   beta_t = matlist[b[ib].imatid].beta_t;
				   // ��������� ���������� �������� ��������� �� ���������� ����� �������
				   // � ��������� ��������� �����.
			   }
			   // �������� ���������� ������������ ������:
			   f[iflow].prop_b[RHO][iG - f[iflow].maxelm] = rho;
			   f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG - f[iflow].maxelm] = mu;
			   f[iflow].prop_b[BETA_T][iG - f[iflow].maxelm] = beta_t;

		   } // G Side
	   }

	   if (f[iflow].neighbors_for_the_internal_node[G][3] != nullptr) {
		   iG = f[iflow].neighbors_for_the_internal_node[G][3][iP];
		   if (iG >= f[iflow].maxelm) {
			   // ��� ��������� ����
			   rho = 1.1614; mu = 1.84e-5; beta_t = 0.003331; // ������������� default dry air 300K 1atm properties
			   if (matlist[b[ib].imatid].blibmat == 1) {
				   doublereal cp, lam;
				   cp = 1005; lam = 0.025;
				   // ������������ ����������� ������ ��������� AliceFlow ��������
				   doublereal temperature = 20.0;
				   if (t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]] >= t.maxelm) {
					   // ��������� ��� ����������� �����
					   temperature = t.potent[t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]]];
				   }
				   else {
					   if (t.neighbors_for_the_internal_node[G][3][f[iflow].ptr[iP]] > -1) {
						   // ����������� �� ����� ����������������� �������� �������������:
						   TOCHKA p1, p2;
						   center_cord3D(f[iflow].ptr[iP], t.nvtx, t.pa, p1, 100); // ���������� ��������� ������ ��.
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
				   // ������������ ����������� ������ ��������� AliceFlow ��������
				   if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
					   my_fluid_properties(temperature, f[iflow].potent[PRESS][iG], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				   }
			   }
			   else if (matlist[b[ib].imatid].blibmat == 0) {
				   // �������� ����������� �������������:
				   // ���������� ��������.

				   // ���������
				   rho = matlist[b[ib].imatid].rho;
				   // ������������ ��������
				   mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[iflow].SInvariantStrainRateTensor[iG]);
				   // ����������� ��������� �������������� ����������.
				   beta_t = matlist[b[ib].imatid].beta_t;
				   // ��������� ���������� �������� ��������� �� ���������� ����� �������
				   // � ��������� ��������� �����.
			   }
			   // �������� ���������� ������������ ������:
			   f[iflow].prop_b[RHO][iG - f[iflow].maxelm] = rho;
			   f[iflow].prop_b[MU_DYNAMIC_VISCOSITY][iG - f[iflow].maxelm] = mu;
			   f[iflow].prop_b[BETA_T][iG - f[iflow].maxelm] = beta_t;

		   } // G Side
	   }
} // gran_prop_flow

// ���������� ������� ���������� ��� �������������
void update_flow_properties(TEMPER &t, FLOW* &f, BLOCK* b, integer lb, integer flow_interior, TPROP* matlist, bool bfirst_start) {

	// ���� ������ ��� ��� �� �������� ���������� ������� ���������� ����������� �� �����.
	for (integer ifi=0; ifi<flow_interior; ifi++) {
		// ������ �� ���� ������ �����.
		
		//TOCHKA p; // ����� ����� ���������������� ��.
		
		doublereal dmin = 1.0e30;
		doublereal dmax = -1.0e30;

		for (integer iP=0; iP<f[ifi].maxelm; iP++) {


			integer ib; // ����� ����� �������� ����������� ��������������� ����������� �����.
			doublereal rho, mu, beta_t;

			// ������ �� ���� ���������� ����������� �������
			//center_cord3D(iP, f[ifi].nvtx, f[ifi].pa, p,100); // ���������� ��������� ������ ��.
			//in_model_flow(p, ib, b, lb); // ���������� ����� ����� ib �������� ����������� ����������� ����� � ������� iP.
			// ����� ������� �������� - �������������� � hash �������.
			ib=t.whot_is_block[f[ifi].ptr[iP]];

			rho=1.1614; mu=1.84e-5; beta_t=0.003331; // ������������� default dry air 300K 1atm properties

			if (ib > -1) {
				if (matlist[b[ib].imatid].blibmat == 1) {
					doublereal cp, lam;
					cp = 1005; lam = 0.025;
					// ������������ ����������� ������ ��������� AliceFlow ��������
					if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						//printf("rho=%e, mu=%e",rho,mu); // debug
						//getchar();
						my_fluid_properties(t.potent[f[ifi].ptr[iP]], f[ifi].potent[PRESS][iP], rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
						//printf("rho=%e, mu=%e",rho,mu); // debug
						//getchar();
					}
				}
				else if (matlist[b[ib].imatid].blibmat == 0) {
					// �������� ����������� �������������:
					// ���������� ��������.
					rho = matlist[b[ib].imatid].rho;
					// ���������� ������������ ��������:
					// � ��� ����� � ��� �������������� ���������.
					mu = return_dynamic_viscosity(matlist, b[ib].imatid, bfirst_start, f[ifi].SInvariantStrainRateTensor[iP]);
					// ����������� ��������� �������������� ����������
					beta_t = matlist[b[ib].imatid].beta_t;
				}
			}
			// �������� ����������� ������������ ������:
			if (mu > dmax) dmax = mu;
			if (mu < dmin) dmin = mu;
			f[ifi].prop[RHO][iP]=rho;
			f[ifi].prop[MU_DYNAMIC_VISCOSITY][iP]=mu;
			f[ifi].prop[BETA_T][iP]=beta_t;
			//printf("rho=%e, mu=%e",f[ifi].prop[RHO][iP],f[ifi].prop[MU][iP]); // debug
			//getchar();

			// ������ ��������� ���������� ��������� ����������� ������,
			// ������� ����������� � ������ ���������� ��.
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

// ������������ ������� ����������� ��������� ������������ � .stl �������.
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
				// �� ������� �������

				// ����� ��� ����������� �����
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
				// ������ ������� �������.

				// ����� ��� ����������� �����
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
				// �������� �� ������� �������.

				// ����� ��� ����������� �����
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
				// �������� ������ ������� �������.

				// ����� ��� ����������� �����
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
				// �������� �� ������� �������.

				// ����� ��� ����������� �����
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
				// �������� ������ ������� �������.

				// ����� ��� ����������� �����
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

// ��������� ����� �������������� ������.
// 12.03.2017
doublereal massa_cabinet(TEMPER &t, FLOW* &f, 
	integer &flow_interior,
	BLOCK* b, integer lb, doublereal temp_ref,
	TPROP* matlist) {

	doublereal massa = 0.0;

	for (integer iP = 0; iP < t.maxelm; iP++) {
		// ���������� �������� �������� ������������ ������:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
		//volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
		volume3D_q(iP, t.nvtx, t.pa, dx, dy, dz);
		integer ib = t.whot_is_block[iP];
		doublereal rho, cp, lam;
		rho = 1.1614; cp = 1005; lam = 0.025; // ������������� default  dry air 300K 1atm properties
		if (matlist[b[ib].imatid].blibmat == 1) {
			// ������������, ����������� ������ ��������� AliceFlow ��������.
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				my_solid_properties(t.potent[iP], rho, cp, lam, matlist[b[ib].imatid].ilibident);
				// �������� �� ������������ ����������.
				diagnostic_critical_temperature(t.potent[iP],f,t,b,lb);
			} // SOLID
			if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				doublereal mu, beta_t; // �������� �� ������������ �� ���������.
				doublereal pressure;
				if (t.ptr[1][iP] == -1) {
					pressure = 0.0; // �������� ������ ������� ���� (����� �� ����� ����, �.�. ����� ����������� ��������).
				}
				else pressure = f[t.ptr[1][iP]].potent[PRESS][t.ptr[0][iP]];
				my_fluid_properties(t.potent[iP], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
			} // FLUID
		}
		else if (matlist[b[ib].imatid].blibmat == 0) {
			// �������� ����������� �������������:
			// ���������� ��������.
			rho = matlist[b[ib].imatid].rho;
			//cp=matlist[b[ib].imatid].cp;
			//lam=matlist[b[ib].imatid].lam;
			cp = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, t.potent[iP]);
			lam = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, t.potent[iP]);

		}
		// �������� ��� ����������� ������������ ������.
		massa += rho*dx*dy*dz;
	}

	//printf("massa=%1.3f kg\n", massa);
	std::cout << "massa=" << massa << " kg" << std::endl;

	return massa;
}


// �������� ������ ����� �������� ������� ������.
// �������� ������� �� ����� �������� ��������� ������� ����� �������� ������� ������.
// 28.10.2019
void report_out_boundary(FLOW &f, TEMPER &t, integer ls, integer lw, WALL* &w, BLOCK* &b, integer lb,
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
				integer ib; // ����� �������� �����
				in_model_temp(p, ib, b, lb);

				doublereal lam= t.prop[LAM][iP]; // �������� �� ������������ �� ���������
				doublereal temperature_i = t.potent[iP]; // �� �� ����� ���� �������� ��������� � ����������� ���������� ����.
				doublereal temperature_ii = temperature_i;
				if (t.border_neighbor[j].iII > -1) {
					temperature_ii = t.potent[t.border_neighbor[j].iII];
				}
				doublereal temperature_w=t.potent[t.border_neighbor[j].iB];

				switch (w[iwall_scan].iPlane) {
				case XY_PLANE: if (t.border_neighbor[j].Norm == T_SIDE) {// ���, ���������� ������.
					//+ �������
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dy) / (0.5 * (dz + dz1))) >
						fabs((lam * (temperature_i-temperature_w) * dx * dy) / (0.5 * (dz)))) {
					wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dy) / (0.5 * (dz + dz1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dy) / (0.5 * (dz));
					}
				}
				if (t.border_neighbor[j].Norm == B_SIDE) {// ����, ���������� ������.
					//+ �������
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
				if (t.border_neighbor[j].Norm == N_SIDE) {// ��, ���������� ������.
					//+ �������
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1)))>
						fabs((lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy)))) {
					    wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy));
					}
				}
				if (t.border_neighbor[j].Norm == S_SIDE) {// �����, ���������� ������.
					//+ �������
					if (fabs((lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1))) >
						fabs((lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy)))) {
					wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dx * dz) / (0.5 * (dy + dy1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dx * dz) / (0.5 * (dy));
					}
				}
				break;
				case YZ_PLANE:  if (t.border_neighbor[j].Norm == E_SIDE) {// �����, ���������� ������.
					//+ �������
					if (fabs((lam * (temperature_ii - temperature_i) * dy * dz) / (0.5 * (dx + dx1)))>
						fabs((lam * (temperature_i-temperature_w) * dy * dz) / (0.5 * (dx)))) {
						wall_power[iwall_scan] += (lam * (temperature_ii - temperature_i) * dy * dz) / (0.5 * (dx + dx1));
					}
					else {
						wall_power[iwall_scan] += (lam * (temperature_i-temperature_w) * dy * dz) / (0.5 * (dx));
					}
				}
				if (t.border_neighbor[j].Norm == W_SIDE) {// ������, ���������� ������.
					//+ �������
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
		// ���������� ����������� �������� �����������.
		if ((!w[iwall_scan].bpressure) && (!w[iwall_scan].bsymmetry)) {
			if (w[iwall_scan].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
				Tamb0 = fmin(Tamb0, w[iwall_scan].Tamb);
			}
		}
	}
	if (Tamb0 > 1.0e10) Tamb0 = Tamb;

	doublereal cp=0.0;

	for (int iwall_scan = 0; iwall_scan < lw; iwall_scan++) {

		doublereal Tmax_wall = -1.0e30;

		//if (iwall_scan == 0) printf("\n");

		if (w[iwall_scan].bpressure) {
			doublereal rashod = 0.0; // m!3/s
			doublereal rashod2 = 0.0; // kg/s
									  // �������� � �� ������� �������� ����� �������� ������� ������.
									  // ������� ����� ����� �.�. ��������� ����� �������� ����������
									  // ���� ������ ������� ������ ������ ������� ������� �� �����������.
			doublereal Qout = 0.0; // ��
			for (integer j = 0; j < f.maxbound; j++) {

				if (f.border_neighbor[j].MCB == (ls + iwall_scan)) {

				integer iP = f.border_neighbor[j].iI;// f.maxelm + j;

				TOCHKA p;
				center_cord3D(iP, f.nvtx, f.pa, p, 100);
				integer ib; // ����� �������� �����
				in_model_flow(p, ib, b, lb);

				doublereal rho, mu, beta_t, lam; // �������� �� ������������ �� ���������
				doublereal pressure = f.potent[PRESS][iP]; // �� �� ����� ���� �������� ��������� � ����������� ���������� ����.
				my_fluid_properties(t.potent[f.ptr[iP]], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				Tmax_wall = fmax(Tmax_wall, t.potent[f.ptr[iP]]);

				// ��������!!! f.prop_b[HEAT_CAPACITY][j] ������������ ������, �.�. ��� fluid ��� �� ����������.

				
					switch (w[iwall_scan].iPlane) {
					case XY_PLANE:
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
									  // �������� � �� ������� �������� ����� �������� ������� ������.
									  // ������� ����� ����� �.�. ��������� ����� �������� ����������
									  // ���� ������ ������� ������ ������ ������� ������� �� �����������.
			doublereal Qout = 0.0; // ��
			for (integer j = 0; j < f.maxbound; j++) {

				if (f.border_neighbor[j].MCB == (ls + iwall_scan)) {

				integer iP = f.border_neighbor[j].iI;// f.maxelm + j;

				

				TOCHKA p;
				center_cord3D(iP, f.nvtx, f.pa, p, 100);
				integer ib; // ����� �������� �����
				in_model_flow(p, ib, b, lb);

				doublereal rho, mu, beta_t, lam; // �������� �� ������������ �� ���������
				doublereal pressure = f.potent[PRESS][iP]; // �� �� ����� ���� �������� ��������� � ����������� ���������� ����.
				my_fluid_properties(t.potent[f.ptr[iP]], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
				Tmax_wall = fmax(Tmax_wall, t.potent[f.ptr[iP]]);

				// 1 - minx, 2-maxx; 3 - miny, 4-maxy; 5 - minz, 6-maxz;
				//if ((iwall_scan == 2)||(iwall_scan == 5)) {
					//printf("iP=%ld rho=%e cp=%e ptr=%d T=%e iPlane=%d Norm=%d\n", iP,rho,cp, f.ptr[iP], t.potent[f.ptr[iP]], w[iwall_scan].iPlane, f.border_neighbor[j].Norm);
					//printf("xE=%e xS=%e\n", w[iwall_scan].g.xE, w[iwall_scan].g.xS);
					//getchar();
				//}

				// ��������!!! f.prop_b[HEAT_CAPACITY][j] ������������ ������, �.�. ��� fluid ��� �� ����������.

				//printf("HEAT_CAPACITY=%e\n", cp); // debug
				//system("pause"); // debug
				
					switch (w[iwall_scan].iPlane) {
					case XY_PLANE:
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
			// ������� ������� ������.

			doublereal rashod = 0.0; // m!3/s
			doublereal rashod2 = 0.0; // kg/s
									  // �������� � �� ������� �������� ����� �������� ������� ������.
									  // ������� ����� ����� �.�. ��������� ����� �������� ����������
									  // ���� ������ ������� ������ ������ ������� ������� �� �����������.
			doublereal Qout = 0.0; // ��
			for (integer j = 0; j < f.maxbound; j++) {

				if (f.border_neighbor[j].MCB == (ls + iwall_scan)) {

					integer iP = f.border_neighbor[j].iI;// f.maxelm + j;

					TOCHKA p;
					center_cord3D(iP, f.nvtx, f.pa, p, 100);
					integer ib; // ����� �������� �����
					in_model_flow(p, ib, b, lb);

					doublereal rho, mu, beta_t, lam; // �������� �� ������������ �� ���������
					doublereal pressure = f.potent[PRESS][iP]; // �� �� ����� ���� �������� ��������� � ����������� ���������� ����.
					my_fluid_properties(t.potent[f.ptr[iP]], pressure, rho, cp, lam, mu, beta_t, matlist[b[ib].imatid].ilibident);
					Tmax_wall = fmax(Tmax_wall, t.potent[f.ptr[iP]]);

					// ��������!!! f.prop_b[HEAT_CAPACITY][j] ������������ ������, �.�. ��� fluid ��� �� ����������.


					switch (w[iwall_scan].iPlane) {
					case XY_PLANE:
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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
						// ������� ����������
						// ������: �� ��� �������� �� � ������.
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